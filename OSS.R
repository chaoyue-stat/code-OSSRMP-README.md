# oss.R
# rm(list = ls())
# gc()
#
# suppressPackageStartupMessages({
#   library(nnet)
# })
# source("I:/softmax/Syn.R", echo = TRUE)

############################################################################
#################### Pure Optimal Subsampling: OSS ##########################
############################################################################
OSS <- function(X, y, n0, n_sub, poisson = FALSE,
                pilot_sampling = "proportional", seed = NULL,
                tol = 1e-8, maxit = 50, add_pilot_variance = TRUE,
                add_full_variance = FALSE, verbose = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  X.all <- X
  N.all <- nrow(X)
  p <- ncol(X)
  if (!is.factor(y)) {
    lev <- sort(unique(y))
    y <- factor(as.character(y), levels = as.character(lev))
  }
  cls <- levels(y)
  K1 <- length(cls)
  K <- K1 - 1
  d <- p + 1
  pilot_sampling <- match.arg(pilot_sampling, c("proportional", "uniform"))
  if (n0 < K1) stop("n0 must be at least the number of classes.")
  if (n0 >= N.all) stop("n0 must be smaller than the training sample size.")

  # Helpers: probabilities, individual gradient/Hessian, and pilot coefficients.
  .softmax_prob <- function(xa, beta) {
    eta <- c(0, as.numeric(xa %*% beta))
    ee <- exp(eta - max(eta))
    ee / sum(ee)
  }
  .score_one <- function(xa, yb, beta) {
    outer(xa, .softmax_prob(xa, beta)[-1] - yb)
  }
  .bread_one <- function(xa, beta) {
    pb <- .softmax_prob(xa, beta)[-1]
    kronecker(diag(pb, nrow = K) - tcrossprod(pb), tcrossprod(xa))
  }
  .onehot <- function(yv, cls) {
    id <- match(as.character(yv), cls)
    if (anyNA(id)) stop("Some class labels are not found in cls.")
    Ymat <- matrix(0, length(yv), length(cls), dimnames = list(NULL, cls))
    Ymat[cbind(seq_along(yv), id)] <- 1
    Ymat
  }
  .get_beta0 <- function(fit0, cls, p, d, K) {
    cf0 <- coef(fit0)
    if (is.vector(cf0)) {
      cf0 <- matrix(cf0, nrow = 1, dimnames = list(cls[2], names(cf0)))
    }
    cf0 <- as.matrix(cf0)
    beta0 <- matrix(0, d, K,
                    dimnames = list(c("(Intercept)", paste0("x", 1:p)), cls[-1]))
    rn <- rownames(cf0)
    if (is.null(rn)) rn <- cls[-1][seq_len(nrow(cf0))]
    for (kk in seq_along(rn)) {
      if (rn[kk] %in% cls[-1]) {
        beta0[, rn[kk]] <- c(cf0[rn[kk], "(Intercept)"],
                             cf0[rn[kk], paste0("x", 1:p)])
      }
    }
    beta0
  }
  .find_H_sub <- function(scr, n_sub) {
    scr <- as.numeric(scr)
    ssum <- sum(scr)
    if (ssum <= 0 || any(!is.finite(scr))) return(Inf)
    if (n_sub * max(scr) <= ssum) return(max(scr))
    uniroot(function(H) sum(pmin(scr, H)) - n_sub * H,
            interval = c(.Machine$double.eps, max(scr)), tol = 1e-12)$root
  }
  .make_poisson_design <- function(scr, n_sub, N) {
    if (n_sub > N)
      stop("For Poisson subsampling, n_sub must be no larger than the remaining sample size.")
    scr <- as.numeric(scr)
    if (sum(scr) <= 0 || any(!is.finite(scr))) {
      H.sub <- Inf
      lam.full <- rep(1, N)
      pi.full <- rep(1 / N, N)
    } else {
      H.sub <- .find_H_sub(scr, n_sub)
      lam.full <- if (is.infinite(H.sub)) rep(1, N) else pmin(scr, H.sub)
      pi.full <- if (is.infinite(H.sub)) rep(1 / N, N) else lam.full / sum(lam.full)
    }
    p.full <- n_sub * pi.full
    if (any(!is.finite(p.full)) || any(p.full < 0) || any(p.full > 1 + 1e-10))
      stop("Invalid Poisson inclusion probabilities.")
    list(H.sub = H.sub, lam.full = lam.full, pi.full = pi.full,
         p.full = pmin(p.full, 1))
  }

  # 1. Draw and fit a pilot sample, disjoint from the remaining data.
  if (pilot_sampling == "uniform") {
    pilot.prob <- rep(1 / N.all, N.all)
  } else {
    n.cls <- table(factor(y, levels = cls))
    if (any(n.cls == 0)) stop("Some classes have zero observations.")
    pilot.prob <- 1 / (K1 * as.numeric(n.cls[as.character(y)]))
    pilot.prob <- pilot.prob / sum(pilot.prob)
  }
  idx0 <- NULL
  for (b in 1:100) {
    idx.try <- sample.int(N.all, n0, replace = FALSE, prob = pilot.prob)
    if (length(unique(y[idx.try])) == K1) {
      idx0 <- idx.try
      break
    }
  }
  if (is.null(idx0))
    stop("Pilot subsample repeatedly missed some classes; increase n0.")
  idxR <- setdiff(seq_len(N.all), idx0)
  X0 <- X.all[idx0, , drop = FALSE]
  XR <- X.all[idxR, , drop = FALSE]
  y0 <- factor(y[idx0], levels = cls)
  yR <- factor(y[idxR], levels = cls)
  N <- nrow(XR)
  if (poisson && n_sub > N)
    stop("For Poisson subsampling, n_sub must be no larger than the remaining sample size.")
  dat0 <- data.frame(y = y0, X0, check.names = FALSE)
  colnames(dat0) <- c("y", paste0("x", 1:p))
  w0.fit <- if (pilot_sampling == "uniform") rep(1, length(idx0)) else 1 / pilot.prob[idx0]
  if (pilot_sampling != "uniform") w0.fit <- w0.fit / mean(w0.fit)
  fit0 <- nnet::multinom(y ~ ., data = dat0, weights = w0.fit, trace = FALSE,
                         MaxNWts = max(10000, 10 * (p + 1) * K))
  beta0 <- .get_beta0(fit0, cls, p, d, K)

  # 2. Evaluate pilot probabilities and L-optimal sampling scores.
  X <- XR
  y <- yR
  Xaug <- cbind(1, X)
  Xaug0 <- cbind(1, X0)
  eta0 <- cbind(0, Xaug %*% beta0)
  eta0 <- eta0 - apply(eta0, 1, max)
  ee0 <- exp(eta0)
  P0 <- ee0 / rowSums(ee0)
  colnames(P0) <- cls
  Y.b <- .onehot(y, cls)[, -1, drop = FALSE]
  Yb0 <- .onehot(y0, cls)[, -1, drop = FALSE]
  R.b <- Y.b - P0[, -1, drop = FALSE]
  scr <- sqrt(rowSums(Xaug^2)) * sqrt(rowSums(R.b^2))

  # 3. Draw the second-stage sample and set its estimating-equation weights.
  if (!poisson) {
    invalid <- sum(scr) <= 0 || any(!is.finite(scr))
    pi.full <- if (invalid) rep(1 / N, N) else scr / sum(scr)
    lam.full <- if (invalid) NULL else scr
    H.sub <- p.full <- nu.full <- sel.pos <- p.sel <- NULL
    idx.sel <- sample.int(N, n_sub, replace = TRUE, prob = pi.full)
    pi.sel <- pi.full[idx.sel]
    n.sel <- n_sub
    cc.sel <- 1 / ((n0 + N) * n_sub * pi.sel)
  } else {
    des.out <- .make_poisson_design(scr, n_sub, N)
    H.sub <- des.out$H.sub
    lam.full <- des.out$lam.full
    pi.full <- des.out$pi.full
    p.full <- des.out$p.full
    nu.full <- rbinom(N, 1, prob = p.full)
    idx.sel <- sel.pos <- which(nu.full == 1L)
    p.sel <- p.full[sel.pos]
    n.sel <- length(sel.pos)
    if (n.sel == 0) stop("Poisson subsample is empty; increase n_sub or change seed.")
    pi.sel <- NULL
    cc.sel <- 1 / ((n0 + N) * p.sel)
  }

  # 4. Fit the combined pilot and IPW estimating equation.
  .bread_pilot <- function(beta) {
    M <- matrix(0, d * K, d * K)
    for (j in 1:n0) M <- M + .bread_one(Xaug0[j, ], beta) / (n0 + N)
    M
  }
  .oss_eval <- function(theta, add_jac = TRUE) {
    beta <- matrix(theta, d, K)
    G <- matrix(0, d, K)
    if (add_jac) J <- matrix(0, d * K, d * K)
    for (j in 1:n0) {
      G <- G + .score_one(Xaug0[j, ], Yb0[j, ], beta) / (n0 + N)
      if (add_jac) J <- J + .bread_one(Xaug0[j, ], beta) / (n0 + N)
    }
    for (t in 1:n.sel) {
      ii <- idx.sel[t]
      G <- G + cc.sel[t] * .score_one(Xaug[ii, ], Y.b[ii, ], beta)
      if (add_jac) J <- J + cc.sel[t] * .bread_one(Xaug[ii, ], beta)
    }
    if (add_jac) list(G = G, J = J) else G
  }
  theta <- as.vector(beta0)
  conv <- FALSE
  for (it in 1:maxit) {
    out.it <- .oss_eval(theta)
    gvec <- as.vector(out.it$G)
    gnorm <- sqrt(sum(gvec^2))
    if (verbose) cat("iter =", it, " |Psi| =", format(gnorm, scientific = TRUE), "\n")
    if (!is.finite(gnorm)) stop("Non-finite gradient norm encountered.")
    if (gnorm < tol) {
      conv <- TRUE
      break
    }
    step <- tryCatch(qr.solve(out.it$J, gvec), error = function(e) solve(out.it$J, gvec))
    alpha <- 1
    ok <- FALSE
    for (ls in 0:20) {
      theta.new <- theta - alpha * step
      G.new <- .oss_eval(theta.new, add_jac = FALSE)
      gnorm.new <- sqrt(sum(as.vector(G.new)^2))
      if (is.finite(gnorm.new) && gnorm.new < gnorm) {
        theta <- theta.new
        ok <- TRUE
        break
      }
      alpha <- alpha / 2
    }
    if (!ok) theta <- theta - step
  }
  beta <- matrix(theta, d, K, dimnames = dimnames(beta0))

  # 5. Estimate the sandwich covariance, retaining both optional contributions.
  Mhat <- .bread_pilot(beta)
  Dhat <- matrix(0, d * K, d * K)
  for (t in 1:n.sel) {
    ii <- idx.sel[t]
    Mhat <- Mhat + cc.sel[t] * .bread_one(Xaug[ii, ], beta)
    psi_t <- matrix(as.vector(.score_one(Xaug[ii, ], Y.b[ii, ], beta)), ncol = 1)
    if (poisson) {
      Dhat <- Dhat + ((1 - p.sel[t]) / p.sel[t]^2) * (psi_t %*% t(psi_t))
    } else {
      Dhat <- Dhat + (psi_t %*% t(psi_t)) / pi.sel[t]^2
    }
  }
  Dhat <- Dhat / if (poisson) (n0 + N)^2 else ((n0 + N)^2 * n_sub^2)
  .meat_data <- function(xa, yb) {
    m <- nrow(xa)
    psi <- matrix(0, m, d * K)
    for (j in 1:m) psi[j, ] <- as.vector(.score_one(xa[j, ], yb[j, ], beta))
    centered <- sweep(psi, 2, colMeans(psi), "-")
    (m / (n0 + N)^2) * crossprod(centered) / max(m - 1, 1)
  }
  Dhat.pilot <- Dhat.full <- matrix(0, d * K, d * K)
  if (add_pilot_variance) {
    Dhat.pilot <- .meat_data(Xaug0, Yb0)
    Dhat <- Dhat + Dhat.pilot
  }
  if (add_full_variance) {
    Dhat.full <- .meat_data(Xaug, Y.b)
    Dhat <- Dhat + Dhat.full
  }
  Mhat <- (Mhat + t(Mhat)) / 2
  Dhat <- (Dhat + t(Dhat)) / 2
  Mhat.inv <- tryCatch(solve(Mhat), error = function(e) qr.solve(Mhat))
  cov.beta <- Mhat.inv %*% Dhat %*% Mhat.inv
  cov.beta <- (cov.beta + t(cov.beta)) / 2
  se.beta <- matrix(sqrt(pmax(diag(cov.beta), 0)), d, K, dimnames = dimnames(beta))
  zvalue <- beta / se.beta
  pvalue <- 2 * pnorm(-abs(zvalue))
  dimnames(zvalue) <- dimnames(pvalue) <- dimnames(beta)

  list(beta = beta, se.beta = se.beta, zvalue = zvalue, pvalue = pvalue,
       cov.beta = cov.beta, beta.pilot = beta0, score = scr,
       lam.full = lam.full, pi.full = pi.full, H.sub = H.sub,
       p.full = p.full, nu.full = nu.full, idx.sel = idx.sel,
       idx.sel.full = idxR[idx.sel], n.sel = n.sel, idx0 = idx0, idxR = idxR,
       N.rem = N, N.all = N.all, pilot_sampling = pilot_sampling,
       pilot.prob = pilot.prob, w0.fit = w0.fit, Dhat.pilot = Dhat.pilot,
       Dhat.full = Dhat.full, poisson = poisson,
       add_pilot_variance = add_pilot_variance,
       add_full_variance = add_full_variance, converged = conv)
}
