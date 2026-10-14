############################################################################
############################ Our Method: OSSRMP #############################
############################################################################
ossrmp <- function(X, y, n0, n_dag, rho = 0.999, RB = TRUE,
                   pilot_sampling = "proportional", seed = NULL,
                   tol = 1e-8, maxit = 50, add_pilot_variance = TRUE,
                   add_full_variance = FALSE, stable_bread = FALSE,
                   verbose = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  X.all <- X
  N.all <- nrow(X.all)
  p <- ncol(X.all)
  if (N.all <= 1) stop("The sample size must be larger than 1.")
  if (n0 >= N.all) stop("n0 must be smaller than the training sample size.")
  if (n_dag <= 0) stop("n_dag must be positive.")
  if (!is.factor(y)) {
    lev <- sort(unique(y))
    y <- factor(as.character(y), levels = as.character(lev))
  } else {
    y <- factor(y, levels = levels(y))
  }
  cls <- levels(y)
  K1 <- length(cls)
  K <- K1 - 1
  d <- p + 1
  if (K < 1) stop("The response must contain at least two classes.")
  if (n0 < K1) stop("n0 must be at least the number of classes.")
  pilot_sampling <- match.arg(pilot_sampling, c("proportional", "uniform"))

  # Helpers: probabilities, gradients, Hessians, and calibration.
  .onehot <- function(yv, cls) {
    id <- match(as.character(yv), cls)
    if (anyNA(id)) stop("Some class labels are not found in cls.")
    Ymat <- matrix(0, length(yv), length(cls), dimnames = list(NULL, cls))
    Ymat[cbind(seq_along(yv), id)] <- 1
    Ymat
  }
  .onehot_id <- function(id, K1) {
    Ymat <- matrix(0, length(id), K1)
    Ymat[cbind(seq_along(id), id)] <- 1
    Ymat
  }
  .safe_solve <- function(A, B = NULL, ridge = TRUE) {
    A <- as.matrix(A)
    if (ridge) {
      eps <- 1e-8 * max(1, mean(abs(diag(A))))
      A <- A + diag(eps, ncol(A))
    }
    if (is.null(B)) {
      tryCatch(solve(A), error = function(e) qr.solve(A))
    } else {
      tryCatch(solve(A, B), error = function(e) qr.solve(A, B))
    }
  }
  .softmax_prob_one <- function(xa, beta) {
    eta <- c(0, as.numeric(xa %*% beta))
    ee <- exp(eta - max(eta))
    ee / sum(ee)
  }
  .softmax_prob_mat <- function(Xa, beta) {
    eta <- cbind(0, Xa %*% beta)
    ee <- exp(eta - apply(eta, 1, max))
    ee / rowSums(ee)
  }
  .score_one <- function(xa, yb, beta) {
    outer(xa, .softmax_prob_one(xa, beta)[-1] - yb)
  }
  .score_mat <- function(Xa, Yb, beta) {
    rr <- .softmax_prob_mat(Xa, beta)[, -1, drop = FALSE] - Yb
    out <- matrix(0, nrow(Xa), d * K)
    for (kk in seq_len(K)) {
      jj <- ((kk - 1) * d + 1):(kk * d)
      out[, jj] <- Xa * as.numeric(rr[, kk])
    }
    out
  }
  .bread_one <- function(xa, beta) {
    pb <- .softmax_prob_one(xa, beta)[-1]
    kronecker(diag(pb, nrow = K) - tcrossprod(pb), tcrossprod(xa))
  }
  .get_beta0 <- function(fit0, cls, p, d, K) {
    cf0 <- coef(fit0)
    if (is.vector(cf0)) {
      cf0 <- matrix(cf0, nrow = 1, dimnames = list(cls[2], names(cf0)))
    }
    cf0 <- as.matrix(cf0)
    beta0 <- matrix(0, d, K,
                    dimnames = list(c("(Intercept)", paste0("x", 1:p)), cls[-1]))
    cn <- colnames(cf0)
    if (is.null(rownames(cf0))) rownames(cf0) <- cls[-1][seq_len(nrow(cf0))]
    for (kk in seq_len(K)) {
      rr <- match(cls[-1][kk], rownames(cf0))
      if (is.na(rr) && K == 1) rr <- 1
      if (is.na(rr)) next
      if ("(Intercept)" %in% cn) beta0[1, kk] <- cf0[rr, "(Intercept)"]
      for (jj in seq_len(p)) {
        nm <- paste0("x", jj)
        if (nm %in% cn) beta0[jj + 1, kk] <- cf0[rr, nm]
      }
    }
    beta0
  }
  .make_g_mat <- function(Xa, Yb, beta_pilot) {
    cbind(1, Yb, .score_mat(Xa, Yb, beta_pilot))
  }
  .rb_weight <- function(gbar, gsel, pi.sel, N.dag, n_dag) {
    gtilde <- colSums(gsel / pi.sel) / (N.dag * n_dag)
    Mg <- crossprod(sqrt(1 / pi.sel) * gsel) / (N.dag * n_dag)
    v.rb <- .safe_solve(Mg, gtilde - gbar, ridge = TRUE)
    1 - as.vector(gsel %*% v.rb)
  }
  .accumulate_score_cov <- function(Xmat, yfac, beta, chunk_size = 50000) {
    n <- nrow(Xmat)
    s1 <- rep(0, d * K)
    s2 <- matrix(0, d * K, d * K)
    for (st in seq(1, n, by = chunk_size)) {
      id <- st:min(st + chunk_size - 1, n)
      Xa <- cbind(1, Xmat[id, , drop = FALSE])
      Yb <- .onehot(yfac[id], cls)[, -1, drop = FALSE]
      psimat <- .score_mat(Xa, Yb, beta)
      s1 <- s1 + colSums(psimat)
      s2 <- s2 + crossprod(psimat)
    }
    if (n <= 1) return(matrix(0, d * K, d * K))
    (s2 - tcrossprod(s1) / n) / (n - 1)
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
  for (bb in 1:100) {
    idx.try <- sample.int(N.all, n0, replace = FALSE, prob = pilot.prob)
    if (length(unique(y[idx.try])) == K1) {
      idx0 <- idx.try
      break
    }
  }
  if (is.null(idx0)) stop("Pilot subsample repeatedly missed some classes; increase n0.")
  idxR <- setdiff(seq_len(N.all), idx0)
  X0 <- X.all[idx0, , drop = FALSE]
  X <- X.all[idxR, , drop = FALSE]
  y0 <- factor(y[idx0], levels = cls)
  y <- factor(y[idxR], levels = cls)
  N <- nrow(X)
  dat0 <- data.frame(y = y0, X0, check.names = FALSE)
  colnames(dat0) <- c("y", paste0("x", 1:p))
  if (pilot_sampling == "uniform") {
    w0.fit <- rep(1, length(idx0))
  } else {
    w0.fit <- 1 / pilot.prob[idx0]
    w0.fit <- w0.fit / mean(w0.fit)
  }
  fit0 <- nnet::multinom(y ~ ., data = dat0, weights = w0.fit, trace = FALSE,
                         MaxNWts = max(10000, 10 * (p + 1) * K),
                         maxit = 300, decay = 1e-4)
  beta0 <- .get_beta0(fit0, cls, p, d, K)
  Xaug0 <- cbind(1, X0)
  Yb0 <- .onehot(y0, cls)[, -1, drop = FALSE]

  # 2. Scan the remaining data: partition, centroids, scores, and auxiliary mean.
  Omega <- setNames(rep(list(integer(0)), K1), cls)
  Nk <- rep(0, K1)
  xsum <- rep(list(rep(0, d)), K1)
  Omega.dag <- integer(0)
  scr.dag <- numeric(0)
  N.dag <- 0
  gbar <- gsum <- NULL
  chunk_size <- 50000
  for (st in seq(1, N, by = chunk_size)) {
    id <- st:min(st + chunk_size - 1, N)
    Xa <- cbind(1, X[id, , drop = FALSE])
    P0.block <- .softmax_prob_mat(Xa, beta0)
    colnames(P0.block) <- cls
    yid <- as.integer(y[id])
    low.flag <- P0.block[cbind(seq_along(id), yid)] >= rho
    if (any(low.flag)) {
      for (kk in seq_len(K1)) {
        pos <- which(low.flag & yid == kk)
        if (length(pos) == 0) next
        Omega[[kk]] <- c(Omega[[kk]], id[pos])
        Nk[kk] <- Nk[kk] + length(pos)
        xsum[[kk]] <- xsum[[kk]] + colSums(Xa[pos, , drop = FALSE])
      }
    }
    high.pos <- which(!low.flag)
    if (length(high.pos) > 0) {
      Omega.dag <- c(Omega.dag, id[high.pos])
      Yb.block <- .onehot_id(yid[high.pos], K1)[, -1, drop = FALSE]
      Rb.block <- Yb.block - P0.block[high.pos, -1, drop = FALSE]
      xnorm <- sqrt(rowSums(Xa[high.pos, , drop = FALSE]^2))
      rnorm <- sqrt(rowSums(Rb.block^2))
      scr.block <- xnorm * rnorm
      scr.block[!is.finite(scr.block) | scr.block < 0] <- 0
      scr.dag <- c(scr.dag, scr.block)
      if (RB == TRUE) {
        g.block <- .make_g_mat(Xa[high.pos, , drop = FALSE], Yb.block, beta0)
        if (is.null(gsum)) gsum <- rep(0, ncol(g.block))
        gsum <- gsum + colSums(g.block)
      }
      N.dag <- N.dag + length(high.pos)
    }
  }
  if (N.dag == 0) stop("Omega.dag is empty. Try a smaller rho.")
  xbar <- vector("list", K1)
  for (kk in seq_len(K1)) {
    xbar[[kk]] <- if (Nk[kk] == 0) rep(NA_real_, d) else xsum[[kk]] / Nk[kk]
  }
  if (RB == TRUE) gbar <- gsum / N.dag

  # 3. Draw a with-replacement L-optimal sample from the high-resolution region.
  pi.dag <- if (sum(scr.dag) <= 0 || any(!is.finite(scr.dag))) {
    rep(1 / N.dag, N.dag)
  } else {
    scr.dag / sum(scr.dag)
  }
  idx.sel.pos <- sample.int(N.dag, n_dag, replace = TRUE, prob = pi.dag)
  idx.sel <- Omega.dag[idx.sel.pos]
  pi.sel <- pi.dag[idx.sel.pos]
  Xaug.sel <- cbind(1, X[idx.sel, , drop = FALSE])
  Yb.sel <- .onehot(y[idx.sel], cls)[, -1, drop = FALSE]

  # 4. Compute RB weights once; use unit weights when RB = FALSE.
  Ahat <- NULL
  if (RB == TRUE) {
    gsel.pilot <- .make_g_mat(Xaug.sel, Yb.sel, beta0)
    w.sel <- .rb_weight(gbar, gsel.pilot, pi.sel, N.dag, n_dag)
  } else {
    gsel.pilot <- NULL
    w.sel <- rep(1, n_dag)
  }
  .bread_base <- function(beta) {
    M <- matrix(0, d * K, d * K)
    for (jj in seq_len(n0)) M <- M + .bread_one(Xaug0[jj, ], beta) / (n0 + N)
    for (kk in seq_len(K1)) {
      if (Nk[kk] > 0) M <- M + (Nk[kk] / (n0 + N)) * .bread_one(xbar[[kk]], beta)
    }
    M
  }
  .mr_eval <- function(theta, add_jac = TRUE) {
    beta <- matrix(theta, d, K)
    G <- matrix(0, d, K)
    if (add_jac) J <- matrix(0, d * K, d * K)
    for (jj in seq_len(n0)) {
      G <- G + .score_one(Xaug0[jj, ], Yb0[jj, ], beta) / (n0 + N)
      if (add_jac) J <- J + .bread_one(Xaug0[jj, ], beta) / (n0 + N)
    }
    for (kk in seq_len(K1)) {
      if (Nk[kk] == 0) next
      pb <- .softmax_prob_one(xbar[[kk]], beta)[-1]
      ek <- rep(0, K)
      if (kk > 1) ek[kk - 1] <- 1
      G <- G + (Nk[kk] / (n0 + N)) * outer(xbar[[kk]], pb - ek)
      if (add_jac) J <- J + (Nk[kk] / (n0 + N)) * .bread_one(xbar[[kk]], beta)
    }
    for (tt in seq_len(n_dag)) {
      cc <- w.sel[tt] / ((n0 + N) * n_dag * pi.sel[tt])
      G <- G + cc * .score_one(Xaug.sel[tt, ], Yb.sel[tt, ], beta)
      if (add_jac) J <- J + cc * .bread_one(Xaug.sel[tt, ], beta)
    }
    if (add_jac) list(G = G, J = J) else G
  }

  # 5. Solve the combined estimating equation by Newton iteration.
  theta <- as.vector(beta0)
  conv <- FALSE
  for (it in seq_len(maxit)) {
    out.it <- .mr_eval(theta, add_jac = TRUE)
    gvec <- as.vector(out.it$G)
    gnorm <- sqrt(sum(gvec^2))
    if (verbose) cat("iter =", it, " |Psi| =", format(gnorm, scientific = TRUE), "\n")
    if (!is.finite(gnorm)) stop("Non-finite gradient norm encountered.")
    if (gnorm < tol) {
      conv <- TRUE
      break
    }
    step <- .safe_solve(out.it$J, gvec, ridge = TRUE)
    if (any(!is.finite(step))) stop("Non-finite Newton step encountered.")
    alpha <- 1
    ok <- FALSE
    theta.best <- theta
    gnorm.best <- gnorm
    for (ls in 0:25) {
      theta.new <- theta - alpha * step
      G.new <- .mr_eval(theta.new, add_jac = FALSE)
      gnorm.new <- sqrt(sum(as.vector(G.new)^2))
      if (is.finite(gnorm.new) && gnorm.new < gnorm.best) {
        theta.best <- theta.new
        gnorm.best <- gnorm.new
        ok <- TRUE
      }
      if (is.finite(gnorm.new) && gnorm.new < gnorm) {
        theta <- theta.new
        ok <- TRUE
        break
      }
      alpha <- alpha / 2
    }
    if (!ok) theta <- theta - 0.25 * step else if (gnorm.best < gnorm) theta <- theta.best
  }
  beta <- matrix(theta, d, K,
                 dimnames = list(c("(Intercept)", paste0("x", 1:p)), cls[-1]))

  # 6. Estimate the sandwich covariance, with optional pilot/full-data terms.
  Mhat <- .bread_base(beta)
  for (tt in seq_len(n_dag)) {
    cc <- if (stable_bread) 1 else w.sel[tt]
    cc <- cc / ((n0 + N) * n_dag * pi.sel[tt])
    Mhat <- Mhat + cc * .bread_one(Xaug.sel[tt, ], beta)
  }
  psi.sel <- .score_mat(Xaug.sel, Yb.sel, beta)
  if (RB == TRUE) {
    S_pg <- t(psi.sel) %*% ((1 / pi.sel) * gsel.pilot)
    S_gg <- crossprod(sqrt(1 / pi.sel) * gsel.pilot)
    Ahat <- S_pg %*% .safe_solve(S_gg, ridge = TRUE)
    Esel <- psi.sel - gsel.pilot %*% t(Ahat)
  } else {
    Esel <- psi.sel
  }
  Dhat <- matrix(0, d * K, d * K)
  for (tt in seq_len(n_dag)) {
    et <- matrix(Esel[tt, ], ncol = 1)
    Dhat <- Dhat + (et %*% t(et)) / (pi.sel[tt]^2)
  }
  Dhat <- Dhat / ((n0 + N)^2 * n_dag^2)
  Dhat.pilot <- Dhat.full <- matrix(0, d * K, d * K)
  if (add_pilot_variance) {
    psi0.mat <- .score_mat(Xaug0, Yb0, beta)
    psi0.center <- sweep(psi0.mat, 2, colMeans(psi0.mat), "-")
    Dhat.pilot <- (n0 / (n0 + N)^2) * crossprod(psi0.center) / max(n0 - 1, 1)
    Dhat <- Dhat + Dhat.pilot
  }
  if (add_full_variance) {
    cov.full <- .accumulate_score_cov(X, y, beta, chunk_size)
    Dhat.full <- (N / (n0 + N)^2) * cov.full
    Dhat <- Dhat + Dhat.full
  }
  Mhat <- (Mhat + t(Mhat)) / 2
  Dhat <- (Dhat + t(Dhat)) / 2
  Mhat.inv <- .safe_solve(Mhat, ridge = TRUE)
  cov.beta <- Mhat.inv %*% Dhat %*% Mhat.inv
  cov.beta <- (cov.beta + t(cov.beta)) / 2
  se.beta <- matrix(sqrt(pmax(diag(cov.beta), 0)), d, K, dimnames = dimnames(beta))
  zvalue <- beta / se.beta
  pvalue <- 2 * pnorm(-abs(zvalue))

  # 7. Return estimates, sampling details, calibration, and convergence status.
  list(beta = beta, se.beta = se.beta, zvalue = zvalue, pvalue = pvalue,
       cov.beta = cov.beta, beta.pilot = beta0, pi.dag = pi.dag,
       idx.sel = idx.sel, idx.sel.full = idxR[idx.sel], idx.sel.pos = idx.sel.pos,
       Omega = Omega, Omega.dag = Omega.dag, Omega.dag.full = idxR[Omega.dag],
       Nk = Nk, N.rem = N, N.all = N.all, N.dag = N.dag,
       idx0 = idx0, idxR = idxR, pilot_sampling = pilot_sampling,
       pilot.prob = pilot.prob, w0.fit = w0.fit, w.sel = w.sel, Ahat = Ahat,
       Dhat.pilot = Dhat.pilot, Dhat.full = Dhat.full, RB = RB,
       add_pilot_variance = add_pilot_variance, add_full_variance = add_full_variance,
       stable_bread = stable_bread, converged = conv)
}
