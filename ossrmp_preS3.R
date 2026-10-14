# Compute four plug-in trace criteria for prediction-oriented OSSRMP.
# T1/T2 use calibrated final equations; T3/T4 use uncalibrated equations.
# T1/T4 use the residual-score design; T2/T3 use the raw-score design.
# This file is loaded with source("ossrmp_preS3.R") from S3SM.R.

suppressPackageStartupMessages({
  library(nnet)
})

ossrmp_preS3 <- function(X, y, n0, n_dag, rho = 0.999,
                         pilot_sampling = "proportional",
                         seed = NULL, tol = 1e-8, maxit = 50,
                         add_pilot_variance = FALSE,
                         add_full_variance = FALSE,
                         stable_bread = FALSE,
                         verbose = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  X <- as.matrix(X)
  storage.mode(X) <- "double"

  X.all <- X
  N.all <- nrow(X.all)
  p <- ncol(X.all)

  pilot_sampling <- match.arg(pilot_sampling, c("proportional", "uniform"))

  if (!is.factor(y)) {
    lev <- sort(unique(y))
    y <- factor(as.character(y), levels = as.character(lev))
  }

  cls <- levels(y)
  K1 <- length(cls)
  K  <- K1 - 1
  d  <- p + 1

  if (length(y) != N.all || anyNA(y) || any(!is.finite(X))) {
    stop("X and y must have matching sample sizes and contain no missing or non-finite values.")
  }
  if (K1 < 2L) stop("At least two classes are required.")
  if (length(n0) != 1L || !is.finite(n0) || n0 != floor(n0) || n0 < K1) {
    stop("n0 must be an integer at least as large as the number of classes.")
  }
  if (n0 >= N.all) stop("n0 must be smaller than the full sample size.")
  stopifnot(n_dag >= 1, n_dag == floor(n_dag), maxit >= 1, maxit == floor(maxit),
            rho > 0, rho < 1, tol > 0)
  counts <- table(factor(y, levels = cls))
  if (any(counts == 0)) stop("Every factor level must have at least one observation.")

  # Probability, score, curvature, and prediction derivatives.
  prob_one <- function(xa, beta) {
    eta <- c(0, as.numeric(xa %*% beta))
    eta <- eta - max(eta)
    ee  <- exp(eta)
    ee / sum(ee)
  }

  score_one <- function(xa, yb, beta) {
    pp <- prob_one(xa, beta)
    outer(xa, pp[-1] - yb)
  }

  bread_one <- function(xa, beta) {
    pp <- prob_one(xa, beta)
    pb <- pp[-1]
    Vp <- diag(pb, nrow = K) - tcrossprod(pb)
    xx <- tcrossprod(xa)
    kronecker(Vp, xx)
  }

  prediction_one <- function(xa, pfull) {
    pb <- pfull[-1]
    S <- -tcrossprod(pfull, pb)
    S[cbind(seq.int(2, K1), seq_len(K))] <- S[cbind(seq.int(2, K1), seq_len(K))] + pfull[-1]
    xx <- tcrossprod(xa)
    kronecker(crossprod(S), xx)
  }

  one_hot <- function(yv, cls) {
    yy <- as.character(yv)
    id <- match(yy, cls)
    if (any(is.na(id))) stop("Some class labels are not found in cls.")
    Ymat <- matrix(0, nrow = length(yy), ncol = length(cls))
    colnames(Ymat) <- cls
    Ymat[cbind(seq_along(yy), id)] <- 1
    Ymat
  }

  safe_solve <- function(A, B = NULL, ridge = TRUE) {
    if (ridge) {
      eps <- 1e-8 * max(1, mean(diag(A)))
      A <- A + diag(eps, ncol(A))
    }
    if (is.null(B)) {
      tryCatch(solve(A), error = function(e) qr.solve(A))
    } else {
      tryCatch(solve(A, B), error = function(e) qr.solve(A, B))
    }
  }

  # nnet returns a coefficient vector for binary responses.
  pilot_coefficients <- function(fit) {
    cf <- coef(fit)
    if (is.null(dim(cf))) {
      cf <- matrix(cf, nrow = 1, dimnames = list(cls[-1], names(cf)))
    }
    columns <- c("(Intercept)", paste0("x", seq_len(p)))
    t(cf[cls[-1], columns, drop = FALSE])
  }

  prediction_score <- function(v, M.inv, SigmaP) {
    vv  <- M.inv %*% v
    val <- drop(t(vv) %*% SigmaP %*% vv)
    sqrt(max(val, 0))
  }

  sampling_scores <- function(scr) {
    scr <- as.numeric(scr)
    scr[!is.finite(scr) | scr < 0] <- 0
    pos <- scr[scr > 0]
    if (length(pos) == 0L) return(rep(1, length(scr)))
    eps <- max(1e-12, min(pos) * 1e-8)
    scr[scr <= 0 | !is.finite(scr)] <- eps
    scr
  }

  # 1. Draw a pilot sample containing every class.
  if (pilot_sampling == "uniform") {
    pilot.prob <- rep(1 / N.all, N.all)
  } else {
    # Proportional mode assigns equal total sampling weight to each class.
    pilot.prob <- 1 / (K1 * as.numeric(counts[as.character(y)]))
    pilot.prob <- pilot.prob / sum(pilot.prob)
  }

  idx0 <- NULL

  for (b in 1:100) {
    idx.try <- sample.int(N.all, size = n0, replace = FALSE, prob = pilot.prob)

    if (length(unique(y[idx.try])) == K1) {
      idx0 <- idx.try
      break
    }
  }

  if (is.null(idx0)) {
    stop("Pilot subsample repeatedly missed some classes; increase n0.")
  }

  idxR <- setdiff(seq_len(N.all), idx0)

  X0 <- X.all[idx0, , drop = FALSE]
  XR <- X.all[idxR, , drop = FALSE]

  y0 <- factor(y[idx0], levels = cls)
  yR <- factor(y[idxR], levels = cls)

  N <- nrow(XR)

  dat0 <- data.frame(y = y0, X0, check.names = FALSE)
  colnames(dat0) <- c("y", paste0("x", seq_len(p)))

  if (pilot_sampling == "uniform") {
    w0.fit <- rep(1, length(idx0))
  } else {
    w0.fit <- 1 / pilot.prob[idx0]
    w0.fit <- w0.fit / mean(w0.fit)
  }

  # Pilot fitting weights retain the original normalization.
  fit0 <- nnet::multinom(y ~ ., data = dat0, weights = w0.fit, trace = FALSE,
                         MaxNWts = max(10000, 10 * (p + 1) * K))

  beta0 <- pilot_coefficients(fit0)

  # 2. Evaluate the pilot probabilities on the remaining observations.
  X <- XR
  y <- yR
  Xaug  <- cbind(1, X)
  Xaug0 <- cbind(1, X0)

  eta0 <- Xaug %*% beta0
  eta0 <- cbind(0, eta0)
  eta0 <- eta0 - apply(eta0, 1, max)
  ee0  <- exp(eta0)
  P0   <- ee0 / rowSums(ee0)
  colnames(P0) <- cls

  Y.full <- one_hot(y, cls)
  Y.b <- Y.full[, -1, drop = FALSE]
  Y0.full <- one_hot(y0, cls)
  Yb0 <- Y0.full[, -1, drop = FALSE]

  eta0.pilot <- Xaug0 %*% beta0
  eta0.pilot <- cbind(0, eta0.pilot)
  eta0.pilot <- eta0.pilot - apply(eta0.pilot, 1, max)
  ee0.pilot  <- exp(eta0.pilot)
  P0.pilot   <- ee0.pilot / rowSums(ee0.pilot)
  colnames(P0.pilot) <- cls

  # 3. Estimate the prediction derivative matrix and pilot curvature.
  SigmaP.hat <- matrix(0, nrow = d * K, ncol = d * K)
  Mhat.plugin <- matrix(0, nrow = d * K, ncol = d * K)

  w0.design <- w0.fit / sum(w0.fit)

  for (j in seq_len(n0)) {
    SigmaP.hat <- SigmaP.hat +
      w0.design[j] * prediction_one(Xaug0[j, ], P0.pilot[j, ])

    Mhat.plugin <- Mhat.plugin +
      w0.design[j] * bread_one(Xaug0[j, ], beta0)
  }

  Mhat.plugin.inv <- safe_solve(Mhat.plugin, ridge = FALSE)

  # 4. Partition the remaining data and summarize low-resolution regions.
  Omega <- setNames(lapply(seq_along(cls), function(k) {
    which(y == cls[k] & P0[, k] >= rho)
  }), cls)
  Omega.dag <- setdiff(seq_len(N), unlist(Omega))
  N.dag <- length(Omega.dag)
  if (N.dag == 0) stop("Omega.dag is empty. Try a smaller rho.")

  Nk <- vapply(Omega, length, numeric(1))
  xbar <- lapply(Omega, function(ii) {
    if (length(ii)) c(1, colMeans(X[ii, , drop = FALSE])) else rep(NA_real_, d)
  })

  Xaug.dag <- Xaug[Omega.dag, , drop = FALSE]
  Yb.dag   <- Y.b[Omega.dag, , drop = FALSE]

  # 5. Construct calibrated and uncalibrated prediction-oriented designs.
  psi.pilot <- matrix(0, nrow = N.dag, ncol = d * K)
  for (i in seq_len(N.dag)) {
    psi.pilot[i, ] <- as.vector(score_one(Xaug.dag[i, ], Yb.dag[i, ], beta0))
  }

  g.dag <- cbind(1, Yb.dag, psi.pilot)
  Gg <- crossprod(g.dag) / N.dag
  Gg.inv <- safe_solve(Gg, ridge = TRUE)
  A.pilot <- (t(psi.pilot) %*% g.dag / N.dag) %*% Gg.inv
  E.pilot <- psi.pilot - g.dag %*% t(A.pilot)

  score.P.raw <- apply(E.pilot, 1, prediction_score,
                       M.inv = Mhat.plugin.inv, SigmaP = SigmaP.hat)
  score.PO.raw <- apply(psi.pilot, 1, prediction_score,
                        M.inv = Mhat.plugin.inv, SigmaP = SigmaP.hat)

  lambda.P  <- sampling_scores(score.P.raw)
  lambda.PO <- sampling_scores(score.PO.raw)
  pi.P  <- lambda.P / sum(lambda.P)
  pi.PO <- lambda.PO / sum(lambda.PO)

  bread_base <- function(beta) {
    M <- matrix(0, nrow = d * K, ncol = d * K)
    for (j in seq_len(n0)) {
      M <- M + bread_one(Xaug0[j, ], beta) / (n0 + N)
    }
    for (k in seq_len(K1)) {
      if (Nk[k] == 0) next
      M <- M + (Nk[k] / (n0 + N)) * bread_one(xbar[[k]], beta)
    }
    M
  }

  # 6. Solve one estimating equation and calculate its plug-in trace.
  run_one <- function(idx.sel.pos, pi.dag, final_calibrated, tag) {
    idx.sel <- Omega.dag[idx.sel.pos]
    pi.sel  <- pi.dag[idx.sel.pos]

    if (final_calibrated) {
      gbar   <- colMeans(g.dag)
      gsel   <- g.dag[idx.sel.pos, , drop = FALSE]
      gtilde <- colSums(gsel / pi.sel) / (N.dag * n_dag)
      Mg     <- crossprod(sqrt(1 / pi.sel) * gsel) / (N.dag * n_dag)
      delta.g <- gtilde - gbar
      v.rb <- safe_solve(Mg, delta.g, ridge = TRUE)
      w.sel <- 1 - as.vector(gsel %*% v.rb)
    } else {
      v.rb <- NULL
      w.sel <- rep(1, n_dag)
    }

    evaluate_equation <- function(theta, add_jac = TRUE) {
      beta <- matrix(theta, nrow = d, ncol = K)
      G <- matrix(0, nrow = d, ncol = K)
      if (add_jac) J <- matrix(0, nrow = d * K, ncol = d * K)

      for (j in seq_len(n0)) {
        Sj <- score_one(Xaug0[j, ], Yb0[j, ], beta)
        G  <- G + Sj / (n0 + N)
        if (add_jac) J <- J + bread_one(Xaug0[j, ], beta) / (n0 + N)
      }

      for (k in seq_len(K1)) {
        if (Nk[k] == 0) next
        pp <- prob_one(xbar[[k]], beta)
        pb <- pp[-1]
        ek <- rep(0, K)
        if (k > 1) ek[k - 1] <- 1
        G <- G + (Nk[k] / (n0 + N)) * outer(xbar[[k]], pb - ek)
        if (add_jac) J <- J + (Nk[k] / (n0 + N)) * bread_one(xbar[[k]], beta)
      }

      for (t in seq_len(n_dag)) {
        St <- score_one(Xaug[idx.sel[t], ], Y.b[idx.sel[t], ], beta)
        cc <- w.sel[t] / ((n0 + N) * n_dag * pi.sel[t])
        G <- G + cc * St
        if (add_jac) J <- J + cc * bread_one(Xaug[idx.sel[t], ], beta)
      }

      if (add_jac) list(G = G, J = J) else G
    }

    # Newton iterations with gradient-norm backtracking.
    theta <- as.vector(beta0)
    for (it in seq_len(maxit)) {
      out.it <- evaluate_equation(theta, add_jac = TRUE)
      G <- out.it$G
      J <- out.it$J
      gvec  <- as.vector(G)
      gnorm <- sqrt(sum(gvec^2))

      if (verbose) cat("[", tag, "] iter =", it,
                       " |Psi| =", format(gnorm, scientific = TRUE), "\n")
      if (!is.finite(gnorm)) stop("Non-finite gradient norm encountered.")
      if (gnorm < tol) {
        break
      }

      step <- tryCatch(qr.solve(J, gvec), error = function(e) solve(J, gvec))
      alpha <- 1
      ok <- FALSE
      for (ls in 0:20) {
        theta.new <- theta - alpha * step
        G.new <- evaluate_equation(theta.new, add_jac = FALSE)
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

    beta <- matrix(theta, nrow = d, ncol = K)

    # Sandwich covariance and optional pilot/full-data contributions.
    Mhat <- bread_base(beta)
    for (t in seq_len(n_dag)) {
      Ht <- bread_one(Xaug[idx.sel[t], ], beta)
      if (stable_bread) {
        cc <- 1 / ((n0 + N) * n_dag * pi.sel[t])
      } else {
        cc <- w.sel[t] / ((n0 + N) * n_dag * pi.sel[t])
      }
      Mhat <- Mhat + cc * Ht
    }

    psi.sel <- matrix(0, nrow = n_dag, ncol = d * K)
    for (t in seq_len(n_dag)) {
      psi.sel[t, ] <- as.vector(score_one(Xaug[idx.sel[t], ], Y.b[idx.sel[t], ], beta))
    }

    if (final_calibrated) {
      gsel <- g.dag[idx.sel.pos, , drop = FALSE]
      S_pg <- t(psi.sel) %*% ((1 / pi.sel) * gsel)
      S_gg <- crossprod(sqrt(1 / pi.sel) * gsel)
      Ahat <- S_pg %*% safe_solve(S_gg, ridge = TRUE)
      Esel <- psi.sel - gsel %*% t(Ahat)
    } else {
      Esel <- psi.sel
    }

    Dhat <- matrix(0, nrow = d * K, ncol = d * K)
    for (t in seq_len(n_dag)) {
      et <- matrix(Esel[t, ], ncol = 1)
      Dhat <- Dhat + (et %*% t(et)) / (pi.sel[t]^2)
    }
    Dhat <- Dhat / ((n0 + N)^2 * n_dag^2)

    if (add_pilot_variance) {
      psi0.mat <- matrix(0, nrow = n0, ncol = d * K)
      for (j in seq_len(n0)) {
        psi0.mat[j, ] <- as.vector(score_one(Xaug0[j, ], Yb0[j, ], beta))
      }
      psi0.center <- sweep(psi0.mat, 2, colMeans(psi0.mat), "-")
      Dhat.pilot <- (n0 / (n0 + N)^2) * crossprod(psi0.center) / max(n0 - 1, 1)
      Dhat <- Dhat + Dhat.pilot
    }

    if (add_full_variance) {
      psiR.mat <- matrix(0, nrow = N, ncol = d * K)
      for (i in seq_len(N)) {
        psiR.mat[i, ] <- as.vector(score_one(Xaug[i, ], Y.b[i, ], beta))
      }
      psiR.center <- sweep(psiR.mat, 2, colMeans(psiR.mat), "-")
      Dhat.full <- (N / (n0 + N)^2) * crossprod(psiR.center) / max(N - 1, 1)
      Dhat <- Dhat + Dhat.full
    }

    Mhat <- (Mhat + t(Mhat)) / 2
    Dhat <- (Dhat + t(Dhat)) / 2
    Mhat.inv <- tryCatch(solve(Mhat), error = function(e) qr.solve(Mhat))
    cov.beta <- Mhat.inv %*% Dhat %*% Mhat.inv
    cov.beta <- (cov.beta + t(cov.beta)) / 2

    sum(t(SigmaP.hat) * cov.beta)
  }

  errors <- setNames(rep(NA_character_, 4), paste0("T", 1:4))
  safe_run <- function(expr, tag) {
    tryCatch(expr, error = function(e) {
      errors[tag] <<- conditionMessage(e)
      NA_real_
    })
  }

  # 7. T1/T4 share one draw; T2/T3 share the other.
  idx.sel.pos.P  <- sample.int(N.dag, size = n_dag, replace = TRUE, prob = pi.P)
  idx.sel.pos.PO <- sample.int(N.dag, size = n_dag, replace = TRUE, prob = pi.PO)

  T1 <- safe_run(run_one(idx.sel.pos.P,  pi.P,  TRUE,  "T1"), "T1")
  T2 <- safe_run(run_one(idx.sel.pos.PO, pi.PO, TRUE,  "T2"), "T2")
  T3 <- safe_run(run_one(idx.sel.pos.PO, pi.PO, FALSE, "T3"), "T3")
  T4 <- safe_run(run_one(idx.sel.pos.P,  pi.P,  FALSE, "T4"), "T4")

  out <- c(T1 = T1, T2 = T2, T3 = T3, T4 = T4)

  attr(out, "pilot_sampling") <- pilot_sampling
  attr(out, "n0") <- n0
  attr(out, "n_dag") <- n_dag
  attr(out, "rho") <- rho
  attr(out, "error_messages") <- errors

  out
}

# Example: run separately after loading Syn.R.
# dat <- DataGenSM(n = 2000, n_class = 4, x_dist = "gaussian",
#                  cov_type = "ar1", rho = 0.5, standardize_X = FALSE, seed = 12345)
# traces <- ossrmp_preS3(dat$X, dat$y, n0 = 200, n_dag = 400,
#                        seed = 12345, maxit = 200, add_pilot_variance = TRUE)
# print(traces)
# attr(traces, "error_messages")
