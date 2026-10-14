############################################################################
######################## Data generating procedure ##########################
############################################################################
DataGenSM <- function(n, n_class = 4,
                      x_dist = c("gaussian", "t", "skew_quadratic", "gaussian_mixture"),
                      cov_type = c("ar1", "exchangeable", "identity"),
                      rho = 0.5, df_t = 4, mix_prob = 0.3, standardize_X = FALSE,
                      coef_scale = 3.0, intercept_scale = 1.0, mr_friendly = TRUE,
                      seed = NULL, use_rfast = TRUE) {
  x_dist <- match.arg(x_dist)
  cov_type <- match.arg(cov_type)
  if (!is.null(seed)) set.seed(seed)
  if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' is required.")
  if (!requireNamespace("mvtnorm", quietly = TRUE)) stop("Package 'mvtnorm' is required.")
  has_rfast <- use_rfast && requireNamespace("Rfast", quietly = TRUE)
  if (n_class < 2) stop("n_class must be at least 2.")
  if (rho < 0 || rho >= 1) stop("rho must satisfy 0 <= rho < 1.")
  if (mix_prob <= 0 || mix_prob >= 1) stop("mix_prob must be in (0, 1).")
  p <- 7
  K <- n_class - 1

  # 1. Covariate covariance/scale matrix and random-generation helpers.
  if (cov_type == "identity") {
    Sigma <- diag(1, p)
  } else if (cov_type == "ar1") {
    Sigma <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
  } else {
    Sigma <- matrix(rho, p, p)
    diag(Sigma) <- 1
  }
  mu0 <- rep(0, p)
  .rmvnorm_fast <- function(nn, mu, Sigma) {
    if (nn <= 0) return(matrix(numeric(0), 0, length(mu)))
    out <- if (has_rfast) {
      Rfast::rmvnorm(nn, mu = mu, sigma = Sigma, seed = NULL)
    } else {
      MASS::mvrnorm(n = nn, mu = mu, Sigma = Sigma)
    }
    out <- as.matrix(out)
    if (ncol(out) != length(mu)) out <- matrix(out, ncol = length(mu))
    out
  }
  .rmvt_fast <- function(nn, mu, Sigma, df_t) {
    if (nn <= 0) return(matrix(numeric(0), 0, length(mu)))
    out <- if (has_rfast) {
      Rfast::rmvt(nn, mu = mu, sigma = Sigma, v = df_t, seed = NULL)
    } else {
      mvtnorm::rmvt(n = nn, sigma = Sigma, df = df_t, delta = mu)
    }
    out <- as.matrix(out)
    if (ncol(out) != length(mu)) out <- matrix(out, ncol = length(mu))
    out
  }

  # 2. Generate covariates, with optional sample standardization.
  if (x_dist == "gaussian") {
    X <- .rmvnorm_fast(n, mu0, Sigma)
  } else if (x_dist == "t") {
    X <- .rmvt_fast(n, mu0, Sigma, df_t)
  } else if (x_dist == "skew_quadratic") {
    Z <- .rmvnorm_fast(n, mu0, Sigma)
    X <- Z + 0.6 * (Z^2 - 1)
  } else {
    # mix_prob N(1_p, Sigma) + (1 - mix_prob) N(-1_p, Sigma).
    comp <- rbinom(n, 1, mix_prob)
    n1 <- sum(comp == 1L)
    n2 <- n - n1
    X <- matrix(0, n, p)
    if (n1 > 0) X[comp == 1L, ] <- .rmvnorm_fast(n1, rep(1, p), Sigma)
    if (n2 > 0) X[comp == 0L, ] <- .rmvnorm_fast(n2, rep(-1, p), Sigma)
  }
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  colnames(X) <- paste0("x", 1:p)
  if (standardize_X) {
    X <- as.matrix(scale(X))
    storage.mode(X) <- "double"
    colnames(X) <- paste0("x", 1:p)
  }

  # 3. True coefficients; the baseline column remains zero.
  B_true <- matrix(0, p + 1, n_class,
                   dimnames = list(c("(Intercept)", paste0("x", 1:p)), as.character(0:K)))
  for (k in 2:n_class) {
    kk <- k - 1
    intercept_k <- 0.35 * cos(0.8 * kk) + 0.15 * sin(1.3 * kk)
    beta_k <- c(0.80 * cos(0.5 * kk) + 0.30,
                -0.70 * sin(0.8 * kk) + 0.25,
                0.60 * cos(1.1 * kk) - 0.35,
                -0.90 * sin(0.6 * kk) + 0.15,
                0.75 * cos(0.9 * kk) + 0.20,
                -0.65 * sin(1.0 * kk) - 0.25,
                0.55 * cos(1.2 * kk) + 0.10)
    B_true[, k] <- c(intercept_scale * intercept_k, coef_scale * beta_k)
  }
  if (mr_friendly) {
    for (k in 2:n_class) {
      kk <- k - 1
      dir_k <- c(cos(0.7 * kk), sin(0.7 * kk), cos(1.1 * kk), sin(1.1 * kk),
                 cos(1.5 * kk), sin(1.5 * kk), (-1)^kk)
      dir_k <- dir_k / sqrt(sum(dir_k^2))
      B_true[-1, k] <- B_true[-1, k] + coef_scale * dir_k
    }
  }

  # 4. Compute softmax probabilities and draw class labels.
  Xaug <- cbind(1, X)
  eta_true <- Xaug %*% B_true
  eta_true <- eta_true - apply(eta_true, 1, max)
  ee <- exp(eta_true)
  P_true <- ee / rowSums(ee)
  colnames(P_true) <- colnames(B_true)
  y <- apply(P_true, 1, function(pr) sample(colnames(P_true), size = 1, prob = pr))
  y <- factor(y, levels = colnames(P_true))

  # 5. Return data and true parameters under both coefficient constraints.
  B_true_nonbase <- B_true[, -1, drop = FALSE]
  B_true_sum <- sweep(B_true, 1, rowMeans(B_true), "-")
  B_true_sum_nonbase <- B_true_sum[, -1, drop = FALSE]
  list(data = data.frame(y = y, X, check.names = FALSE), X = X, Xaug = Xaug,
       y = y, eta_true = eta_true, P_true = P_true,
       B_true = B_true, B_true_nonbase = B_true_nonbase,
       B_true_sum = B_true_sum, B_true_sum_nonbase = B_true_sum_nonbase,
       theta_true = as.vector(B_true_nonbase), theta_true_sum = as.vector(B_true_sum_nonbase),
       Sigma = Sigma,
       setting = list(n = n, p = p, K = K, n_class = n_class, x_dist = x_dist,
                      cov_type = cov_type, rho = rho, df_t = df_t, mix_prob = mix_prob,
                      standardize_X = standardize_X, use_rfast = has_rfast))
}
