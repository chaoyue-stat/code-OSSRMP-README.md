rm(list = ls())
gc()

# Setting 1 compares eight estimators over four covariate distributions.
# Within each replication, all methods and subsample sizes use the same data
# and train-test split. OSSRMP-Sum is obtained by transforming OSSRMP.
# The original MSE and elapsed fitting-time definitions are retained.

suppressPackageStartupMessages({
  library(MASS)
  library(mvtnorm)
  library(nnet)
  library(dplyr)
  library(parallel)
})

# Run this script from the directory containing the companion files.
source("Syn.R")
source("ossrmp_pre.R")
source("ossrmp.R")
source("ossrmp_pos.R")
source("OSS.R")

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

make_seed_s1 <- function(seed, offset = 0) {
  if (is.null(seed)) return(NULL)
  value <- as.double(seed) + offset
  if (length(value) != 1L || !is.finite(value) || value < 0 ||
      value > .Machine$integer.max || value != floor(value)) stop("Invalid seed.")
  as.integer(value)
}

distribution_args_s1 <- function(Dist) {
  if (!Dist %in% 1:4) stop("Dist must be one of 1, 2, 3, 4.")
  list(
    list(DistLabel = "D1-Gaussian-AR1", x_dist = "gaussian", cov_type = "ar1", rho = 0.5),
    list(DistLabel = "D2-t4-EXCH", x_dist = "t", cov_type = "exchangeable", rho = 0.5, df_t = 4),
    list(DistLabel = "D3-SkewQuad-AR1", x_dist = "skew_quadratic", cov_type = "ar1", rho = 0.5),
    list(DistLabel = "D4-GMix-EXCH", x_dist = "gaussian_mixture", cov_type = "exchangeable", rho = 0.5, mix_prob = 0.3)
  )[[as.integer(Dist)]]
}

# Convert coefficients and probabilities between baseline and sum constraints.
to_sum_constraint_s1 <- function(beta) {
  B <- cbind(0, beta)
  sweep(B, 1, rowMeans(B), "-")[, -1, drop = FALSE]
}

sum_constraint_map_s1 <- function(d, K) {
  kronecker(diag(K) - matrix(1 / (K + 1), K, K), diag(d))
}

softmax_prob_s1 <- function(X, beta, constraint = "baseline") {
  b0 <- if (constraint == "sum") -rowSums(beta) else rep(0, nrow(beta))
  eta <- cbind(1, X) %*% cbind(b0, beta)
  E <- exp(sweep(eta, 1, apply(eta, 1, max), "-"))
  E / rowSums(E)
}

# CP and interval length refer to the first coefficient in column-major order.
evaluate_fit_s1 <- function(beta, covariance, truth, X_test, P_test,
                            constraint = "baseline") {
  theta <- as.vector(beta)
  if (is.null(dim(beta)) || length(theta) != length(truth) ||
      nrow(beta) != ncol(X_test) + 1L || ncol(beta) + 1L != ncol(P_test) ||
      any(!is.finite(theta))) stop("Invalid fitted coefficient matrix.")
  if (!is.matrix(covariance) || !identical(dim(covariance), c(length(theta), length(theta))) ||
      any(!is.finite(covariance))) stop("Invalid fitted covariance matrix.")
  se <- sqrt(pmax(diag(covariance), 0))
  err <- theta[1] - truth[1]
  P_hat <- softmax_prob_s1(X_test, beta, constraint)
  c(MSE = sum((theta - truth)^2), CP = as.numeric(abs(err) <= 1.96 * se[1]),
    LEN = 2 * 1.96 * se[1], ERR = err, SE = se[1],
    MSPE = mean(rowSums((P_hat - P_test)^2)))
}

fit_once_s1 <- function(fun) {
  start <- proc.time()[["elapsed"]]
  fit <- tryCatch(fun(), error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - start
  failed <- inherits(fit, "error")
  list(ok = !failed, fit = if (failed) NULL else fit, cpu = elapsed,
       err = if (failed) conditionMessage(fit) else NA_character_)
}

# Pack one method result; evaluation failures affect only that method.
result_row_s1 <- function(context, method, result, truth = NULL, X_test = NULL,
                          P_test = NULL, constraint = "baseline", fit_source = method) {
  fit <- result$fit
  met <- setNames(rep(NA_real_, 6), c("MSE", "CP", "LEN", "ERR", "SE", "MSPE"))
  if (result$ok) {
    value <- tryCatch(evaluate_fit_s1(fit$beta, fit$cov.beta, truth, X_test, P_test, constraint),
                      error = function(e) e)
    if (inherits(value, "error")) {
      result$ok <- FALSE
      result$err <- conditionMessage(value)
    } else met <- value
  }
  n_sel <- if (!is.null(fit$n.sel)) as.integer(fit$n.sel) else
    if (!is.null(fit$idx.sel)) length(fit$idx.sel) else NA_integer_
  n_hr <- if (!is.null(fit$Omega.dag)) length(fit$Omega.dag) else NA_integer_
  data.frame(context, Method = method, FitSource = fit_source, N_hr = n_hr,
             n_sel = n_sel, as.list(met), CPU = result$cpu, Success = result$ok,
             Converged = result$ok && isTRUE(fit$converged), ErrorMsg = result$err,
             stringsAsFactors = FALSE, row.names = NULL)
}

replication_summary_s1 <- function(raw, ok = TRUE) {
  data.frame(rep = raw$rep[1], Dist = raw$Dist[1], DistLabel = raw$DistLabel[1],
             K = raw$K[1], N = raw$N[1], n_rows = nrow(raw),
             n_success = sum(raw$Success), n_fail = sum(!raw$Success), ok = ok)
}

# ------------------------------------------------------------
# One replication
# ------------------------------------------------------------
# A common dataset and split are reused throughout each replication.
run_replication <- function(i, N, Dist, K, n0, n_dag_grid,
                            train_ratio = 0.7, rho_mr = 0.999, tol = 1e-8,
                            maxit = 200L, verbose_fit = FALSE,
                            seed_base = 12345L) {
  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  if (!is.null(seed_base)) set.seed(make_seed_s1(seed_base, i - 1))
  args <- distribution_args_s1(Dist)
  label <- args$DistLabel
  args$DistLabel <- NULL
  dat <- do.call(DataGenSM, c(list(n = N, n_class = K + 1,
                                 standardize_X = FALSE, seed = NULL), args))
  X <- as.matrix(dat$X)
  y <- dat$y
  P_true <- as.matrix(dat$P_true)
  B_true <- as.matrix(dat$B_true)
  if (nrow(X) != N || length(y) != N || nrow(P_true) != N ||
      ncol(P_true) != K + 1L || ncol(B_true) != K + 1L ||
      nrow(B_true) != ncol(X) + 1L || any(!is.finite(B_true)) ||
      any(B_true[, 1] != 0)) stop("Invalid baseline-constraint DataGenSM output.")
  base_truth <- as.vector(B_true[, -1, drop = FALSE])
  sum_truth <- as.vector(sweep(B_true, 1, rowMeans(B_true), "-")[, -1, drop = FALSE])
  Tmat <- sum_constraint_map_s1(nrow(B_true), K)
  idx_train <- sample.int(N, floor(train_ratio * N))
  is_train <- logical(N)
  is_train[idx_train] <- TRUE
  idx_test <- which(!is_train)
  X_train <- X[idx_train, , drop = FALSE]
  y_train <- y[idx_train]
  X_test <- X[idx_test, , drop = FALSE]
  P_test <- P_true[idx_test, , drop = FALSE]
  rows <- vector("list", 8L * length(n_dag_grid))
  at <- 0L

  for (jj in seq_along(n_dag_grid)) {
    n_dag <- n_dag_grid[jj]
    seed_now <- make_seed_s1(seed_base, 100000 * i + 1000 * jj)
    context <- list(rep = i, Dist = Dist, DistLabel = label, K = K, N = N,
                    n_train = length(idx_train), n_test = length(idx_test),
                    n0 = n0, n_dag = n_dag)
    fits <- list(
      OSSRMP = fit_once_s1(function() ossrmp(
        X = X_train, y = y_train, n0 = n0, n_dag = n_dag,
        rho = rho_mr, RB = TRUE, seed = make_seed_s1(seed_now, 11),
        tol = tol, maxit = maxit, verbose = verbose_fit
      )),
      `MR-O` = fit_once_s1(function() ossrmp(
        X = X_train, y = y_train, n0 = n0, n_dag = n_dag,
        rho = rho_mr, RB = FALSE, seed = make_seed_s1(seed_now, 21),
        tol = tol, maxit = maxit, verbose = verbose_fit
      )),
      OSS = fit_once_s1(function() OSS(
        X = X_train, y = y_train, n0 = n0, n_sub = n_dag,
        poisson = FALSE, seed = make_seed_s1(seed_now, 31),
        tol = tol, maxit = maxit, verbose = verbose_fit
      )),
      `OSS-Pos` = fit_once_s1(function() OSS(
        X = X_train, y = y_train, n0 = n0, n_sub = n_dag,
        poisson = TRUE, seed = make_seed_s1(seed_now, 41),
        tol = tol, maxit = maxit, verbose = verbose_fit
      )),
      POSSRMP = fit_once_s1(function() ossrmp_pre(
        X = X_train, y = y_train, n0 = n0, n_dag = n_dag,
        rho = rho_mr, calibrated = TRUE, seed = make_seed_s1(seed_now, 51),
        tol = tol, maxit = maxit, verbose = verbose_fit
      )),
      `POSSRMP-O` = fit_once_s1(function() ossrmp_pre(
        X = X_train, y = y_train, n0 = n0, n_dag = n_dag,
        rho = rho_mr, calibrated = FALSE, seed = make_seed_s1(seed_now, 61),
        tol = tol, maxit = maxit, verbose = verbose_fit
      )),
      `OSSRMP-Pos` = fit_once_s1(function() ossrmp_poisson(
        X = X_train, y = y_train, n0 = n0, n_dag = n_dag,
        rho = rho_mr, RB = TRUE, seed = make_seed_s1(seed_now, 71),
        tol = tol, maxit = maxit, verbose = verbose_fit
      ))
    )

    for (method in names(fits)) {
      result <- fits[[method]]
      at <- at + 1L
      rows[[at]] <- result_row_s1(context, method, result, base_truth, X_test, P_test)
      if (method == "OSSRMP") {
        derived <- result
        if (result$ok) {
          transformed <- tryCatch(list(
            beta = to_sum_constraint_s1(result$fit$beta),
            covariance = Tmat %*% result$fit$cov.beta %*% t(Tmat)
          ), error = function(e) e)
          if (inherits(transformed, "error")) {
            derived$ok <- FALSE
            derived$err <- conditionMessage(transformed)
          } else {
            derived$fit$beta <- transformed$beta
            derived$fit$cov.beta <- transformed$covariance
          }
        }
        at <- at + 1L
        rows[[at]] <- result_row_s1(context, "OSSRMP-Sum", derived, sum_truth,
                             X_test, P_test, "sum", "OSSRMP")
      }
    }
  }
  raw <- dplyr::bind_rows(rows)
  list(raw = raw, rep_summary = replication_summary_s1(raw))
}

# Retain failed replications in the raw table so failure counts remain complete.
failed_replication_s1 <- function(i, cfg, message) {
  methods <- c("OSSRMP", "OSSRMP-Sum", "MR-O", "OSS", "OSS-Pos",
               "POSSRMP", "POSSRMP-O", "OSSRMP-Pos")
  grid <- expand.grid(Method = methods, n_dag = cfg$n_dag_grid, stringsAsFactors = FALSE)
  rows <- lapply(seq_len(nrow(grid)), function(j) {
    context <- list(rep = i, Dist = cfg$Dist, DistLabel = distribution_args_s1(cfg$Dist)$DistLabel,
                    K = cfg$K, N = cfg$N, n_train = floor(cfg$train_ratio * cfg$N),
                    n_test = cfg$N - floor(cfg$train_ratio * cfg$N),
                    n0 = cfg$n0, n_dag = grid$n_dag[j])
    source <- if (grid$Method[j] == "OSSRMP-Sum") "OSSRMP" else grid$Method[j]
    result_row_s1(context, grid$Method[j], list(ok = FALSE, fit = NULL, cpu = NA_real_, err = message),
            fit_source = source)
  })
  raw <- dplyr::bind_rows(rows)
  list(raw = raw, rep_summary = replication_summary_s1(raw, FALSE))
}

safe_mean_s1 <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
safe_sd_s1 <- function(x) if (sum(is.finite(x)) > 1L) sd(x[is.finite(x)]) else NA_real_
                       
                       # As in the original script, summaries use successful fits, including returned
                       # fits flagged as nonconverged. Convergence rates are reported separately.
                       summarize_results_s1 <- function(raw) {
  keys <- c("Method", "FitSource", "Dist", "DistLabel", "K", "N", "n_train", "n_test", "n0", "n_dag")
  metrics <- c("MSE", "CP", "LEN", "ERR", "SE", "MSPE", "CPU")
  prepared <- raw
  prepared[!prepared$Success, metrics] <- NA_real_
  out <- dplyr::summarise(
    dplyr::group_by(prepared, dplyr::across(dplyr::all_of(keys))),
    B = sum(Success), Fail = sum(!Success),
    ConvRate = safe_mean_s1(as.numeric(Converged[Success])),
    Bias = safe_mean_s1(ERR), AbsBias = safe_mean_s1(abs(ERR)), EmpSD = safe_sd_s1(ERR),
    dplyr::across(dplyr::all_of(metrics), list(mean = safe_mean_s1, sd = safe_sd_s1), .names = "{.col}_{.fn}"),
    N_hr = safe_mean_s1(N_hr), n_sel = safe_mean_s1(n_sel), .groups = "drop"
  )
  for (metric in metrics) out[[paste0(metric, "_se")]] <- out[[paste0(metric, "_sd")]] / sqrt(pmax(out$B, 1))
  out$SERatio <- ifelse(out$EmpSD > 0, out$SE_mean / out$EmpSD, NA_real_)
  out <- dplyr::rename(out, MSE = MSE_mean, CP = CP_mean, LEN = LEN_mean,
                       SE = SE_mean, MSPE = MSPE_mean, CPU = CPU_mean)
  out$ERR_mean <- NULL
  dplyr::arrange(out, Dist, K, n_dag, Method)
}

# ------------------------------------------------------------
# Independent replications
# ------------------------------------------------------------
# Use PSOCK on Windows; either backend is available on Unix.
run_repetitions <- function(m, N, Dist, K, n0, n_dag_grid,
                            train_ratio = 0.7, rho_mr = 0.999, tol = 1e-8,
                            maxit = 200L, verbose_fit = FALSE, cores = 4L,
                            backend = c("auto", "psock", "multicore"),
                            project_dir = getwd(),
                            inner_out_dir = file.path(project_dir, "S1SMK3_inner_rds"),
                            seed_base = 12345L) {
  backend <- match.arg(backend)
  if (backend == "auto") backend <- if (.Platform$OS.type == "windows") "psock" else "multicore"
  if (backend == "multicore" && .Platform$OS.type == "windows") stop("Use backend = 'psock' on Windows.")
  cores <- min(as.integer(cores), as.integer(m))
  if (is.na(cores) || cores < 1L || m < 1L) stop("m and cores must be positive integers.")
  dir.create(inner_out_dir, recursive = TRUE, showWarnings = FALSE)
  cfg <- list(N = N, Dist = Dist, K = K, n0 = n0, n_dag_grid = n_dag_grid,
              train_ratio = train_ratio, rho_mr = rho_mr, tol = tol, maxit = maxit,
              verbose_fit = verbose_fit, seed_base = seed_base)
  one_rep <- function(i) {
    cat(sprintf("[Setting 1] rep %d/%d | Dist=%d | K=%d\n", i, m, Dist, K))
    out <- tryCatch(do.call(run_replication, c(list(i = i), cfg)), error = function(e) e)
    ok <- !inherits(out, "error")
    error <- if (ok) NULL else conditionMessage(out)
    if (!ok) out <- failed_replication_s1(i, cfg, error)
    tag <- sprintf("rep%03d_D%d_K%d", i, Dist, K)
    stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
    record <- list(ok = ok, context = list(rep = i, Dist = Dist, K = K),
                   timestamp = Sys.time(), result = out, error = error)
    saveRDS(record, file.path(inner_out_dir, sprintf("S1SMK3_%s_%s_%s.rds", tag, if (ok) "ok" else "err", stamp)))
    failures <- out$raw[!out$raw$Success, c("Method", "n_dag", "ErrorMsg"), drop = FALSE]
    if (nrow(failures)) write.table(failures, file.path(inner_out_dir, paste0("errors_", tag, ".tsv")),
                                    sep = "\t", quote = TRUE, row.names = FALSE)
    out
  }
  if (cores == 1L) {
    results <- lapply(seq_len(m), one_rep)
  } else if (backend == "psock") {
    cl <- parallel::makeCluster(cores, type = "PSOCK", outfile = "")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    exports <- c("make_seed_s1", "distribution_args_s1", "to_sum_constraint_s1", "sum_constraint_map_s1", "softmax_prob_s1",
                 "evaluate_fit_s1", "fit_once_s1", "result_row_s1", "replication_summary_s1", "failed_replication_s1", "run_replication")
    parallel::clusterExport(cl, exports, envir = environment(run_replication))
    parallel::clusterCall(cl, function(path) {
      Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
      setwd(path)
      suppressPackageStartupMessages({
        library(MASS)
        library(mvtnorm)
        library(nnet)
        library(dplyr)
      })
      source("Syn.R")
      source("ossrmp_pre.R")
      source("ossrmp.R")
      source("ossrmp_pos.R")
      source("OSS.R")
      NULL
    }, project_dir)
    results <- parallel::parLapply(cl, seq_len(m), one_rep)
  } else {
    results <- parallel::mclapply(seq_len(m), one_rep, mc.cores = cores,
                                  mc.preschedule = TRUE, mc.set.seed = is.null(seed_base))
  }
  results <- lapply(seq_len(m), function(i) {
    out <- results[[i]]
    if (inherits(out, "try-error") || is.null(out)) {
      failed_replication_s1(i, cfg, if (is.null(out)) "Worker returned NULL." else as.character(out))
    } else out
  })
  raw <- dplyr::bind_rows(lapply(results, `[[`, "raw"))
  reps <- dplyr::bind_rows(lapply(results, `[[`, "rep_summary"))
  stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  saveRDS(raw, file.path(inner_out_dir, sprintf("S1SMK3_raw_D%d_K%d_m%d_%s.rds", Dist, K, m, stamp)))
  saveRDS(reps, file.path(inner_out_dir, sprintf("S1SMK3_rep_summary_D%d_K%d_m%d_%s.rds", Dist, K, m, stamp)))
  summarize_results_s1(raw)
}

# ------------------------------------------------------------
# Setting 1
# ------------------------------------------------------------
# Run each (Dist, K) combination and save the cumulative summary.
run_setting1 <- function(N = 1e6, Dist = 1:4, K = c(3, 5, 7), m = 200L,
                         n0 = 1000L, n_dag_grid = seq(2000L, 6000L, 1000L),
                         train_ratio = 0.7, rho_mr = 0.999, tol = 1e-8,
                         maxit = 200L, verbose_fit = FALSE, cores = 4L,
                         backend = "auto", project_dir = getwd(),
                         out_dir = file.path(project_dir, "S1SMK3_runs"),
                         log_to_console = TRUE, log_file = file.path(out_dir, "run.log"),
                         seed_base = 12345L) {
  stopifnot(N == floor(N), N > 2, m == floor(m), m >= 1,
            all(Dist %in% 1:4), all(K >= 1 & K == floor(K)),
            train_ratio > 0, train_ratio < 1, rho_mr > 0, rho_mr < 1,
            n0 >= 1, n0 == floor(n0), n0 <= floor(train_ratio * N),
            length(n_dag_grid) > 0, all(n_dag_grid >= 1 & n_dag_grid == floor(n_dag_grid)),
            cores >= 1, cores == floor(cores))
  backend <- match.arg(backend, c("auto", "psock", "multicore"))
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(log_file) && !is.na(log_file)) dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  log <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%F %T"), "] ", paste0(..., collapse = ""))
    if (log_to_console) cat(msg, "\n")
    if (!is.null(log_file) && !is.na(log_file)) cat(msg, "\n", file = log_file, append = TRUE)
  }
  params <- expand.grid(Dist = Dist, K = K, KEEP.OUT.ATTRS = FALSE)
  if (nrow(params) == 0L) stop("The parameter grid is empty.")
  metadata <- list(N = N, Dist = Dist, K = K, m = m, n0 = n0, n_dag_grid = n_dag_grid,
                   train_ratio = train_ratio, rho_mr = rho_mr, tol = tol, maxit = maxit,
                   cores = cores, backend = backend, seed_base = seed_base,
                   MSE_definition = "sum((theta_hat - theta_true)^2)",
                   CPU_definition = "elapsed fitting time in seconds", session = sessionInfo())
  saveRDS(metadata, file.path(out_dir, "S1SMK3_metadata.rds"))
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  results <- vector("list", nrow(params))
  start <- Sys.time()
  log("Starting ", nrow(params), " tasks; output: ", out_dir)
  for (idx in seq_len(nrow(params))) {
    tag <- sprintf("task-%03d_D%d_K%d", idx, params$Dist[idx], params$K[idx])
    log("START ", tag)
    task_start <- Sys.time()
    result <- tryCatch(run_repetitions(
      m = m, N = N, Dist = params$Dist[idx], K = params$K[idx], n0 = n0,
      n_dag_grid = n_dag_grid, train_ratio = train_ratio, rho_mr = rho_mr,
      tol = tol, maxit = maxit, verbose_fit = verbose_fit, cores = cores,
      backend = backend, project_dir = project_dir,
      inner_out_dir = file.path(out_dir, "inner_rds"), seed_base = make_seed_s1(seed_base, 1000 * idx)
    ), error = function(e) {
      log("ERROR ", tag, ": ", conditionMessage(e))
      structure(list(error = conditionMessage(e)), class = "sim_error")
    })
    saveRDS(result, file.path(out_dir, paste0("S1SMK3_", tag, ".rds")))
    results[[idx]] <- result
    valid <- vapply(results, is.data.frame, logical(1))
    if (any(valid)) saveRDS(dplyr::bind_rows(results[valid]), file.path(out_dir, "S1SMK3_cumulative.rds"))
    log("DONE ", tag, " (", round(as.numeric(difftime(Sys.time(), task_start, units = "mins")), 2), " min)")
  }
  valid <- vapply(results, is.data.frame, logical(1))
  if (!any(valid)) stop("All tasks failed; see run.log and the task files.")
  final <- dplyr::bind_rows(results[valid])
  attr(final, "elapsed_hours") <- as.numeric(difftime(Sys.time(), start, units = "hours"))
  attr(final, "failed_tasks") <- params[!valid, , drop = FALSE]
  saveRDS(final, file.path(out_dir, "S1SMK3_final.rds"))
  write.csv(final, file.path(out_dir, "S1SMK3_final.csv"), row.names = FALSE)
  log("ALL DONE: ", sum(valid), "/", nrow(params), " tasks completed.")
  final
}

# ------------------------------------------------------------
# Parameters and simulation call
# ------------------------------------------------------------

N <- 1e6
Dist <- 1:4
K <- 3                       # Use c(3, 5, 7) for the full study.
m <- 200
n0 <- 1000
n_dag_grid <- c(2000, 3000, 4000, 5000, 6000)
train_ratio <- 0.7
rho_mr <- 0.999
tol <- 1e-8
maxit <- 200
cores <- 20
backend <- "auto"            # PSOCK on Windows, multicore on Unix.
seed_base <- 12345
project_dir <- getwd()
out_dir <- file.path(project_dir, "S1SMK3_runs")

# For a short trial, use N = 2e4, Dist = 1, m = 2, n0 = 100,
# n_dag_grid = c(100, 200), cores = 2, and a separate output directory.

start_time <- Sys.time()
res_s1 <- run_setting1(
  N = N, Dist = Dist, K = K, m = m,
  n0 = n0, n_dag_grid = n_dag_grid,
  train_ratio = train_ratio, rho_mr = rho_mr,
  tol = tol, maxit = maxit, cores = cores, backend = backend,
  project_dir = project_dir, out_dir = out_dir, seed_base = seed_base
)
print(difftime(Sys.time(), start_time, units = "hours"))
