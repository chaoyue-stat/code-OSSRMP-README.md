# Setting 3: assess sensitivity to rho_N under four covariate distributions.
# Each replicate uses one dataset for all rho_N values and five methods.
# Save replicate results, failure records, and Monte Carlo summaries.
# MSE is the sum of squared coefficient errors, as in the original script.
# The legacy CPU column records elapsed fitting time in seconds.

rm(list = ls())
gc()

library(MASS)
library(mvtnorm)
library(nnet)
library(dplyr)
library(parallel)

# Run from the directory containing the companion files.
source("Syn.R")
source("ossrmp_pre.R")
source("ossrmp.R")
source("ossrmp_pos.R")

# Helper functions

.a1_methods <- function() c("OSSRMP", "MR-O", "POSSRMP", "POSSRMP-O", "OSSRMP-Pos")
.a1_seed <- function(seed, offset = 0) {
  if (is.null(seed)) return(NULL)
  value <- as.double(seed) + offset
  if (length(value) != 1L || !is.finite(value) || value < 0 ||
      value > .Machine$integer.max || value != floor(value)) stop("Invalid seed.")
  as.integer(value)
}
.a1_threads <- function() {
  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  invisible(NULL)
}
.a1_dist <- function(Dist) {
  switch(as.character(Dist),
    "1" = list(DistLabel = "D1-Gaussian-AR1", x_dist = "gaussian", cov_type = "ar1", rho = 0.5),
    "2" = list(DistLabel = "D2-t2-EXCH", x_dist = "t", cov_type = "exchangeable", rho = 0.5, df_t = 2),
    "3" = list(DistLabel = "D3-SkewQuad-AR1", x_dist = "skew_quadratic", cov_type = "ar1", rho = 0.5),
    "4" = list(DistLabel = "D4-GMix-EXCH", x_dist = "gaussian_mixture", cov_type = "exchangeable", rho = 0.5, mix_prob = 0.3),
    stop("Dist must be one of 1, 2, 3, 4."))
}
.a1_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}
.a1_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) > 1L) sd(x) else NA_real_
}
.a1_mse <- function(beta, theta_true) {
  theta <- as.vector(beta)
  if (!is.numeric(theta) || length(theta) != length(theta_true) || any(!is.finite(theta)))
    stop("Invalid coefficient estimates or coefficient dimensions.")
  sum((theta - theta_true)^2)
}
.a1_nsel <- function(fit) {
  if (!is.null(fit$n.sel)) return(as.integer(fit$n.sel))
  if (!is.null(fit$idx.sel)) return(as.integer(length(fit$idx.sel)))
  NA_integer_
}
.a1_nhr <- function(fit) {
  if (!is.null(fit$Omega.dag)) return(as.integer(length(fit$Omega.dag)))
  if (!is.null(fit$N.dag)) return(as.integer(fit$N.dag))
  NA_integer_
}
.a1_row <- function(i, method, rho, cfg, error = NA_character_) {
  data.frame(rep = i, Method = method, FitSource = method, Dist = cfg$Dist,
             DistLabel = .a1_dist(cfg$Dist)$DistLabel, K = cfg$K, N = cfg$N,
             n0 = cfg$n0, n_dag = cfg$n_dag, rho_N = rho,
             N_hr = NA_integer_, n_sel = NA_integer_, MSE = NA_real_, CPU = NA_real_,
             Success = FALSE, Converged = FALSE, ErrorMsg = error, stringsAsFactors = FALSE)
}
.a1_rep_summary <- function(raw, ok = TRUE) {
  raw %>%
    group_by(rep, Dist, DistLabel, K, N) %>%
    summarise(n_rows = n(), n_success = sum(Success), n_fail = sum(!Success),
              ok = ok, .groups = "drop")
}
.a1_failed_rep <- function(i, cfg, message) {
  rows <- lapply(cfg$rho_grid, function(rho) bind_rows(lapply(.a1_methods(), function(method)
    .a1_row(i, method, rho, cfg, error = message))))
  raw <- bind_rows(rows)
  list(raw = raw, rep_summary = .a1_rep_summary(raw, ok = FALSE))
}
.a1_try_fit <- function(fun) {
  start <- proc.time()["elapsed"]
  fit <- tryCatch(fun(), error = function(e) e)
  elapsed <- unname(proc.time()["elapsed"] - start)
  list(fit = fit, cpu = elapsed)
}
.a1_call <- function(method, args) {
  switch(method,
    OSSRMP = do.call(ossrmp, c(args, list(RB = TRUE))),
    `MR-O` = do.call(ossrmp, c(args, list(RB = FALSE))),
    POSSRMP = do.call(ossrmp_pre, c(args, list(calibrated = TRUE))),
    `POSSRMP-O` = do.call(ossrmp_pre, c(args, list(calibrated = FALSE))),
    `OSSRMP-Pos` = do.call(ossrmp_poisson, c(args, list(RB = TRUE))))
}
.a1_log_error <- function(e, i, cfg) {
  dir.create(cfg$error_dir, recursive = TRUE, showWarnings = FALSE)
  file <- file.path(cfg$error_dir, sprintf("err_D%s_K%s_rep%04d.txt", cfg$Dist, cfg$K, i))
  text <- c(sprintf("rep=%d, Dist=%s, K=%s", i, cfg$Dist, cfg$K),
            paste("message:", conditionMessage(e)),
            paste("call:", paste(deparse(conditionCall(e)), collapse = " ")))
  writeLines(text, file)
  file
}
.a1_run_rep <- function(i, cfg) {
  args <- cfg[c("N", "Dist", "K", "n0", "n_dag", "rho_grid", "tol", "maxit",
                "verbose_fit", "project_dir", "seed_base")]
  error_file <- NA_character_
  out <- tryCatch(do.call(A1RhoK3_inner_grid, c(list(i = i), args)), error = function(e) {
    error_file <<- .a1_log_error(e, i, cfg)
    .a1_failed_rep(i, cfg, conditionMessage(e))
  })
  file <- file.path(cfg$inner_out_dir, sprintf("A1RhoK3_rep%04d_D%s_K%s.rds", i, cfg$Dist, cfg$K))
  record <- list(ok = isTRUE(out$rep_summary$ok[1L]),
                 context = list(rep = i, Dist = cfg$Dist, K = cfg$K), timestamp = Sys.time(),
                 result = out, error_file = error_file)
  saveRDS(record, file)
  out
}
.a1_summarize <- function(raw) {
  # Retain all successful fits; report convergence separately.
  raw %>%
    group_by(Method, FitSource, Dist, DistLabel, K, N, n0, n_dag, rho_N) %>%
    summarise(B = sum(Success), Fail = sum(!Success),
              ConvRate = .a1_mean(Converged[Success]),
              MSE_mean = .a1_mean(MSE[Success]), MSE_sd = .a1_sd(MSE[Success]),
              MSE_se = if (B > 0L) MSE_sd / sqrt(B) else NA_real_,
              CPU_mean = .a1_mean(CPU[Success]), N_hr = .a1_mean(N_hr),
              n_sel = .a1_mean(n_sel), .groups = "drop") %>%
    rename(MSE = MSE_mean, CPU = CPU_mean) %>%
    arrange(rho_N, Method)
}
.a1_validate <- function(N, Dist, K, m, n0, n_dag, rho_grid, tol, maxit, cores) {
  is_count <- function(x) length(x) == 1L && is.finite(x) && x >= 1 &&
    x == floor(x) && x <= .Machine$integer.max
  stopifnot(is_count(N), is_count(K), is_count(m), is_count(n0), is_count(n_dag),
            is_count(maxit), is_count(cores), n0 < N, n_dag <= N - n0,
            length(Dist) > 0L, all(Dist %in% 1:4), !anyDuplicated(Dist),
            length(rho_grid) > 0L, all(is.finite(rho_grid)), !anyDuplicated(rho_grid),
            all(rho_grid > 0 & rho_grid < 1), length(tol) == 1L, is.finite(tol), tol > 0)
}

# One replicate: generate data once, then fit every method over rho_N.

A1RhoK3_inner_grid <- function(
    i, N = 1e6, Dist = 1L, K = 5L, n0 = 3000L, n_dag = 3000L,
    rho_grid = c(0.90, 0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999),
    tol = 1e-8, maxit = 200L, verbose_fit = FALSE,
    project_dir = getwd(), seed_base = NULL) {
  .a1_threads()
  seed <- .a1_seed(seed_base, i - 1L)
  if (!is.null(seed)) set.seed(seed)
  dist_args <- .a1_dist(Dist)
  data <- do.call(DataGenSM, c(list(n = N, n_class = K + 1L, standardize_X = FALSE, seed = NULL),
                              dist_args[setdiff(names(dist_args), "DistLabel")]))
  if (!is.matrix(data$B_true) || ncol(data$B_true) != K + 1L || any(!is.finite(data$B_true)))
    stop("DataGenSM must return a finite B_true matrix with K + 1 columns.")
  theta_true <- as.vector(data$B_true[, -1L, drop = FALSE])
  cfg <- list(N = N, Dist = Dist, K = K, n0 = n0, n_dag = n_dag)
  methods <- .a1_methods()
  seed_offsets <- c(11L, 21L, 31L, 41L, 51L)
  rows <- vector("list", length(rho_grid) * length(methods))
  row_id <- 0L
  for (j in seq_along(rho_grid)) {
    rho <- rho_grid[j]
    seed_now <- .a1_seed(seed_base, 100000 * i + 1000 * j)
    args <- list(X = data$X, y = data$y, n0 = n0, n_dag = n_dag, rho = rho,
                 tol = tol, maxit = maxit, verbose = verbose_fit)
    for (k in seq_along(methods)) {
      method <- methods[k]
      row <- .a1_row(i, method, rho, cfg)
      # Keep the original method-specific seed offsets.
      fit_args <- c(args, list(seed = .a1_seed(seed_now, seed_offsets[k])))
      result <- .a1_try_fit(function() .a1_call(method, fit_args))
      row$CPU <- result$cpu
      fit <- result$fit
      if (inherits(fit, "error")) row$ErrorMsg <- conditionMessage(fit) else {
        metrics <- tryCatch(list(MSE = .a1_mse(fit$beta, theta_true),
                                 N_hr = .a1_nhr(fit), n_sel = .a1_nsel(fit),
                                 Converged = isTRUE(fit$converged)), error = function(e) e)
        if (inherits(metrics, "error")) row$ErrorMsg <- conditionMessage(metrics) else {
          row[names(metrics)] <- metrics
          row$Success <- TRUE
        }
      }
      row_id <- row_id + 1L
      rows[[row_id]] <- row
      rm(fit, result)
    }
  }
  raw <- bind_rows(rows)
  list(raw = raw, rep_summary = .a1_rep_summary(raw))
}

# Parallel repetitions: use PSOCK on Windows and multicore on Unix by default.

A1RhoK3_mclapply_grid <- function(
    m, N = 1e6, Dist = 1L, K = 5L, n0 = 3000L, n_dag = 3000L,
    rho_grid = c(0.90, 0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999),
    tol = 1e-8, maxit = 200L, verbose_fit = FALSE, cores = 4L,
    backend = c("auto", "psock", "multicore"), project_dir = getwd(),
    inner_out_dir = file.path(project_dir, "A1RhoK3_inner_rds"), seed_base = 12345) {
  .a1_validate(N, Dist, K, m, n0, n_dag, rho_grid, tol, maxit, cores)
  stopifnot(length(Dist) == 1L)
  .a1_seed(seed_base)
  .a1_threads()
  backend <- match.arg(backend)
  if (backend == "auto") backend <- if (.Platform$OS.type == "windows") "psock" else "multicore"
  if (.Platform$OS.type == "windows") backend <- "psock"
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  dir.create(inner_out_dir, recursive = TRUE, showWarnings = FALSE)
  inner_out_dir <- normalizePath(inner_out_dir, mustWork = TRUE)
  cfg <- list(N = N, Dist = Dist, K = K, n0 = n0, n_dag = n_dag, rho_grid = rho_grid,
              tol = tol, maxit = maxit, verbose_fit = verbose_fit, project_dir = project_dir,
              inner_out_dir = inner_out_dir, error_dir = file.path(inner_out_dir, "errors"), seed_base = seed_base)
  workers <- min(as.integer(cores), as.integer(m))
  if (workers == 1L) {
    results <- lapply(seq_len(m), .a1_run_rep, cfg = cfg)
  } else if (backend == "psock") {
    cl <- makePSOCKcluster(workers)
    on.exit(stopCluster(cl), add = TRUE)
    clusterCall(cl, function(project_dir) {
      Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
      setwd(project_dir)
      suppressPackageStartupMessages({
        library(MASS)
        library(mvtnorm)
        library(nnet)
        library(dplyr)
      })
      for (file in c("Syn.R", "ossrmp_pre.R", "ossrmp.R", "ossrmp_pos.R")) sys.source(file, .GlobalEnv)
      TRUE
    }, project_dir)
    export_env <- environment(A1RhoK3_mclapply_grid)
    exports <- c("A1RhoK3_inner_grid", ls(envir = export_env, all.names = TRUE, pattern = "^\\.a1_"))
    clusterExport(cl, exports, envir = export_env)
    results <- parLapply(cl, seq_len(m), .a1_run_rep, cfg = cfg)
  } else {
    results <- mclapply(seq_len(m), .a1_run_rep, cfg = cfg, mc.cores = workers, mc.preschedule = TRUE)
  }
  raw <- bind_rows(lapply(results, `[[`, "raw"))
  rep_summary <- bind_rows(lapply(results, `[[`, "rep_summary"))
  tag <- sprintf("D%s_K%s_m%d", Dist, K, m)
  saveRDS(raw, file.path(inner_out_dir, paste0("A1RhoK3_raw_", tag, ".rds")))
  saveRDS(rep_summary, file.path(inner_out_dir, paste0("A1RhoK3_rep_summary_", tag, ".rds")))
  cat(sprintf("[A1RhoK3] Dist=%s | K=%s | %d/%d successful method fits\n", Dist, K, sum(raw$Success), nrow(raw)))
  .a1_summarize(raw)
}

# Run the four distribution settings and save cumulative results.

RunA1RhoK3Study <- function(
    N = 1e6, Dist = 1:4, K = 5L, m = 200L, n0 = 3000L, n_dag = 3000L,
    rho_grid = c(0.90, 0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999),
    tol = 1e-8, maxit = 200L, verbose_fit = FALSE, cores = 4L, backend = "auto",
    project_dir = getwd(), out_dir = file.path(project_dir, "A1RhoK3_runs"),
    log_to_console = TRUE, log_file = file.path(out_dir, "run.log"), seed_base = 12345) {
  .a1_validate(N, Dist, K, m, n0, n_dag, rho_grid, tol, maxit, cores)
  .a1_seed(seed_base)
  backend <- match.arg(backend, c("auto", "psock", "multicore"))
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, mustWork = TRUE)
  if (!is.null(log_file) && !is.na(log_file)) dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  log <- function(...) {
    message <- paste0("[", format(Sys.time(), "%F %T"), "] ", paste0(..., collapse = ""))
    if (log_to_console) cat(message, "\n")
    if (!is.null(log_file) && !is.na(log_file)) cat(message, "\n", file = log_file, append = TRUE)
  }
  files <- file.path(project_dir, c("Syn.R", "ossrmp_pre.R", "ossrmp.R", "ossrmp_pos.R"))
  metadata <- list(N = N, Dist = Dist, K = K, m = m, n0 = n0, n_dag = n_dag,
                   rho_grid = rho_grid, tol = tol, maxit = maxit, cores = cores, backend = backend,
                   seed_base = seed_base, source_md5 = tools::md5sum(files), timestamp = Sys.time(),
                   session = capture.output(sessionInfo()), timing = "CPU stores elapsed seconds")
  saveRDS(metadata, file.path(out_dir, "metadata.rds"))
  start <- Sys.time()
  results <- vector("list", length(Dist))
  log(sprintf("RunA1RhoK3Study: %d distributions | K=%s | N=%s | m=%s", length(Dist), K, N, m))
  for (j in seq_along(Dist)) {
    tag <- sprintf("task-%03d_D%s_K%s", j, Dist[j], K)
    log("START ", tag)
    results[[j]] <- tryCatch(A1RhoK3_mclapply_grid(
      m = m, N = N, Dist = Dist[j], K = K, n0 = n0, n_dag = n_dag,
      rho_grid = rho_grid, tol = tol, maxit = maxit, verbose_fit = verbose_fit,
      cores = cores, backend = backend, project_dir = project_dir,
      inner_out_dir = file.path(out_dir, "replicates"),
      seed_base = .a1_seed(seed_base, 1000 * j)), error = function(e) {
        log("ERROR ", tag, ": ", conditionMessage(e))
        structure(list(error = conditionMessage(e)), class = "sim_error")
      })
    saveRDS(results[[j]], file.path(out_dir, paste0("A1RhoK3_", tag, ".rds")))
    ok <- vapply(results, function(x) !is.null(x) && !inherits(x, "sim_error"), logical(1))
    if (any(ok)) saveRDS(bind_rows(results[ok]), file.path(out_dir, "A1RhoK3_cumulative.rds"))
    log("DONE ", tag)
  }
  ok <- vapply(results, function(x) !is.null(x) && !inherits(x, "sim_error"), logical(1))
  if (!any(ok)) stop("All distribution tasks failed; see run.log and task records.")
  final <- bind_rows(results[ok])
  attr(final, "elapsed_hours") <- as.numeric(difftime(Sys.time(), start, units = "hours"))
  saveRDS(final, file.path(out_dir, "A1RhoK3_final.rds"))
  write.csv(final, file.path(out_dir, "A1RhoK3_final.csv"), row.names = FALSE)
  log(sprintf("ALL DONE: %.2f hours.", attr(final, "elapsed_hours")))
  final
}

# Optional MSE plot; each distribution has its own panel.

PlotA1RhoK3MSE <- function(res_a1, log_y = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install package 'ggplot2' to plot results.")
  methods <- .a1_methods()
  labels <- expression(hat(vartheta)[MR], hat(vartheta)[MR * ",O"],
                       hat(vartheta)[MR * ",P"], hat(vartheta)[MR * ",P,O"], hat(vartheta)[MR]^Pos)
  rho_values <- sort(unique(res_a1$rho_N))
  rho_labels <- format(rho_values, digits = 10, trim = TRUE, scientific = FALSE)
  data <- res_a1
  data$Method <- factor(data$Method, levels = methods)
  data$rho_f <- factor(data$rho_N, levels = rho_values, labels = rho_labels)
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = rho_f, y = MSE, color = Method,
                                          shape = Method, linetype = Method, group = Method)) +
    ggplot2::geom_line(linewidth = 0.9) + ggplot2::geom_point(size = 2.4) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = pmax(MSE - 1.96 * MSE_se, 0),
                                     ymax = MSE + 1.96 * MSE_se), width = 0.08, linewidth = 0.5, na.rm = TRUE) +
    ggplot2::facet_wrap(~DistLabel, ncol = 2) +
    ggplot2::scale_color_discrete(breaks = methods, labels = labels) +
    ggplot2::scale_shape_discrete(breaks = methods, labels = labels) +
    ggplot2::scale_linetype_discrete(breaks = methods, labels = labels) +
    ggplot2::labs(x = expression(rho[N]), y = if (log_y) "MSE (log scale)" else "MSE") +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom",
                   legend.box = "horizontal", legend.title = ggplot2::element_blank()) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1), shape = ggplot2::guide_legend(nrow = 1),
                    linetype = ggplot2::guide_legend(nrow = 1))
  if (log_y) plot <- plot + ggplot2::scale_y_log10()
  plot
}

# Simulation settings and calls
# Set seed_base to an integer for reproducible repetitions.
# Set options(A1RhoK3.autorun = FALSE) before sourcing to skip the formal run.

N <- 1e6
Dist <- 1:4
K <- 3L
m <- 200L
n0 <- 3000L
n_dag <- 3000L
rho_grid <- c(0.90, 0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999)
cores <- 20L
seed_base <- NULL
project_dir <- getwd()
out_dir <- file.path(project_dir, "A1RhoK3_runs")

if (isTRUE(getOption("A1RhoK3.autorun", TRUE))) {
  res_a1 <- RunA1RhoK3Study(
    N = N, Dist = Dist, K = K, m = m, n0 = n0, n_dag = n_dag,
    rho_grid = rho_grid, tol = 1e-8, maxit = 200L, verbose_fit = FALSE,
    cores = cores, backend = "auto", project_dir = project_dir,
    out_dir = out_dir, seed_base = seed_base)
  print(attr(res_a1, "elapsed_hours"))
}

# Optional plot (requires ggplot2).
# p_a1 <- PlotA1RhoK3MSE(res_a1, log_y = TRUE)
# print(p_a1)
# ggplot2::ggsave(file.path(out_dir, "A1RhoK3_MSE_plot.pdf"), p_a1, width = 10, height = 7)
