# Study A2: examine how pilot size affects the high-resolution set size.
# Generate one Gaussian-AR1 dataset per replicate and vary n0 at fixed rho_N.
# Pilot observations are excluded from the remaining sample before partitioning.
# Report N_dagger and N_dagger / N for K = 3, 5, 7, where K excludes the baseline class.
# The legacy CPU column records elapsed time for pilot fitting and partitioning.

rm(list = ls())
gc()

library(MASS)
library(mvtnorm)
library(nnet)
library(dplyr)
library(parallel)
library(ggplot2)

# Run from the directory containing the companion generator.
source("Syn.R")

# Helper functions

.a2_seed <- function(seed, offset = 0) {
  if (is.null(seed)) return(NULL)
  value <- as.double(seed) + offset
  if (length(value) != 1L || !is.finite(value) || value < 0 ||
      value > .Machine$integer.max || value != floor(value)) stop("Invalid seed.")
  as.integer(value)
}
.a2_threads <- function() {
  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  invisible(NULL)
}
.a2_dist <- function(Dist) {
  if (length(Dist) != 1L || is.na(Dist) || Dist != 1L)
    stop("This study uses only Example 1 (Gaussian-AR1).")
  list(DistLabel = "D1-Gaussian-AR1", x_dist = "gaussian", cov_type = "ar1", rho = 0.5)
}
.a2_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}
.a2_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) > 1L) sd(x) else NA_real_
}
.a2_beta <- function(fit0, cls, p) {
  coef0 <- coef(fit0)
  if (is.vector(coef0))
    coef0 <- matrix(coef0, nrow = 1L, dimnames = list(cls[2L], names(coef0)))
  coef0 <- as.matrix(coef0)
  if (is.null(rownames(coef0))) rownames(coef0) <- cls[-1L][seq_len(nrow(coef0))]
  predictors <- c("(Intercept)", paste0("x", seq_len(p)))
  if (!all(predictors %in% colnames(coef0)) || any(!is.finite(coef0)))
    stop("Pilot coefficients are missing or non-finite.")
  beta0 <- matrix(0, p + 1L, length(cls) - 1L, dimnames = list(predictors, cls[-1L]))
  fitted_classes <- intersect(rownames(coef0), cls[-1L])
  beta0[, fitted_classes] <- t(coef0[fitted_classes, predictors, drop = FALSE])
  beta0
}
.a2_partition <- function(y, cls, P0, rho) {
  setNames(lapply(seq_along(cls), function(k) which(y == cls[k] & P0[, k] >= rho)), cls)
}
.a2_row <- function(i, n0, cfg, error = NA_character_) {
  data.frame(rep = i, Dist = cfg$Dist, DistLabel = .a2_dist(cfg$Dist)$DistLabel,
             K = cfg$K, N = cfg$N, n0 = n0, rho_N = cfg$rho_mr, N_rem = cfg$N - n0,
             N_dag = NA_real_, N_dag_ratio = NA_real_, CPU = NA_real_, Success = FALSE,
             ErrorMsg = error, PilotSampling = cfg$pilot_sampling, stringsAsFactors = FALSE)
}
.a2_rep_summary <- function(raw, ok = TRUE) {
  raw %>%
    group_by(rep, Dist, DistLabel, K, N) %>%
    summarise(n_rows = n(), n_success = sum(Success), n_fail = sum(!Success),
              ok = ok, .groups = "drop")
}
.a2_failed_rep <- function(i, cfg, message) {
  raw <- bind_rows(lapply(cfg$n0_grid, function(n0) .a2_row(i, n0, cfg, error = message)))
  list(raw = raw, rep_summary = .a2_rep_summary(raw, ok = FALSE))
}
.a2_try_fit <- function(fun) {
  start <- proc.time()["elapsed"]
  fit <- tryCatch(fun(), error = function(e) e)
  list(fit = fit, cpu = unname(proc.time()["elapsed"] - start))
}
.a2_log_error <- function(e, i, cfg) {
  dir.create(cfg$error_dir, recursive = TRUE, showWarnings = FALSE)
  file <- file.path(cfg$error_dir, sprintf("err_D%s_K%s_rep%04d.txt", cfg$Dist, cfg$K, i))
  writeLines(c(sprintf("rep=%d, Dist=%s, K=%s", i, cfg$Dist, cfg$K),
               paste("message:", conditionMessage(e)),
               paste("call:", paste(deparse(conditionCall(e)), collapse = " "))), file)
  file
}
.a2_run_rep <- function(i, cfg) {
  args <- cfg[c("N", "Dist", "K", "n0_grid", "rho_mr", "pilot_sampling", "project_dir", "seed_base")]
  error_file <- NA_character_
  out <- tryCatch(do.call(A2PilotNdag_inner_grid, c(list(i = i), args)), error = function(e) {
    error_file <<- .a2_log_error(e, i, cfg)
    .a2_failed_rep(i, cfg, conditionMessage(e))
  })
  file <- file.path(cfg$inner_out_dir, sprintf("A2PilotNdag_rep%04d_D%s_K%s.rds", i, cfg$Dist, cfg$K))
  saveRDS(list(ok = isTRUE(out$rep_summary$ok[1L]), context = list(rep = i, Dist = cfg$Dist, K = cfg$K),
               timestamp = Sys.time(), result = out, error_file = error_file), file)
  out
}
.a2_summarize <- function(raw) {
  raw %>%
    group_by(Dist, DistLabel, PilotSampling, K, N, rho_N, n0) %>%
    summarise(B = sum(Success), Fail = sum(!Success), N_rem_mean = .a2_mean(N_rem[Success]),
              N_dag_mean = .a2_mean(N_dag[Success]), N_dag_sd = .a2_sd(N_dag[Success]),
              N_dag_se = if (B > 0L) N_dag_sd / sqrt(B) else NA_real_,
              N_dag_ratio_mean = .a2_mean(N_dag_ratio[Success]), N_dag_ratio_sd = .a2_sd(N_dag_ratio[Success]),
              N_dag_ratio_se = if (B > 0L) N_dag_ratio_sd / sqrt(B) else NA_real_,
              CPU_mean = .a2_mean(CPU[Success]), CPU_sd = .a2_sd(CPU[Success]),
              CPU_se = if (B > 0L) CPU_sd / sqrt(B) else NA_real_, .groups = "drop") %>%
    rename(N_rem = N_rem_mean, N_dag = N_dag_mean, N_dag_ratio = N_dag_ratio_mean, CPU = CPU_mean) %>%
    arrange(K, n0)
}
.a2_validate <- function(N, Dist, K, m, n0_grid, rho_mr, cores) {
  is_count <- function(x) length(x) == 1L && is.finite(x) && x >= 1 &&
    x == floor(x) && x <= .Machine$integer.max
  stopifnot(is_count(N), is_count(m), is_count(cores), length(K) > 0L,
            all(vapply(as.list(K), is_count, logical(1))), !anyDuplicated(K),
            length(n0_grid) > 0L, all(is.finite(n0_grid)), !anyDuplicated(n0_grid),
            all(n0_grid == floor(n0_grid)), all(n0_grid >= max(K) + 1L), all(n0_grid < N),
            length(rho_mr) == 1L, is.finite(rho_mr), rho_mr > 0, rho_mr < 1)
  .a2_dist(Dist)
  invisible(NULL)
}

# Fit a weighted pilot multinomial model and partition the remaining data.

.a2_fit_pilot <- function(X, y, n0, rho, pilot_sampling) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  if (!is.factor(y)) y <- factor(as.character(y), levels = as.character(sort(unique(y))))
  cls <- levels(y)
  N <- nrow(X)
  p <- ncol(X)
  n_class <- length(cls)
  if (n_class < 2L || n0 < n_class || n0 >= N) stop("Invalid pilot size or class count.")
  if (pilot_sampling == "uniform") {
    pilot_prob <- rep(1 / N, N)
  } else {
    counts <- table(factor(y, levels = cls))
    if (any(counts == 0L)) stop("Some classes have zero observations.")
    # Keep pi_i = 1 / (n_class * N_{y_i}) and the original normalization.
    pilot_prob <- 1 / (n_class * as.numeric(counts[as.character(y)]))
    pilot_prob <- pilot_prob / sum(pilot_prob)
  }
  idx0 <- NULL
  for (attempt in seq_len(100L)) {
    candidate <- sample.int(N, n0, replace = FALSE, prob = pilot_prob)
    if (length(unique(y[candidate])) == n_class) {
      idx0 <- candidate
      break
    }
  }
  if (is.null(idx0)) stop("Pilot subsample repeatedly missed some classes; increase n0.")
  idxR <- setdiff(seq_len(N), idx0)
  X0 <- X[idx0, , drop = FALSE]
  XR <- X[idxR, , drop = FALSE]
  y0 <- factor(y[idx0], levels = cls)
  yR <- factor(y[idxR], levels = cls)
  data0 <- data.frame(y = y0, X0, check.names = FALSE)
  colnames(data0) <- c("y", paste0("x", seq_len(p)))
  if (pilot_sampling == "uniform") {
    weights <- rep(1, length(idx0))
  } else {
    weights <- 1 / pilot_prob[idx0]
    weights <- weights / mean(weights)
  }
  # Retain multinom's original fitting defaults and weight limit.
  fit0 <- nnet::multinom(y ~ ., data = data0, weights = weights, trace = FALSE,
                        MaxNWts = max(10000, 10 * (p + 1L) * (n_class - 1L)))
  beta0 <- .a2_beta(fit0, cls, p)
  eta0 <- cbind(0, cbind(1, XR) %*% beta0)
  eta0 <- eta0 - apply(eta0, 1L, max)
  P0 <- exp(eta0)
  P0 <- P0 / rowSums(P0)
  colnames(P0) <- cls
  Omega <- .a2_partition(yR, cls, P0, rho)
  Omega_dag <- setdiff(seq_len(nrow(XR)), unlist(Omega, use.names = FALSE))
  list(idx0 = idx0, idxR = idxR, N.rem = nrow(XR), N.dag = length(Omega_dag),
       N.dag.ratio = length(Omega_dag) / N, pilot_sampling = pilot_sampling)
}

# One replicate: use the same dataset for all pilot sizes.

A2PilotNdag_inner_grid <- function(
    i, N = 1e6, Dist = 1L, K = 3L,
    n0_grid = c(100, 300, 500, 700, 900, 1200, 1500, 2000, 3000),
    rho_mr = 0.999, pilot_sampling = "proportional", project_dir = getwd(), seed_base = NULL) {
  pilot_sampling <- match.arg(pilot_sampling, c("proportional", "uniform"))
  .a2_threads()
  seed <- .a2_seed(seed_base, i - 1L)
  if (!is.null(seed)) set.seed(seed)
  dist_args <- .a2_dist(Dist)
  data <- do.call(DataGenSM, c(list(n = N, n_class = K + 1L, standardize_X = FALSE, seed = NULL),
                              dist_args[setdiff(names(dist_args), "DistLabel")]))
  if (nrow(as.matrix(data$X)) != N || length(data$y) != N) stop("Unexpected generated data dimensions.")
  cfg <- list(N = N, Dist = Dist, K = K, rho_mr = rho_mr, pilot_sampling = pilot_sampling)
  rows <- vector("list", length(n0_grid))
  for (j in seq_along(n0_grid)) {
    n0 <- as.integer(n0_grid[j])
    row <- .a2_row(i, n0, cfg)
    # One replication seed controls data generation and all successive pilot draws.
    result <- .a2_try_fit(function() .a2_fit_pilot(data$X, data$y, n0, rho_mr, pilot_sampling))
    row$CPU <- result$cpu
    if (inherits(result$fit, "error")) row$ErrorMsg <- conditionMessage(result$fit) else {
      row$N_rem <- result$fit$N.rem
      row$N_dag <- result$fit$N.dag
      row$N_dag_ratio <- result$fit$N.dag.ratio
      row$Success <- TRUE
    }
    rows[[j]] <- row
    rm(result)
  }
  raw <- bind_rows(rows)
  list(raw = raw, rep_summary = .a2_rep_summary(raw))
}

# Parallel repetitions: PSOCK on Windows and multicore on Unix by default.

A2PilotNdag_mclapply_grid <- function(
    m, N = 1e6, Dist = 1L, K = 3L,
    n0_grid = c(100, 300, 500, 700, 900, 1200, 1500, 2000, 3000),
    rho_mr = 0.999, pilot_sampling = "proportional", cores = 4L,
    backend = c("auto", "psock", "multicore"), project_dir = getwd(),
    inner_out_dir = file.path(project_dir, "A2PilotNdag_inner_rds"), seed_base = 12345) {
  .a2_validate(N, Dist, K, m, n0_grid, rho_mr, cores)
  stopifnot(length(K) == 1L)
  .a2_seed(seed_base)
  .a2_threads()
  pilot_sampling <- match.arg(pilot_sampling, c("proportional", "uniform"))
  backend <- match.arg(backend)
  if (backend == "auto") backend <- if (.Platform$OS.type == "windows") "psock" else "multicore"
  if (.Platform$OS.type == "windows") backend <- "psock"
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  dir.create(inner_out_dir, recursive = TRUE, showWarnings = FALSE)
  inner_out_dir <- normalizePath(inner_out_dir, mustWork = TRUE)
  cfg <- list(N = N, Dist = Dist, K = K, n0_grid = n0_grid, rho_mr = rho_mr,
              pilot_sampling = pilot_sampling, project_dir = project_dir, inner_out_dir = inner_out_dir,
              error_dir = file.path(inner_out_dir, "errors"), seed_base = seed_base)
  workers <- min(as.integer(cores), as.integer(m))
  if (workers == 1L) {
    results <- lapply(seq_len(m), .a2_run_rep, cfg = cfg)
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
      sys.source("Syn.R", .GlobalEnv)
      TRUE
    }, project_dir)
    export_env <- environment(A2PilotNdag_mclapply_grid)
    exports <- c("A2PilotNdag_inner_grid", ls(envir = export_env, all.names = TRUE, pattern = "^\\.a2_"))
    clusterExport(cl, exports, envir = export_env)
    results <- parLapply(cl, seq_len(m), .a2_run_rep, cfg = cfg)
  } else {
    results <- mclapply(seq_len(m), .a2_run_rep, cfg = cfg, mc.cores = workers, mc.preschedule = TRUE)
  }
  raw <- bind_rows(lapply(results, `[[`, "raw"))
  rep_summary <- bind_rows(lapply(results, `[[`, "rep_summary"))
  tag <- sprintf("D%s_K%s_m%d", Dist, K, m)
  saveRDS(raw, file.path(inner_out_dir, paste0("A2PilotNdag_raw_", tag, ".rds")))
  saveRDS(rep_summary, file.path(inner_out_dir, paste0("A2PilotNdag_rep_summary_", tag, ".rds")))
  cat(sprintf("[A2PilotNdag] K=%s | %d/%d successful pilot fits\n", K, sum(raw$Success), nrow(raw)))
  .a2_summarize(raw)
}

# Run the K grid and save cumulative summaries and run metadata.

RunA2PilotNdagStudy <- function(
    N = 1e6, Dist = 1L, K = c(3, 5, 7), m = 200L,
    n0_grid = c(100, 300, 500, 700, 900, 1200, 1500, 2000, 3000),
    rho_mr = 0.999, pilot_sampling = "proportional", cores = 4L, backend = "auto",
    project_dir = getwd(), out_dir = file.path(project_dir, "A2PilotNdag_runs"),
    log_to_console = TRUE, log_file = file.path(out_dir, "run.log"), seed_base = 12345) {
  .a2_validate(N, Dist, K, m, n0_grid, rho_mr, cores)
  .a2_seed(seed_base)
  pilot_sampling <- match.arg(pilot_sampling, c("proportional", "uniform"))
  backend <- match.arg(backend, c("auto", "psock", "multicore"))
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, mustWork = TRUE)
  if (!is.null(log_file) && !is.na(log_file)) dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  log <- function(...) {
    text <- paste0("[", format(Sys.time(), "%F %T"), "] ", paste0(..., collapse = ""))
    if (log_to_console) cat(text, "\n")
    if (!is.null(log_file) && !is.na(log_file)) cat(text, "\n", file = log_file, append = TRUE)
  }
  metadata <- list(N = N, Dist = Dist, K = K, m = m, n0_grid = n0_grid, rho_mr = rho_mr,
                   pilot_sampling = pilot_sampling, cores = cores, backend = backend, seed_base = seed_base,
                   source_md5 = tools::md5sum(file.path(project_dir, "Syn.R")), timestamp = Sys.time(),
                   session = capture.output(sessionInfo()), timing = "CPU stores elapsed seconds", ratio = "N_dag / N")
  saveRDS(metadata, file.path(out_dir, "metadata.rds"))
  start <- Sys.time()
  results <- vector("list", length(K))
  log(sprintf("RunA2PilotNdagStudy: %d tasks | Dist=%s | N=%s | m=%s", length(K), Dist, N, m))
  for (j in seq_along(K)) {
    tag <- sprintf("task-%03d_D%s_K%s", j, Dist, K[j])
    log("START ", tag)
    results[[j]] <- tryCatch(A2PilotNdag_mclapply_grid(
      m = m, N = N, Dist = Dist, K = K[j], n0_grid = n0_grid, rho_mr = rho_mr,
      pilot_sampling = pilot_sampling, cores = cores, backend = backend, project_dir = project_dir,
      inner_out_dir = file.path(out_dir, "replicates"), seed_base = .a2_seed(seed_base, 1000 * j)),
      error = function(e) {
        log("ERROR ", tag, ": ", conditionMessage(e))
        structure(list(error = conditionMessage(e)), class = "sim_error")
      })
    saveRDS(results[[j]], file.path(out_dir, paste0("A2PilotNdag_", tag, ".rds")))
    ok <- vapply(results, function(x) !is.null(x) && !inherits(x, "sim_error"), logical(1))
    if (any(ok)) saveRDS(bind_rows(results[ok]), file.path(out_dir, "A2PilotNdag_cumulative.rds"))
    log("DONE ", tag)
  }
  ok <- vapply(results, function(x) !is.null(x) && !inherits(x, "sim_error"), logical(1))
  if (!any(ok)) stop("All K tasks failed; see run.log and task records.")
  final <- bind_rows(results[ok])
  attr(final, "elapsed_hours") <- as.numeric(difftime(Sys.time(), start, units = "hours"))
  saveRDS(final, file.path(out_dir, "A2PilotNdag_final.rds"))
  write.csv(final, file.path(out_dir, "A2PilotNdag_final.csv"), row.names = FALSE)
  log(sprintf("ALL DONE: %.2f hours.", attr(final, "elapsed_hours")))
  final
}

# Plot mean counts or ratios; shaded bands show one Monte Carlo SD.

PlotA2PilotNdag <- function(res_a2, plot_ratio = FALSE, show_band = TRUE, log_y = FALSE) {
  yvar <- if (plot_ratio) "N_dag_ratio" else "N_dag"
  ylabel <- if (plot_ratio) expression(N["\u2020"] / N) else expression(N["\u2020"])
  data <- res_a2
  levels <- sort(unique(data$K))
  data$K <- factor(data$K, levels = levels, labels = paste("K =", levels))
  if (show_band) {
    data$y_lo <- pmax(data[[yvar]] - data[[paste0(yvar, "_sd")]], 0)
    data$y_hi <- data[[yvar]] + data[[paste0(yvar, "_sd")]]
  }
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = n0, y = .data[[yvar]],
                                           group = K, color = K, shape = K, linetype = K))
  if (show_band) plot <- plot + ggplot2::geom_ribbon(
    ggplot2::aes(ymin = y_lo, ymax = y_hi, fill = K), alpha = 0.12, colour = NA, na.rm = TRUE)
  reference <- if (plot_ratio) 1 else unique(data$N)
  plot <- plot + ggplot2::geom_line(linewidth = 0.9) + ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_x_continuous(breaks = sort(unique(data$n0))) +
    ggplot2::geom_hline(yintercept = reference, linetype = 2, linewidth = 0.5) +
    ggplot2::labs(x = expression(n[0]), y = ylabel, color = "Class", fill = "Class", shape = "Class", linetype = "Class") +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom", legend.box = "horizontal") +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1), fill = ggplot2::guide_legend(nrow = 1),
                    shape = ggplot2::guide_legend(nrow = 1), linetype = ggplot2::guide_legend(nrow = 1))
  if (log_y) plot <- plot + ggplot2::scale_y_log10()
  plot
}

# Simulation settings and calls
# Set seed_base to an integer for reproducible repetitions.
# Set options(A2PilotNdag.autorun = FALSE) before sourcing to skip the formal run.

N <- 1e6
Dist <- 1L
K <- c(3L, 5L, 7L)
m <- 200L
n0_grid <- c(100, 300, 500, 700, 900, 1200, 1500, 2000, 3000)
rho_mr <- 0.999
pilot_sampling <- "proportional"
cores <- 20L
seed_base <- NULL
project_dir <- getwd()
out_dir <- file.path(project_dir, "A2PilotNdag_runs")

if (isTRUE(getOption("A2PilotNdag.autorun", TRUE))) {
  res_a2 <- RunA2PilotNdagStudy(
    N = N, Dist = Dist, K = K, m = m, n0_grid = n0_grid, rho_mr = rho_mr,
    pilot_sampling = pilot_sampling, cores = cores, backend = "auto",
    project_dir = project_dir, out_dir = out_dir, seed_base = seed_base)
  print(attr(res_a2, "elapsed_hours"))

  p_a2_ndag <- PlotA2PilotNdag(res_a2, plot_ratio = FALSE, show_band = TRUE)
  print(p_a2_ndag)
  ggplot2::ggsave(file.path(out_dir, "A2PilotNdag_plot.pdf"), p_a2_ndag, width = 8, height = 5, encoding = "WinAnsi.enc")

  p_a2_ratio <- PlotA2PilotNdag(res_a2, plot_ratio = TRUE, show_band = TRUE)
  print(p_a2_ratio)
  ggplot2::ggsave(file.path(out_dir, "A2PilotNdag_ratio_plot.pdf"), p_a2_ratio, width = 8, height = 5, encoding = "WinAnsi.enc")
}
