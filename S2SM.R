rm(list = ls())
gc()

# Setting 2 studies pilot allocation under a fixed subsampling budget.


suppressPackageStartupMessages({
  library(MASS)
  library(mvtnorm)
  library(nnet)
  library(dplyr)
  library(parallel)
  library(ggplot2)
  library(ggh4x)
  library(latex2exp)
})

# Run this script from the directory containing the companion files.
source("Syn.R")
source("ossrmp.R")
source("ossrmp_pre.R")
source("ossrmp_pos.R")

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

make_seed_s2 <- function(seed, offset = 0) {
  if (is.null(seed)) return(NULL)
  value <- as.double(seed) + offset
  if (length(value) != 1L || !is.finite(value) || value < 0 ||
      value > .Machine$integer.max || value != floor(value)) stop("Invalid seed.")
  as.integer(value)
}

make_budget_grid_s2 <- function(total_budget, alpha_grid) {
  stopifnot(length(total_budget) == 1L, is.finite(total_budget),
            total_budget == floor(total_budget), total_budget >= 2,
            length(alpha_grid) > 0L, all(is.finite(alpha_grid)),
            all(alpha_grid > 0 & alpha_grid < 1))
  alpha <- round(alpha_grid, 4)
  n0 <- round(total_budget * alpha)
  grid <- data.frame(alpha = as.numeric(alpha), n0 = as.numeric(n0),
                     n_dag = as.numeric(total_budget - n0))
  if (any(grid$n0 <= 0 | grid$n_dag <= 0)) stop("Both subsample sizes must be positive.")
  grid
}

to_sum_constraint_s2 <- function(beta) {
  B <- cbind(0, beta)
  sweep(B, 1, rowMeans(B), "-")[, -1, drop = FALSE]
}

evaluate_mse_s2 <- function(beta, truth, use_mean = TRUE) {
  theta <- as.vector(beta)
  if (is.null(dim(beta)) || length(theta) != length(truth) ||
      any(!is.finite(theta)) || any(!is.finite(truth))) stop("Invalid coefficients or target vector.")
  error <- (theta - truth)^2
  if (use_mean) mean(error) else sum(error)
}

safe_mean_s2 <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

safe_sd_s2 <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) > 1L) sd(x) else NA_real_
}

# CPU records elapsed fitting time in seconds, as in the original script.
fit_once_s2 <- function(fun) {
  start <- proc.time()[["elapsed"]]
  fit <- tryCatch(fun(), error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - start
  if (inherits(fit, "error")) {
    return(list(ok = FALSE, fit = NULL, cpu = elapsed, err = conditionMessage(fit)))
  }
  if (is.null(fit$beta) || any(!is.finite(as.numeric(fit$beta)))) {
    return(list(ok = FALSE, fit = fit, cpu = elapsed, err = "Non-finite beta returned."))
  }
  list(ok = TRUE, fit = fit, cpu = elapsed, err = NA_character_)
}

result_row_s2 <- function(context, method, result, truth = NULL,
                          use_mean = TRUE, fit_source = method) {
  fit <- result$fit
  mse <- NA_real_
  if (result$ok) {
    value <- tryCatch(evaluate_mse_s2(fit$beta, truth, use_mean), error = function(e) e)
    if (inherits(value, "error")) {
      result$ok <- FALSE
      result$err <- conditionMessage(value)
    } else mse <- value
  }
  n_hr <- n_sel <- NA_real_
  if (result$ok) {
    if (!is.null(fit$Omega.dag)) n_hr <- as.numeric(length(fit$Omega.dag))
    n_sel <- if (!is.null(fit$n.sel)) as.numeric(fit$n.sel) else
      if (!is.null(fit$idx.sel)) as.numeric(length(fit$idx.sel)) else NA_real_
  }
  data.frame(context, Method = method, FitSource = fit_source, N_hr = n_hr,
             n_sel = n_sel, MSE = mse, CPU = result$cpu, Success = result$ok,
             Converged = result$ok && isTRUE(fit$converged), ErrorMsg = result$err,
             stringsAsFactors = FALSE, row.names = NULL)
}

replication_summary_s2 <- function(raw, ok = TRUE) {
  data.frame(rep = raw$rep[1], K = raw$K[1], N = raw$N[1], n_rows = nrow(raw),
             n_success = sum(raw$Success), n_fail = sum(!raw$Success), ok = ok)
}

# ------------------------------------------------------------
# One replication
# ------------------------------------------------------------

run_replication_s2 <- function(i, N, K, total_budget = 6000,
                               alpha_grid = seq(0.1, 0.9, 0.1), rho_mr = 0.999,
                               tol = 1e-8, maxit = 200, verbose_fit = FALSE,
                               seed_base = 12345, use_mean_mse = TRUE) {
  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  if (!is.null(seed_base)) set.seed(make_seed_s2(seed_base, i - 1))
  grid <- make_budget_grid_s2(total_budget, alpha_grid)
  dat <- DataGenSM(n = N, n_class = K + 1, x_dist = "gaussian", cov_type = "ar1",
                   rho = 0.5, standardize_X = FALSE, seed = NULL)
  X <- as.matrix(dat$X)
  y <- dat$y
  truth <- as.vector(dat$theta_true)
  truth_sum <- as.vector(dat$theta_true_sum)
  if (nrow(X) != N || length(y) != N ||
      length(truth) != (ncol(X) + 1L) * K || length(truth_sum) != length(truth) ||
      any(!is.finite(truth)) || any(!is.finite(truth_sum))) stop("Invalid DataGenSM output.")
  rows <- vector("list", 4L * nrow(grid))
  at <- 0L

  for (jj in seq_len(nrow(grid))) {
    n0 <- grid$n0[jj]
    n_dag <- grid$n_dag[jj]
    seed_now <- make_seed_s2(seed_base, 100000 * i + 1000 * jj)
    context <- list(rep = i, Dist = 1L, DistLabel = "D1-Gaussian-AR1", K = K,
                    N = as.numeric(N), total_budget = as.numeric(total_budget),
                    alpha = grid$alpha[jj], n0 = n0, n_dag = n_dag)

    fits <- list(
      OSSRMP = fit_once_s2(function() ossrmp(
        X = X, y = y, n0 = n0, n_dag = n_dag, rho = rho_mr,
        RB = TRUE, seed = make_seed_s2(seed_now, 11),
        tol = tol, maxit = maxit, verbose = verbose_fit
      )),
      `OSSRMP-P` = fit_once_s2(function() ossrmp_pre(
        X = X, y = y, n0 = n0, n_dag = n_dag, rho = rho_mr,
        calibrated = TRUE, seed = make_seed_s2(seed_now, 21),
        tol = tol, maxit = maxit, verbose = verbose_fit
      )),
      `OSSRMP-Pos` = fit_once_s2(function() ossrmp_poisson(
        X = X, y = y, n0 = n0, n_dag = n_dag, rho = rho_mr,
        RB = TRUE, seed = make_seed_s2(seed_now, 31),
        tol = tol, maxit = maxit, verbose = verbose_fit
      ))
    )

    for (method in names(fits)) {
      result <- fits[[method]]
      at <- at + 1L
      rows[[at]] <- result_row_s2(context, method, result, truth, use_mean_mse)
      if (method == "OSSRMP") {
        derived <- result
        if (result$ok) {
          beta_sum <- tryCatch(to_sum_constraint_s2(result$fit$beta), error = function(e) e)
          if (inherits(beta_sum, "error")) {
            derived$ok <- FALSE
            derived$err <- conditionMessage(beta_sum)
          } else derived$fit$beta <- beta_sum
        }
        at <- at + 1L
        rows[[at]] <- result_row_s2(context, "OSSRMP-Sum", derived,
                                    truth_sum, use_mean_mse, "OSSRMP")
      }
    }
  }
  raw <- bind_rows(rows)
  list(raw = raw, rep_summary = replication_summary_s2(raw))
}

# Keep one failed row per method and allocation when a replication fails.
failed_replication_s2 <- function(i, cfg, message) {
  grid <- make_budget_grid_s2(cfg$total_budget, cfg$alpha_grid)
  methods <- c("OSSRMP", "OSSRMP-Sum", "OSSRMP-P", "OSSRMP-Pos")
  rows <- list()
  at <- 0L
  for (j in seq_len(nrow(grid))) {
    context <- list(rep = i, Dist = 1L, DistLabel = "D1-Gaussian-AR1", K = cfg$K,
                    N = as.numeric(cfg$N), total_budget = as.numeric(cfg$total_budget),
                    alpha = grid$alpha[j], n0 = grid$n0[j], n_dag = grid$n_dag[j])
    for (method in methods) {
      at <- at + 1L
      source <- if (method == "OSSRMP-Sum") "OSSRMP" else method
      rows[[at]] <- result_row_s2(context, method,
        list(ok = FALSE, fit = NULL, cpu = NA_real_, err = message), fit_source = source)
    }
  }
  raw <- bind_rows(rows)
  list(raw = raw, rep_summary = replication_summary_s2(raw, FALSE))
}

# Summaries retain successful fits; convergence is reported separately.
summarize_results_s2 <- function(raw) {
  raw %>%
    group_by(Method, FitSource, Dist, DistLabel, K, N, total_budget, alpha, n0, n_dag) %>%
    summarise(
      B = sum(Success & is.finite(MSE)), Fail = sum(!Success | !is.finite(MSE)),
      ConvRate = safe_mean_s2(as.numeric(Converged[Success])),
      MSE_m = safe_mean_s2(MSE[Success]), CPU_m = safe_mean_s2(CPU[Success]),
      MSE_sd = safe_sd_s2(MSE[Success]), CPU_sd = safe_sd_s2(CPU[Success]),
      MSE_se = if (B > 0) MSE_sd / sqrt(B) else NA_real_,
      CPU_se = if (B > 0) CPU_sd / sqrt(B) else NA_real_,
      N_hr = safe_mean_s2(N_hr[Success]), n_sel = safe_mean_s2(n_sel[Success]),
      .groups = "drop"
    ) %>%
    rename(MSE = MSE_m, CPU = CPU_m) %>%
    arrange(K, alpha, Method)
}

# ------------------------------------------------------------
# Independent replications
# ------------------------------------------------------------

run_repetitions_s2 <- function(m, N, K, total_budget, alpha_grid,
                               rho_mr = 0.999, tol = 1e-8, maxit = 200,
                               verbose_fit = FALSE, cores = 4, backend = "auto",
                               project_dir = getwd(), inner_out_dir = "S2SM_inner_rds",
                               seed_base = 12345, use_mean_mse = TRUE) {
  backend <- match.arg(backend, c("auto", "psock", "multicore"))
  if (backend == "auto") backend <- if (.Platform$OS.type == "windows") "psock" else "multicore"
  if (backend == "multicore" && .Platform$OS.type == "windows") stop("Use backend = 'psock' on Windows.")
  cores <- min(as.integer(cores), as.integer(m))
  if (is.na(cores) || cores < 1L || m < 1L) stop("m and cores must be positive integers.")
  dir.create(inner_out_dir, recursive = TRUE, showWarnings = FALSE)
  cfg <- list(N = N, K = K, total_budget = total_budget, alpha_grid = alpha_grid,
              rho_mr = rho_mr, tol = tol, maxit = maxit, verbose_fit = verbose_fit,
              seed_base = seed_base, use_mean_mse = use_mean_mse)

  one_rep <- function(i) {
    cat(sprintf("[Setting 2] rep %d/%d | K=%d\n", i, m, K))
    out <- tryCatch(do.call(run_replication_s2, c(list(i = i), cfg)), error = function(e) e)
    ok <- !inherits(out, "error")
    error <- if (ok) NULL else conditionMessage(out)
    if (!ok) out <- failed_replication_s2(i, cfg, error)
    tag <- sprintf("rep%03d_K%d", i, K)
    stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
    record <- list(ok = ok, context = list(rep = i, K = K), timestamp = Sys.time(),
                   result = out, error = error)
    saveRDS(record, file.path(inner_out_dir,
      sprintf("S2SM_%s_%s_%s.rds", tag, if (ok) "ok" else "err", stamp)))
    failures <- out$raw[!out$raw$Success, c("Method", "alpha", "ErrorMsg"), drop = FALSE]
    if (nrow(failures)) write.table(failures, file.path(inner_out_dir, paste0("errors_", tag, ".tsv")),
                                    sep = "\t", quote = TRUE, row.names = FALSE)
    out
  }

  if (cores == 1L) {
    results <- lapply(seq_len(m), one_rep)
  } else if (backend == "psock") {
    cl <- makeCluster(cores, type = "PSOCK", outfile = "")
    on.exit(stopCluster(cl), add = TRUE)
    exports <- c("make_seed_s2", "make_budget_grid_s2", "to_sum_constraint_s2", "evaluate_mse_s2",
                 "fit_once_s2", "result_row_s2", "replication_summary_s2",
                 "failed_replication_s2", "run_replication_s2")
    clusterExport(cl, exports, envir = environment(run_replication_s2))
    clusterCall(cl, function(path) {
      Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
      setwd(path)
      suppressPackageStartupMessages({
        library(MASS)
        library(mvtnorm)
        library(nnet)
        library(dplyr)
      })
      source("Syn.R")
      source("ossrmp.R")
      source("ossrmp_pre.R")
      source("ossrmp_pos.R")
      NULL
    }, project_dir)
    results <- parLapply(cl, seq_len(m), one_rep)
  } else {
    results <- mclapply(seq_len(m), one_rep, mc.cores = cores,
                        mc.preschedule = TRUE, mc.set.seed = is.null(seed_base))
  }
  results <- lapply(seq_len(m), function(i) {
    out <- results[[i]]
    if (inherits(out, "try-error") || is.null(out)) {
      failed_replication_s2(i, cfg, if (is.null(out)) "Worker returned NULL." else as.character(out))
    } else out
  })
  raw <- bind_rows(lapply(results, `[[`, "raw"))
  reps <- bind_rows(lapply(results, `[[`, "rep_summary"))
  stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  saveRDS(raw, file.path(inner_out_dir, sprintf("S2SM_raw_K%d_m%d_%s.rds", K, m, stamp)))
  saveRDS(reps, file.path(inner_out_dir, sprintf("S2SM_rep_summary_K%d_m%d_%s.rds", K, m, stamp)))
  summarize_results_s2(raw)
}

# ------------------------------------------------------------
# Setting 2
# ------------------------------------------------------------

run_setting2 <- function(N = 1e6, K = c(3, 5, 7), total_budget = 6000,
                         alpha_grid = seq(0.1, 0.9, 0.1), m = 200,
                         rho_mr = 0.999, tol = 1e-8, maxit = 200,
                         verbose_fit = FALSE, cores = 4, backend = "auto",
                         project_dir = getwd(), out_dir = file.path(project_dir, "S2SM_runs"),
                         log_to_console = TRUE, log_file = file.path(out_dir, "run.log"),
                         seed_base = 12345, use_mean_mse = TRUE) {
  stopifnot(N == floor(N), N > 1, length(K) > 0, all(K >= 1 & K == floor(K)),
            m >= 1, m == floor(m), cores >= 1, cores == floor(cores),
            rho_mr > 0, rho_mr < 1)
  grid <- make_budget_grid_s2(total_budget, alpha_grid)
  if (any(grid$n0 > N)) stop("The pilot sample size exceeds N.")
  backend <- match.arg(backend, c("auto", "psock", "multicore"))
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(log_file) && !is.na(log_file)) dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  log <- function(...) {
    message <- paste0("[", format(Sys.time(), "%F %T"), "] ", paste0(..., collapse = ""))
    if (log_to_console) cat(message, "\n")
    if (!is.null(log_file) && !is.na(log_file)) cat(message, "\n", file = log_file, append = TRUE)
  }
  metadata <- list(N = N, K = K, total_budget = total_budget, budget_grid = grid,
                   m = m, rho_mr = rho_mr, tol = tol, maxit = maxit, cores = cores,
                   backend = backend, seed_base = seed_base, use_mean_mse = use_mean_mse,
                   CPU_definition = "elapsed fitting time in seconds", session = sessionInfo())
  saveRDS(metadata, file.path(out_dir, "S2SM_metadata.rds"))
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  results <- vector("list", length(K))
  start <- Sys.time()
  log("Starting ", length(K), " tasks; output: ", out_dir)

  for (idx in seq_along(K)) {
    tag <- sprintf("task-%03d_K%d", idx, K[idx])
    log("START ", tag)
    task_start <- Sys.time()
    result <- tryCatch(run_repetitions_s2(
      m = m, N = N, K = K[idx], total_budget = total_budget, alpha_grid = alpha_grid,
      rho_mr = rho_mr, tol = tol, maxit = maxit, verbose_fit = verbose_fit,
      cores = cores, backend = backend, project_dir = project_dir,
      inner_out_dir = file.path(out_dir, "inner_rds"),
      seed_base = make_seed_s2(seed_base, 1000 * idx), use_mean_mse = use_mean_mse
    ), error = function(e) {
      log("ERROR ", tag, ": ", conditionMessage(e))
      structure(list(error = conditionMessage(e)), class = "sim_error")
    })
    results[[idx]] <- result
    saveRDS(result, file.path(out_dir, paste0("S2SM_", tag, ".rds")))
    valid <- vapply(results, is.data.frame, logical(1))
    if (any(valid)) saveRDS(bind_rows(results[valid]), file.path(out_dir, "S2SM_cumulative.rds"))
    log("DONE ", tag, " (", round(as.numeric(difftime(Sys.time(), task_start, units = "mins")), 2), " min)")
  }
  valid <- vapply(results, is.data.frame, logical(1))
  if (!any(valid)) stop("All tasks failed; see run.log and the task files.")
  final <- bind_rows(results[valid])
  attr(final, "elapsed_hours") <- as.numeric(difftime(Sys.time(), start, units = "hours"))
  attr(final, "failed_K") <- K[!valid]
  saveRDS(final, file.path(out_dir, "S2SM_final.rds"))
  write.csv(final, file.path(out_dir, "S2SM_final.csv"), row.names = FALSE)
  log("ALL DONE: ", sum(valid), "/", length(K), " tasks completed.")
  final
}

# ------------------------------------------------------------
# Plotting
# ------------------------------------------------------------

plot_setting2 <- function(results, show_band = FALSE, conf_level = 0.95, eps = 1e-20) {
  stopifnot(conf_level > 0, conf_level < 1, eps > 0)
  methods <- c("OSSRMP", "OSSRMP-P", "OSSRMP-Pos")
  labels <- c(
    latex2exp::TeX(r"($\widehat{\vartheta}_{MR}$)"),
    latex2exp::TeX(r"($\widehat{\vartheta}_{MR,P}$)"),
    latex2exp::TeX(r"($\widehat{\vartheta}_{MR}^{Pos}$)")
  )
  colors <- c("OSSRMP" = "#1B9E77", "OSSRMP-P" = "#7570B3", "OSSRMP-Pos" = "#E7298A")
  lines <- c("OSSRMP" = "solid", "OSSRMP-P" = "solid", "OSSRMP-Pos" = "dotdash")
  shapes <- c("OSSRMP" = 16, "OSSRMP-P" = 3, "OSSRMP-Pos" = 4)
  z <- qnorm(1 - (1 - conf_level) / 2)
  data <- results %>%
    filter(Method %in% methods, is.finite(MSE)) %>%
    mutate(Method = factor(Method, levels = methods),
           log_MSE = log(pmax(MSE, eps)),
           log_MSE_lo = log(pmax(MSE - z * MSE_se, eps)),
           log_MSE_hi = log(pmax(MSE + z * MSE_se, eps)))
  if (!nrow(data)) stop("No finite MSE values are available for plotting.")
  k_values <- sort(unique(data$K))
  data$KLabel <- factor(data$K, levels = k_values, labels = paste("K =", k_values))
  p <- ggplot(data, aes(alpha, log_MSE, colour = Method, fill = Method,
                         shape = Method, linetype = Method, group = Method))
  if (show_band) {
    bands <- filter(data, is.finite(log_MSE_lo), is.finite(log_MSE_hi))
    p <- p + geom_ribbon(data = bands, aes(ymin = log_MSE_lo, ymax = log_MSE_hi),
                          alpha = 0.12, colour = NA)
  }
  p +
    geom_line(linewidth = 0.85, na.rm = TRUE) +
    geom_point(size = 2.2, na.rm = TRUE) +
    ggh4x::facet_grid2(cols = vars(KLabel), scales = "free_y", independent = "y") +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
    scale_colour_manual(values = colors, breaks = methods, labels = labels) +
    scale_fill_manual(values = colors, breaks = methods, labels = labels) +
    scale_linetype_manual(values = lines, breaks = methods, labels = labels) +
    scale_shape_manual(values = shapes, breaks = methods, labels = labels) +
    labs(x = expression(alpha == n[0] / (n[0] + n[dagger])), y = "Log(MSE)",
         colour = "Method", fill = "Method", shape = "Method", linetype = "Method") +
    theme_bw(base_size = 14) +
    theme(panel.grid.major = element_line(linewidth = 0.3, linetype = "dashed", colour = "grey85"),
          panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(face = "bold", size = 13), legend.position = "bottom",
          legend.box = "horizontal", legend.title = element_text(size = 13),
          legend.text = element_text(size = 14), plot.margin = margin(5.5, 10, 5.5, 10)) +
    guides(colour = guide_legend(nrow = 1, byrow = TRUE), fill = guide_legend(nrow = 1, byrow = TRUE),
           shape = guide_legend(nrow = 1, byrow = TRUE), linetype = guide_legend(nrow = 1, byrow = TRUE))
}

# ------------------------------------------------------------
# Parameters and simulation call
# ------------------------------------------------------------

N <- 1e6
K <- c(3, 5, 7)
total_budget <- 6000
alpha_grid <- seq(0.1, 0.9, 0.1)
m <- 200
rho_mr <- 0.999
tol <- 1e-8
maxit <- 200
cores <- 20
backend <- "auto"
seed_base <- 12345
use_mean_mse <- TRUE
show_band <- FALSE
project_dir <- getwd()
out_dir <- file.path(project_dir, "S2SM_runs")

# For a short trial, use N = 2e4, K = 3, total_budget = 600,
# alpha_grid = c(0.2, 0.5, 0.8), m = 2, cores = 2, and a separate output directory.

start_time <- Sys.time()
res_s2 <- run_setting2(
  N = N, K = K, total_budget = total_budget, alpha_grid = alpha_grid, m = m,
  rho_mr = rho_mr, tol = tol, maxit = maxit, cores = cores, backend = backend,
  project_dir = project_dir, out_dir = out_dir, seed_base = seed_base,
  use_mean_mse = use_mean_mse
)
print(difftime(Sys.time(), start_time, units = "hours"))

p_s2 <- plot_setting2(res_s2, show_band = show_band)
print(p_s2)
ggsave(file.path(out_dir, "S2_logMSE_bigpanel_noSum.png"), p_s2,
       width = 13, height = 4.8, dpi = 300)
