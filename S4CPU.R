# Setting 4: compare total fitting time and computational cost by stage.
# Run repeated simulations over N, with pilot and high-resolution budgets
# proportional to N. Save replicate results, summaries, and timing plots.
# CPU is user.self + sys.self; Elapsed is wall-clock time, both in seconds.
# Primary stages are disjoint. Diagnostic timers overlap these stages.

rm(list = ls())
gc()

library(MASS)
library(mvtnorm)
library(nnet)
library(parallel)

# Run from the directory containing these exact companion method files.
project_dir <- getwd()
method_files <- c("Syn.R", "OSS.R", "ossrmp.R", "ossrmp_pre.R", "ossrmp_pos.R")
method_env <- new.env(parent = globalenv())
for (file in method_files) source(file.path(project_dir, file), local = method_env)

# Timing helpers

.S4CPU_WORKER <- NULL
.s4_primary <- c("Setup", "Pilot", "PredictionDesign", "FullProb", "Partition",
                 "SamplingScores", "AuxiliaryMean", "Sampling", "RBCalibration",
                 "Iteration", "Inference", "Output")
.s4_diagnostic <- c("ScanTotal", "HRGradientOnly", "HRScoreJacobian")

.s4_text <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
.s4_block <- function(x) as.call(c(list(as.name("{")), x))
.s4_slice <- function(x, from, to) if (from <= to) x[seq.int(from, to)] else list()
.s4_lhs <- function(x) {
  if (is.call(x) && identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]))
    as.character(x[[2L]]) else ""
}
.s4_is_if <- function(x, condition) {
  is.call(x) && identical(x[[1L]], as.name("if")) &&
    identical(.s4_text(x[[2L]]), condition)
}
.s4_is_for <- function(x, variable) {
  is.call(x) && identical(x[[1L]], as.name("for")) &&
    identical(x[[2L]], as.name(variable))
}
.s4_wrap <- function(stage, x) {
  if (length(x)) as.call(list(as.name(".s4_time"), stage, .s4_block(x))) else NULL
}
.s4_group <- function(x, stages) {
  if (!length(x)) return(list())
  stopifnot(length(x) == length(stages))
  starts <- c(1L, which(stages[-1L] != stages[-length(stages)]) + 1L)
  ends <- c(starts[-1L] - 1L, length(x))
  lapply(seq_along(starts), function(i) .s4_wrap(stages[starts[i]], x[starts[i]:ends[i]]))
}
.s4_anchor <- function(hits, label, after = 0L) {
  hits <- hits[hits > after]
  if (!length(hits)) stop("Timing anchor not found: ", label,
                         ". Use the companion method versions supplied with this script.")
  hits[1L]
}
.s4_find <- function(x, target, after = 0L) {
  .s4_anchor(which(vapply(x, .s4_lhs, character(1)) == target), target, after)
}
.s4_cond <- function(x, condition, after = 0L) {
  .s4_anchor(which(vapply(x, .s4_is_if, logical(1), condition)), condition, after)
}

# Wrap existing expressions without changing formulas or random draws.
# Evaluation in the caller's frame preserves the method's local assignments.
.s4_instrument_scan <- function(loop, prediction = FALSE) {
  scan <- as.list(loop[[4L]])[-1L]
  high <- .s4_cond(scan, "length(high.pos) > 0")
  branch <- as.list(scan[[high]][[3L]])[-1L]
  lhs <- vapply(branch, .s4_lhs, character(1))
  labels <- rep("SamplingScores", length(branch))
  labels[lhs %in% c("rows.high", "Omega.dag", "N.dag")] <- "Partition"
  if (prediction) {
    labels[lhs %in% c("psi.block", "g.block", "gsum")] <- "AuxiliaryMean"
    labels[lhs %in% c("S_pg_all", "S_gg_all")] <- "PredictionDesign"
    labels[vapply(branch, .s4_is_if, logical(1), "is.null(gsum)")] <- "AuxiliaryMean"
  } else {
    labels[vapply(branch, .s4_is_if, logical(1), "RB == TRUE")] <- "AuxiliaryMean"
  }
  scan[[high]][[3L]] <- .s4_block(.s4_group(branch, labels))
  prob_end <- which(vapply(scan, function(x) grepl("^colnames\\(P", .s4_text(x)), logical(1)))
  if (length(prob_end) != 1L || prob_end >= high)
    stop("Cannot identify the scan probability block.")
  loop[[4L]] <- .s4_block(list(
    .s4_wrap("FullProb", scan[seq_len(prob_end)]),
    .s4_wrap("Partition", .s4_slice(scan, prob_end + 1L, high - 1L)),
    scan[[high]],
    .s4_wrap("Partition", .s4_slice(scan, high + 1L, length(scan)))))
  loop
}
.s4_instrument_eval <- function(expr) {
  fun <- expr[[3L]]
  code <- as.list(fun[[3L]])[-1L]
  hits <- if (.s4_lhs(expr) == ".oss_eval")
    which(vapply(code, .s4_is_if, logical(1), "!poisson")) else
      which(vapply(code, .s4_is_for, logical(1), "tt"))
  if (length(hits) != 1L) stop("Cannot identify the high-resolution gradient loop.")
  stage <- quote(if (add_jac) "HRScoreJacobian" else "HRGradientOnly")
  code[[hits]] <- .s4_wrap(stage, list(code[[hits]]))
  fun[[3L]] <- .s4_block(code)
  expr[[3L]] <- fun
  expr
}
.s4_instrument <- function(fun, kind = c("ordinary", "mr", "prediction")) {
  kind <- match.arg(kind)
  code <- as.list(body(fun))[-1L]
  pilot <- .s4_cond(code, 'pilot_sampling == "uniform"')
  beta <- .s4_find(code, "beta0")
  inference <- .s4_find(code, "Mhat", beta)
  output <- .s4_find(code, "out", inference)
  if (!identical(.s4_text(code[[length(code)]]), "return(out)"))
    stop("Unexpected method output; timing instrumentation stopped.")
  wrap <- function(stage, from, to) .s4_wrap(stage, .s4_slice(code, from, to))
  pieces <- list(wrap("Setup", 1L, pilot - 1L))
  if (kind == "ordinary") {
    scan <- .s4_find(code, "X", beta)
    scores <- .s4_find(code, "xnorm", scan)
    sampling <- .s4_cond(code, "!poisson", scores)
    iteration <- .s4_find(code, ".bread_pilot", sampling)
    pieces <- c(pieces, list(wrap("Pilot", pilot, scan - 1L),
                            wrap("FullProb", scan, scores - 1L),
                            wrap("SamplingScores", scores, sampling - 1L),
                            wrap("Sampling", sampling, iteration - 1L)))
    gradient <- .s4_find(code, ".oss_eval", iteration - 1L)
  } else {
    partition <- .s4_find(code, "Omega", beta)
    if (kind == "prediction") {
      design <- .s4_find(code, "w0.design", beta)
      pieces <- c(pieces, list(wrap("Pilot", pilot, design - 1L),
                              wrap("PredictionDesign", design, partition - 1L)))
    } else pieces <- c(pieces, list(wrap("Pilot", pilot, partition - 1L)))
    hits <- which(vapply(code, .s4_is_for, logical(1), "st"))
    scan <- .s4_anchor(hits, "full-data scan loop", partition)
    pieces <- c(pieces, list(wrap("Partition", partition, scan - 1L),
                            .s4_wrap("ScanTotal", list(.s4_instrument_scan(
                              code[[scan]], kind == "prediction")))))
    if (kind == "prediction") {
      sampling <- .s4_cond(code, "calibrated_p", scan)
    } else if (any(vapply(code, .s4_lhs, character(1))[-seq_len(scan)] == "des.out")) {
      sampling <- .s4_find(code, "des.out", scan)
    } else {
      hits <- which(vapply(code, function(x) is.call(x) &&
        identical(x[[1L]], as.name("if")) &&
        grepl("sum\\(scr.dag\\)", .s4_text(x[[2L]])), logical(1)))
      sampling <- .s4_anchor(hits, "sampling", scan)
    }
    tail_scan <- .s4_slice(code, scan + 1L, sampling - 1L)
    labels <- rep("Partition", length(tail_scan))
    if (kind == "prediction") {
      lhs <- vapply(tail_scan, .s4_lhs, character(1))
      labels[lhs == "gbar"] <- "AuxiliaryMean"
      labels[lhs == "A.pilot"] <- "PredictionDesign"
    } else {
      labels[vapply(tail_scan, .s4_is_if, logical(1), "RB == TRUE")] <- "AuxiliaryMean"
    }
    pieces <- c(pieces, .s4_group(tail_scan, labels))
    response <- .s4_find(code, "Yb.sel", sampling)
    iteration <- .s4_find(code, ".bread_base", response)
    if (kind == "prediction") {
      # A calibrated design may require a second high-resolution pass.
      pieces <- c(pieces, list(wrap("PredictionDesign", sampling, sampling)))
      sampling <- sampling + 1L
    }
    calibration <- if (kind == "prediction")
      quote(if (calibrated) "RBCalibration" else "Sampling") else
        quote(if (RB) "RBCalibration" else "Sampling")
    pieces <- c(pieces, list(wrap("Sampling", sampling, response),
                            wrap(calibration, response + 1L, iteration - 1L)))
    gradient <- .s4_find(code, ".mr_eval", iteration - 1L)
  }
  code[[gradient]] <- .s4_instrument_eval(code[[gradient]])
  init <- as.list(substitute({
    .s4_clock <- new.env(parent = emptyenv())
    .s4_clock$times <- matrix(0, 2L, length(STAGES), dimnames = list(c("CPU", "Elapsed"), STAGES))
    .s4_clock$calls <- setNames(integer(ncol(.s4_clock$times)), colnames(.s4_clock$times))
    .s4_time <- function(stage, code) {
      force(stage)
      start <- proc.time()
      on.exit({
        delta <- proc.time() - start
        .s4_clock$times[, stage] <- .s4_clock$times[, stage] +
          c(unname(delta["user.self"] + delta["sys.self"]), unname(delta["elapsed"]))
        .s4_clock$calls[stage] <- .s4_clock$calls[stage] + 1L
      })
      eval(substitute(code), parent.frame())
    }
  }, list(STAGES = c(.s4_primary, .s4_diagnostic))))[-1L]
  finish <- as.list(quote({
    out$timing <- .s4_clock$times
    out$timing_calls <- .s4_clock$calls
    out$iterations <- it
    out$pilot_converged <- identical(fit0$convergence, 0L)
    out$q <- d * K
    out$d_g <- 1L + K + d * K
    return(out)
  }))[-1L]
  body(fun) <- .s4_block(c(init, pieces, list(
    wrap("Iteration", iteration, inference - 1L),
    wrap("Inference", inference, output - 1L),
    wrap("Output", output, length(code) - 1L)), finish))
  fun
}

LoadS4CPU <- function(project_dir = getwd(), method_env = NULL) {
  for (package in c("MASS", "mvtnorm", "nnet"))
    suppressPackageStartupMessages(library(package, character.only = TRUE))
  if (is.null(method_env)) {
    files <- c("Syn.R", "OSS.R", "ossrmp.R", "ossrmp_pre.R", "ossrmp_pos.R")
    paths <- file.path(normalizePath(project_dir, mustWork = TRUE), files)
    if (any(!file.exists(paths))) stop("Missing files: ", paste(files[!file.exists(paths)], collapse = ", "))
    method_env <- new.env(parent = globalenv())
    for (file in paths) sys.source(file, envir = method_env)
  }
  method_env$OSS_timed <- .s4_instrument(method_env$OSS, "ordinary")
  method_env$ossrmp_timed <- .s4_instrument(method_env$ossrmp, "mr")
  method_env$ossrmp_pre_timed <- .s4_instrument(method_env$ossrmp_pre, "prediction")
  method_env$ossrmp_poisson_timed <- .s4_instrument(method_env$ossrmp_poisson, "mr")
  method_env$warmed <- FALSE
  method_env
}

# Simulation helpers

.s4_methods <- function(include_uncalibrated) {
  methods <- list(
    O = list(fun = "OSS_timed", args = list(poisson = FALSE)),
    `O-Pos` = list(fun = "OSS_timed", args = list(poisson = TRUE)),
    MR = list(fun = "ossrmp_timed", args = list(RB = TRUE)),
    `MR-P` = list(fun = "ossrmp_pre_timed", args = list(calibrated = TRUE, calibrated_p = FALSE)),
    `MR-Pos` = list(fun = "ossrmp_poisson_timed", args = list(RB = TRUE)))
  if (include_uncalibrated) methods$`MR-NoRB` <- list(fun = "ossrmp_timed", args = list(RB = FALSE))
  methods
}
.s4_seed <- function(seed, grid_id, rep_id, offset = 0L) {
  as.integer((as.double(seed) + grid_id * 1000003 + rep_id * 1009 + offset) %% 2147483646 + 1)
}
.s4_threads <- function() {
  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
             VECLIB_MAXIMUM_THREADS = "1", BLIS_NUM_THREADS = "1")
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
    return("RhpcBLASctl + environment")
  }
  "environment only; set BLAS threads before starting R when applicable"
}
.s4_clock_delta <- function(start) {
  delta <- proc.time() - start
  c(CPU = unname(delta["user.self"] + delta["sys.self"]), Elapsed = unname(delta["elapsed"]))
}
.s4_call <- function(env, method, X, y, n0, n_dag, cfg, seed) {
  spec <- .s4_methods(cfg$include_uncalibrated)[[method]]
  args <- list(X = X, y = y, n0 = n0, seed = seed, tol = cfg$tol, maxit = cfg$maxit,
               pilot_sampling = cfg$pilot_sampling, verbose = FALSE,
               add_pilot_variance = cfg$add_pilot_variance, add_full_variance = FALSE)
  if (method %in% c("O", "O-Pos")) args$n_sub <- n_dag else {
    args$n_dag <- n_dag
    args$rho <- cfg$rho_mr
  }
  do.call(env[[spec$fun]], c(args, spec$args))
}
.s4_empty_row <- function(job, method, cfg, n0, n_dag) {
  row <- data.frame(rep = job$rep, N = job$N, K = cfg$K, Method = method,
                    n0 = n0, n_dag = n_dag, PilotRatio = n0 / job$N,
                    HighResRatio = n_dag / job$N, BudgetRatio = (n0 + n_dag) / job$N,
                    N_hr = NA_integer_, n_sel = NA_integer_, q = NA_integer_, d_g = NA_integer_,
                    CPU = NA_real_, Elapsed = NA_real_, EstimationCPU = NA_real_,
                    EstimationElapsed = NA_real_, SamplingPreparationCPU = NA_real_,
                    SamplingPreparationElapsed = NA_real_, RBSetupCPU = NA_real_,
                    RBSetupElapsed = NA_real_, Success = FALSE, Converged = FALSE,
                    PilotConverged = FALSE, Iterations = NA_integer_, MSE = NA_real_,
                    FitSeed = .s4_seed(cfg$seed_base, job$grid, job$rep, 200L),
                    Order = NA_integer_, ErrorMsg = NA_character_, stringsAsFactors = FALSE)
  for (stage in c(.s4_primary, "Other", .s4_diagnostic))
    for (metric in c("CPU", "Elapsed")) row[[paste0(stage, metric)]] <- NA_real_
  for (stage in .s4_diagnostic) row[[paste0(stage, "Calls")]] <- NA_integer_
  row
}
.s4_record_fit <- function(row, fit, data, total) {
  if (!is.list(fit) || !is.matrix(fit$timing) ||
      !all(c("CPU", "Elapsed") %in% rownames(fit$timing)) ||
      !all(c(.s4_primary, .s4_diagnostic) %in% colnames(fit$timing)) ||
      any(!is.finite(fit$timing)) || any(fit$timing < 0) ||
      any(lengths(fit[c("iterations", "q", "d_g")]) != 1L) ||
      !all(.s4_diagnostic %in% names(fit$timing_calls))) stop("Invalid timed method result.")
  row$Converged <- isTRUE(fit$converged)
  row$PilotConverged <- isTRUE(fit$pilot_converged)
  row$Iterations <- fit$iterations
  row$q <- fit$q
  row$d_g <- fit$d_g
  row$n_sel <- if (!is.null(fit$n.sel)) fit$n.sel else length(fit$idx.sel)
  row$N_hr <- if (!is.null(fit$N.dag)) fit$N.dag else NA_integer_
  # Preserve the original sum of squared coefficient errors.
  if (!is.null(data$theta_true) && length(data$theta_true) == length(fit$beta))
    row$MSE <- sum((as.vector(fit$beta) - data$theta_true)^2)
  for (stage in c(.s4_primary, .s4_diagnostic))
    for (metric in c("CPU", "Elapsed")) row[[paste0(stage, metric)]] <- fit$timing[metric, stage]
  for (metric in c("CPU", "Elapsed")) {
    remainder <- total[metric] - sum(fit$timing[metric, .s4_primary])
    if (remainder < -1e-6) stop("Overlapping primary timers detected.")
    # Other includes bookkeeping and uninstrumented loop overhead.
    row[[paste0("Other", metric)]] <- max(0, remainder)
    row[[paste0("Estimation", metric)]] <- total[metric] -
      sum(fit$timing[metric, c("Inference", "Output")])
    row[[paste0("SamplingPreparation", metric)]] <- sum(fit$timing[metric,
      c("FullProb", "Partition", "SamplingScores", "Sampling")])
    row[[paste0("RBSetup", metric)]] <- sum(fit$timing[metric, c("AuxiliaryMean", "RBCalibration")])
  }
  for (stage in .s4_diagnostic) row[[paste0(stage, "Calls")]] <- fit$timing_calls[stage]
  row$Success <- TRUE
  row
}
.s4_task <- function(job) {
  cfg <- job$cfg
  methods <- names(.s4_methods(cfg$include_uncalibrated))
  n0 <- max(as.integer(round(cfg$pilot_ratio * job$N)), cfg$K + 1L)
  n_dag <- max(as.integer(round(cfg$highres_ratio * job$N)), 1L)
  rows <- setNames(lapply(methods, function(method)
    .s4_empty_row(job, method, cfg, n0, n_dag)), methods)
  set.seed(.s4_seed(cfg$seed_base, job$grid, job$rep))
  data <- tryCatch({
    if (is.null(cfg$data_generator))
      .S4CPU_WORKER$DataGenSM(n = job$N, n_class = cfg$K + 1L,
        x_dist = "gaussian", cov_type = "ar1", rho = 0.5,
        standardize_X = FALSE, seed = NULL, use_rfast = FALSE) else
          cfg$data_generator(job$N, cfg$K)
  }, error = function(e) e)
  if (inherits(data, "error")) {
    for (method in methods) rows[[method]]$ErrorMsg <- paste("Data generation:", conditionMessage(data))
  } else {
    # Warm-up fits and garbage collection are outside the measured fit time.
    if (cfg$warmup && !isTRUE(.S4CPU_WORKER$warmed)) {
      nw <- min(2000L, nrow(data$X))
      for (method in methods) invisible(tryCatch(.s4_call(.S4CPU_WORKER, method,
        data$X[seq_len(nw), , drop = FALSE], data$y[seq_len(nw)],
        min(200L, nw - 1L), min(200L, nw - 1L), cfg, 4321L), error = function(e) NULL))
      .S4CPU_WORKER$warmed <- TRUE
    }
    set.seed(.s4_seed(cfg$seed_base, job$grid, job$rep, 100L))
    method_order <- sample(methods)
    for (i in seq_along(method_order)) {
      method <- method_order[i]
      row <- rows[[method]]
      row$Order <- i
      invisible(gc())
      start <- proc.time()
      fit <- tryCatch(.s4_call(.S4CPU_WORKER, method, data$X, data$y,
                              n0, n_dag, cfg, row$FitSeed), error = function(e) e)
      total <- .s4_clock_delta(start)
      row$CPU <- total["CPU"]
      row$Elapsed <- total["Elapsed"]
      if (inherits(fit, "error")) row$ErrorMsg <- conditionMessage(fit) else {
        recorded <- tryCatch(.s4_record_fit(row, fit, data, total), error = function(e) e)
        if (inherits(recorded, "error")) row$ErrorMsg <- conditionMessage(recorded) else row <- recorded
      }
      rows[[method]] <- row
      rm(fit)
    }
    # MR and MR-NoRB share fit seeds, pilots, and sampling scores.
    # The companion methods retain their original pilot fitting controls.
  }
  raw <- do.call(rbind, rows)
  rownames(raw) <- NULL
  file <- file.path(cfg$out_dir, "replicates", sprintf("N%.0f_rep%04d.rds", job$N, job$rep))
  saveRDS(raw, file)
  raw
}

# Summaries and plots

.s4_stats <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  c(mean = if (n) mean(x) else NA_real_, sd = if (n > 1L) sd(x) else NA_real_,
    se = if (n > 1L) sd(x) / sqrt(n) else NA_real_, median = if (n) median(x) else NA_real_)
}
.s4_add_stats <- function(row, x, name, median = TRUE) {
  stats <- .s4_stats(x)
  row[[name]] <- stats["mean"]
  for (stat in c("sd", "se", if (median) "median")) row[[paste0(name, "_", stat)]] <- stats[stat]
  row
}
SummarizeS4CPU <- function(raw, require_converged = TRUE) {
  if (!nrow(raw)) return(list(summary = data.frame(), stages = data.frame(),
                             diagnostics = data.frame(), paired_RB = data.frame()))
  groups <- split(seq_len(nrow(raw)), interaction(raw$N, raw$Method, drop = TRUE))
  totals <- stages <- diagnostics <- list()
  for (ix in groups) {
    all <- raw[ix, , drop = FALSE]
    used <- all[all$Success & (!require_converged | all$Converged), , drop = FALSE]
    base <- data.frame(N = all$N[1L], Method = all$Method[1L], n_rep = nrow(all),
                       n_success = sum(all$Success), n_converged = sum(all$Success & all$Converged),
                       n_pilot_converged = sum(all$Success & all$PilotConverged), n_used = nrow(used),
                       FailureRate = mean(!all$Success), NonConvergenceRate = mean(all$Success & !all$Converged))
    row <- base
    for (v in c("CPU", "Elapsed", "EstimationCPU", "EstimationElapsed", "SamplingPreparationCPU",
                "SamplingPreparationElapsed", "RBSetupCPU", "RBSetupElapsed", "Iterations", "N_hr", "n_sel", "MSE"))
      row <- .s4_add_stats(row, used[[v]], v)
    totals[[length(totals) + 1L]] <- row
    for (stage in c(.s4_primary, "Other")) {
      row <- base
      row$Stage <- stage
      for (metric in c("CPU", "Elapsed")) {
        values <- used[[paste0(stage, metric)]]
        row <- .s4_add_stats(row, values, metric, median = FALSE)
        positive <- is.finite(used[[metric]]) & used[[metric]] > 0
        row[[paste0(metric, "Share")]] <- if (any(positive))
          mean(values[positive] / used[[metric]][positive]) else NA_real_
      }
      stages[[length(stages) + 1L]] <- row
    }
    for (stage in .s4_diagnostic) {
      row <- base
      row$Diagnostic <- stage
      calls <- used[[paste0(stage, "Calls")]]
      row$MeanCalls <- if (nrow(used)) mean(calls) else NA_real_
      for (metric in c("CPU", "Elapsed")) {
        values <- used[[paste0(stage, metric)]]
        row[[metric]] <- .s4_stats(values)["mean"]
        present <- is.finite(calls) & calls > 0L
        row[[paste0(metric, "PerCall")]] <- if (any(present)) mean(values[present] / calls[present]) else NA_real_
      }
      diagnostics[[length(diagnostics) + 1L]] <- row
    }
  }
  summary <- do.call(rbind, totals)
  summary <- summary[order(summary$N, summary$Method), , drop = FALSE]
  paired <- list()
  for (N in sort(unique(raw$N))) {
    rb <- raw[raw$N == N & raw$Method == "MR", , drop = FALSE]
    no_rb <- raw[raw$N == N & raw$Method == "MR-NoRB", , drop = FALSE]
    if (!nrow(rb) || !nrow(no_rb)) next
    pairs <- merge(rb, no_rb, by = "rep", suffixes = c(".RB", ".NoRB"))
    ok <- pairs$Success.RB & pairs$Success.NoRB
    if (require_converged) ok <- ok & pairs$Converged.RB & pairs$Converged.NoRB
    pairs <- pairs[ok, , drop = FALSE]
    row <- data.frame(N = N, n_pairs = nrow(pairs))
    # Paired differences include changes in iteration and inference time.
    for (v in c("CPU", "Elapsed", "EstimationCPU", "EstimationElapsed", "MSE")) {
      stats <- .s4_stats(pairs[[paste0(v, ".RB")]] - pairs[[paste0(v, ".NoRB")]])
      row[[paste0(v, "Difference")]] <- stats["mean"]
      row[[paste0(v, "DifferenceSE")]] <- stats["se"]
    }
    paired[[length(paired) + 1L]] <- row
  }
  list(summary = summary, stages = do.call(rbind, stages), diagnostics = do.call(rbind, diagnostics),
       paired_RB = if (length(paired)) do.call(rbind, paired) else data.frame())
}
.s4_colors <- c("#D55E00", "#999900", "#009E73", "#0072B2", "#CC79A7", "#666666")
PlotS4CPU <- function(res_s4, metric = c("CPU", "Elapsed"), file = NULL, show_se = TRUE) {
  metric <- match.arg(metric)
  data <- if (is.list(res_s4) && !is.data.frame(res_s4)) res_s4$summary else res_s4
  if (!is.null(file)) {
    pdf(file, width = 7, height = 4.8)
    on.exit(dev.off())
  }
  methods <- unique(data$Method)
  upper <- data[[metric]]
  if (show_se) upper <- upper + ifelse(is.finite(data[[paste0(metric, "_se")]]), data[[paste0(metric, "_se")]], 0)
  if (!any(is.finite(upper))) {
    plot.new()
    title("No eligible fits")
    return(invisible(NULL))
  }
  plot(range(data$N), c(0, max(upper, na.rm = TRUE) * 1.1 + 1e-3), type = "n",
       xlab = "Full-data sample size N", ylab = paste(metric, "time (seconds)"))
  for (i in seq_along(methods)) {
    points <- data[data$Method == methods[i], , drop = FALSE]
    points <- points[order(points$N), , drop = FALSE]
    lines(points$N, points[[metric]], type = "b", pch = i, col = .s4_colors[i])
    if (show_se) {
      se <- points[[paste0(metric, "_se")]]
      ok <- is.finite(se) & se > 1e-8
      if (any(ok)) arrows(points$N[ok], pmax(0, points[[metric]][ok] - se[ok]),
        points$N[ok], points[[metric]][ok] + se[ok], code = 3, angle = 90, length = 0.035, col = .s4_colors[i])
    }
  }
  legend("topleft", methods, col = .s4_colors[seq_along(methods)], pch = seq_along(methods), lty = 1, bty = "n", cex = 0.85)
  invisible(data)
}
.s4_stage_plot <- function(result, file, metric = "CPU") {
  pdf(file, width = 10, height = 7)
  on.exit(dev.off())
  old <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  on.exit(par(old), add = TRUE, after = FALSE)
  data <- result$stages
  stages <- c("Partition", "AuxiliaryMean", "RBCalibration", "PredictionDesign", "Iteration", "Inference")
  methods <- unique(data$Method)
  for (stage in stages) {
    panel <- data[data$Stage == stage, , drop = FALSE]
    finite <- panel[[metric]][is.finite(panel[[metric]])]
    ymax <- if (length(finite)) max(finite) * 1.12 + 1e-3 else 1
    plot(range(data$N), c(0, ymax), type = "n", xlab = "N", ylab = paste(metric, "seconds"), main = stage)
    for (i in seq_along(methods)) {
      points <- panel[panel$Method == methods[i], , drop = FALSE]
      points <- points[order(points$N), , drop = FALSE]
      lines(points$N, points[[metric]], type = "b", col = .s4_colors[i], pch = i)
    }
    if (stage == stages[1L]) legend("topleft", methods, col = .s4_colors[seq_along(methods)],
                                   lty = 1, pch = seq_along(methods), cex = 0.65, bty = "n")
  }
  invisible(NULL)
}

# Main simulation

RunS4CPU_ParamGrid <- function(
    N_grid = c(2e5, 4e5, 6e5, 8e5, 1e6), K = 5L, m = 200L,
    pilot_ratio = 0.0005, highres_ratio = 0.0015, rho_mr = 0.999,
    tol = 1e-8, maxit = 200L, cores = 1L, seed_base = 20261012L,
    project_dir = getwd(), out_dir = file.path(project_dir, "S4CPU_stage_runs"),
    include_uncalibrated = TRUE, pilot_sampling = "proportional",
    add_pilot_variance = TRUE, require_converged = TRUE, warmup = TRUE,
    make_plots = TRUE, data_generator = NULL, method_env = NULL) {
  is_count <- function(x) length(x) == 1L && is.finite(x) && x >= 1 && x == floor(x)
  is_ratio <- function(x) length(x) == 1L && is.finite(x) && x > 0 && x < 1
  stopifnot(length(N_grid) > 0L, all(is.finite(N_grid)), all(N_grid > 1),
            !anyDuplicated(N_grid), all(N_grid == floor(N_grid)),
            is_count(K), is_count(m), is_count(cores), is_count(maxit),
            length(seed_base) == 1L, is.finite(seed_base), seed_base >= 0,
            is_ratio(pilot_ratio), is_ratio(highres_ratio), is_ratio(rho_mr),
            length(tol) == 1L, is.finite(tol), tol > 0,
            length(pilot_sampling) == 1L, pilot_sampling %in% c("uniform", "proportional"))
  n0 <- pmax(round(pilot_ratio * N_grid), K + 1L)
  n_dag <- pmax(round(highres_ratio * N_grid), 1L)
  if (any(n0 >= N_grid)) stop("Pilot size must be smaller than N.")
  if (any(n_dag > N_grid - n0)) stop("Requested Poisson budget exceeds remaining observations.")
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, mustWork = TRUE)
  dir.create(file.path(out_dir, "replicates"), showWarnings = FALSE)
  thread_control <- .s4_threads()
  cfg <- list(K = as.integer(K), m = as.integer(m), pilot_ratio = pilot_ratio,
              highres_ratio = highres_ratio, rho_mr = rho_mr, tol = tol, maxit = as.integer(maxit),
              seed_base = seed_base, include_uncalibrated = include_uncalibrated,
              pilot_sampling = pilot_sampling, add_pilot_variance = add_pilot_variance,
              require_converged = require_converged, warmup = warmup,
              project_dir = project_dir, out_dir = out_dir, data_generator = data_generator)
  .S4CPU_WORKER <<- LoadS4CPU(project_dir, method_env)
  workers <- min(as.integer(cores), as.integer(m))
  cl <- NULL
  if (workers > 1L) {
    cl <- makePSOCKcluster(workers)
    on.exit(stopCluster(cl), add = TRUE)
    # Workers load methods, but never source this driver or its run block.
    export_env <- environment(RunS4CPU_ParamGrid)
    exports <- c("LoadS4CPU", ls(envir = export_env, all.names = TRUE, pattern = "^\\.s4_"))
    clusterExport(cl, exports, envir = export_env)
    clusterCall(cl, function(project_dir) {
      .s4_threads()
      assign(".S4CPU_WORKER", LoadS4CPU(project_dir), envir = .GlobalEnv)
      TRUE
    }, project_dir)
  }
  files <- file.path(project_dir, c("Syn.R", "OSS.R", "ossrmp.R", "ossrmp_pre.R", "ossrmp_pos.R"))
  metadata <- list(config = cfg, N_grid = N_grid, cores = cores, workers = workers,
                   metric = "CPU=user.self+sys.self", thread_control = thread_control,
                   source_md5 = tools::md5sum(files), timestamp = Sys.time(), R = R.version.string,
                   session = capture.output(sessionInfo()),
                   note = "Diagnostic timers overlap primary stages. MR-NoRB pairs share fit seeds.")
  saveRDS(metadata, file.path(out_dir, "metadata.rds"))
  raw <- data.frame()
  for (g in seq_along(N_grid)) {
    cat(sprintf("[S4CPU] N=%.0f | %d repetitions | %d worker(s)\n", N_grid[g], m, workers))
    jobs <- lapply(seq_len(m), function(i) list(N = N_grid[g], rep = i, grid = g, cfg = cfg))
    results <- if (is.null(cl)) lapply(jobs, .s4_task) else parLapplyLB(cl, jobs, .s4_task)
    raw <- rbind(raw, do.call(rbind, results))
    saveRDS(raw, file.path(out_dir, "S4CPU_raw.rds"))
    summaries <- SummarizeS4CPU(raw, require_converged)
    cat(sprintf("[S4CPU] completed N=%.0f: %d/%d converged method fits\n", N_grid[g],
                sum(raw$N == N_grid[g] & raw$Converged), sum(raw$N == N_grid[g])))
    for (name in names(summaries)) write.csv(summaries[[name]], file.path(out_dir, paste0(name, ".csv")),
                                            row.names = FALSE, fileEncoding = "UTF-8")
    saveRDS(summaries$summary, file.path(out_dir, "S4CPU_final.rds"))
  }
  result <- c(list(raw = raw, metadata = metadata), summaries)
  saveRDS(result, file.path(out_dir, "S4CPU_results.rds"))
  if (make_plots) {
    PlotS4CPU(result, "CPU", file.path(out_dir, "S4CPU_total_CPU.pdf"))
    PlotS4CPU(result, "Elapsed", file.path(out_dir, "S4CPU_total_elapsed.pdf"))
    .s4_stage_plot(result, file.path(out_dir, "S4CPU_stages_CPU.pdf"))
  }
  cat("[S4CPU] Results: ", out_dir, "\n", sep = "")
  invisible(result)
}

# Simulation settings and calls
# Set options(S4CPU.autorun = FALSE) before sourcing to skip the formal run.

N_grid <- c(2e5, 4e5, 6e5, 8e5, 1e6)
K <- 5L
m <- 200L
cores <- 30L
pilot_ratio <- 0.0005
highres_ratio <- 0.0015
rho_mr <- 0.999
seed_base <- 20261012L
out_dir <- file.path(project_dir, "S4CPU_stage_runs")

if (isTRUE(getOption("S4CPU.autorun", TRUE))) {
  res_s4 <- RunS4CPU_ParamGrid(
    N_grid = N_grid, K = K, m = m, cores = cores,
    pilot_ratio = pilot_ratio, highres_ratio = highres_ratio,
    rho_mr = rho_mr, tol = 1e-8, maxit = 200L, seed_base = seed_base,
    project_dir = project_dir, out_dir = out_dir,
    include_uncalibrated = TRUE, method_env = method_env)
  PlotS4CPU(res_s4, metric = "CPU")
}
