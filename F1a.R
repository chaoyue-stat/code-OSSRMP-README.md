# ============================================================
# Visualize L-optimal sampling probabilities for softmax regression.
# Simulate three-class data, fit a uniform pilot estimator, and color
# observations by their log sampling probabilities within each class.
# Solid lines show pilot-based decision boundaries. Optional dashed
# lines show logit-margin contours at C_r0 = log(rho / (1 - rho));
# optional crosses mark the centroids of the margin-based partitions.
# ============================================================

rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(nnet)
  library(ggplot2)
  library(ggnewscale)
  library(scales)
})

# Parameters
seed <- 20261013L
N <- 50000L; n0 <- 3000L; n_show <- 50000L
rho <- 0.999; grid_n <- 250L; amp <- 2.5
show_margin <- FALSE; show_centroids <- FALSE
cls <- c("0", "1", "2")
C_r0 <- qlogis(rho)
set.seed(seed)

softmax <- function(M) {
  E <- exp(sweep(M, 1, apply(M, 1, max), "-"))
  E / rowSums(E)
}

logit_margin <- function(E, k) {
  other <- setdiff(cls, k)
  E[, k] - pmax(E[, other[1]], E[, other[2]])
}

# Simulate three-class softmax data
X <- matrix(rnorm(N * 2), N, 2, dimnames = list(NULL, c("x1", "x2")))
B <- amp * cbind(c(1.8, -0.9), c(-1, 1.7), c(-0.8, -0.8))
P <- softmax(X %*% B)
Y <- factor(apply(P, 1, function(pr) sample(0:2, 1, prob = pr)), levels = 0:2)
dat <- data.frame(y = Y, x1 = X[, 1], x2 = X[, 2])

# Uniform pilot subsample and pilot estimator
idx0 <- sample.int(N, n0)
if (any(table(Y[idx0]) == 0L)) stop("The pilot sample must contain all classes.")
fit <- multinom(y ~ x1 + x2, dat[idx0, ], trace = FALSE, maxit = 1000)
if (fit$convergence != 0L) stop("The pilot fit did not converge.")
Phat <- as.matrix(predict(fit, dat, type = "probs"))[, cls, drop = FALSE]

# L-optimal sampling probabilities; class 0 is the reference class
Ycirc <- cbind(as.numeric(Y == "1"), as.numeric(Y == "2"))
score <- sqrt(rowSums(cbind(1, X)^2)) *
  sqrt(rowSums((Ycirc - Phat[, c("1", "2")])^2))
pi_L <- score / sum(score)
df <- data.frame(
  x1 = X[, 1], x2 = X[, 2], class = Y,
  pi = pi_L, logpi = log(pmax(pi_L, .Machine$double.xmin))
)

# Margin-based partitions: Omega_k = {Y = k, D_k(X) >= C_r0}
Bhat <- cbind("0" = c(0, 0, 0), t(coef(fit)[c("1", "2"), ]))
Eta_pt <- cbind(1, X) %*% Bhat
margins <- vapply(cls, function(k) logit_margin(Eta_pt, k), numeric(N))
Delta_true <- margins[cbind(seq_len(N), as.integer(Y))]
region <- ifelse(Delta_true >= C_r0, paste0("Omega", Y), "Omega_dag")
Omega <- setNames(
  lapply(cls, function(k) which(region == paste0("Omega", k))), cls
)
OmegaD <- which(region == "Omega_dag")

centroids <- do.call(rbind, lapply(cls, function(k) {
  ii <- Omega[[k]]
  data.frame(
    class = k,
    x1 = if (length(ii)) mean(X[ii, 1]) else NA_real_,
    x2 = if (length(ii)) mean(X[ii, 2]) else NA_real_
  )
}))
centroids <- centroids[complete.cases(centroids), ]

cat("n0 =", n0, "; rho =", rho, "; C_r0 =", round(C_r0, 4), "\n")
print(table(factor(region, levels = c(paste0("Omega", cls), "Omega_dag"))))

# Plotting sample and prediction grid
set.seed(seed + 1L)
df_plot <- df[sample.int(N, min(n_show, N)), ]
xlim <- range(df$x1); ylim <- range(df$x2)
g <- expand.grid(
  x1 = seq(xlim[1], xlim[2], length.out = grid_n),
  x2 = seq(ylim[1], ylim[2], length.out = grid_n)
)
eta <- cbind(1, g$x1, g$x2) %*% Bhat
pred <- cls[max.col(eta, ties.method = "first")]

# Add a black contour with a white halo
add_boundary <- function(plot, z, level, lty = "solid", lw = 1.15) {
  zz <- z[is.finite(z)]
  if (!length(zz) || level <= min(zz) || level >= max(zz)) return(plot)
  d <- data.frame(g, z = z)
  plot +
    geom_contour(
      data = d, aes(x1, x2, z = z), breaks = level,
      colour = "white", linewidth = 2.6,
      na.rm = TRUE, inherit.aes = FALSE
    ) +
    geom_contour(
      data = d, aes(x1, x2, z = z), breaks = level,
      colour = "black", linewidth = lw, linetype = lty,
      na.rm = TRUE, inherit.aes = FALSE
    )
}

# Class-specific gradients with a common numerical scale
lo <- c("#FEE5D9", "#DEEBF7", "#E5F5E0")
hi <- c("#B10026", "#08519C", "#006D2C")
lty <- c("longdash", "twodash", "dashed")
rng <- range(df$logpi)
brks <- pretty(rng, n = 6)
brks <- brks[brks >= rng[1] & brks <= rng[2]]

p3 <- ggplot()
for (j in seq_along(cls)) {
  k <- cls[j]
  if (j > 1L) p3 <- p3 + new_scale_color()
  p3 <- p3 +
    geom_point(
      data = df_plot[df_plot$class == k, ],
      aes(x1, x2, colour = logpi), size = 0.9, alpha = 0.9
    ) +
    scale_color_gradient(
      low = lo[j], high = hi[j], limits = rng, breaks = brks,
      oob = squish, name = bquote("Class " * .(k) * ":  " * log(pi[i]^L)),
      guide = guide_colorbar(
        order = j, title.position = "top", title.hjust = 0.5,
        barwidth = grid::unit(0.5, "cm"),
        barheight = grid::unit(2.2, "cm"),
        ticks.colour = "grey40", frame.colour = "grey70"
      )
    )
}

# Decision boundaries: retain ties between the two largest logits
for (ab in combn(cls, 2, simplify = FALSE)) {
  other <- setdiff(cls, ab)
  active <- pmax(eta[, ab[1]], eta[, ab[2]]) >= eta[, other] - 1e-8
  z <- ifelse(active, eta[, ab[1]] - eta[, ab[2]], NA_real_)
  p3 <- add_boundary(p3, z, 0)
}

# Optional piecewise-linear logit-margin boundaries
if (show_margin) {
  for (j in seq_along(cls)) {
    z <- ifelse(pred == cls[j], logit_margin(eta, cls[j]), NA_real_)
    p3 <- add_boundary(p3, z, C_r0, lty[j], 1.25)
  }
}

# Optional partition centroids
if (show_centroids) {
  p3 <- p3 + geom_point(
    data = centroids, aes(x1, x2), inherit.aes = FALSE,
    shape = 3, size = 5, stroke = 1.5, colour = "black"
  )
}

p3 <- p3 +
  coord_equal(xlim = xlim, ylim = ylim) +
  theme_bw() +
  theme(
    legend.box = "vertical",
    legend.spacing.y = grid::unit(0.15, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  labs(x = expression(X[1]), y = expression(X[2]))

print(p3)