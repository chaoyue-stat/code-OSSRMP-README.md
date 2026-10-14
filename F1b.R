# ============================================================
# Visualize L-optimal sampling probabilities for softmax regression.
#
# This script simulates three-class data, fits a multinomial logistic
# model on a uniform pilot subsample, and computes L-optimal sampling
# probabilities from the norms of the individual pilot score vectors.
# Points are colored by their log sampling probabilities within each
# observed class. Solid lines show pilot-based decision boundaries,
# dashed lines show predicted-probability contours at level rho, and
# crosses mark the centroids of high-confidence sample partitions.
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
cls <- c("0", "1", "2")
set.seed(seed)

softmax <- function(M) {
  E <- exp(sweep(M, 1, apply(M, 1, max), "-"))
  E / rowSums(E)
}

# Simulate three-class softmax data
X <- matrix(rnorm(N * 2), N, 2, dimnames = list(NULL, c("x1", "x2")))
B <- amp * cbind(c(1.8, -0.9), c(-1, 1.7), c(-0.8, -0.8))
P <- softmax(X %*% B)
u <- runif(N)
Y <- factor(
  as.integer(u > P[, 1]) + as.integer(u > P[, 1] + P[, 2]),
  levels = 0:2
)
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

# High-confidence partitions and their sample centroids
p_y <- Phat[cbind(seq_len(N), as.integer(Y))]
region <- ifelse(p_y >= rho, paste0("Omega", Y), "Omega_dag")
Omega <- lapply(cls, function(k) which(region == paste0("Omega", k)))
names(Omega) <- cls

centroids <- do.call(rbind, lapply(cls, function(k) {
  ii <- Omega[[k]]
  data.frame(
    class = k,
    x1 = if (length(ii)) mean(X[ii, 1]) else NA_real_,
    x2 = if (length(ii)) mean(X[ii, 2]) else NA_real_
  )
}))
centroids <- centroids[complete.cases(centroids), ]

cat("n0 =", n0, "; rho =", rho, "\n")
print(table(factor(region, levels = c(paste0("Omega", cls), "Omega_dag"))))

# Plotting sample and prediction grid
set.seed(seed + 1L)
df_plot <- df[sample.int(N, min(n_show, N)), ]
xlim <- range(df$x1); ylim <- range(df$x2)
g <- expand.grid(
  x1 = seq(xlim[1], xlim[2], length.out = grid_n),
  x2 = seq(ylim[1], ylim[2], length.out = grid_n)
)

Bhat <- cbind("0" = c(0, 0, 0), t(coef(fit)[c("1", "2"), ]))
eta <- cbind(1, g$x1, g$x2) %*% Bhat
Pg <- softmax(eta)
pred <- cls[max.col(Pg, ties.method = "first")]

# Add a black contour with a white halo
add_boundary <- function(plot, z, level, lty = "solid", lw = 1.2) {
  zz <- z[is.finite(z)]
  if (!length(zz) || level <= min(zz) || level >= max(zz)) return(plot)
  d <- data.frame(g, z = z)
  plot +
    geom_contour(
      data = d, aes(x1, x2, z = z), breaks = level,
      colour = "white", linewidth = 2.7,
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
      aes(x1, x2, colour = logpi), size = 0.9, alpha = 0.95
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

# High-confidence boundaries within the predicted class regions
for (j in seq_along(cls)) {
  z <- ifelse(pred == cls[j], Pg[, j], NA_real_)
  p3 <- add_boundary(p3, z, rho, lty[j], 1.3)
}

# Add centroids and finalize the plot
p3 <- p3 +
  geom_point(
    data = centroids, aes(x1, x2), inherit.aes = FALSE,
    shape = 3, size = 4.2, stroke = 1.4, colour = "black"
  ) +
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