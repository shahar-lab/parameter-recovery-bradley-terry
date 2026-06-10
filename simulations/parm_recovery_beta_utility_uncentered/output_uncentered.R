rm(list = ls())

#### SETUP ####
library(tidyverse)
library(posterior)
library(ggdist)
library(patchwork)

data_dir <- "./simulations/parm_recovery_beta_utility_uncentered/data/"
figs_dir <- "./simulations/parm_recovery_beta_utility_uncentered/figs/"

load(paste0(data_dir, "cfg.rdata"))
load(paste0(data_dir, "true_parameters.rdata"))
fit <- readRDS(paste0(data_dir, "fit.rds"))

#### EXTRACT POSTERIOR SAMPLES ####

draws <- as_draws_df(fit)

mu_log_beta_samples    <- draws$mu_log_beta
sigma_log_beta_samples <- draws$sigma_log_beta

beta_cols <- grep("^beta\\[", names(draws), value = TRUE)
u_cols    <- grep("^u_matrix\\[", names(draws), value = TRUE)

median_beta <- sapply(beta_cols, function(col) median(draws[[col]]))

u_medians    <- sapply(u_cols, function(col) median(draws[[col]]))
u_matrix_est <- matrix(
  u_medians,
  nrow  = cfg$Nsubjects,
  ncol  = cfg$Noffer,
  byrow = FALSE
)

estimates <- list(
  median_beta  = median_beta,
  u_matrix_est = u_matrix_est
)

#### FIGURE 1: POSTERIOR HYPERPARAMETERS ####

true_mu    <- true_parameters$beta_mu
true_sigma <- true_parameters$beta_sigma

med_mu    <- median(mu_log_beta_samples)
med_sigma <- median(sigma_log_beta_samples)

p_mu <- ggplot(data.frame(theta = mu_log_beta_samples), aes(x = theta, y = 0)) +
  stat_slab(fill = "gray80") +
  stat_pointinterval(
    .width     = c(0.80, 0.90),
    point_size = 3,
    linewidth  = c(2, 1)
  ) +
  geom_vline(xintercept = true_mu, linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  geom_vline(xintercept = med_mu,  linetype = "dashed", colour = "grey65", linewidth = 0.4) +
  annotate(
    "text", x = med_mu, y = Inf,
    label = sprintf("[median = %.2f]", med_mu),
    hjust = -0.05, vjust = 1.4, size = 3.2, colour = "grey40"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid   = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    axis.line.x  = element_line(colour = "grey30")
  ) +
  labs(x = expression(mu[log~beta])) +
  coord_cartesian(ylim = c(0, 1.3), clip = "off")

p_sigma <- ggplot(data.frame(theta = sigma_log_beta_samples), aes(x = theta, y = 0)) +
  stat_slab(fill = "gray80") +
  stat_pointinterval(
    .width     = c(0.80, 0.90),
    point_size = 3,
    linewidth  = c(2, 1)
  ) +
  geom_vline(xintercept = true_sigma, linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  geom_vline(xintercept = med_sigma,  linetype = "dashed", colour = "grey65", linewidth = 0.4) +
  annotate(
    "text", x = med_sigma, y = Inf,
    label = sprintf("[median = %.2f]", med_sigma),
    hjust = -0.05, vjust = 1.4, size = 3.2, colour = "grey40"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid   = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    axis.line.x  = element_line(colour = "grey30")
  ) +
  labs(x = expression(sigma[log~beta])) +
  coord_cartesian(ylim = c(0, 1.3), clip = "off")

fig1 <- p_mu + p_sigma + plot_annotation(tag_levels = "A")

ggsave(
  paste0(figs_dir, "posterior_hyperparams.pdf"),
  fig1,
  width = 10, height = 3
)

#### FIGURE 2: U_MATRIX PARAMETER RECOVERY ####

df_u <- data.frame(
  true      = as.vector(true_parameters$u_matrix),
  estimated = as.vector(u_matrix_est),
  offer     = rep(1:cfg$Noffer, each = cfg$Nsubjects)
) |>
  mutate(offer = factor(offer, labels = paste0("Offer ", 1:cfg$Noffer)))

shared_limits <- range(c(df_u$true, df_u$estimated), na.rm = TRUE)
axis_breaks   <- seq(shared_limits[1], shared_limits[2], length.out = 4)

pearson_labels <- df_u |>
  group_by(offer) |>
  summarise(r = cor(true, estimated), .groups = "drop") |>
  mutate(label = sprintf("[Pearson r = %.2f]", r))

fig2 <- ggplot(df_u, aes(x = true, y = estimated)) +
  geom_point(colour = "#4477AA", alpha = 0.75, size = 2) +
  geom_smooth(method = "lm", se = FALSE, colour = "#EE6677", linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_text(
    data = pearson_labels,
    aes(label = label),
    x = Inf, y = Inf,
    hjust = 1.05, vjust = 1.4,
    size = 3.5, colour = "grey30",
    inherit.aes = FALSE
  ) +
  scale_x_continuous(breaks = round(axis_breaks, 1)) +
  scale_y_continuous(breaks = round(axis_breaks, 1)) +
  coord_equal(xlim = shared_limits, ylim = shared_limits, clip = "off") +
  facet_wrap(~offer, nrow = 2) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "True utility", y = "Estimated utility")

ggsave(
  paste0(figs_dir, "u_matrix_recovery.pdf"),
  fig2,
  width = 12, height = 6
)

#### FIGURE 3: BETA PARAMETER RECOVERY ####

df_beta <- data.frame(
  true      = true_parameters$beta,
  estimated = unname(median_beta)
)
df_beta <- df_beta[complete.cases(df_beta$true, df_beta$estimated), ]
stopifnot(nrow(df_beta) == length(true_parameters$beta))

shared_limits_beta <- range(c(df_beta$true, df_beta$estimated), na.rm = TRUE)
axis_breaks_beta   <- seq(shared_limits_beta[1], shared_limits_beta[2], length.out = 4)
pearson_r_beta     <- cor(df_beta$true, df_beta$estimated, method = "pearson")

fig3 <- ggplot(df_beta, aes(x = true, y = estimated)) +
  geom_point(colour = "#4477AA", alpha = 0.75, size = 2) +
  geom_smooth(method = "lm", se = FALSE, colour = "#EE6677", linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  annotate(
    "text",
    x     = Inf,
    y     = Inf,
    label = sprintf("[Pearson r = %.2f]", pearson_r_beta),
    hjust = 1.05, vjust = 1.4,
    size  = 3.5, colour = "grey30"
  ) +
  scale_x_continuous(breaks = round(axis_breaks_beta, 2)) +
  scale_y_continuous(breaks = round(axis_breaks_beta, 2)) +
  coord_equal(xlim = shared_limits_beta, ylim = shared_limits_beta, clip = "off") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = expression("True " * beta), y = expression("Estimated " * beta * " (median)"))

ggsave(
  paste0(figs_dir, "beta_recovery.pdf"),
  fig3,
  width = 5, height = 5
)
