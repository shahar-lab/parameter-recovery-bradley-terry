rm(list = ls())

#### SETUP ####
library(tidyverse)
library(posterior)
library(ggdist)
library(patchwork)

data_dir <- "./simulations/parm_recovery_beta_utility_compared_5_for_real_data/data/"
figs_dir <- "./simulations/parm_recovery_beta_utility_compared_5_for_real_data/figs/"

load(paste0(data_dir, "cfg.rdata"))
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


med_mu    <- median(mu_log_beta_samples)
med_sigma <- median(sigma_log_beta_samples)

p_mu <- ggplot(data.frame(theta = mu_log_beta_samples), aes(x = theta, y = 0)) +
  stat_slab(fill = "gray80") +
  stat_pointinterval(
    .width     = c(0.80, 0.90),
    point_size = 3,
    linewidth  = c(2, 1)
  ) +
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
  xlim(0, 1) +
  labs(x = expression(mu[log~beta])) +
  coord_cartesian(ylim = c(0, 1.3), clip = "off")

p_sigma <- ggplot(data.frame(theta = sigma_log_beta_samples), aes(x = theta, y = 0)) +
  stat_slab(fill = "gray80") +
  stat_pointinterval(
    .width     = c(0.80, 0.90),
    point_size = 3,
    linewidth  = c(2, 1)
  ) +
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
  xlim(0, 1) +
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
  estimated = as.vector(u_matrix_est),
  subject   = factor(rep(1:cfg$Nsubjects, times = cfg$Noffer)),
  offer     = rep(1:cfg$Noffer, each = cfg$Nsubjects)
) |>
  mutate(offer_label = factor(offer, labels = c("2", "3", "4", "5", "6", "7", "8")))

fig2 <- ggplot(df_u, aes(x = offer_label, y = estimated)) +
  geom_boxplot(outlier.shape = NA, fill = "#4477AA", alpha = 0.3, width = 0.5) +
  geom_jitter(aes(color = subject), width = 0.15, alpha = 0.7, size = 2) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none" 
  ) +
  labs(
    x = "Offer", 
    y = "Estimated Utility",
    title = "Estimated Utility by Offer"
  )

ggsave(
  paste0(figs_dir, "u_matrix_estimated.pdf"),
  fig2,
  width = 8, height = 6
)

#### FIGURE 3: BETA PARAMETER RECOVERY ####

df_beta <- data.frame(
  estimated = unname(median_beta)
)

fig3 <- ggplot(df_beta, aes(x = estimated)) +
  geom_histogram(fill = "#4477AA", color = "white", bins = 15, alpha = 0.8) +
  geom_vline(aes(xintercept = median(estimated)), color = "#EE6677", linetype = "dashed", linewidth = 1) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    x = expression("Estimated " * beta * " (median)"), 
    y = "Count (Number of Subjects)",
    title = expression("Distribution of Estimated " * beta)
  )

ggsave(
  paste0(figs_dir, "beta_estimated.pdf"),
  fig3,
  width = 6, height = 5
)
