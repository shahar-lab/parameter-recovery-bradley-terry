rm(list = ls())

#### SETUP ####
library(tidyverse)
library(posterior)
library(ggdist)
library(patchwork)

data_dir <- "./simulations/parm_recovery_beta_utility_compared_5_for_real_data/data_local/"
figs_dir <- "./simulations/parm_recovery_beta_utility_compared_5_for_real_data/figs_local/"
data_path = './data/data_filtered/'
raw_data <- load(paste0(data_dir, "df.rdata"))

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


### ggplot for individual u values

samples_rvar <- as_draws_rvars(fit)

library(tidybayes)

draws_for_plot <- samples_rvar |> 
  spread_draws(u_matrix[subject, offer]) |> 
  mutate(offer_label = as.factor(offer + 1))
#head(grid_for_plot)

# fig4 <- ggplot(grid_for_plot, aes(x = offer_label, y = u_matrix)) +
#   geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0.2, color = "darkgray") +
#   geom_point(color = "navy", size = 2) + 
#   facet_wrap( ~ subject) +
#   theme_minimal()


fig5 <- ggplot(draws_for_plot, aes(x = offer_label, y = u_matrix)) +
  stat_halfeye(fill = "lightblue", alpha = 0.7, .width = 0.95, linewidth = 0.8) + 
  facet_wrap( ~ subject) +
  theme_minimal() +
  labs(y = "Utility (u)", x = "Valence")

ggsave(
  paste0(figs_dir, "individual_u_values.pdf"),
  fig5,
  width = 12, height = 12
)

### RT as a Function of Utility Difference

u_medians <- samples_rvar |> 
  spread_draws(u_matrix[subject, offer]) |> 
  median_qi(u_matrix) |>   
  mutate(offer_label = as.factor(offer + 1)) |> 
  select(subject, offer, u_median = u_matrix, offer_label)

raw_data_with_delta <- raw_data |> 
  left_join(u_medians, by = c("subject" = "subject", "offer_A" = "offer")) |> 
  rename(u_first = u_median) |> 
  
  # מיזוג עבור ההצעה השנייה בטרייל
  left_join(u_medians, by = c("subject" = "subject", "offer2" = "offer")) |> 
  rename(u_second = u_median) |> 
  
  # חישוב ההפרש בערך מוחלט
  mutate(delta_u = abs(u_first - u_second))

# 3. בניית הגרף: זמן תגובה כפונקציה של דלתא u
ggplot(raw_data_with_delta, aes(x = delta_u, y = subject_rt)) + # החלף ל-rt או לשם עמודת ה-RT שלך
  geom_point(alpha = 0.2, color = "black") + # נקודות חצי שקופות בגלל עומס הנתונים
  geom_smooth(method = "lm", color = "blue", se = TRUE) + # קו מגמה לינארי לכל נבדק
  facet_wrap(~ subject) + 
  theme_minimal() +
  labs(
    x = "Absolute Utility Difference (|Δu|)",
    y = "Reaction Time (RT)",
    title = "Reaction Time as a Function of Utility Difference per Subject"
  )

