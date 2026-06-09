rm(list = ls())

#### SETUP ####
library(cmdstanr)
library(posterior)

model_name <- "beta_utility"
model_dir  <- paste0("./models/", model_name, "/")
data_dir   <- "./simulations/parm_recovery_beta_utility/data/"
figs_dir   <- "./simulations/parm_recovery_beta_utility/figs/"

load(paste0(data_dir, "cfg.rdata"))
load(paste0(data_dir, "df.rdata"))

# 1. Compile the model
model <- cmdstan_model(paste0(model_dir, model_name, ".stan"))

# 2. Prepare data list
stan_data <- list(
  N_trials   = nrow(df),
  N_subjects = cfg$Nsubjects,
  N_options  = cfg$Noffer,
  subject    = df$subject,
  offer_A    = df$offer_A,
  offer_B    = df$offer_B,
  is_choice_A = df$is_choice_A
)

#### VB ESTIMATION ####

# 3. Run variational inference
fit <- model$variational(
  data = stan_data,
  algorithm = "meanfield",
  iter = 10000,
  draws = 4000
)

save(fit, file = paste0(data_dir, "fit.rdata"))

#### MCMC SAMPLING ####

# # 3. Run MCMC sampling
# fit <- model$sample(
#   data = stan_data,
#   chains = 4,
#   parallel_chains = 4,
#   iter_warmup = 1000,
#   iter_sampling = 1000
# )

#pars <- as_draws_rvars(fit) # Add all relevant parameters
#pars

#summarize_draws(pars, posterior::rhat, ess_mean, ess_tail)


# AI code 

# # ==========================================
# # PART A: RECOVERY FOR UTILITIES (u_matrix)
# # ==========================================
# u_summary <- fit$summary("u_matrix", "mean")
# estimated_u_matrix <- matrix(u_summary$mean, nrow = 10, ncol = 6)
# 
# u_correlation <- cor(as.vector(estimated_u_matrix), as.vector(true_parameters$u_matrix))
# print(paste("Utility (U) Recovery Correlation:", round(u_correlation, 3)))
# 
# # ==========================================
# # PART B: RECOVERY FOR BETA
# # ==========================================
# beta_summary <- fit$summary("beta", "mean")
# estimated_betas <- beta_summary$mean
# 
# # Check correlation for beta (true_parameters$beta contains the 10 true betas)
# beta_correlation <- cor(estimated_betas, true_parameters$beta)
# print(paste("Beta Recovery Correlation:", round(beta_correlation, 3)))
# 
# # Quick diagnostic plot for Beta
# plot(true_parameters$beta, estimated_betas,
#      xlab = "True Beta", ylab = "Estimated Beta",
#      main = "Beta Parameter Recovery", pch = 19, col = "darkgreen")
# abline(0, 1, col = "red", lty = 2)