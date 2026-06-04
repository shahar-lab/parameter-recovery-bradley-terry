library(cmdstanr)
library(posterior)

# 1. Compile the model
setwd("~/GitHub/parameter-recovery-bradley-terry/models/beta_utility")
model <- cmdstan_model("beta_utility.stan")
df <- read.csv("C://Users//lihin//OneDrive//מסמכים//GitHub//parameter-recovery-bradley-terry//simulations//parm_recovery_beta_utility//data//df.csv") 

# 2. Prepare data list (assuming 'df' has the 'choice_bin' column)
stan_data <- list(
  N_trials   = nrow(df),
  N_subjects = max(df$subject), 
  N_options  = 6,               
  subject    = df$subject,
  offer_1    = df$offer.1,
  offer_2    = df$offer.2,
  choice     = df$choice_bin
)

# 3. Run MCMC sampling
fit <- model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

pars <- as_draws_rvars(fit) # Add all relevant parameters
pars

summarize_draws(pars, posterior::rhat, ess_mean, ess_tail)


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