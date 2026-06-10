rm(list = ls())

#### SETUP ####
library(cmdstanr)
library(posterior)

model_name <- "beta_utility_uncentered"
model_dir  <- paste0("./models/", model_name, "/")
data_dir   <- "./simulations/parm_recovery_beta_utility_uncentered/data/"
figs_dir   <- "./simulations/parm_recovery_beta_utility_uncentered/figs/"

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
fit$save_object(file = paste0(data_dir, "fit.rds"))
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


