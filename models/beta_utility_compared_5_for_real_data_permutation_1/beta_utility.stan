data {
  int<lower=1> N_trials;     // Total number of trials
  int<lower=1> N_subjects;   // Number of unique subjects
  int<lower=1> N_options;    // Number of available options

  array[N_trials] int<lower=1, upper=N_subjects> subject; // Subject ID for each trial
  array[N_trials] int<lower=1, upper=N_options> offer_A;  // Option A ID
  array[N_trials] int<lower=1, upper=N_options> offer_B;  // Option B ID
  array[N_trials] int<lower=0, upper=1> is_choice_A;      // Choice indicator (1 if A, 0 if B)
}

parameters {
  // Raw, unconstrained utility parameters
  matrix[N_subjects, N_options] u_raw;

  // Group-level beta parameters
  real mu_log_beta;
  real<lower=0> sigma_log_beta;

  // Subject-level betas
  vector<lower=0>[N_subjects] beta;
}

transformed parameters {
  // Constrained utilities
  matrix[N_subjects, N_options] u_matrix;
  
  // Hard constraint: Force mean = 0 and variance = 1 for each subject
  for (i in 1:N_subjects) {
    real mu_u = mean(to_vector(u_raw[i, ]));
    real sd_u = sd(to_vector(u_raw[i, ]));
    
    // Standardize the raw parameters
    u_matrix[i, ] = (u_raw[i, ] - mu_u) / sd_u; 
  }
}

model {
  // Hierarchical priors for beta
  mu_log_beta ~ normal(0, 2);
  sigma_log_beta ~ exponential(2);
  beta ~ lognormal(mu_log_beta, sigma_log_beta);

  // Weak prior on u_raw to provide initial geometry before transformation
  to_vector(u_raw) ~ normal(0, 1);
  
  // Likelihood function
  for (t in 1:N_trials) {
    is_choice_A[t] ~ bernoulli_logit(
      beta[subject[t]] * (u_matrix[subject[t], offer_A[t]] - 
       u_matrix[subject[t], offer_B[t]])
    );
  }
}
