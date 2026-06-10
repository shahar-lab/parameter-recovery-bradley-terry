data {
  int<lower=1> N_trials;
  int<lower=1> N_subjects;
  int<lower=1> N_options;

  array[N_trials] int<lower=1, upper=N_subjects> subject;
  array[N_trials] int<lower=1, upper=N_options> offer_A;
  array[N_trials] int<lower=1, upper=N_options> offer_B;
  array[N_trials] int<lower=0, upper=1> is_choice_A;
}

parameters {
  // Utility parameters - no mean/SD constraint
  matrix[N_subjects, N_options] u_matrix;

  // Group-level beta parameters
  real mu_log_beta;
  real<lower=0> sigma_log_beta;

  // Subject-level betas
  vector<lower=0>[N_subjects] beta;
}

model {
  // Hierarchical priors for beta
  mu_log_beta ~ normal(0, 1);
  sigma_log_beta ~ exponential(1);

  beta ~ lognormal(mu_log_beta, sigma_log_beta);

  // Prior on utilities
  to_vector(u_matrix) ~ normal(0, 1);

  // Likelihood
  for (t in 1:N_trials) {
    is_choice_A[t] ~ bernoulli_logit(
      beta[subject[t]] *
      (
        u_matrix[subject[t], offer_A[t]] -
        u_matrix[subject[t], offer_B[t]]
      )
    );
  }
}
