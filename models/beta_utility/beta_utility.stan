data {
  int<lower=1> N_trials;     // Total number of trials
  int<lower=1> N_subjects;   // Number of unique subjects (10)
  int<lower=1> N_options;    // Number of available options (6)

  array[N_trials] int<lower=1, upper=N_subjects> subject; // tells who is the subject in each trial
  array[N_trials] int<lower=1, upper=N_options> offer_A; // tells what is offer 1 in each trial
  array[N_trials] int<lower=1, upper=N_options> offer_B; // tells what is offer 2 in each trial
  array[N_trials] int<lower=0, upper=1> is_choice_A;  // tells us the choice in each trial
}

parameters {

  // subject-level utilities
  matrix<lower=-1, upper=1>[N_subjects, N_options] u_matrix;

  // group-level beta parameters
  real mu_log_beta;
  real<lower=0> sigma_log_beta;

  // subject-level betas
  vector<lower=0>[N_subjects] beta;

}

model {
  mu_log_beta ~ normal(0, 1);
  sigma_log_beta ~ exponential(1);

  beta ~ lognormal(mu_log_beta, sigma_log_beta);

  for (t in 1:N_trials) {
    is_choice_A[t] ~ bernoulli_logit(
      beta[subject[t]] *
      (u_matrix[subject[t], offer_A[t]] -
        u_matrix[subject[t], offer_B[t]]));
  }
}
