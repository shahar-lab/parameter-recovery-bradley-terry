data {
  int<lower=1> N_trials;     // Total number of trials
  int<lower=1> N_subjects;   // Number of unique subjects (10)
  int<lower=1> N_options;    // Number of available options (6)

  array[N_trials] int<lower=1, upper=N_subjects> subject; // tells who is the subject in each trial
  array[N_trials] int<lower=1, upper=N_options> offer_1; // tells what is offer 1 in each trial
  array[N_trials] int<lower=1, upper=N_options> offer_2; // tells what is offer 2 in each trial
  array[N_trials] int<lower=0, upper=1> choice;  // tells us the choice in each trial
}

parameters {
  // Matrix containing the estimated utility for each subject and option (10 x 6)
  matrix[N_subjects, N_options] u_matrix;
  
  // Vector containing the estimated sensitivity parameter (beta) for each subject (10)
  vector<lower=0>[N_subjects] beta; 
}

model {
  // PRIORS
  // Standard normal prior for utilities
  target += normal_lpdf(to_vector(u_matrix) | 0, 2);
  
  // Lognormal prior for beta (ensures beta is strictly positive, as forced by lower=0)
  target += lognormal_lpdf(beta | 0, 1); 
  
  // LIKELIHOOD
  for (t in 1:N_trials) {
    int s = subject[t];
    int opt1 = offer_1[t];
    int opt2 = offer_2[t];
    
    // Extract the utilities for the specific subject and presented options
    real u1 = u_matrix[s, opt1];
    real u2 = u_matrix[s, opt2];
    
    // The choice probability is based on: beta * (u1 - u2)
    // This perfectly mirrors 'plogis(beta * (u1 - u2))' from your R simulation
    target += bernoulli_logit_lpmf(choice[t] | beta[s] * (u1 - u2));
  }
}
