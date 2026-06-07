rm(list = ls())
library(tidyverse)

source('./models/beta_utility/beta_utility.R')

data_path = './simulations/parm_recovery_beta_utility/data/'

#### SETUP STUDY CONFIG ####

cfg = data.frame(
  Nsubjects = 500,
  Ntrials = 100,
  Noffer  = 6
)
beta_mu = 0.4
beta_sigma = 0.3
beta = rlnorm(cfg$Nsubjects, 0.4, 0.3) #500 values of beta, one per subject. 

# This is equivalent to sample from N(0.4, 0.3) |> exp()
# log(beta) ~ Normal(mu, sigma) but we need beta. So, beta ~ exp(Normal(mu, sigma)) which the same as beta ~ rlnorm(0.4, 0.3)

hist(beta)

u_matrix = matrix(0, nrow = cfg$Nsubjects, ncol = cfg$Noffer)
for (i in 1:cfg$Nsubjects) {
u_matrix[i,]  = rnorm(cfg$Noffer, 0, 1) #for each raw (subject), sample 6 values from N(0, 1)
}


true_parameters = list( beta = beta,
                              u_matrix = u_matrix 
                              )
#### SIMULATE DATA ####

df = df_temp = data.frame()

for (subject in 1:cfg$Nsubjects) {
  
  df_temp = sim.block(subject, 
                      true_parameters$u[subject,], 
                      true_parameters$beta[subject], 
                      cfg)
  df = rbind(df, df_temp)
  
}


  
save(df, file = paste0(data_path,"df.rdata"))
save(true_parameters, file = paste0(data_path,"true_parameters.rdata"))
write.csv(df, paste0(data_path,"df.csv"))
  