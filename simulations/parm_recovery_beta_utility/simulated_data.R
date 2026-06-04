rm(list = ls())
library(tidyverse)

source('./models/beta_utility/beta_utility.R')

data_path = './simulations/parm_recovery_beta_utility/data/'

#### SETUP STUDY CONFIG ####

cfg = data.frame(
  Nsubjects = 10,
  Ntrials = 100,
  Noffer  = 6
)

beta = runif(cfg$Nsubjects, 0.5, 2)

u_matrix = matrix(0, nrow = cfg$Nsubjects, ncol = cfg$Noffer)
for (i in 1:cfg$Nsubjects) {
u_matrix[i,]  = rnorm(cfg$Noffer, 0, 1)
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


df$choice_bin <- ifelse(df$choice == 1, 1, 0)
  
save(df, file = paste0(data_path,"df.rdata"))
save(true_parameters, file = paste0(data_path,"true_parameters.rdata"))
write.csv(df, paste0(data_path,"df.csv"))
  