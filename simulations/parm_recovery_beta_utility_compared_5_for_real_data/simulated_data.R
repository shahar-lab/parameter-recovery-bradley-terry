rm(list = ls())
library(tidyverse)

data_path = './simulations/parm_recovery_beta_utility_compared_5_for_real_data/data/'

#### SETUP STUDY CONFIG ####

raw_data <- read.csv(paste0(data_path, "lihis_data_20_subs_for_bradley_terry.csv"))

cfg = data.frame( 
  Nsubjects = length(unique(raw_data$subject)),
  Noffer  = 7
)

#### PROCESS REAL DATA ####

# We initialize df exactly like the original script
df = df_temp = data.frame()

# Loop through each unique subject to extract their real empirical data
unique_subjects <- unique(raw_data$subject)

for (subject_idx in 1:cfg$Nsubjects) {
  print(subject_idx)  
  
  # Get the actual subject ID from the data
  current_sub <- unique_subjects[subject_idx]
  
  # Subset the real data for the current subject
  df_temp = raw_data[raw_data$subject == current_sub, ]
  
  # Ensure the subject column contains the 1:Nsubjects index for the Stan model
  df_temp$subject = subject_idx
  
  df = rbind(df, df_temp)
}

#### CONVERT OFFER VALUES TO STAN INDICES (1 TO 7) ####

df$offer_A = df$offer_A - 1
df$offer_B = df$offer_B - 1
df$choice  = df$choice - 1

#### SAVE ARTIFACTS ####

save(cfg, file = paste0(data_path,"cfg.rdata"))
save(df, file = paste0(data_path,"df.rdata"))
