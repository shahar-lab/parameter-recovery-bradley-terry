sim.block <- function(subject, u, beta, cfg){
  
  Noffer  <- cfg$Noffer
  Ntrials <- cfg$Ntrials
  
  df <- data.frame()
  
  for (trial in 1:Ntrials) {
    offer <- sample(1:Noffer, 2, replace = FALSE)
    
    u1 <- u[offer[1]]
    u2 <- u[offer[2]]
    
    choice_prob <- plogis(beta * (u1 - u2))
    
    choice_num <- sample(c(1, 2), size = 1, prob = c(choice_prob, 1 - choice_prob))
    
    choice <- factor(
      ifelse(choice_num == 1, "A", "B"),
      levels = c("A", "B"))
    
    is_choice_A <- as.numeric(choice_num == 1)
    
    df <- rbind(df, data.frame(
      subject = subject,
      trial = trial,
      offer_A = offer[1],
      offer_B = offer[2],
      u1 = u1,
      u2 = u2,
      choice_prob = choice_prob,
      choice = choice,
      is_choice_A = is_choice_A
    ))
  }
  
  return(df)
}