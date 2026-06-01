sim.block <- function(subject, u, beta, cfg){

  subject = subject
  
  Noffer = cfg$Noffer
  Ntrials = cfg$Ntrials
  
  df = data.frame()


  
  
  for (trial in 1:Ntrials)  { 
    offer = sample(1:Noffer, 2)
    u1 = u[offer[1]]
    u2 = u[offer[2]]
    choice_prob = plogis(beta *(u1 - u2))
    choice = sample(1:2, 1, prob = c(choice_prob, 1 - choice_prob))
    
    df = rbind(df, data.frame(
      subject = subject,
      trial = trial,
      offer   = t(offer),
      u1      = u1,
      u2      = u2,
      choice_prob = choice_prob,
      choice = choice
      
    ))
  }
  return(df)
}

