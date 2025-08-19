
f0 <- rep(0,r)

rM <- rankMatrix(M)[1]


  # average Km per substrate
  aKm <- colSums(t(K[-n,])) + 2*colSums(t(KA[-n,]))  
  
  aKm[aKm == 0] <- 10
  
  lowerp <- 0.4 # for E. coli ~ 55% of DW is protein [Cayley]
  
  lowerm <- aKm[-p]/1000
  
  lowerb <- c(lowerm,lowerp)
  
  # for models with water
  if ("H2O" %in% i_reactant) {
    
    # normalizes again
    lowerb <- lowerb/(10*sum(lowerb))
    
    lowerb[i_reactant == "H2O"] <- 0.7
    
    # now these refer to the total bouyant density
    lowerb[p] <- lowerp*340/1100
    
  }
  
  # defining parameters
  objective.fn <- 1/kcatf    # to minimize the sum f/kcat
  
  const.mat <- rbind(M,sM,diag(r)) # flux balance, density, non-negative f
  
  const.dir <- c("=",rep(">=",p-2),"=","=",rep(">=",r) )
  
  const.rhs <- c(lowerb, 1, min_f)
  
  # solving model
  
  lp.solution <- lp("min", objective.fn, const.mat, const.dir, const.rhs)
  
  # solution to the linear optimization
  
  f0_linear <- lp.solution$solution
  
  f0 <- f0_linear
  
  #f0
  
  #mu(f0)








