# Singular Growth Modes ########################################################
f0 <- NA  # initiate f0

rM <- rankMatrix(M)[1]

#H <-  t(M)%*%M

#G <- M%*%t(M)

# eigenvalues (= singular values^2)
#lambda <- eigen(H)$values

# Singular Value Decomposition of M
U   <- svd(M,nu=p,nv=r)$u 

V <- svd(M,nu=p,nv=r)$v 


U[abs(U) < 1e-10] <- 0

V[abs(V) < 1e-10] <- 0

sigma <- svd(M,nu=p,nv=r)$d

# adds zeros if necessary
#if (rM < max(p,r)) sigma <- c(sigma,rep(0,max(p,r)-rM))

# redefine U,V such that sign of V^r are always positive (and then SVD is also #
# unique)

correction <- sign(V[r,1:rM])

for (i in 1:rM) U[,i] <- U[,i]*correction[i]

for (i in 1:rM) V[,i] <- V[,i]*correction[i]

################################################################################

# definition of mode values ####################################################

alpha <- 1:r
beta  <- 1:p

# rho at first condition
rho <- rho_cond[1]

# external concentrations x at first condition
x  <- x_cond[,1]

#h <- (1/rho)*colSums(U)*sigma

#ci <-  U[,4]/sum(U[,4])*340

#orderh <- order(h)

# cost-benefit balance for protein production ##################################
# (assuming inactive null modes x^nu = 0)

if (r>rM) nu <- (rM+1):r

#V[r,beta] - U[p,beta]*sigma[beta]

# Matrix \Sigma
#S <-  diag(sigma)

#S  <- S[beta,alpha]

#maxh <- orderh[p]


# test whether there is SGM with all non-negative c

signU <- 1*(U > 0)

is.BGS <- as.numeric(colSums(signU) == p)

pBGS <- grep(p,colSums(signU))

# estimate f0 based on SGM

# if mode with highest value is a valid
if (length(pBGS) == 1) {
  print("Estimating f0 based on SGM")

  # for a pure state (only x[maxh] != 0), then

  f0_sgm <- V[,pBGS]/(sum(U[,pBGS])*sigma[pBGS])


  b0_sgm <- M%*%f0_sgm

  f0 <- f0_sgm

  #print("Initial f0 (BGM):")
  #print(f0)
  #print("mu(f0)")
  #print(mu(f0))
  #print("tau(ci(f0_sgm))")
  #print(tau(ci(f0_sgm)))

}


# if there is no BGM, we guess initial f0 by minimizing total protein at given mu_data 
if((length(pBGS) == 0) | (sum(f0 < 0) > 0)) {

  print("Estimating f0 by minimizing total protein")
  # average Km per substrate
  aKm <- colSums(t(K[-n,])) + 2*colSums(t(KA[-n,]))  
  
  aKm[aKm == 0] <- 10
  
  lowerp <- 0.55 # for E. coli ~ 55% of DW is protein [Cayley]
  
  lowerm <- aKm[-p]/1000
  
  lowerb <- c(lowerm,lowerp)
  
  # for models with water
  if ("H2O" %in% i_reactant) {
    
    # normalizes again
    lowerb <- lowerb/(10*sum(lowerb))
    
    lowerb[i_reactant == "H2O"] <- 0.7
    
    # now these refer to the total bouyant density
    lowerb[p] <- 0.55*340/1100
    
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

}






