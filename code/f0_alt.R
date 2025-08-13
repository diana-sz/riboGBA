# f0 ########################################################
# estimates f0 by minimizing total protein at given mu_data 
# lower bound for metabolite concentrations is based on experimental composition
# if no solution is found, lower bounds are reduced

f0 <- rep(0,r)

number.tests <- 0

while(sum(f0) == 0 & number.tests < 30 ) {
  
  # defining parameters
  #objective.fn <- 1/kcatf    # to minimize the sum f/kcat
  objective.fn <- rep(1, length(kcatf))

  const.mat <- rbind(M, sM, diag(r)) # flux balance, density, non-negative f
  
  const.dir <- c(rep(">=",p),"=",rep(">=",r))
  
  min_f0 <- min_f
  const.rhs <- c(biomass_lb, 1, min_f0*1.001) # to avoid: Error in is.nloptr(ret) : at least one element in x0 < lb
  
  # solving model
  
  lp.solution <- lp("min", objective.fn, const.mat, const.dir, const.rhs)
  
  # solution to the linear optimization
  print(lp.solution)
  f0 <- lp.solution$solution

  number.tests <- number.tests + 1
  
  if (sum(f0) == 0) biomass_lb <- biomass_lb/1.1
  
}

print("Initial f0:")
print(f0)

# set maintenance energy
maint_index <- which(reaction=="Maint")
min_f[maint_index] <- f0[maint_index]
max_f[maint_index] <- f0[maint_index]+1e-5

rho <- rho_cond[1]
x  <- x_cond[,1]
mu <- function(f) as.numeric(M[p,r]*f[r]/(tau(ci(f))%*%f)) 
ci <- function(f) rho*M%*%f
mu_f0 <- mu(f0)






