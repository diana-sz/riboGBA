# f0 ########################################################
# estimates f0 by minimizing total protein at given mu_data 

# min protein in DW parameter 
min_cp <- 0.6

# average Km per substrate
aKm <- colSums(t(K[-n,])) + 2*colSums(t(KA[-n,]))  

aKm[aKm == 0] <- 10

lowerp <- min_cp # for E. coli ~ 55% of DW is protein [Cayley]

# Here we estimate the lower bound for metabolite concentrations factions := c_m/rho
# for a given protein density fraction min_cp as: the density fraction left for metabolites (1 - min_cp) x 
# the estimated c_m via their average Km, normalized by their total (aKm[-p]/sum(aKm[-p])) divided by a factor of 2
# to give some "slack" for the minimal metabolite bounds (if factor = 1, there would be no slack). Then
# later in the while loop below we always give more "slack" for the lower bounds if the linear optimization
# finds no solution
lowerm <- (1 - min_cp)*(aKm[-p]/sum(aKm[-p]))/2
           
lowerb <- c(lowerm,lowerp)

# for models with water
if ("H2O" %in% i_reactant) {
  
  # normalizes again
  lowerb <- lowerb/(10*sum(lowerb))
  
  lowerb[i_reactant == "H2O"] <- 0.7
  
  # now these refer to the total bouyant density
  lowerb[p] <- min_cp*340/1100
  
}

# now tries to find f0 for the given min_cp, if not successful, tries lower min_cp

f0 <- rep(0,r)

number.tests <- 0

while(sum(f0) == 0 & number.tests < 30 ) {

# defining parameters
objective.fn <- 1/kcatf    # to minimize the sum f/kcat

const.mat <- rbind(M,sM,diag(r)) # flux balance, density, non-negative f

const.dir <- c(rep(">=",p),"=",rep(">=",r) )

const.rhs <- c(lowerb, 1, min_f+0.001)

# solving model

lp.solution <- lp("min", objective.fn, const.mat, const.dir, const.rhs)

# solution to the linear optimization

f0_linear <- lp.solution$solution

f0 <- f0_linear

number.tests <- number.tests + 1

if (sum(f0) == 0) lowerb <- lowerb/1.1

}

print("Initial f0:")
print(f0)






