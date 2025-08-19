rm(list=ls(all=TRUE))

library(viridis)
library(here)
require('rstudioapi') 
require('readODS')
require('nloptr')
require('lpSolve')

directory <- paste0(here(), "/code")
setwd(directory) 

is.reversible <- 0
predict.parameters <- 0
#fer_res_factor <- 5
maintenance_fun <- "constant"
#rescale_kcats <- TRUE


modelname <- "A10inh"


suppressMessages(source("Readmodelods_v2.R"))

rho_cond <- rho_cond[1]
n_conditions <- 1

# get optimal solution
source("GBA_Kinetics.R")
source('f0_alt.R')
#source("Parameter_prediction.R")
source("GBA_solver.R")

f0_wt <- f0
p_opt <- prot(f_opt[1,])
opt_phis <- (p_opt/rho_cond)/sum(p_opt/rho_cond)
mu_orig <- mu_opt

results <- data.frame(matrix(ncol=4, nrow=0))
colnames(results) <- c("protein", "phi", "mu_opt", "convergence")
phis_to_test <- c(0.001, 0.01,0.02, 0.05, 0.1, 0.2, 0.4)
ferm_phi <- opt_phis[reaction == "FERM"]
resp_phi <- opt_phis[reaction == "RESP"]

inhibitor <- which(reactant == "INH")
inhibition <- which(reaction == "INH")

for (KI_inh in c(100,10,1,0.1)){
  KI[inhibitor, r] <- KI_inh
  
  for (fraction in phis_to_test){
    min_phi[inhibition] <- fraction
    max_phi[inhibition] <- fraction+1e-5

    error_check <- try({source("GBA_solver.R")}, silent = TRUE)
    if(class(error_check) == "try-error"){
      print(paste("Solver error with fraction =", fraction))
      next
    }
    
    f0 <- f0_wt
    
    i <- nrow(results)+1
    results[i, "KI"] <- KI_inh
    results[i, "phi"] <- fraction
    results[i, "mu"] <- mu_opt
    results[i, "convergence"] <- res$convergence
  }
  
  # reset phis
  min_phi <- rep(0, length(min_phi))
  max_phi <- rep(1, length(max_phi))
}

results$mu_norm <- results$mu/mu_orig
results$shape <- 19
results$shape[results$convergence == -1] <- 21

pdf(paste0("../figures/",modelname,"_KI_cost_test_", as.integer(is.reversible), ".pdf"),
    title=paste("Protein cost testing", modelname), width=5, height=5)

par(mar=c(5,5,5,1))
plot(mu_norm ~ phi, data=results,
     xlab = "phi",
     ylab = "mu/mu_orig", #bquote("Growth rate" ~ "[" * h^-1 * "]"),
     ylim = c(0,1), xlim = c(0,0.5),
     pch = results$shape, 
     cex = 1.3,
     cex.lab = 1.3,
     col = as.factor(results$KI))


dev.off()


