rm(list=ls(all=TRUE))

library(viridis)
library(here)
library(RColorBrewer)
library('rstudioapi') 
library('readODS')
library('nloptr')
library('lpSolve')

directory <- paste0(here(), "/code")
setwd(directory) 

is.reversible <- 0
predict.parameters <- 0

modelname <- "A9dekel2"

suppressMessages(source("Readmodelods_v2.R"))

rho_cond <- rho_cond[1]
n_conditions <- 1

phis_to_test <- c(0.001, seq(0.01, 0.2, 0.01))
kcats_to_test <- c(30, 300)
alt_concentrations <- c(1e-5, 0.05, 0.2)


# get optimal solution
source("GBA_Kinetics.R")
source('f0_alt.R')
#source("Parameter_prediction.R")
source("GBA_solver.R")

mu_orig <- mu_opt
f0_wt <- f0
last_feasible_f0 <- f0

results <- data.frame()


for(kcat in kcats_to_test){
  
  for (x_C2 in alt_concentrations){
    
    x_cond[1,1] <- 0.02
    x_cond[2,1] <- x_C2
    
    alt_ind <- which(reaction == "BGAL")
    kcatf[alt_ind] <- kcat
    
    source("GBA_solver.R")
    mu_orig <- mu_opt
    last_feasible_f0 <- f0
    p_opt <- prot(f_opt[1,])
    opt_phis <- (p_opt/rho_cond)/sum(p_opt/rho_cond)
    opt_phi <- opt_phis[alt_ind]
    if(opt_phi < 1e-8){
      opt_phi_nonzero <- FALSE
      opt_phi <- 0.01
    }else{
      opt_phi_nonzero <- TRUE
    }
    
    for (fraction in phis_to_test){
      min_phi[alt_ind] <- fraction
      max_phi[alt_ind] <- fraction+1e-5
      
      f0 <- last_feasible_f0
      error_check <- try({source("GBA_solver.R")}, silent = TRUE)
      if(class(error_check) == "try-error"){
        print(paste("Solver error with fraction =", fraction))
        next
      }
      
      # if(res$convergence == -1){
      #   # try again with the initial solution
      #   f0 <- f0_wt
      #   error_check <- try({source("GBA_solver.R")}, silent = TRUE)
      #   if(class(error_check) == "try-error"){
      #     print(paste("Solver error with fraction =", fraction))
      #     next
      #   }
      # }
      # if there is a converged solution, save the latest f0
      if(res$convergence == 4){
        last_feasible_f0 <- f0
      }
      
      i <- nrow(results)+1
      results[i, "x_C2"] <- x_C2
      results[i, "kcat"] <- kcat
      results[i, "phi"] <- fraction
      results[i, "rel_phi"] <- fraction/opt_phi
      results[i, "mu"] <- mu_opt
      results[i, "mu_norm"] <- mu_opt/mu_orig
      results[i, "opt_phi_nonzero"] <- opt_phi_nonzero
      results[i, "convergence"] <- res$convergence
    }
    
    # reset phis
    min_phi <- rep(0, length(min_phi))
    max_phi <- rep(1, length(max_phi))
  }
  #}
  
  
  results$shape <- 19
  results$shape[results$convergence == -1] <- 21

  par(mar=c(5,2,5,1), mfrow=c(1,2))
  #for(r in unique(results$reaction)){
  #  one_rxn <- results[results$reaction == r,]


  }


write.csv(results, paste0("../data/", modelname, is.reversible, "_dekel.csv"))



