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

modelname <- "A9dekel3"

suppressMessages(source("Readmodelods_v2.R"))

rho_cond <- rho_cond[1]
n_conditions <- 1

phis_to_test <- c(0.001, seq(0.01, 0.5, 0.02))
kcats_to_test <- c(34, 150)
alt_concentrations <- c(0.0001, 0.001, 0.01, 0.1)


# get optimal solution
source("GBA_Kinetics.R")
source('f0_alt.R')

if(is.reversible == 1){
  predict.parameters <- 1
  fer_res_factor <- 5
  rescale_kcats <- TRUE
  
  # get optimal solution
  source("Parameter_prediction.R")

}

# KP[3,2] <- 3
source("GBA_solver.R")


kcat_orig <- kcatf
mu_orig <- mu_opt
f0_wt <- f0
last_feasible_f0 <- f0

results <- data.frame()

for(protein in c("tC2", "tI")){
  for(kcat in kcats_to_test){
    for(x_C2 in alt_concentrations){
      
      x_cond[1,1] <- 0.02
      if(protein == "tC2"){
        x_cond[2,1] <- x_C2
        x_cond[3,1] <- 1e-5
      }else{
        x_cond[3,1] <- x_C2
        x_cond[2,1] <- 1e-5
      }
      
      alt_ind <- which(reaction == protein)
      kcatf <- kcat_orig  # reset kcats
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
        results[i, "protein"] <- protein
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

    results$shape <- 19
    results$shape[results$convergence == -1] <- 21
    
    
  }
}

write.csv(results, paste0("../data/", modelname, is.reversible, "_dekel.csv"))



