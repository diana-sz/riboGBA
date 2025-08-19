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

models <- c("A10fill", "A10cons", "A10inh")
phis_to_test <- c(0.001, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1, 2)
kcats_to_test <- c(0, 0.001,0.01,0.1,0.5)

for (modelname in models){
  
  suppressMessages(source("Readmodelods_v2.R"))
  
  rho_cond <- rho_cond[1]
  n_conditions <- 1
  
  # get optimal solution
  source("GBA_Kinetics.R")
  source('f0_alt.R')
  #source("Parameter_prediction.R")
  source("GBA_solver.R")
  f0_wt <- f0
  last_feasible_f0 <- f0
  
  p_opt <- prot(f_opt[1,])
  opt_phis <- (p_opt/rho_cond)/sum(p_opt/rho_cond)
  mu_orig <- mu_opt
  
  results <- data.frame(matrix(ncol = 2+nx+2+p+r+r+r+r, nrow = 0))
  
  colnames(results) <- c("kcat", "phi", "convergence","mu",reactant,paste("tau",reaction),
                           paste("v",reaction),paste("p",reaction),paste("f",reaction))
  

  protein <- which(reaction == "INH")
  
  for (fraction in phis_to_test){
    min_phi[protein] <- fraction
    max_phi[protein] <- fraction+1e-5
    
    for (kcat in kcats_to_test){
      kcatf[protein] <- kcat
    
      f0 <- last_feasible_f0
      error_check <- try({source("GBA_solver.R")}, silent = TRUE)
      if(class(error_check) == "try-error"){
        print(paste("Solver error with fraction =", fraction))
        next
      }
      
      if(res$convergence == -1){
        # try again with the initial solution
        f0 <- f0_wt
        error_check <- try({source("GBA_solver.R")}, silent = TRUE)
        if(class(error_check) == "try-error"){
          print(paste("Solver error with fraction =", fraction))
          next
        }
      }
      # if there is a converged solution, save the latest f0
      if(res$convergence == 4){
        last_feasible_f0 <- f0
      }

      i <- nrow(results)+1
      results[i, "kcat"] <- kcat
      results[i, "phi"] <- fraction
      # results[i, "mu"] <- mu_opt
      # results[i, "convergence"] <- res$convergence
      rho <- rho_cond[1]
      x  <- x_cond[,1]
      f <- f_opt[1,]
      results[i, 3:ncol(results)] <- c(conv[1],mu(f),x,ci(f),tau(ci(f)),v(f),prot(f),f)

    }
    
    # reset phis
    min_phi <- rep(0, length(min_phi))
    max_phi <- rep(1, length(max_phi))
  }
  
  results$mu_norm <- results$mu/mu_orig
  results$shape <- 19
  results$shape[results$convergence == -1] <- 21
  
  pdf(paste0("../figures/",modelname,"_INH_cost_test_", as.integer(is.reversible), ".pdf"),
      title=paste("Protein cost testing", modelname), width=5, height=5)
  
  par(mar=c(5,5,5,1))
  plot(mu_norm ~ phi, data=results,
       xlab = "phi",
       ylab = "mu/mu_orig", #bquote("Growth rate" ~ "[" * h^-1 * "]"),
       main = modelname,
       ylim = c(0,1), xlim = c(0,1),
       pch = results$shape, 
       cex = 1.3,
       cex.lab = 1.3,
       col = as.factor(results$kcat))
  
  legend("topright", col = unique(as.factor(results$kcat)), legend = unique(results$kcat), pch=19)
  dev.off()
  
  write.csv(results, paste0("../data/", modelname, "_useless_protein_test.csv"))
  
}


