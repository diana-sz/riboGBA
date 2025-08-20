rm(list=ls(all=TRUE))

library(viridis)
library(here)
library(RColorBrewer)
require('rstudioapi') 
require('readODS')
require('nloptr')
require('lpSolve')

directory <- paste0(here(), "/code")
setwd(directory) 

is.reversible <- 0
# predict.parameters <- 3
# fer_res_factor <- 5
# maintenance_fun <- "constant"
# rescale_kcats <- TRUE

modelname <- "A7simple_p"


suppressMessages(source("Readmodelods_v2.R"))

rho_cond <- rho_cond[1]
n_conditions <- 1
# x_cond[1,1] <- x_Glc

phis_to_test <- c(0, 0.001, seq(0.01, 0.5, 0.01))


# get optimal solution
source("GBA_Kinetics.R")
source('f0_alt.R')
source("GBA_solver.R")

results <- data.frame()

f0_wt <- f0
last_feasible_f0 <- f0
p_opt <- prot(f_opt[1,])
opt_phis <- (p_opt/rho_cond)/sum(p_opt/rho_cond)
mu_orig <- mu_opt

for(protein in 1:length(reaction)){
  opt_phi <- opt_phis[protein]
  if(opt_phi < 1e-8){
    opt_phi_nonzero <- FALSE
    opt_phi <- 0.01
  }else{
    opt_phi_nonzero <- TRUE
    }
  
  for (fraction in phis_to_test){
    min_phi[protein] <- fraction
    max_phi[protein] <- fraction+1e-5
    
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
    results[i, "x_Glc"] <- x_cond[1,1]
    results[i, "reaction"] <- reaction[protein]
    results[i, "phi"] <- fraction
    results[i, "rel_phi"] <- fraction/opt_phi
    results[i, "mu"] <- mu_opt
    results[i, "mu_norm"] <- mu_opt/mu_orig
    results[i, "opt_phi_nonzero"] <- opt_phi_nonzero
    results[i, "convergence"] <- res$convergence
    fs <- f_opt[1,]
    results[i, paste0("f.", reaction)] <- fs
    results[i, paste0("v.", reaction)] <- v(fs)
    results[i, paste0("p.", reaction)] <- prot(fs)
    results[i, paste0("c.", reactant)] <- c(x, ci(fs))
  }
  
  # reset phis
  min_phi <- rep(0, length(min_phi))
  max_phi <- rep(1, length(max_phi))
}


results$shape <- 19
results$shape[results$convergence == -1] <- 21
write.csv(results, paste0("../data/", modelname, is.reversible, "_protein_cost.csv"))

pdf(paste0("../figures/",modelname,"_conditions_cost_test_", as.integer(is.reversible), ".pdf"),
    title=paste("Protein cost testing", modelname), width=10, height=5)

par(mar=c(5,5,5,1), mfrow=c(1,2))
for(r in unique(results$reaction)){
  one_rxn <- results[results$reaction == r,]
  colors <- as.factor(one_rxn$x_Glc)
  plot(mu_norm ~ rel_phi, data=one_rxn[one_rxn$opt_phi_nonzero == TRUE, ],
       xlab = "Titrated level / opt. level",
       ylab = "mu", #bquote("Growth rate" ~ "[" * h^-1 * "]"),
       ylim = c(0,1.1), xlim = c(0,4),
       pch = results$shape,
       cex = 1.3,
       cex.lab = 1.3,
       main = r,
       col = colors)
  plot(mu ~ phi, data=one_rxn,
       xlab = "Proteome fraction",
       ylab = bquote("Growth rate" ~ "[" * h^-1 * "]"),
       ylim = c(0,1.1), xlim = c(0,0.15),
       pch = results$shape, 
       cex = 1.3,
       cex.lab = 1.3,
       main = r,
       col = colors)
  #abline(v = one_rxn$phi[which.max(one_rxn$mu_norm)])
}

dev.off()


