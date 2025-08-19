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
predict.parameters <- 0
fer_res_factor <- 5
maintenance_fun <- "constant"
rescale_kcats <- TRUE


modelname <- "A10cons"


suppressMessages(source("Readmodelods_v2.R"))

rho_cond <- rho_cond[1]
n_conditions <- 1

phis_to_test <- c(0.001, seq(0.01, 0.30, 0.01))
kcats_to_test <- c(1e-5, 5, 20, 2000)
glucose_concentration <- 0.003

# get optimal solution
source("GBA_Kinetics.R")
source('f0_alt.R')
#source("Parameter_prediction.R")
source("GBA_solver.R")
mu_orig <- mu_opt
f0_wt <- f0
last_feasible_f0 <- f0

results <- data.frame()

#for (x_Glc in glucose_concentrations){
for(kcat in kcats_to_test){
  
  x_cond[1,1] <- glucose_concentration
  
  
  ferm_ind <- which(reaction == "FERM")
  resp_ind <- which(reaction == "RESP")
  
  kcatf[ferm_ind] <- kcat

  source("GBA_solver.R")
  mu_orig <- mu_opt
  last_feasible_f0 <- f0
  p_opt <- prot(f_opt[1,])
  opt_phis <- (p_opt/rho_cond)/sum(p_opt/rho_cond)
  opt_phi <- opt_phis[resp_ind]
  if(opt_phi < 1e-8){
    opt_phi_nonzero <- FALSE
    opt_phi <- 0.01
  }else{
    opt_phi_nonzero <- TRUE
    }
    
    for (fraction in phis_to_test){
      min_phi[resp_ind] <- fraction
      max_phi[resp_ind] <- fraction+1e-5
      
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
      #results[i, "x_Glc"] <- x_Glc
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

pdf(paste0("../figures/",modelname,"_conditions2_cost_test_", as.integer(is.reversible), ".pdf"),
    title=paste("Protein cost testing", modelname), width=10, height=5)

par(mar=c(5,2,5,1), mfrow=c(1,2))
#for(r in unique(results$reaction)){
#  one_rxn <- results[results$reaction == r,]
  

results$kcat <- factor(results$kcat, levels = kcats_to_test)
colors <- results$kcat
palette_colors <- brewer.pal(length(levels(colors)), "Dark2")
color_map <- setNames(palette_colors, levels(colors))
point_colors <- color_map[as.character(colors)]


plot(mu_norm ~ rel_phi, data=results,
     xlab = "Titrated level / opt. level",
     ylab = "mu/mu_max", #bquote("Growth rate" ~ "[" * h^-1 * "]"),
     ylim = c(0,1.1), xlim = c(0, 4), #max(results$rel_phi)),
     pch = results$shape,
     cex = 1.3,
     cex.lab = 1.3,
     col = point_colors)
legend("bottomright",
       legend = kcats_to_test,
       col = palette_colors,
       pch = 19, title="FERM kcat")

plot(mu ~ phi, data=results,
     xlab = "Proteome fraction",
     ylab = bquote("Growth rate" ~ "[" * h^-1 * "]"),
     ylim = c(0, max(results$mu)*1.1), xlim = c(0, max(phis_to_test)),
     pch = results$shape, 
     cex = 1.3,
     cex.lab = 1.3,
     col = point_colors)
  #abline(v = one_rxn$phi[which.max(one_rxn$mu_norm)])
#}

legend("bottomright",
       legend = kcats_to_test,
       col = palette_colors,
       pch = 19, title="FERM kcat")


dev.off()


