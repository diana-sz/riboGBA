rm(list=ls(all=TRUE))

library(here)
require('rstudioapi') 
require('readODS')
require('nloptr')
require('lpSolve')
library("RColorBrewer")

directory <- paste0(here(), "/code")
setwd(directory) 

is.reversible <- 0
predict.parameters <- 0
maintenance_fun <- "constant"

modelname <- "A8alt_ATPS"

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

results <- data.frame(matrix(ncol=6, nrow=0))

phis_to_test <- c(0, 0.0001, seq(0.005, 0.15, 0.01))
kms_to_test <- c(0.1)
kcats_to_test <- c(200, 10, 2, 0.01)
  
transporter <- which(reaction == "ATPS")
transporter2 <- which(reaction == "ATPS2")
opt_phi <- opt_phis[transporter]

for (km in kms_to_test){
  K[1,transporter2] <- km
  
  for (kcat_t in kcats_to_test){
    kcatf[transporter2] <- kcat_t
    
    # get optimal phis
    source("GBA_solver.R")
    f0_wt <- f0
    last_feasible_f0 <- f0
    # p_opt <- prot(f_opt[1,])
    # opt_phis <- (p_opt/rho_cond)/sum(p_opt/rho_cond)
    # opt_phi <- opt_phis[transporter]
    # mu_orig <- mu_opt
    
    for (fraction in phis_to_test){
      min_phi[transporter] <- fraction
      max_phi[transporter] <- fraction+1e-5
      
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
      results[i, "kcat"] <- kcat_t
      results[i, "km"] <- km
      results[i, "phi"] <- fraction
      results[i, "rel_phi"] <- fraction/opt_phi
      results[i, "mu"] <- mu_opt
      results[i, "mu_norm"] <- mu_opt/mu_orig
      results[i, "convergence"] <- res$convergence
    }
    
    # reset phis
    min_phi <- rep(0, length(min_phi))
    max_phi <- rep(1, length(max_phi))
  }
  
}


pdf(paste0("../figures/",modelname,"_cost_test_", as.integer(is.reversible), ".pdf"),
    title=paste("Protein cost testing", modelname), width=5, height=5)

par(mar=c(5,5,2,1))

#results <- results[results$convergence == 4,]

# Set factor levels
results$kcat <- factor(results$kcat, levels = kcats_to_test)
results$km <- factor(results$km, levels = kms_to_test)

# Set colors based on kcat
palette_colors <- brewer.pal(length(kcats_to_test), "Dark2")
color_map <- setNames(palette_colors, kcats_to_test)
point_colors <- color_map[as.character(results$kcat)]

# Set shapes based on km
shape_list <- 19#:(18 - length(kms_to_test) + 1)
shape_map <- setNames(shape_list, kms_to_test)
point_shapes <- shape_map[as.character(results$km)]

# Create the plot
plot(
  mu_norm ~ rel_phi, data = results,
  xlab = "Titrated level / opt. level",
  ylab = "Fraction of optimal growth rate",
  xlim = c(0, 2), ylim = c(0, 1.1),
  pch = point_shapes,
  col = point_colors,
  cex = 1.3, cex.lab = 1.5
)

# Add legend
legend(
  "bottomright",
  legend = c(kcats_to_test),#, kms_to_test),
  col = c(palette_colors),#, rep("black", length(kms_to_test))),
  pch = c(rep(19, length(kcats_to_test))),#, shape_list),
  cex = 1.3,
  title = "kcat"
)

dev.off()

print(results)

