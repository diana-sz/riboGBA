rm(list=ls(all=TRUE))

library(viridis)
library(here)
require('rstudioapi') 
require('readODS')
require('nloptr')
require('lpSolve')

directory <- paste0(here(), "/code")
setwd(directory) 

modelname <- "A9ferm4"
is.reversible <- 0
predict.parameters <- 3

suppressMessages(source("Readmodelods_v2.R"))

rho_cond <- rho_cond[1]
n_conditions <- 1

# get optimal solution
source("GBA_Kinetics.R")
source('f0_alt.R')
source("Parameter_prediction.R")
source("GBA_solver.R")

opt_kcats <- kcatf
mu_orig <- mu_opt

results <- data.frame(matrix(ncol=4, nrow=0))
colnames(results) <- c("reaction", "kcat", "mu_opt", "convergence")
kcats_to_test <- c(0.2,0.5,1,2,5)

# ferm_phi <- opt_phis[reaction == "FERM"]
# resp_phi <- opt_phis[reaction == "RESP"]

for (protein in 1:length(kcatf)){
  opt_kcatf <- opt_kcats[protein]
  prot_name <- reaction[protein]

  for (fraction in kcats_to_test){
    kcatf[protein] <- fraction*opt_kcatf

    error_check <- try({source("GBA_solver.R")}, silent = TRUE)
    if(class(error_check) == "try-error"){
      print(paste("Solver error with fraction =", fraction))
      next
    }
    
    i <- nrow(results)+1
    results[i, "protein"] <- prot_name
    results[i, "kcat"] <- kcatf[protein]
    results[i, "mu"] <- mu_opt
    results[i, "convergence"] <- res$convergence
  }
  
  # reset parameters
  kcatf <- opt_kcats
}

results$mu_norm <- results$mu/mu_orig
results$shape <- 19
results$shape[results$convergence == -1] <- 21

pdf(paste0("../figures/",modelname,"_kcat_sensitivity_test_", as.integer(is.reversible), ".pdf"),
    title=paste("Protein cost testing", modelname), width=5, height=5)

par(mar=c(5,5,5,1))
for(prot_name in unique(results$protein)){
  one_prot <- results[results$protein == prot_name,]
  plot(mu_norm ~ kcat, data=one_prot,
       xlab = "kcat",
       ylab = "mu/mu_max", #bquote("Growth rate" ~ "[" * h^-1 * "]"),
       main = paste(prot_name),
       ylim = c(0,2), #xlim = c(0,4),
       pch = one_prot$shape, 
       cex = 1.3,
       cex.lab = 1.3)
  abline(h=1, col="grey70", lty=2)
}

dev.off()


