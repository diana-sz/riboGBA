rm(list=ls(all=TRUE))

library(viridis)
library(here)

require('rstudioapi') 
require('readODS')
require('nloptr')
require('lpSolve')

directory <- paste0(here(), "/code")
setwd(directory) 

modelname <- "D17rrn"
rna_id <- "rRNA"  # rRNA ID for calculation of ribosome composition
is.reversible <- 1

suppressMessages(source("Readmodelods.R"))
source("GBA_Kinetics.R")
source('f0.R')
source("GBA_solver.R")


results <- data.frame(matrix(ncol=3, nrow=0))
colnames(results) <- c("copies", "mu_opt", "convergence")
original_KA <- KA["DNA", "rRNAp"]

KA["DNA", "rRNAp"] <- 1*10^(-3)
mass_one_rrn <- 2*10^(-5)

copies_to_test <- seq(from=10, to=1, length.out=10)

for(copies in copies_to_test){
  
  min_c[i_reactant=="DNA"] <- mass_one_rrn*copies
  max_c[i_reactant=="DNA"] <- mass_one_rrn*copies+10^(-7)

  #source('f0.R')
  
  error_check <- try({source("GBA_solver.R")}, silent = TRUE)
  if(class(error_check) == "try-error"){
    print(paste("Solver error with xrp =", xrp))
    next
  }
  
  
  row <- c(copies, mu_opt, res$convergence)
  
  results[nrow(results)+1,] <- row
  
}


pdf(paste0("../results/GBA_Model_",modelname,"DNA_test.pdf"),
    title=paste("Parameter testing", modelname), width=5, height=5)

plot(mu_opt ~ copies, data=results,
     xlab = "rRNA DNA",
     ylab = "Growth rate",
     ylim = c(0, max(results$mu_opt)*1.1),
     cex.lab = 1.3)

dev.off()
