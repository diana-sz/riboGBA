rm(list=ls(all=TRUE))

library(viridis)
library(here)
require('rstudioapi') 
require('readODS')
require('nloptr')
require('lpSolve')

directory <- paste0(here(), "/code")
setwd(directory) 

modelname <- "P17_1"
rna_id <- "rRNA"  # rRNA ID for calculation of ribosome composition
is.reversible <- 1

suppressMessages(source("Readmodelods.R"))
source("GBA_Kinetics.R")
source('f0_alt.R')
source("GBA_solver.R")


results <- data.frame(matrix(ncol=7, nrow=0))
colnames(results) <- c("xrp", "kdegmax", "Kval", "mu_opt", "free_rRNA", "total_r_mass", "deg_ratio", "convergence")
original_r_kcat <- kcatf[r]
original_rRNase_kcat <- kcatf[which(reaction=="rRNase")]
original_K_rRNase <- K["rRNA", "rRNase"]

kdegmax_to_test <- c(5000,10000,20000)
Ks_to_test <- c(1,10,30)
ribocomp_to_test <- seq(from=0.1, to=0.8, length.out=10)


pdf(paste0("../results/GBA_Model_",modelname,"ribocomp_test.pdf"),
    title=paste("Parameter testing", modelname), width=5, height=5)

for (kdegmax in kdegmax_to_test){
  for (Kval in Ks_to_test){
    for (xrp in ribocomp_to_test){
      
      ribcomp <- xrp
      kcatf[r] <- original_r_kcat/(xrp/0.36)
      K["rRNA", "rRNase"] <- Kval
      hill <- 4
      kcatf[which(reaction=="rRNase")] <- kdegmax*(1-(xrp^hill/(xrp^hill+0.2^hill)))
      
      #suppressMessages(source('f0.R'))
      
      error_check <- try({source("GBA_solver.R")}, silent = TRUE)
      if(class(error_check) == "try-error"){
        print(paste("Solver error with xrp =", xrp))
        next
      }
      
      p_opt <- prot(f_opt[cond,]) #opt_state[,paste("p",reaction)]
      rRNA <- c_opt[,grep(rna_id, i_reactant)] 
      
      brRNA <- rRNA*(rRNA/(rRNA + KA[grep(rna_id, reactant),r]) )
      xrp_pred <- p_opt[r]/(p_opt[r] + brRNA)
      
      free_rRNA <- (rRNA-brRNA)/(rRNA)
      
      # Fraction of degraded RNA #####################################################
      
      RNAP   <- f_opt[,grep("rRNAp", reaction)]
      RNase  <- f_opt[,grep("rRNase", reaction)]
      deg_ratio <- RNase/RNAP
      
      
      row <- c(xrp, kdegmax, Kval, mu_opt, free_rRNA, 
               p_opt[r] + brRNA, deg_ratio, res$convergence)
      
      results[nrow(results)+1,] <- row
      
    }
    
    shapes <- rep(19, nrow(results))
    non_converged <- which(results$convergence == -1)
    
    shapes[non_converged] <- 21
    plot(mu_opt ~ xrp, data=results[(results$kdegmax == kdegmax) | (results$Kval == Kval),],
         xlab = "Protein mass fraction in ribosome",
         ylab = bquote("Growth rate" ~ "[" * h^-1 * "]"),
         main = bquote(k[deg]^max ~ .(kdegmax) ~ K ~ .(Kval)),
         ylim = c(0,1), xlim = c(0,1),
         pch = shapes, cex = 1.3,
         cex.lab = 1.3)
  }
  
  }

dev.off()


