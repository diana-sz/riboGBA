rm(list=ls(all=TRUE))

library(viridis)
library(here)
library('rstudioapi') 
library('readODS')
library('nloptr')
library('lpSolve')

directory <- paste0(here(), "/code")
setwd(directory) 

modelname <- "A15"
rna_id <- "rRNA"  # rRNA ID for calculation of ribosome composition
is.reversible <- 0
rescale_kcats <- TRUE
predict.parameters <- TRUE
fer_res_factor <- 5

suppressMessages(source("Readmodelods_v2.R"))

rho_cond <- rho_cond[1]
n_conditions <- 1


source("GBA_Kinetics.R")
source('f0_alt.R')
source("Parameter_prediction.R")
source("GBA_solver.R")
f0_wt <- f0
last_feasible_f0 <- f0


original_r_kcat <- kcatf[r]
original_rRNase_kcat <- kcatf[which(reaction=="rRNase")]
original_K_rRNase <- K["rRNA", "rRNase"]

kdegmax_to_test <- c(1000, 2000, 5000)
Ks_to_test <- c(10, 100)
ribocomp_to_test <- seq(from=0.15, to=0.8, length.out=10)
glucose_concentrations <- c(0.005, 0.01, 0.02, 0.1)

pdf(paste0("../figures/GBA_Model_", modelname, "_ribocomp_test.pdf"),
    title=paste("Parameter testing", modelname), width=5, height=5)

for (kdegmax in kdegmax_to_test){
  for (Kval in Ks_to_test){
    for(glc in glucose_concentrations){
      print(paste0("Running kdegmax=", kdegmax, ", K=", Kval, ", glc=", glc))
      
      x_cond[1,1] <- glc
      results <- data.frame(matrix(ncol=9, nrow=0))
      colnames(results) <- c("glc", "xrp", "kdegmax", "Kval", "mu_opt",
                             "free_rRNA", "total_r_mass", "deg_ratio", "convergence")
      
      for (xrp in ribocomp_to_test){
        ribcomp <- xrp
        kcatf[r] <- original_r_kcat/(xrp/0.36)
        K["rRNA", "rRNase"] <- Kval
        hill <- 3
        kcatf[which(reaction=="rRNase")] <- kdegmax*(1-(xrp^hill/(xrp^hill+0.2^hill)))
        
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
        
        p_opt <- prot(f_opt[cond,]) #opt_state[,paste("p",reaction)]
        rRNA <- c_opt[,grep(rna_id, i_reactant)] 
        
        brRNA <- rRNA*(rRNA/(rRNA + KA[grep(rna_id, reactant),r]) )
        xrp_pred <- p_opt[r]/(p_opt[r] + brRNA)
        
        free_rRNA <- (rRNA-brRNA)/(rRNA)
        
        # Fraction of degraded RNA #####################################################
        
        RNAP   <- f_opt[,grep("rRNAp", reaction)]
        RNase  <- f_opt[,grep("rRNase", reaction)]
        deg_ratio <- RNase/RNAP
        
        
        row <- c(glc, xrp, kdegmax, Kval, mu_opt, free_rRNA, 
                 p_opt[r] + brRNA, deg_ratio, res$convergence)
        
        results[nrow(results)+1,] <- row
        
      }
      
      shapes <- rep(19, nrow(results))
      shapes[results$convergence == -1] <- 21

      results$glc <- factor(results$glc, levels = glucose_concentrations)
      colors <- results$glc
      palette_colors <- viridis(length(levels(colors)))
      color_map <- setNames(palette_colors, levels(colors))
      point_colors <- color_map[as.character(colors)]

      
      plot(mu_opt ~ xrp, data=results,
           xlab = "Protein mass fraction in ribosome",
           ylab = bquote("Growth rate" ~ "[" * h^-1 * "]"),
           main = bquote(k[deg]^max ~ .(kdegmax) ~ K ~ .(Kval)),
           ylim = c(0,1.2), xlim = c(0,1),
           pch = shapes, cex = 1.3,
           col = point_colors,
           cex.lab = 1.3)
      par(new=TRUE)
    }
    
    legend("topleft", legend=glucose_concentrations,
           col=palette_colors)
    par(new=FALSE)
  } 
}
dev.off()
  


