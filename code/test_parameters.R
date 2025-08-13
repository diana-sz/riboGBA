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

# to do - make this into a csv and import it / make it automatic to test all parameters
tested_parameters <- read.csv("../data/parameters_to_test.csv")

# store original parameters so they can be reset after each round
original_K <- K
original_KA <- KA
original_kcatf <- kcatf

pdf(paste0("../figures/GBA_Model_",modelname,"_parameter_test.pdf"),
    title=paste("Parameter testing", modelname), width=15, height=5)


for(par_i in 1:nrow(tested_parameters)){
  K <- original_K
  KA <- original_KA
  kcatf <- original_kcatf
  
  results <- data.frame(matrix(ncol=11, nrow=0))
  colnames(results) <- c("parameter_type", "Ki" , "Kj" , "k_val" , "kcat_val", 
                         "mu_opt", "xRP", "free_rRNA", "total_r_mass", "r_mass", "convergence")
  
  Ki <- tested_parameters$Ki[par_i]
  Kj <- tested_parameters$Kj[par_i]

  if(grepl("KA", Ki) > 0){
    parameter_type <- "KA"
    Ki <- gsub("KA ", "", Ki)
    Kj <- gsub("KA ", "", Kj)
    wt_K <- KA[Ki, Kj]
    
  }else{
    parameter_type <- "KM"
    Ki <- gsub("KM ", "", Ki)
    Kj <- gsub("KM ", "", Kj)
    wt_K <- K[Ki, Kj]
  }
  
  wt_kcat <- kcatf[which(reaction == Kj)]
  
  
  # generate range to test
  tested_kcats <- exp(seq(log(wt_kcat/50), log(wt_kcat*50), length.out=8))
  tested_Ks <- exp(seq(log(wt_K/50), log(wt_K*50), length.out=8))
  
  for(kcat_val in tested_kcats){
    r_kcat <- kcatf[r]
    rp_mass <- 21*109*3600/r_kcat
    kcatf[which(reaction == Kj)] <- kcat_val
    
    for(k_val in tested_Ks){
      
      if(parameter_type == "KA"){
        KA[Ki, Kj] <- k_val
      }else{
        K[Ki, Kj] <- k_val
      }
      
      # run optimization
      #source("GBA_Kinetics.R")
      source('f0_alt.R')
      
      error_check <- try({source("GBA_solver.R")}, silent = TRUE)
      if(class(error_check) == "try-error"){
        print(paste("Solver error with kcat =", kcat_val, " KA =", k_val))
        next
      }
      if(res$convergence != 4){
        print(paste("Model did not converge with kcat =", kcat_val, " KA =", k_val))
        #next
      }
    
      # opt_state <- matrix(rep(0,(nx+2+p+r+r+r+r)*n_conditions),nrow = n_conditions)
      # for (cond in 1:n_conditions) {
      #   
      #   rho <- rho_cond[cond]
      #   x  <- x_cond[,cond]
      #   f <- f_opt[cond,]
      #   opt_state[cond,] <- c(conv[cond],mu(f),x,ci(f),tau(ci(f)),v(f),prot(f),f)
      # }
      # colnames(opt_state) <- c("convergence","mu",reactant,paste("tau",reaction),
      #                          paste("v",reaction),paste("p",reaction),paste("f",reaction))
  
      p_opt <- prot(f_opt[cond,]) #opt_state[,paste("p",reaction)]
      rRNA <- c_opt[,grep(rna_id, i_reactant)] 
      
      brRNA <- rRNA*(rRNA/(rRNA + KA[grep(rna_id, reactant),r]) )
      xRP <- p_opt[r]/(p_opt[r] + brRNA)
      
      free_rRNA <- (rRNA-brRNA)/(rRNA)

      row <- c(parameter_type, Ki, Kj, k_val, kcat_val, mu_opt, xRP, free_rRNA, 
               p_opt[r] + brRNA, rp_mass / xRP, res$convergence)

      results[nrow(results)+1,] <- row
    }
  }
  
  cols_to_convert <- c( "k_val" , "kcat_val", "mu_opt", "xRP", "free_rRNA", "total_r_mass", "r_mass")
  results[cols_to_convert] <- lapply(results[cols_to_convert], as.numeric)
  

  par(mfrow=c(1,3))
  mu_normalized <- results$mu_opt/max(results$mu_opt)
  color_palette <- magma(100, end=0.9)
  color.legend <- rev(colorRampPalette(color_palette)(10))
  mus_legend <-  rev(seq(min(results$mu_opt), max(results$mu_opt), length.out = 10))

  palette_mapped <- color_palette[as.numeric(cut(mu_normalized, breaks = 100))]
  shapes <- rep(19, nrow(results))
  non_converged <- which(results$convergence == -1)
  #shapes[non_converged] <- 21
  
  plot(results$kcat_val, results$xRP,
       col = palette_mapped,
       pch = shapes, cex = 1.5, ylim = c(0,1), log = "x",
       xlab=paste("kcat", Kj),
       ylab="Protein mass fraction ribosome")
  legend("bottomright", legend = round(mus_legend, 2), 
         fill = color.legend, title = "Growth rate", cex = 0.8)
  title(main=paste("Parameters of", Kj), outer=TRUE, cex.main = 2, line=-2)
  
  plot(results$k_val, results$xRP,
       col = palette_mapped,
       pch = shapes, cex = 1.5,  ylim = c(0,1), log = "x",
       xlab=paste(parameter_type, Ki),
       ylab="Protein mass fraction ribosome")

  plot(results$kcat_val/results$k_val, results$xRP,
       col = palette_mapped,
       pch = shapes, cex = 1.5,  ylim = c(0,1), log = "x",
       xlab=paste("kcat/",parameter_type, Ki),
       ylab="Protein mass fraction ribosome")

  # plot(results$mu_opt, results$xRP, ylim = c(0,1), pch=19,
  #      xlab = "Growth rate", ylab = "%")
  # par(new=TRUE)
  # plot(results$mu_opt, results$free_rRNA, ylim = c(0,1), pch=18, 
  #      xlab = NA, ylab = NA)
  # par(new=FALSE)
  par(mfrow=c(1,1))
  #legend("bottomright", legend = c("% protein in r", "% free RNA"), pch = c(19,18))

}

dev.off()