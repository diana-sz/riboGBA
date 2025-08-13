rm(list=ls(all=TRUE))

library(here)
require('rstudioapi') 
require('readODS')
require('nloptr')
require('lpSolve')

directory <- paste0(here(), "/code")
setwd(directory) 

is.reversible <- 0


modelname <- "mmsyn_fcr_v2_test2"
print(paste0("Running ", modelname))


suppressMessages(source("Readmodelods_v2.R"))

f0_data <- read.csv("Models/f_opt.csv", sep = ",")

#rho_cond <- rho_cond[1]
#n_conditions <- 1

source("GBA_Kinetics.R")
f0 <- f0_data$f0

#source('f0_alt.R')

source("GBA_solver.R")

source("Exportcsv.R")

source("GBA_Plots_mmsyn.R")
