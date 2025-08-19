rm(list=ls(all=TRUE))

library(here)
require('rstudioapi') 
require('readODS')
require('nloptr')
require('Matrix')
require('MASS')
require('lpSolve')
require("latex2exp")
require("lattice")
require("deSolve")
require("rgl")
require("plot3D")
require("car")
require('gplots')
require('Rlibeemd')
#require('sfsmisc')
require('plotly')

directory <- paste0(here(), "/code")
setwd(directory) 

orig_name <- "B17b"
modelname <- orig_name
predict.parameters <- 1
is.reversible <- 1

suppressMessages(source("Readmodelods.R"))
source("GBA_Kinetics.R")
source('SGM.R')
source("Parameter_prediction.R")

original_kcat <- kcatb

for(i in 1:length(kcatb)){
  # reset kcats to wt
  kcatb <- original_kcat
  
  # skip if it is always irreversible
  if (kcatb[i] == 0){next}
  
  # make reaction i irreversible
  kcatb[i] <- 0
  current_reaction <- reaction[i]
  print(paste("*** Running model with", current_reaction, "irreversible ***"))
  
  source("GBA_Kinetics.R")
  source('SGM.R')
  source("GBA_solver.R")

  modelname <- paste0(orig_name, current_reaction, "i")
  source("Exportcsv.R")
  source("GBA_Plots_diana.R")
  
  
  # Make everything irreversible except reaction i
  kcatb <- rep(0, length(original_kcat))
  kcatb[i] <- original_kcat[i]
  
  print(paste("*** Running model with", current_reaction, "reversible ***"))
  
  source("GBA_Kinetics.R")
  source('SGM.R')
  source("GBA_solver.R")

  modelname <- paste0(orig_name, current_reaction, "r")
  source("Exportcsv.R")
  source("GBA_Plots_diana.R")
  
}