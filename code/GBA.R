# GBA solver ###################################################################

# # Clear variables
# rm(list=ls(all=TRUE))
# 
# require('rstudioapi') 
# require('readODS')
# require('nloptr')
# require('Matrix')
# require('MASS')
# require('lpSolve')
# require('gplots')
# 
# # Sets working directory to source file location ###############################
# 
# directory <- dirname(getActiveDocumentContext()$path)
# 
# setwd(directory) 
# 
# # Model name here ##############################################################
# 
# modelname <- "L3"
# 
# # Forces all reactions to be irreversible if desired
# 
# is.reversible <- 0
# 
# predict.parameters <- 0

# Reads model saved as .ods file ###############################################

suppressMessages(source("Readmodelods.R"))

if (is.reversible == 0) modelname <- paste(modelname,"i",sep="")

# kinetics #####################################################################

source("GBA_Kinetics.R")

# Singular Growth Modes ########################################################

#source('f0.R')
source('f0_alt.R')

# Predicts kinetic parameters based on mu and phi data #########################

if (predict.parameters > 0) source("Parameter_prediction.R")


# Optimization on f ############################################################

source("GBA_solver.R") 


# Exporting results ############################################################

# Exporting csv file with results #######

source("Exportcsv.R")


# Plots #################################

source("GBA_Plots_diana.R")

# dev.off()
