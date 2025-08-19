rm(list=ls(all=TRUE))

library(here)

directory <- paste0(here(), "/code")

#setwd(directory) 
#setwd("C:/Users/dajas/RiboComp_GBA/RiboComp_GBA/code")

setwd("~/RiboComp_GBA/code")

source("Interface.R")

modelnames <- c("A9dekel2")#, "A5", "A9ferm4", "A17") #"A9ferm4" A11 "A17", A9_r_inh c("A5", "A7", "A8", "A10")

for(m in modelnames){
  print(paste0("Running ", m))
  modelname_orig <- m

  keep_ribosome_kcat <- FALSE
  keep_transport_kcat <- FALSE
  rescale_kcats <- TRUE

  fer_res_factor <- 10
  maintenance_fun <- "constant"
  modelname <- modelname_orig
  GBA(modelname, 0, 0)

}



# # # # ################################################################################
# modelname="A16"
# suppressMessages(source("Readmodelods_v2.R"))
# source("GBA_Kinetics.R")
# #
# #source('f0.R')
# source('f0_alt.R')
# 
# rho <- rho_cond[1]
# x  <- x_cond[,1]
# 
# mu <- function(f) as.numeric(M[p,r]*f[r]/(tau(ci(f))%*%f))
# ci <- function(f) rho*M%*%f
# print(ci(f0))
# 
# predict.parameters <- 4
# is.reversible <- 0
# source("Parameter_prediction.R")
# "A11ferm"
# kcatf
# print(mu(f0))


# # # source("GBA_Plots_diana.R")


modelname=modelnames[1]
suppressMessages(source("Readmodelods_v2.R"))
# source("GBA_Kinetics.R")

# Assume df is your data frame
X <- as.matrix(M)  # Convert to numeric matrix if not already

# Perform QR decomposition
qr_decomp <- qr(X)

# Get linearly independent columns
independent_cols <- qr_decomp$pivot[1:qr_decomp$rank]

# Identify dependent columns
all_cols <- seq_len(ncol(X))
dependent_cols <- setdiff(all_cols, independent_cols)

# Print results
cat("Linearly independent columns:", independent_cols, "\n")
cat("Linearly dependent columns:", dependent_cols, "\n")

