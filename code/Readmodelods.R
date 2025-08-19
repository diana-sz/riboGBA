# Loads model saved as .ods file ###############################################

setwd(paste(directory,"/Models",sep=""))

odsfile <- paste(modelname,".ods", sep = "")

sheets <- list_ods_sheets(odsfile)

nsheets    <- length(sheets)

# Position of parameters
posM          <- (1:nsheets)[sheets == "M"]
posK          <- (1:nsheets)[sheets == "K"]
posKI         <- (1:nsheets)[sheets == "KI"]
posKA         <- (1:nsheets)[sheets == "KA"]
poskcat       <- (1:nsheets)[sheets == "kcat"]
posconditions <- (1:nsheets)[sheets == "conditions"]
posx          <- (1:nsheets)[sheets == "x"]
posq          <- (1:nsheets)[sheets == "q"]
posphi        <- (1:nsheets)[sheets == "phi"]
posrmw        <- (1:nsheets)[sheets == "rmw"]
posminphi     <- (1:nsheets)[sheets == "min_phi"]
posmaxphi     <- (1:nsheets)[sheets == "max_phi"]
posminf       <- (1:nsheets)[sheets == "min_f"]
posmaxf       <- (1:nsheets)[sheets == "max_f"]
posminc       <- (1:nsheets)[sheets == "min_c"]
posmaxc       <- (1:nsheets)[sheets == "max_c"]
posribcomp    <- (1:nsheets)[sheets == "ribcomp"]

# Reads q0, qT, T ##############################################################

#q0  <- as.numeric(read_ods(odsfile, sheet= posq)[1,-1])
#dq0 <- as.numeric(read_ods(odsfile, sheet= posq)[2,-1]) 

# Reads x(t) function ##########################################################

eval(parse(text=noquote(read_ods(odsfile, sheet= posx)[1,2]))) 

# Reads dxdt(t) function #######################################################

eval(parse(text=noquote(read_ods(odsfile, sheet= posx)[2,2])))

# Getting data from sheets

# Mass fraction matrix Mtotal including external reactants #####################

Mtotal <- as.matrix(read_ods(odsfile, sheet= posM)[,-1])

# reaction and reactant names ##################################################

reaction <- colnames(Mtotal)
reactant <- unlist(read_ods(odsfile, sheet= posM)[,1])

rownames(Mtotal) <- reactant

# Michaelis constant matrix K ##################################################
if (length(posK) > 0) {
  
  K <- as.matrix(read_ods(odsfile, sheet= posK)[,-1])
  
} else K <- 0.1*(Mtotal<0)

rownames(K) <- reactant

# inhibition constant matrix KI ################################################
if (length(posKI) > 0) {
  
  KI <- as.matrix(read_ods(odsfile, sheet= posKI)[,-1])
  
} else KI <- 0*K

rownames(KI) <- reactant

# activation constant matrix KA ################################################
if (length(posKA) > 0) {
  
  KA <- as.matrix(read_ods(odsfile, sheet= posKA)[,-1])
  
} else KA <- 0*K

rownames(KA) <- reactant

# kcat
kcatf <- as.numeric(read_ods(odsfile, sheet= poskcat)[1,-1])
kcatb <- as.numeric(read_ods(odsfile, sheet= poskcat)[2,-1])

# Growth condition names #######################################################

condition <- colnames(read_ods(odsfile, sheet= posconditions))[-1]

# Mass density rho at each condition ###########################################

rho_cond <- as.numeric(read_ods(odsfile, sheet= posconditions)[1,-1])

# external concentrations at each condition ####################################

x_cond  <- as.matrix(read_ods(odsfile, sheet= posconditions)[-1,-1])

# Experimental data ############################################################

exp_data <- as.matrix(read_ods(odsfile, sheet= posphi))

# reactant molecular weights ###################################################

if (length(posrmw) > 0) rmw <- as.matrix(read_ods(odsfile, sheet= posrmw)[,-1])

# reads ribosome composition constraint ########################################

if (length(posribcomp) > 0) {
  
  ribcomp <- as.numeric(read_ods(odsfile, sheet= posribcomp)[1,1])
  
} else  ribcomp <- 0

# Definitions ##################################################################

# index for external reactants
n <- 1:dim(x_cond)[1]

# internal matrix M
M <- Mtotal[-n,]

# number of external reactants
nx <- dim(x_cond)[1]

# number of growth conditions
n_conditions <- dim(x_cond)[2]

# names of internal reactants
i_reactant <- reactant[-n]

# number of internal reactants
p <- dim(M)[1]

# number of reactions
r <- dim(M)[2]

# the sum of each M column 
sM <- colSums(M)

# delete numerical artifacts
sM[abs(sM) < 1e-10] <- 0

# indexes for reactions: s (transport), e (enzymatic), and ribosome r 

e <- c(1:(r-1))[sM[1:(r-1)] == 0]  

s <- c(1:(r-1))[sM[1:(r-1)] != 0] 

# indexes: m (metabolite), a (all proteins)

m <- 1:(p-1)

# number of transport reactions
ns <- length(s)

# bounds for minimal and maximal c #############################################

if (length(posminc) > 0) {
  
  min_c <- as.numeric(as.matrix(read_ods(odsfile, sheet= posminc)[,2]))
  
} else min_c <- rep(0,p)

if (length(posmaxc) > 0) {
  
  max_c <- as.numeric(as.matrix(read_ods(odsfile, sheet= posmaxc)[,2]))
  
} else max_c <- rep(rho_cond[1],p)

# bounds for minimal and maximal f #############################################
if (length(posminf) > 0) {
  
  min_f <- as.numeric(read_ods(odsfile, sheet= posminf)[1,])
  
} else min_f <- rep(0,r)

if (length(posmaxf) > 0) {
  
  max_f <- as.numeric(read_ods(odsfile, sheet= posmaxf)[1,])
  
} else max_f <- rep(100,r)

# bounds for minimal phi #######################################################
if (length(posminphi) > 0) {
  
  min_phi <- as.numeric(read_ods(odsfile, sheet= posminphi)[1,])
  
} else min_phi <- rep(0,r)

# bounds for maximal phi
if (length(posmaxphi) > 0) {
  
  max_phi <- as.numeric(read_ods(odsfile, sheet= posmaxphi)[1,])
  
} else max_phi <- rep(1,r)


setwd(directory)

