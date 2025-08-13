# Kinetic parameter prediction

mu_data  <- exp_data[1]

phi_data <- exp_data[-1]

cm0 <-  as.numeric(rho_cond[1]*M%*%f0)[-p]

esat <- predict.parameters

# Estimates Km values assuming typical ratio Km = cm/esat for substrates and 

# forces some reactions to be irreversible if kcatb = 0 in the file
rev <- as.vector(1*(kcatb > 0))

K[(nx+1):(p+nx-1),] <- ceiling(diag(cm0/esat)%*%(1*(M[-p,] < 0)) + is.reversible*diag(cm0)%*%(1*(M[-p,] > 0))%*%diag(rev) )

### test KM ###
# offset <- ci(f0)["p",]*exp_data0/10
# offset_matrix <- matrix(offset, nrow = nrow(K), ncol = ncol(K), byrow = TRUE)
# result <- (K != 0) * 1 * offset_matrix
# 
# cols_to_adjust <- which(!reaction %in% c("Maint", "tC"))
# K[, cols_to_adjust] <- (result[, cols_to_adjust]+K[, cols_to_adjust])/2



# Estimates total saturation of each reaction assuming irrev. with cm = esat*Km

sat <- (esat/(esat + 1))^(colSums(1*(M[-p,] < 0)))


# assumes transporters are completely saturated

#sat[s] <- rep(1,ns)
sat[1] <- 1

# Estimates now kcatf based on ...

b0 <- M%*%f0

kcatf <- round(f0/(b0[p]*phi_data*sat))


# if some kcat was rounded to zero
kcatf[kcatf == 0] <- 10
kcatf[kcatf == Inf] <- 10
kcatf[is.nan(kcatf)] <- 10


fer_ind <- grep("FERM",reaction)
res_ind <- grep("RESP",reaction)
avg_kcat <- (kcatf[fer_ind] + kcatf[res_ind])/2

if ("RESP" %in% reaction & "FERM" %in% reaction) {
  print(paste0("Original kcats: FER=", kcatf[fer_ind], ", RES=", kcatf[res_ind] ))
  
  #ferm_predicted <- f0[fer_ind] > 0
  
  # if(ferm_predicted){
  if(fer_res_factor > 0){
    kcatf[res_ind] <- (1/fer_res_factor)*kcatf[fer_ind]
    new_avg_kcat <- (kcatf[fer_ind] + kcatf[res_ind])/2
    kcatf[c(res_ind, fer_ind)] <- kcatf[c(res_ind, fer_ind)]*(avg_kcat/new_avg_kcat)
  }

  #   
  # }else{
  #   kcatf[fer_ind] <- factor*kcatf[res_ind]
  #   kcatf[c(res_ind, fer_ind)] <- kcatf[c(res_ind, fer_ind)]/(factor/2)
  # }
  
  K[, res_ind] <-  0.2*K[, fer_ind]

  print(paste0("New kcats: FER=", kcatf[fer_ind], ", RES=", kcatf[res_ind] ))
  
}


# estimating kcatb
kcatb <- round(kcatf/10)
kcatb[kcatb == 0] <- 1
kcatb <- is.reversible*rev*kcatb
# ribosome is irreversible
kcatb[r] <- 0

# first separate Km matrices for substrates and for products
KS <- K*(Mtotal<0)
KP <- K*(Mtotal>0)

KS[KS == 0] <- Inf
KP[KP == 0] <- Inf
KI[KI == 0] <- Inf


if(rescale_kcats){
  tol <- 0.05
  all_conds <- rho_cond
  rho_cond <- rho_cond[1]
  n_conditions <- 1
  source("GBA_solver.R") 

  # scale kcats so that mu=1
  print(paste0("Initial mu = ", mu_opt))
  if(abs(mu_opt - mu_data) > 0.05){
    print(paste0(" -> scaling kcats to mu = ", mu_data))
    
    n = 0
    while(abs(mu_opt - mu_data) > 0.05){
      n = n+1
      
      if(mu_opt < mu_data){
        kcatf <- round(kcatf*1.05,1)
      }else{
        kcatf <- round(kcatf*0.95,1)
      }
      source("GBA_solver.R") 
      print(paste0("mu after ", n, " rounds = ", mu_opt))
      
      if(n > 20){
        break
      }
    }
    
    kcatb <- round(kcatf/10)
    kcatb[kcatb == 0] <- 1
    kcatb <- is.reversible*rev*kcatb
    # ribosome is irreversible
    kcatb[r] <- 0
    
  }
  
  source("GBA_solver.R") 
  print(paste0("final mu after ", n, " rounds = ", mu_opt))
  
  rho_cond <- all_conds
  n_conditions <- length(rho_cond)  
}




