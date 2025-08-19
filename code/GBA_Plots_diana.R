# plot results #######################################################################################################
library(RColorBrewer)

setwd(paste(directory,"/Results GBA",sep=""))


add_axes <- function(){
  axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
  axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
  axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
  axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
}

# extract feasible solutions and relevant results
# feasible <- conv == 4
# mu_opt <- mu_opt[feasible]
# opt_state <- opt_state[feasible,]
# c_opt <- c_opt[feasible,]

f_opt <- opt_state[,paste("f",reaction)]
p_opt <- opt_state[,paste("p",reaction)]
phi_opt <- p_opt/c_opt[,p]
cm_opt <- c_opt[,-p]
v_opt <- opt_state[,paste("v",reaction)]
tau_opt <- opt_state[,paste("tau",reaction)]


plots_per_row <- min(5,r)

pdf(paste("GBA Model ",modelname," ",solver,", mean time (",mean_time,"s) results.pdf",sep=""),
    title="Optimization results",  width=2.5*plots_per_row,height=5)
par(mfrow=c(2,plots_per_row))


# graphic
coloraxis <- "grey10"
colorlab <- "grey10"

maxmu <- max(mu_opt)+0.15
lowmu  <-  max(mu_opt[conv == 4])/4 # min(mu_opt) 
highmu <- max(mu_opt[conv == 4]) #  4*max(mu_opt)/5
mediummu <- lowmu < mu_opt & mu_opt < highmu


# reaction colors
if (r <= 12) colorj <- brewer.pal(r, "Paired")
if (r > 12) colorj <- c(brewer.pal(12, "Paired"), rainbow(r-12))

# metabolite colors
if (r <= 12) colori <- brewer.pal(p, "Paired")
if (r > 12) colori <- c(brewer.pal(12, "Paired"), rainbow(p-12))


# Predicting growth laws #####################################################
# 
# # average v at medium mu
# av <- rep(0,r)
# for (j in 1:r) av[j] <- mean(v_opt[mediummu,j])
# 
# # estimate dtau/d ln mu
# dtaudlnmu <- rep(0,r)
# for (j in 1:r) dtaudlnmu[j] <- lm(tau_opt[mediummu,j] ~ log(mu_opt[mediummu]) )$coefficients[2]
# 
# # acp = average cp at medium growth
# 
# acp <- mean(c_opt[mediummu,p])
# 
# # estimation phi0
# ephi0 <- - signif(av*dtaudlnmu/acp, digits=3)

# Individual phi ###############################################################

phi0 <- 0
slope <- 0

# Growth rate vs. log10(first external concentration) ##############################################################

cx1    <- opt_state[,reactant[1]]

minx1 <- 0
maxx1 <- max(cx1)*1.1

plot(cx1,mu_opt,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= 1,xlim=c(minx1,maxx1),ylim=c(0,maxmu), xaxt = "n",
     yaxt = "n", ylab=bquote("Growth rate"~ mu ~ (h^-1)), xlab=paste("Concentration",reactant[1],"(g/L)"), yaxs="i", xaxs="i")
add_axes()

# Internal reactant concentrations vs.mu #############################################################################

maxci <- 1.05*max(c_opt) 

matplot(mu_opt,c_opt,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= colori,xlim=c(0,maxmu),ylim=c(0,maxci), xaxt = "n",
        yaxt = "n", xlab=bquote("Growth rate"~ mu ~ (h^-1)), ylab=bquote("Internal concentrations " ~ c^i ~ " (g/L)"), yaxs="i", xaxs="i")
add_axes()
legend("topright",i_reactant,col=colori,pch=16,cex=0.5)

# Fluxes vs. mu ######################################################################################################

maxv <- ceiling(max(v_opt)*1.1)

matplot(mu_opt,v_opt,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= colorj,xlim=c(0,maxmu),ylim=c(0,maxv), xaxt = "n",
        yaxt = "n", xlab=bquote("Growth rate"~ mu ~ (h^-1)), ylab=bquote("Fluxes " ~  v ~ " (g/L/h)"), yaxs="i", xaxs="i")
add_axes()
legend("topleft",reaction,col=colorj,pch=16,cex=0.5)


# "Growth law" plots

for (j in 1:r) {
  
  maxphi <- 1.1*(max(phi_opt[,j]))
  
  # regression on the points where phi_opt[,j] > 1e-6
  if (sum(phi_opt[,j]) > 1e-5) philm <- lm(phi_opt[phi_opt[,j]/max(phi_opt[,j]) > 1e-3,j] ~ mu_opt[phi_opt[,j]/max(phi_opt[,j]) > 1e-3] )
   
  #phi0[j] <- signif(philm$coefficients[1],3)
  
  #slope[j] <- signif(philm$coefficients[2],3)
  
  # plot #######################################################################
  
  maxphi <- 1.1*max(phi_opt[,j])
  
  #, main=paste(reaction[j],": offset =",phi0[j],", slope = ",slope[j])
  matplot(mu_opt,phi_opt[,j],pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, 
          main=reaction[j], col= colorj[j],xlim=c(0,maxmu),ylim=c(0,maxphi), xaxt = "n", yaxt = "n", 
          xlab=bquote("Growth rate"~ mu ~ (h^-1)), ylab=bquote("Proteome fraction " ~  phi ), 
          yaxs="i", xaxs="i")
  add_axes()
  #legend("top",reaction,col=colorj,pch=16,cex=0.7)
  if (sum(phi_opt[,j]) > 1e-4 & (sum(is.na(philm$coefficients)) == 0)) abline(philm)
  #abline(v=lowmu,lty=2)
  #text(0.5*lowmu,0.95*maxphi,labels=bquote(0.2*mu[max]))
  #abline(v=highmu,lty=2)
  #text(1.13*highmu,0.95*maxphi,labels=bquote(0.8*mu[max]))
  #abline(a=ephi0[j],b=0,lty=2,col=3)
  #text(1.13*highmu,1.1*ephi0[j],labels=ephi0[j])
  
}


# # Growth rate vs. log10(first external concentration) ################################################################
# 
# mu_opt <- opt_state[,"mu"]
# cx1    <- log10(opt_state[,reactant[1]])
# 
# minx1 <- floor(min(cx1))
# maxx1 <- ceiling(max(cx1)*1.1)
# maxmu <- max(mu_opt)*1.1
# 
# plot(cx1,mu_opt,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= 1,xlim=c(minx1,maxx1),ylim=c(0,maxmu), xaxt = "n",
#      yaxt = "n", ylab=bquote("Growth rate"~ mu ~ (h^-1)), xlab=paste("Concentration",reactant[1],"log10(g/L)"), yaxs="i", xaxs="i")
# axis(1, at=seq(minx1,maxx1,1), mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
# axis(1, at=seq(minx1,maxx1,1), mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
# axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
# axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
# 
# if (nx > 1) {
#   # Growth rate vs. second external concentration ####################################################################
#   
#   cx2    <- log10(opt_state[,reactant[2]])
#   
#   minx2 <- floor(min(cx2))
#   maxx2 <- ceiling(max(cx2)*1.1)
#   
#   plot(cx2,mu_opt,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= 1,ylim=c(0,maxmu), xaxt = "n",
#        yaxt = "n", ylab=bquote("Growth rate"~ mu ~ (h^-1)), xlab=paste("Concentration",reactant[2],"log10(g/L)"), yaxs="i", xaxs="i")
#   axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#   axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
#   axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#   axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
#   
# } else { 
  
  
#}


# Metabolite concentrations vs. mu ###################################################################################

# removes water if it's the model, to make other more visible

if ("H2O" %in% i_reactant ) {
  
  wp <- grep("H2O",i_reactant)
  
  cm_opt <- cm_opt[,-wp]
  
  maxcm <- ceiling(max(cm_opt))
  
  matplot(mu_opt,cm_opt,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= colori,xlim=c(0,maxmu),ylim=c(0,maxcm), xaxt = "n",
          yaxt = "n", xlab=bquote("Growth rate"~ mu ~ (h^-1)), ylab=bquote("Metabolite concentrations " ~ c^m ~ " (g/L)"), yaxs="i", xaxs="i")
  add_axes()
  legend("topleft",i_reactant[-c(wp,p)],col=colori[1:(p-2)],pch=16,cex=0.5)

} else {
  
  maxcm <- ceiling(max(cm_opt))
  
  matplot(mu_opt,cm_opt,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= colori,xlim=c(0,maxmu),ylim=c(0,maxcm), xaxt = "n",
          yaxt = "n", xlab=bquote("Growth rate"~ mu ~ (h^-1)), ylab=bquote("Metabolite concentrations " ~ c^m ~ " (g/L)"), yaxs="i", xaxs="i")
  add_axes()
  legend("topleft",i_reactant[-p],col=colori,pch=16,cex=0.5)
  
  
  
}

# Protein concentrations #############################################################################################


maxp <- ceiling(max(p_opt)*1.1)

matplot(mu_opt,p_opt,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= colorj,xlim=c(0,maxmu),ylim=c(0,maxp), xaxt = "n",
        yaxt = "n", xlab=bquote("Growth rate"~ mu ~ (h^-1)), ylab=bquote("Protein concentrations " ~  p ~ " (g/L)"), yaxs="i", xaxs="i")
add_axes()
legend("topright",reaction,col=colorj,pch=16,cex=0.5)

# Protein fractions ##################################################################################################

maxphi <- max(phi_opt)*1.1

matplot(mu_opt,phi_opt,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= colorj,xlim=c(0,maxmu),ylim=c(0,maxphi), xaxt = "n",
        yaxt = "n", xlab=bquote("Growth rate"~ mu ~ (h^-1)), ylab=bquote("Proteome fractions " ~  phi ), yaxs="i", xaxs="i")
add_axes()
#legend("top",reaction,col=colorj,pch=16,cex=0.7)


# flux fractions vs. mu ##############################################################################################

#f_opt <- v_opt/(mu_opt*rho_cond)

maxw <- ceiling(max(f_opt))

matplot(mu_opt,f_opt,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= colorj,xlim=c(0,maxmu),ylim=c(0,maxw), xaxt = "n",
        yaxt = "n", xlab=bquote("Growth rate"~ mu ~ (h^-1)), ylab=bquote("Flux fractions " ~  f), yaxs="i", xaxs="i")
add_axes()
#legend("topright",reaction,col=colorj,pch=16,cex=0.7)

# Turnover times ####################################################################################################

maxtau <- 5

matplot(mu_opt,tau_opt,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= colorj,xlim=c(0,maxmu),ylim=c(0,maxtau), xaxt = "n",
        yaxt = "n", xlab=bquote("Growth rate"~ mu ~ (h^-1)), ylab=bquote("Turnover times " ~  tau ~ (h)), yaxs="i", xaxs="i")
add_axes()
legend("topright",reaction,col=colorj,pch=16,cex=0.5)

# Turnover frequencies ##############################################################################################

k_opt <- 1/tau_opt

# deletes some k_opt if it's too high, so we can see the others
pk <- kcatf < 1000

maxf <- ceiling(max(k_opt[,pk]))


matplot(mu_opt,k_opt[,pk],pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= colorj,xlim=c(0,maxmu),ylim=c(0,maxf), xaxt = "n",
        yaxt = "n", xlab=bquote("Growth rate"~ mu ~ (h^-1)), ylab=bquote("Apparent turnover numbers " ~  k[app] ~ (h^-1)), yaxs="i", xaxs="i")
add_axes()
#legend("topleft",reaction,col=colorj,pch=16,cex=0.7)

# if the model has tRNA, mRNA in it:
if (length(grep("RNA",i_reactant)) > 0 ) {

# RNA/P ########################################################################

# the total concentrations of every reactant that has "RNA" in its name
RNA <- rowSums(c_opt[,grep("RNA",i_reactant), drop=FALSE] )

RNAP <- RNA/c_opt[,p]

RNAPr <- RNA/p_opt[,r]

RNAPlm <- lm(RNAP[mediummu] ~ mu_opt[mediummu] )

slopeRNAP <- signif(RNAPlm$coefficients[2],3)

#, main=paste(reaction[j],": offset =",phi0[j],", slope = ",slope[j])
matplot(mu_opt,RNAP,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, 
        main="RNA/P",
        col= 1,xlim=c(0,maxmu),ylim=c(0,1.1*max(RNAP)), xaxt = "n", yaxt = "n", 
        xlab=bquote("Growth rate"~ mu ~ (h^-1)), ylab="RNA/P", 
        yaxs="i", xaxs="i")
add_axes()
#legend("top",reaction,col=colorj,pch=16,cex=0.7)
abline(RNAPlm)

}

# Ribosomal RNA and protein
rrna_row <- grep("rRNA",i_reactant)
rna_row <- grep("RNA",i_reactant)
if (length(rrna_row) | (length(rna_row)==1)) {
  
  if(length(rrna_row)){
    rRNA <- c_opt[, rrna_row]
    KA_RNA <- KA[rrna_row+nx, r]
  }else{
    rRNA <- c_opt[, rna_row] 
    KA_RNA <- KA[rna_row+nx, r]
  }
  rRNArP <- rRNA/p_opt[,r]
  rRNArPlm <- lm(rRNArP[mediummu] ~ mu_opt[mediummu] )
  
  #, main=paste(reaction[j],": offset =",phi0[j],", slope = ",slope[j])
  matplot(mu_opt,rRNArP,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, 
          main="rRNA/rP",
          col= 1,xlim=c(0,maxmu),ylim=c(0,1.1*max(rRNArP)), xaxt = "n", yaxt = "n", 
          xlab=bquote("Growth rate"~ mu ~ (h^-1)), ylab="rRNA/rP", 
          yaxs="i", xaxs="i")
  add_axes()
  #legend("top",reaction,col=colorj,pch=16,cex=0.7)
  abline(rRNArPlm)
  
  
  # Ribosome "composition"
  rPR <- p_opt[,r]/(p_opt[,r] + rRNA)
  
  rPRlm <- lm(rPR[mediummu] ~ mu_opt[mediummu] )
  
  matplot(mu_opt,rPR,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, 
          main="rP/(rP + rRNA)",
          col= 1,xlim=c(0,maxmu),ylim=c(0,1.1*max(rPR)), xaxt = "n", yaxt = "n", 
          xlab=bquote("Growth rate"~ mu ~ (h^-1)), ylab="rP/(rP + rRNA)", 
          yaxs="i", xaxs="i")
  add_axes()
  #legend("top",reaction,col=colorj,pch=16,cex=0.7)
  abline(rPRlm)
  
  # Now only for the "bounded" rRNA
  brRNA <- rRNA*(rRNA/(rRNA + KA_RNA) )
  
  brPR <- p_opt[,r]/(p_opt[,r] + brRNA)
  
  brPRlm <- lm(brPR[mediummu] ~ mu_opt[mediummu] )
  
  matplot(mu_opt,brPR,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, 
          main="Protein mass fraction in ribosome",
          col= 1,xlim=c(0,maxmu),ylim=c(0,1.1*max(brPR)), xaxt = "n", yaxt = "n", 
          xlab=bquote("Growth rate"~ mu ~ (h^-1)), ylab="rP/(rP + brRNA)", 
          yaxs="i", xaxs="i")
  add_axes()
  #legend("top",reaction,col=colorj,pch=16,cex=0.7)
  abline(brPRlm)
  
  
}



# Fraction of degraded RNA #####################################################
if (length(grep("rRNAp", reaction))){
  RNAP   <- v_opt[,grep("rRNAp", reaction)]
}else{
  RNAP   <- v_opt[,grep("RNAp", reaction)]
}

if(length(grep("rRNase", reaction))){
  RNase  <- v_opt[,grep("rRNase", reaction)]*abs(M["rRNA", "rRNase"])
  deg_ratio <- RNase/RNAP
  
  maxv <- max(deg_ratio)*1.1
  
  plot(mu_opt,deg_ratio,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= "black",
       xlim=c(0,maxmu),ylim=c(0,1), xaxt = "n",
       yaxt = "n", xlab=bquote("Growth rate"~ mu ~ (h^-1)), 
       ylab="Fraction of degraded rRNA", yaxs="i", xaxs="i")
  add_axes()
}



# plot transcribed fractions of RNA ############################################
rrnap_cols <- grep("rRNAp", reaction)
mrnap_cols <- grep("mRNAp", reaction)
trnap_cols <- grep("tRNAp", reaction)

if(length(rrnap_cols) & length(mrnap_cols) & length(trnap_cols)){
  
  v_rrna <- v_opt[, rrnap_cols]
  v_mrna <- v_opt[, mrnap_cols]
  v_trna <- v_opt[, trnap_cols]
  sum_rnas <- v_rrna + v_mrna + v_trna
  rna_percent <- data.frame("rRNA"= v_rrna/sum_rnas,
                            "tRNA"= v_trna/sum_rnas,
                            "mRNA"= v_mrna/sum_rnas)
  
  matplot(mu_opt, rna_percent,
          ylim = c(0,1), pch = 19, type="b",
          xlab=bquote("Growth rate"~ mu ~ (h^-1)),
          ylab = "% transcribed RNA",
          xaxt = "n", yaxt = "n")
  legend("topleft", colnames(rna_percent), col = 1:3, pch=19)
  add_axes()
}

# Relative composition #########################################################
par(xpd=TRUE)

# order by growth rate
ordered_mu_ind <- order(mu_opt)
ordered_mu <- mu_opt[ordered_mu_ind]

water_index <- grep("H2O", i_reactant)
if (length(water_index) > 0){
  wo_water <- c_opt[,-water_index]
  reactants <- i_reactant[-water_index]
  
}else{
  wo_water <- c_opt
  reactants <- i_reactant
}

colors <- colori[1:length(reactants)]

wo_water <- wo_water[ordered_mu_ind, ]

rel_comp <- wo_water/rowSums(wo_water)
rel_comp <- rel_comp[, ncol(wo_water):1]
cum_comp <- t(apply(rel_comp, 1, cumsum))

plot(NA, xlim = c(0, max(mu_opt)), ylim=c(0,1),
     yaxs="i", xaxs="i", 
     xaxt = "n", yaxt = "n",
     xlab=bquote("Growth rate"~ mu ~ (h^-1)), 
     ylab="Relative biomass composition")
for(i in ncol(rel_comp):1){
  col <- colors[i]
  polygon(c(ordered_mu, rev(ordered_mu)), c(rep(0, nrow(cum_comp)), rev(cum_comp[,i])),
          col=col, border=col)
  par(new=TRUE)
}
par(new=FALSE)
add_axes()
legend('bottomleft', reactants, col = rev(colors), pch = 15, cex=0.5,  bg="white")


# Macro composition
#rel_comp <- wo_water/rowSums(wo_water)
prot <- ncol(wo_water)
lip <- which(reactants == "LIP")
RNA <- c(grep("RNA", reactants), grep("TC", reactants))
DNA <- grep("DNA", reactants)

prot_sum <- wo_water[, prot, drop = FALSE]
lip_sum <- wo_water[, lip, drop = FALSE]
RNA_sum <- rowSums(wo_water[, unique(RNA), drop = FALSE])
DNA_sum <- wo_water[, DNA, drop = FALSE]
rest_sum <- rowSums(wo_water[, -c(DNA, RNA, lip, prot), drop = FALSE])

sectors <- cbind(rest_sum, lip_sum, DNA_sum, RNA_sum, prot_sum)
sectors <-  t(apply(sectors[,ncol(sectors):1], 1, cumsum))

if(ncol(lip_sum) == 1){
  sector_names <- c("Rest", "Lipids", "DNA", "RNA", "Protein")
}else{
  sector_names <- c("Rest", "DNA", "RNA", "Protein")
}
colors <- brewer.pal(ncol(sectors), "Paired")

# plot(NA, xlim = c(0, max(ordered_mu)), ylim=c(0, max(sectors)),
#      yaxs="i", xaxs="i",
#      xaxt = "n", yaxt = "n",
#      xlab=bquote("Growth rate"~ mu ~ (h^-1)), 
#      ylab="Macro composition")
# for(i in ncol(sectors):1){
#   col <- colors[i]
#   polygon(c(ordered_mu, rev(ordered_mu)), c(rep(0, nrow(sectors)), rev(sectors[,i])),
#           col=col, border=col)
#   par(new=TRUE)
# }
# par(new=FALSE)
# add_axes()
# legend('bottomleft', sector_names, col = rev(colors), pch = 15, cex=0.5)


rel_comp <- sectors/sectors[, ncol(sectors)]
plot(NA,
     xlim=c(0, 1), #max(ordered_mu)),
     ylim=c(0, 1),#max(rel_comp)),
     yaxs="i", xaxs="i",
     xaxt = "n", yaxt = "n",
     xlab=bquote("Growth rate"~ mu ~ (h^-1)), 
     ylab="Predicted biomass")
for(i in ncol(rel_comp):1){
  col <- colors[i]
  polygon(c(ordered_mu, rev(ordered_mu)), c(rep(0, nrow(rel_comp)), rev(rel_comp[,i])),
          col=col, border=col)
  par(new=TRUE)
}
par(new=FALSE)
add_axes()
legend('bottomleft', sector_names, col = rev(colors), pch = 15, cex=0.8, bg="white")

# source("../validation_plots.R")
# setwd("../code/Results GBA")


# Proteome composition #########################################################
cum_phis <- t(apply(phi_opt[,ncol(phi_opt):1], 1, cumsum))

cum_phis <- cum_phis[ordered_mu_ind, ]
plot(NA, xlim = c(0, 1), ylim=c(0, 1),
     yaxs="i", xaxs="i",
     xaxt = "n", yaxt = "n",
     xlab=bquote("Growth rate"~ mu ~ (h^-1)), 
     ylab="Proteome composition")
for(i in ncol(cum_phis):1){
  col <- rev(colorj)[i]
  polygon(c(ordered_mu, rev(ordered_mu)), c(rep(0, nrow(cum_phis)), rev(cum_phis[,i])),
          col=col, border=col)
  par(new=TRUE)
}
par(new=FALSE)
add_axes()
legend('bottomleft', colnames(phi_opt), col = colorj, pch = 15, cex=0.5)



# Proteome sectors #############################################################

sector_names <- c( "Other", "Metabolism", "Translation and trascription")

phi_opt_ordered <- phi_opt[ordered_mu_ind,]

# extract sectors
maint_index <- which(colnames(phi_opt_ordered) == "p Maint")
other <- phi_opt_ordered[, maint_index]
trans_indices <- c(grep("RNA", colnames(phi_opt_ordered)),
                 grep("DNA", colnames(phi_opt_ordered)),
                 grep("Deg", colnames(phi_opt_ordered)),
                 which(colnames(phi_opt_ordered) == "p r"))
translation_transcription <- rowSums(phi_opt_ordered[, unique(trans_indices)])
metabolism <- rowSums(phi_opt_ordered[, -c(trans_indices, maint_index)])

sectors <- cbind(other, metabolism, translation_transcription)
sectors <-  t(apply(sectors[,ncol(sectors):1], 1, cumsum))
  
colors <- brewer.pal(ncol(sectors), "Paired")

plot(NA, xlim = c(0, 1), ylim=c(0, 1),
     yaxs="i", xaxs="i",
     xaxt = "n", yaxt = "n",
     xlab=bquote("Growth rate"~ mu ~ (h^-1)), 
     ylab="Proteome sectors")
for(i in ncol(sectors):1){
  col <- colors[i]
  polygon(c(ordered_mu, rev(ordered_mu)), c(rep(0, nrow(sectors)), rev(sectors[,i])),
          col=col, border=col)
  par(new=TRUE)
}
par(new=FALSE)
add_axes()
legend('bottomleft', sector_names, col = rev(colors), pch = 15, cex=0.5)

par(xpd=FALSE)


#### finer sectors #############################################################

# # extract sectors
# r_indices <- c(which(colnames(phi_opt) == "p r"),
#               which(colnames(phi_opt) == "p tRNAc"))
# r_sector <- rowSums(phi_opt[, unique(r_indices)])
# 
# a_indices <- c(which(colnames(phi_opt) == "p EAA"),
#               which(colnames(phi_opt) == "p ENT"))
# a_sector <- rowSums(phi_opt[, unique(a_indices)])
# 
# c_indices <- c(which(colnames(phi_opt) == "p tC"),
#                which(colnames(phi_opt) == "p tW"),
#                which(colnames(phi_opt) == "p CCM"),
#                which(colnames(phi_opt) == "p FERM"),
#                which(colnames(phi_opt) == "p ATPS"))
# c_sector <- rowSums(phi_opt[, unique(c_indices)])
# 
# sectors <- cbind(c_sector, a_sector, r_sector)
# sector_names <- c("Catabolism (CS-sector)", "Anabolism (UA-sector)", "Translation (R-sector)")
# colors <- brewer.pal(length(sector_names), "Paired")
# 
# matplot(mu_opt,sectors,pch=16, mgp=c(2, 0.5, 0), col = colors, frame.plot=F, 
#         main="Proteome sectors",
#         xlim=c(0,maxmu),ylim=c(0,0.5), 
#         xaxt = "n", yaxt = "n", 
#         xlab=bquote("Growth rate"~ mu ~ (h^-1)), ylab="Fraction of proteome", 
#         yaxs="i", xaxs="i")
# add_axes()
# legend("topleft",sector_names,col=colors,pch=16,cex=0.7)


# Plots Model ##################################################################

par(mfrow=c(1,1))

size <- 5*r^(-0.8)

plot.new()
text(0.2, 0.9, paste("keep_ribosome_kcat", keep_ribosome_kcat, "\n",
                     "keep_transport_kcat",  keep_transport_kcat, "\n",
                     "maintenance_fun", maintenance_fun, "\n"))

textplot(M,cex=size)
title("M")

textplot(K,cex=size)
title("K")

textplot(KA,cex=size)
title("KA")

textplot(rbind(kcatf,kcatb),cex=size)
title("kcat")

# Calculates the corresponding equilibrium constants Keq
Keq <- 0
for (j in 1:r) Keq[j] <- kcatf[j]*prod(KP[KP[,j] < Inf,j])/(kcatb[j]*prod(KS[KS[,j] < Inf,j]))

textplot(Keq,cex=size)
title("Keq")

if(predict.parameters > 0){
  textplot(phi_data,cex=size)
  title("phi input")
  
  textplot(esat,cex=2)
  title("average saturation input")
}

textplot(min_phi,cex=size)
title("minimal phi constraint")

textplot(min_f,cex=size)
title("minimal f constraint")

dev.off()

setwd(directory)

