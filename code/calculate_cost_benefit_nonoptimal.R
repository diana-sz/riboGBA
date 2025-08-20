library(here)
library(readODS)
library(rstudioapi)
library(xtable)
library(RColorBrewer)

directory <- paste0(here(), "/code")
setwd(directory)


modelname="A7simple"
suppressMessages(source("Readmodelods_v2.R"))

source("GBA_Kinetics.R")

opt_data <- read.csv(paste0("../data/", modelname, "0_protein_cost.csv"), row.names = 1)

row <- 1
rho <- rho_cond[1]
x  <- x_cond[,1]


results <- data.frame(
  growth_rate = numeric(),
  tested_protein = character(),
  phi = numeric(),
  variable = character(),
  protein_benefit = numeric(),
  local_cost = numeric(),
  local_benefit = numeric(),
  transport_benefit = numeric(),
  sum = numeric(),
  stringsAsFactors = FALSE
)

for (tested_protein in unique(opt_data$reaction)) {
  one_prot <- opt_data[opt_data$reaction == tested_protein,]
  for (phi in one_prot$phi){
    one_phi <- one_prot[one_prot$phi == phi,]
    
    vs <- one_phi[row, grep("v\\.", colnames(one_phi))]
    fs <- one_phi[row, grep("f\\.", colnames(one_phi))]
    
    cint <- one_phi[row, grep("c\\.", colnames(one_phi))]
    cint <- cint[, (1+length(x)):length(cint)]
    fint <- cint/rho
    taus <- tau(t(cint))
    dtaus <- dtau(t(cint))
    growth_rate <- one_phi[row, "mu"]
    
    for(j in 1:length(vs)){
      Mjp <- M["p",j]
      local_cost <- growth_rate*taus[j]
      local_benefit <- unlist(vs) %*% (dtaus %*% M[, j])
      transport_benefit <- colSums(M)[j] * unlist(vs) %*% dtaus %*% unlist(fint)

      # print(paste(colnames(vs)[j],
      #             "prot", Mjp,
      #             "local cost", round(local_cost, 4),
      #             "local benefit", round(local_benefit, 4),
      #             "transport benefit", round(transport_benefit, 4),
      #             "sum", round(Mjp  - local_cost - local_benefit + transport_benefit, 4)))  
      
      results <- rbind(
        results,
        data.frame(
          growth_rate = growth_rate,
          tested_protein = tested_protein,
          phi = phi,
          variable = colnames(vs)[j],
          protein_benefit = as.numeric(Mjp),
          local_cost = -as.numeric(local_cost),
          local_benefit = -as.numeric(local_benefit),
          transport_benefit = as.numeric(transport_benefit),
          sum = as.numeric(Mjp - local_cost - local_benefit + transport_benefit)
        )
      )
    }
  }
}
  

xlim <- c(0,0.5)
ylim <- c(-1.5,1.5)
plots <- c("local_cost", "local_benefit", "transport_benefit", "sum")
colors <- brewer.pal(length(plots), "Set1")

n_prot <- length(unique(results$tested_protein))
n_rxn  <- length(unique(results$variable))

png(paste0("../figures/", modelname, "_cost_benefit.png"), 
    type="cairo", units="cm",
    width=22, height=20, res=300)

par(mfrow = c(n_rxn, n_prot), mar = c(0,0,0.5,0.5), oma = c(5,5.5,1,1), xpd = FALSE)

rxn_idx <- 0
for(reaction in rev(unique(results$variable))){
  rxn_idx <- rxn_idx + 1
  one_rxn <- results[results$variable == reaction, ]
  
  prot_idx <- 0
  for(prot in unique(one_rxn$tested_protein)){
    prot_idx <- prot_idx + 1
    one_prot <- one_rxn[one_rxn$tested_protein == prot, ]
    
    # plot background
    plot(NA, xlim = xlim, ylim = ylim,
         axes = FALSE, xlab = "", ylab = "")
    abline(h=0)
    abline(v=one_prot$phi[which.max(one_prot$growth_rate)], col="grey70", lty=2)
    box()
    
    # add all curves
    for(plot_idx in seq_along(plots)){
      points(one_prot$phi, one_prot[[plots[plot_idx]]],
             pch = 19, col = colors[plot_idx], cex=0.5)
    }
    
    # only y-axis on first column
    if(prot_idx == 1){
      axis(2, at = c(-1,0,1))
      mtext(gsub("v\\.", "", reaction), side = 2, line = 2.5, cex = 0.8)
    }
    
    # only x-axis on last row
    if(rxn_idx == n_prot){
      axis(1)
      mtext(prot,
            side = 1, line = 2.5, cex = 0.8)
    }
    
  }
}

# Add one common x and y label
mtext("Proteome fraction of tested protein", side = 1, line = 4, outer = TRUE, cex = 1.1)
mtext("Cost / Benefit", side = 2, line = 4, outer = TRUE, cex = 1.1)

# add legend once in outer margin
par(xpd = NA)  # allow legend outside plot region
legend("bottomleft", inset = c(-0.05,-0.05),
       legend = plots, col = colors, pch = 19, bty = "n")

dev.off()
