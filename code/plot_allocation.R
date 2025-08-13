library(RColorBrewer)
library(scales)
library(here)

setwd(paste0(here(), "/code/Results GBA"))

opt_data <- read.csv("GBA Model A11, mean time (25.5s) results.csv", row.names = 1)
reactants <- colnames(opt_data)[5:17]
p_opt <- opt_data[,grep("p.",colnames(opt_data))]
phi_opt <- p_opt/opt_data[,"p"]

plot_and_fit <- function(data, growth_rates, ylim, xlim, 
                         xlab=bquote("Growth rate"~ mu ~ (h^-1)),
                         ylab="Proteome fraction", legend=TRUE){
  if(length(data) > 8){
    print("Too many groups")
    return()
  }
  
  colors <- brewer.pal(8, "Dark2")
  for(i in seq_along(data)){
    df <- data.frame(x=growth_rates, y=data[[names(data)[i]]])
    plot(y~x, data=df,
         col = colors[i],
         pch=20,
         ylim=ylim, 
         xlim=xlim,
         ylab=ylab,
         xlab=xlab,
         cex=1.5,
         cex.lab = 1.3)
    
    fit <- lm(y~x, data=df)
    
    x_seq <- seq(-0.1,1.3,0.01)
    pred <- predict(fit, newdata = data.frame(x = x_seq), interval = "confidence")
    
    if(length(growth_rates)>30){
      print(paste(names(data)[i], "phi at mu=1:", round(pred[which(x_seq==1)],3)))
    }
    
    # Add shaded confidence interval
    polygon(c(x_seq, rev(x_seq)), 
            c(pred[, "lwr"], rev(pred[, "upr"])), 
            col = alpha(colors[i],0.3), border = NA)  # Transparent shading
    abline(fit, col = colors[i])
    
    par(new=TRUE)
    
  }
  if(legend){
    legend("topleft", 
           bty = "n",
           legend = names(data),
           col = colors[1:length(data)], pch=20, cex=1.1)
  }
  par(new=FALSE)
}


#### Make plots ################################################################
legend <- TRUE
xlim <- c(0,1.15)
#+ fig.height = 6.6, fig.width = 12
par(mfrow=c(2, 4))


sectors <- list("Translation (w/o tRNA synthetases)" = phi_opt$p.r,
                "tRNA synthetases" = phi_opt$p.tRNAc)
plot_and_fit(sectors, 
             opt_data$mu, 
             ylim = c(0,0.45), 
             xlim = xlim,
             legend = legend)
title("Predicted")

transcription <- rowSums(phi_opt[, c("p.rRNAp", "p.mRNAp", "p.tRNAp",
                                     "p.rRNase", "p.mRNase", "p.tRNAse")])
sectors <- list("Transcription" = transcription,
                "DNA replication" = phi_opt$p.DNAp)

plot_and_fit(sectors, 
             opt_data$mu, 
             ylim = c(0,0.1), 
             xlim = xlim,
             legend = legend)
title("Predicted")

sectors <- list("Central Carbon metabolism" = phi_opt$p.tC,
                "ETC" = phi_opt$p.ATPS)
                #"All carbohydrate/energy metabolism" = phi_opt$p.tC+phi_opt$p.ATPS)
plot_and_fit(sectors, 
             opt_data$mu, 
             ylim = c(0,0.5), 
             xlim = xlim,
             legend = legend)
title("Predicted")


sectors <- list("Nucleotide metabolism" = phi_opt$p.ENT + phi_opt$p.ADPS,
                "Amino acid metabolism" = phi_opt$p.EAA,
                "Lipid metabolism" = phi_opt$p.LIPS)

plot_and_fit(sectors, 
             opt_data$mu, 
             ylim = c(0,0.2), 
             xlim = xlim,
             legend = legend)
title("Predicted")

# Proteome composition #########################################################
sector_names <- c( "Other", "Metabolism", "Info. processing")


ordered_mu_ind <- order(opt_data$mu)
ordered_mu <- opt_data$mu[ordered_mu_ind]
phi_opt_ordered <- phi_opt[ordered_mu_ind,]

maint_index <- which(colnames(phi_opt_ordered) == "p.Maint")
other <- phi_opt_ordered[, maint_index]
trans_indices <- c(grep("RNA", colnames(phi_opt_ordered)),
                   grep("DNA", colnames(phi_opt_ordered)),
                   grep("RNase", colnames(phi_opt_ordered)),
                   which(colnames(phi_opt_ordered) == "p.r"))
translation_transcription <- rowSums(phi_opt_ordered[, unique(trans_indices)])
metabolism <- rowSums(phi_opt_ordered[, -c(trans_indices, maint_index)])

sectors <- cbind(other, metabolism, translation_transcription)
cum_comp <-  t(apply(sectors[,ncol(sectors):1], 1, cumsum))

colors <- brewer.pal(ncol(sectors), "Paired")

plot(NA, xlim = xlim, ylim=c(0,1),
     yaxs="i", xaxs="i", 
     # xaxt = "n", yaxt = "n",
     xlab=bquote("Growth rate"~ mu ~ (h^-1)), 
     ylab="Relative proteome composition",
     cex.lab = 1.3)
for(i in ncol(cum_comp):1){
  col <- colors[i]
  polygon(c(ordered_mu, rev(ordered_mu)), c(rep(0, nrow(cum_comp)), rev(cum_comp[,i])),
          col=col, border=col)
  par(new=TRUE)
}
par(new=FALSE)
if(legend){
  legend('bottomright', sector_names, #bty = "n",
         col = rev(colors), pch = 15, cex=1.1, border = NA)
}
title("Predicted")




