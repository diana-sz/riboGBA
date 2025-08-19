library(RColorBrewer)
library(scales)
library(here)

setwd(paste0(here(), "/code/Results GBA"))

modelname <- "A9ferm4i"# "A10ferm_decoupled"
filename <- "GBA Model A10ferm_decoupledi mean time (7.85s) results.csv"
filename <- "GBA Model A9ferm4i mean time (8.14s) results.csv"

opt_data <- read.csv(filename, row.names = 1)
reactants <- colnames(opt_data)[5:13]
p_opt <- opt_data[,grep("p.",colnames(opt_data))]
f_opt <- opt_data[,grep("f.",colnames(opt_data))]

phi_opt <- p_opt/opt_data[,"p"]

plot_and_fit <- function(data, growth_rates, ylim, xlim, 
                         xlab=bquote("Growth rate"~ mu ~ (h^-1)),
                         ylab="Proteome fraction", legend=TRUE, fit=TRUE){
  if(length(data) > 8){
    print("Too many groups")
    return()
  }
  
  colors <- brewer.pal(8, "Set1")
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
    
    if(fit){
      fitted <- lm(y~x, data=df)
      
      x_seq <- seq(-0.1,1.3,0.01)
      pred <- predict(fitted, newdata = data.frame(x = x_seq), interval = "confidence")
      
      
      # Add shaded confidence interval
      polygon(c(x_seq, rev(x_seq)), 
              c(pred[, "lwr"], rev(pred[, "upr"])), 
              col = alpha(colors[i],0.3), border = NA)  # Transparent shading
      abline(fitted, col = colors[i])
    }
    
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
par(mfrow=c(1, 2))


sectors <- list("Fermentation" = phi_opt$p.FERM,
                "Respiration" = phi_opt$p.RESP)
                #"All energy metabolism" = phi_opt$p.FERM+phi_opt$p.RESP)
plot_and_fit(sectors, 
             opt_data$mu, 
             ylim = c(0,0.07), 
             xlim = xlim,
             legend = legend,
             fit=FALSE)
title("Proteome fractions")



png(paste0("../../figures/", modelname, "_ferm_resp_fluxes.png"), 
    type="cairo", units="cm", pointsize=10,
    width=10, height=10, res=300)
sectors <- list("Fermentation" = opt_data$f.FERM,
                "Respiration" = opt_data$f.RESP)
#"All energy metabolism" = phi_opt$p.FERM+phi_opt$p.RESP)
plot_and_fit(sectors, 
             opt_data$mu, 
             ylim = c(0,7), 
             xlim = xlim,
             legend = legend,
             fit=FALSE)
title("Flux fractions")
dev.off()

# sectors <- list("Translation" = phi_opt$p.r,
#                 "Transcription" = phi_opt$p.RNAp)
# plot_and_fit(sectors, 
#              opt_data$mu, 
#              ylim = c(0,0.45), 
#              xlim = xlim,
#              legend = legend)
# title("Gene expression")


# sectors <- list("Nucleotide metabolism" = phi_opt$p.ENT + phi_opt$p.ADPS,
#             "Amino acid metabolism" = phi_opt$p.EAA
#             #"Lipid metabolism" = phi_opt$p.LIPS
#             )
# 
# plot_and_fit(sectors, 
#          opt_data$mu, 
#          ylim = c(0,0.25), 
#          xlim = xlim,
#          legend = legend)
# title("Predicted")

# Proteome composition #########################################################
# sector_names <- c( "Other", "Metabolism", "Info. processing")
# 
# 
# ordered_mu_ind <- order(opt_data$mu)
# ordered_mu <- opt_data$mu[ordered_mu_ind]
# phi_opt_ordered <- phi_opt[ordered_mu_ind,]
# 
# trans_indices <- c(grep("RNA", colnames(phi_opt_ordered)),
#                    grep("DNA", colnames(phi_opt_ordered)),
#                    grep("RNase", colnames(phi_opt_ordered)),
#                    which(colnames(phi_opt_ordered) == "p.r"))
# translation_transcription <- rowSums(phi_opt_ordered[, unique(trans_indices)])
# metabolism <- rowSums(phi_opt_ordered[, -c(trans_indices)])
# 
# sectors <- cbind(metabolism, translation_transcription)
# cum_comp <-  t(apply(sectors[,ncol(sectors):1], 1, cumsum))
# 
# colors <- brewer.pal(ncol(sectors), "Paired")
# 
# plot(NA, xlim = xlim, ylim=c(0,1),
#      yaxs="i", xaxs="i", 
#      # xaxt = "n", yaxt = "n",
#      xlab=bquote("Growth rate"~ mu ~ (h^-1)), 
#      ylab="Relative proteome composition",
#      cex.lab = 1.3)
# for(i in ncol(cum_comp):1){
#   col <- colors[i]
#   polygon(c(ordered_mu, rev(ordered_mu)), c(rep(0, nrow(cum_comp)), rev(cum_comp[,i])),
#           col=col, border=col)
#   par(new=TRUE)
# }
# par(new=FALSE)
# if(legend){
#   legend('bottomright', sector_names, #bty = "n",
#          col = rev(colors), pch = 15, cex=1.1, border = NA)
# }
# title("Predicted")
# 
# 
# 
# 
