library(here)
library(RColorBrewer)


modelname <- "A7simple1"
data <- read.csv(paste0("../data/", modelname, "_protein_cost_test.csv"), row.names = 1)



plot_composition <- function(proteome, target_phi, colors = NULL,
                             main = "Proteome Composition",
                             ylab = "Composition",
                             xlab = "Proteome fraction",
                             legend=TRUE) {
  # Dimensions
  n <- nrow(proteome)
  m <- ncol(proteome)

  # Cumulative sum by row (used for stacking)
  proteome <- proteome/rowSums(proteome)
  y_cum <- t(apply(proteome[ncol(proteome):1], 1, cumsum))

  # Base plot (empty)
  plot(NA, xlim = c(0, max(target_phi)), ylim = c(0, 1), xlab = xlab, ylab = ylab, main = main)
  
  # Bottom line (start from 0)
  y_prev <- rep(0, length(target_phi))
  
  # Draw polygons for each protein group
  for (i in m:1) {
    y_top <- y_cum[, i]
    polygon(
      c(target_phi, rev(target_phi)),
      c(y_prev, rev(y_top)),
      col = colors[i],
      border = NA
    )
  }
  
  # Add legend
  legend("topright", legend = rev(colnames(proteome)), fill = colors, bty = "n", cex = 0.8)
}



for(protein in unique(data$protein)){
  
  png(paste0("../figures/", modelname, "_", protein, ".png"), 
      type="cairo", units="cm",
      width=22, height=5, res=300)
  par(mfcol=c(1,4), mar = c(4,4,2,0.5))
  
  one_prot <- data[data$protein == protein, ]
  target_phi <- one_prot$phi
  proteome <- one_prot[, grep("p\\.", colnames(one_prot))]
  biomass <- one_prot[, 8:which(colnames(one_prot)=="p")]
  
  mu <- one_prot$mu
  plot(one_prot$phi, mu, xlim = c(0,  max(target_phi)), ylim = c(0,max(mu)),
       pch=one_prot$shape,
       main = protein,
       xlab = "Proteome fraction", ylab="Growth rate")
  
  plot(one_prot$phi/one_prot$phi[which.max(one_prot$mu)], mu, 
       xlim = c(0,  4), ylim = c(0,max(mu)),
       pch=one_prot$shape,
       xlab = "Expression relative to optimum", ylab="Growth rate")
  
  colors <- brewer.pal(ncol(proteome), "Paired")
  plot_composition(proteome, target_phi, colors, main="", ylab = "Proteome composition")
  plot_composition(biomass, target_phi, colors, main="", ylab = "Biomass composition")

  
  dev.off()
  
  
}
