library(here)
library(RColorBrewer)

directory <- paste0(here(), "/data")
setwd(directory) 

gausing <- read.csv("gausing_RNA_deg.csv")
biomass <- read.csv("biomass.csv")
bremer <- read.csv("bremer_percent_transcription.csv")
 
colnames(biomass)[6] <- 'Rest'
biomass[,2:6] <- biomass[,2:6]/100
size <- 1.3


coloraxis <- "black"
add_axes <- function(){
  axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, 
       tck=0.02,las=1, cex.axis=1.5, padj=0.5)
  axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, 
       tck=-0.02,las=1, cex.axis=1.5, padj=0.5)
  axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, 
       tck=0.02,las=1, cex.axis=1.5, hadj=1.1)
  axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, 
       tck=-0.02,las=1, cex.axis=1.5, hadj=1.1)
}

#+ fig.height = 3, fig.width = 12
par(mar=c(6,6,2,2), mfrow=c(1,3))
plot(percent_degraded/100 ~ growth_rate, data=gausing,
     ylim = c(0,1), xlim=c(0,1.4), pch = 19,
     xaxt = "n",
     yaxt = "n", 
     xlab=bquote("Growth rate"~ mu ~ (h^-1)), 
     ylab="Fraction of degraded rRNA",
     cex.lab = 1.8,
     yaxs="i", xaxs="i")
add_axes()


matplot(bremer$mu, bremer[,2:4],
     ylim = c(0,100), pch = 19, type="b",
     xlab=bquote("Growth rate"~ mu ~ (h^-1)),
     ylab = "% transcribed RNA",
     cex.lab = size, cex.axis = size, cex=size)
legend("topleft", colnames(bremer)[2:4], col = 1:3, pch=19)



# Absolute composition #########################################################
par(new=FALSE, xpd=FALSE)
#coloraxis <- "darkgrey"
#par(mar=c(6,6,2,2), mfrow = c(1,1))

# add_axes <- function(){
#   axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1, cex.axis=2)
#   axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1, cex.axis=2)
#   axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1, cex.axis=2)
#   axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1, cex.axis=2)
# }

colors <- rev(brewer.pal(5, "Paired"))
cum_comp <- t(apply(biomass[,2:ncol(biomass)], 1, cumsum))
plot(NA, xlim = c(0, max(biomass$mu)), ylim=c(0, 1),
     yaxs="i", xaxs="i",
     xaxt = "n", yaxt = "n",
     #cex.lab = size*1.5, cex.axis = size*1.3,
     xlab=bquote("Growth rate"~ mu ~ (h^-1)), 
     ylab="Experimental biomass")
for(i in ncol(cum_comp):1){
  col <- rev(colors)[i]
  polygon(c(biomass$mu, rev(biomass$mu)), c(rep(0, nrow(cum_comp)), rev(cum_comp[,i])),
          col=col, border=col)
  par(new=TRUE)
}
add_axes()
#legend('bottomleft', colnames(biomass)[ncol(biomass):2], col = colors, pch = 15, cex=size)


par(xpd=TRUE)
for(mu in biomass$mu){
  arrows(x0=mu, y0=1.05, x1=mu, y1=1.00, length=0.05, code=2, lwd=1.5)
}
par(xpd=FALSE)

