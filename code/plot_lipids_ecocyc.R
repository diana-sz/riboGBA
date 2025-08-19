library(RColorBrewer)
library(scales)
library(here)

setwd(here())

source("code/process_ecocyc_annotation2.R")
source("code/process_proteomics.R")

scaled <- FALSE

all_row_names <- unique(c(rownames(mori_proteomics), rownames(mori_proteomics2), rownames(schmidt_proteomics)))

# Align each data frame to the combined row names, filling missing rows with NA
df1_aligned <- data.frame(row.names = all_row_names, mori_proteomics[all_row_names, , drop = FALSE])
df2_aligned <- data.frame(row.names = all_row_names, mori_proteomics2[all_row_names, , drop = FALSE])
df3_aligned <- data.frame(row.names = all_row_names, schmidt_proteomics[all_row_names, , drop = FALSE])

# Combine the data frames
merged <- cbind(df1_aligned, df2_aligned, df3_aligned)

proteomics_datasets <- list(Mori_titrated_glucose=mori_proteomics,
                            Mori_c_sources=mori_proteomics2,
                            Schmidt=schmidt_proteomics,
                            Merged=merged)


if(scaled){
  scaled_1 <- mori_growth_rates/max(mori_growth_rates)
  scaled_2 <- mori_growth_rates2/max(mori_growth_rates2)
  scaled_3 <- schmidt_growth_rates/max(schmidt_growth_rates)
  all_growth_rates <- c(scaled_1, scaled_2, scaled_3)
  growth_datasets <- list(Mori_titrated_glucose=scaled_1,
                          Mori_c_sources=scaled_2,
                          Schmidt=scaled_3,
                          Merged=all_growth_rates)
}else{
  all_growth_rates <- c(mori_growth_rates, mori_growth_rates2, schmidt_growth_rates)
  growth_datasets <- list(Mori_titrated_glucose=mori_growth_rates,
                          Mori_c_sources=mori_growth_rates2,
                          Schmidt=schmidt_growth_rates,
                          Merged=all_growth_rates)
}

#### Function definitions ######################################################
get_sector <- function(target, proteomics, gene_annot){
  sector_sum <- data.frame()
  
  target_genes <- c()
  for (gene in names(gene_annot)){
    if (target %in% gene_annot[[gene]]){
      if("amino acid biosynthetic process" %in% gene_annot[[gene]]){
        next
      }
      if("protein maturation" %in% gene_annot[[gene]]){
        next
      }
      if(length(gene_annot[[gene]]) > 1){
        print(gene_annot[[gene]])
      }
      target_genes <- c(target_genes, gene)

    }
  }
  if(is.null(target_genes)){
    print(paste("no genes in", target, "group"))
    next
  }
  target_rows <- proteomics[target_genes, ]

  return(target_rows)
}



plot_and_fit <- function(ydata, growth_rates, ylim, xlim, 
                         xlab=bquote("Growth rate"~ mu ~ (h^-1)),
                         ylab="Proteome fraction", legend=TRUE){

  colors <- brewer.pal(8, "Dark2")
  df <- data.frame(x=growth_rates, y=unlist(ydata))
  plot(y~x, data=df,
       col = colors[1],
       pch=20,
       ylim=ylim, 
       xlim=xlim,
       ylab=ylab,
       xlab=xlab,
       cex=1.5,
       cex.lab = 1.3,
       main=rownames(data))
    
    fit <- lm(y~x, data=df)
    
    x_seq <- seq(-0.1,1.3,0.01)
    pred <- predict(fit, newdata = data.frame(x = x_seq), interval = "confidence")
    
    # Add shaded confidence interval
    polygon(c(x_seq, rev(x_seq)), 
            c(pred[, "lwr"], rev(pred[, "upr"])), 
            col = alpha(colors[1],0.3), border = NA)  # Transparent shading
    abline(fit, col = colors[1])
}


#### Make plots ################################################################
#+ fig.height = 29.7, fig.width = 12
par(mfrow=c(4,3), mar = c(4,4,1,1))
xlim <- c(0, 1.1)

for(dataset in names(proteomics_datasets)){
  legend <- TRUE
  if(dataset != "Merged"){legend <- FALSE}
  current_dataset <- proteomics_datasets[[dataset]]
  current_growth <- growth_datasets[[dataset]]
  lipids <- get_sector("lipid biosynthetic process", current_dataset, gene_annot) 
  summed <- colSums(lipids, na.rm = TRUE)
  fractions <- summed/colSums(current_dataset, na.rm=TRUE)
  plot_and_fit(fractions, current_growth, c(0,0.04), xlim)
  title(paste(dataset, "Lipids"))
  
  transporters <- current_dataset[transporter_genes, ]
  summed_transporters <- colSums(transporters, na.rm = TRUE)
  fractions_transp <- summed_transporters/colSums(current_dataset, na.rm=TRUE)
  plot_and_fit(fractions_transp, current_growth, c(0,0.18), xlim)
  title(paste(dataset, "Transporters"))
  
  
  plot_and_fit(fractions_transp/fractions, current_growth, c(0,10), xlim)
  title(paste(dataset,  "Ratio"))
  
  
  
  # for(row in 1:nrow(lipids)){
  #   one_protein <- lipids[row, ]
  #   max_frac <- max(one_protein, na.rm=TRUE)
  #   if(max_frac < 0.0001){
  #     next
  #   }
  #   if(all(is.na(one_protein))){
  #     next
  #   }
  #   ylim <- c(0, max_frac*1.2)
  #   plot_and_fit(one_protein, current_growth, ylim, xlim)
  #   title(paste(dataset, rownames(one_protein)))
  #}
}

biomass <- read.csv("data/biomass.csv")

prot_corr <- lm(Protein ~ mu, data=biomass)
pred_prot <- predict(prot_corr, newdata = data.frame(mu=current_growth))/100

transporter_mass_frac <- fractions_transp*pred_prot
plot_and_fit(transporter_mass_frac, current_growth, xlim=xlim, ylim=c(0,0.05), ylab = "Transporter mass fraction")


exp_model <- nls(Lipids ~ a * exp(b * mu), start = list(a = 1, b = -0.1), data=biomass)
a_fit <- coef(exp_model)["a"]
b_fit <- coef(exp_model)["b"]

# Generate predicted values
pred_lip <- a_fit * exp(b_fit * current_growth)

# Plot
plot(Lipids~mu, data=biomass, main = "Lipid mass fraction", ylim = c(0,20))
points(current_growth, pred_lip, col = "red", lwd = 1)


plot_and_fit(transporter_mass_frac/pred_lip, current_growth, xlim=xlim, ylim=c(0,0.004), ylab = "Transporters/lipids")



metadata <- read.csv("data/Schmidt_proteomics_table_s23.csv", skip=2, nrows = 26)
metadata <- metadata[metadata$Strain == "BW25113",]
rownames(metadata) <- metadata[,1]
volumes <- metadata$Single.cell.volume..fl.1
growth_rates <- metadata$Growth.rate..h.1.
plot(volumes~growth_rates)

vol_corr <- lm(volumes~growth_rates)
pred_volume <- predict(vol_corr, newdata = data.frame(growth_rates=current_growth))


plot_and_fit((transporter_mass_frac*pred_volume)/pred_volume^(2/3), current_growth, xlim=xlim, ylim=c(0,0.04), ylab = "Transporters/lipids")

sav <- -2.8*current_growth+9.3
plot_and_fit(transporter_mass_frac/sav, current_growth, xlim=xlim, ylim=c(0,0.004), ylab = "Transporters/lipids")
