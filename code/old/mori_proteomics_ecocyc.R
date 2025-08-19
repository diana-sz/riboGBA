library(RColorBrewer)
library(here)

directory <- paste0(here(), "/code")
setwd(directory)

metadata <- read.csv("../data/EV3-Samples-2.csv")
c_lim_samples <- metadata[8:22,]

proteomics <- read.csv("../data/EV9-AbsoluteMassFractions-2.csv")
proteomics <- proteomics[proteomics$Gene.locus != "", ]
rownames(proteomics) <- proteomics$Gene.locus
# drops <- c("Gene.name", "Gene.locus", "Protein.ID")
# proteomics <- proteomics[ , !(names(proteomics) %in% drops)]
proteomics <- proteomics[, c_lim_samples$Sample.ID]

source("process_ecocyc_annotation.R")


get_sector_fraction <- function(targets, proteomics, gene_annot){
  sector_sum <- data.frame()
  
  for (target in targets){
    target_genes <- c()
    for (gene in names(gene_annot)){
      if (target %in% gene_annot[[gene]]){
        target_genes <- c(target_genes, gene)
        sharing_factor <- 1/length(gene_annot[[gene]]) # if gene belongs multiple groups
          }
    }
    if(is.null(target_genes)){
      print(paste("no genes in", target, "group"))
      next
    }
    target_rows <- proteomics[target_genes, ]*sharing_factor
    target_sums <- colSums(target_rows, na.rm = TRUE)
    sector_sum <- rbind(sector_sum, target_sums) 
  }
  
  summed_groups <- colSums(sector_sum, na.rm=TRUE)
  fractions <- summed_groups/colSums(proteomics, na.rm=TRUE)
  names(fractions) <- colnames(proteomics)
  return(fractions)
}


plot_and_fit <- function(data, growth_rates, ylim, ylab, xlab){
  if(length(data) > 8){
    print("Too many groups")
    return()
  }
  
  colors <- brewer.pal(8, "Dark2")
  for(i in seq_along(data)){
    group_data <- data[[names(data)[i]]]
    plot(growth_rates, group_data, col = colors[i], pch=20, ylim=ylim, cex=1.5,
         ylab = ylab, xlab = xlab, cex.lab = 1.3)
    abline(lm(group_data~growth_rates), col = colors[i])
    par(new=TRUE)
  }
  legend("topleft", 
         legend = names(data),
         col = colors[1:length(data)], pch=20, cex=1.3)
  par(new=FALSE)
}


met <- get_sector_fraction(c("lipid biosynthetic process",            
                  "amino acid biosynthetic process", "nucleobase biosynthetic process",
                  "generation of precursor metabolites and energy"), proteomics, gene_annot)
info <- get_sector_fraction(c("DNA-templated transcription", "DNA biosynthetic process",
                              "RNA processing", "translation", "protein maturation",
                              "tRNA aminoacylation"), proteomics, gene_annot)

sectors <- list("Metabolism" = met,
                "Genetic Information Processing" = info)
plot_and_fit(sectors, 
             c_lim_samples$Growth.rate..1.h., 
             ylim = c(0,0.5), 
             ylab = "Proteome fraction",
             bquote("Growth rate"~ mu ~ (h^-1)))


translation <- get_sector_fraction(c("translation", "protein maturation"), proteomics, gene_annot)
ribosomes <- colSums(proteomics[ribo_genes, ], na.rm = TRUE)
trna_synthetases <- get_sector_fraction(c("tRNA aminoacylation"), proteomics, gene_annot)

sectors <- list("Translation (w/o tRNA synthetases)" = translation-trna_synthetases,
                "Ribosomal proteins" = ribosomes,
                "tRNA synthetases" = trna_synthetases)
plot_and_fit(sectors, 
             c_lim_samples$Growth.rate..1.h., 
             ylim = c(0,0.4), 
             ylab = "Proteome fraction",
             bquote("Growth rate"~ mu ~ (h^-1)))


carbs <- get_sector_fraction(c("generation of precursor metabolites and energy"), proteomics, gene_annot)

sectors <- list("Carbohydrate metabolism" = carbs)
plot_and_fit(sectors, 
             c_lim_samples$Growth.rate..1.h., 
             ylim = c(0,0.3), 
             ylab = "Proteome fraction",
             bquote("Growth rate"~ mu ~ (h^-1)))



aas <- get_sector_fraction(c("amino acid biosynthetic process"), proteomics, gene_annot)
nts <- get_sector_fraction(c("nucleobase biosynthetic process"), proteomics, gene_annot)
lip <- get_sector_fraction(c("lipid biosynthetic process"), proteomics, gene_annot)

sectors <- list("Amino acid metabolism" = aas,
                "Nucleotide metabolism" = nts,
                "Lipid metabolism" = lip)
plot_and_fit(sectors, 
             c_lim_samples$Growth.rate..1.h., 
             ylim = c(0,0.2), 
             ylab = "Proteome fraction",
             bquote("Growth rate"~ mu ~ (h^-1)))




