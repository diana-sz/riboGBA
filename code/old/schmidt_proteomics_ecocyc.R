rm(list=ls(all=TRUE))

library(RColorBrewer)
library(here)
require(dplyr)

directory <- paste0(here(), "/code")
setwd(directory)

metadata <- read.csv("../data/Schmidt_proteomics_table_s23.csv", skip=2, nrows = 26)
metadata <- metadata[metadata$Strain == "BW25113",]
rownames(metadata) <- metadata[,1]

proteomics_raw <- read.csv("../data/Schmidt_proteomics_table_s6.csv", skip=1)
start_index <- grep("Mass", colnames(proteomics_raw))
end_index <- grep("Coefficient", colnames(proteomics_raw)) - 1 
colnames(proteomics_raw) <- proteomics_raw[1,]
proteomics_all <- proteomics_raw[2:nrow(proteomics_raw),]
gene_ids <- proteomics_raw$Bnumber

proteomics_masses <- proteomics_all[, start_index:end_index]
proteomics_masses[] <- lapply(proteomics_masses, as.numeric)
proteomics_masses$Bnumber <- proteomics_all$Bnumber


selected_conditions <- c("Glucose", "Acetate", "Fumarate", "Glucosamine", 
                         "Glycerol", "Pyruvate", "Chemostat 0.5", "Chemostat 0.35",
                         "Chemostat 0.20", "Chemostat 0.12", "Xylose", "Mannose",
                         "Galactose ", "Succinate", "Fructose")

proteomics <- proteomics_masses[,selected_conditions]
proteomics <- proteomics %>% 
  group_by(proteomics_masses$Bnumber) %>% 
  summarise_all(list(sum=sum))

proteomics <- data.frame(proteomics)
rownames(proteomics) <- proteomics$`proteomics_masses.Bnumber`
proteomics <- proteomics[,2:ncol(proteomics)]
colnames(proteomics) <- selected_conditions

growth_rates <-  metadata[colnames(proteomics), "Growth.rate..h.1."]


source("process_ecocyc_annotation.R")


#### Function definitions ######################################################
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


get_sector_fraction_gene_list <- function(target_genes, proteomics){
  target_rows <- proteomics[target_genes, ]
  target_sums <- colSums(target_rows, na.rm = TRUE)
  fractions <- target_sums/colSums(proteomics, na.rm=TRUE)
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


#### Make plots ################################################################

met <- get_sector_fraction(c("lipid biosynthetic process",            
                  "amino acid biosynthetic process", "nucleobase biosynthetic process",
                  "generation of precursor metabolites and energy"), proteomics, gene_annot)
info <- get_sector_fraction(c("DNA-templated transcription", "DNA biosynthetic process",
                              "RNA processing", "translation", "protein maturation",
                              "tRNA aminoacylation"), proteomics, gene_annot)

sectors <- list("Metabolism" = met,
                "Genetic Information Processing" = info)
plot_and_fit(sectors, 
             growth_rates, 
             ylim = c(0,0.5), 
             ylab = "Proteome fraction",
             bquote("Growth rate"~ mu ~ (h^-1)))


translation <- get_sector_fraction(c("translation", "protein maturation"), proteomics, gene_annot)
ribosomes <- get_sector_fraction_gene_list(ribo_genes, proteomics)
trna_synthetases <- get_sector_fraction_gene_list(trna_syn, proteomics)

sectors <- list("Translation (w/o tRNA synthetases)" = translation-trna_synthetases,
                "Ribosomal proteins" = ribosomes,
                "tRNA synthetases" = trna_synthetases)
plot_and_fit(sectors, 
             growth_rates, 
             ylim = c(0,0.4), 
             ylab = "Proteome fraction",
             bquote("Growth rate"~ mu ~ (h^-1)))


carbs <- get_sector_fraction(c("generation of precursor metabolites and energy"), proteomics, gene_annot)

sectors <- list("Carbohydrate metabolism" = carbs)
plot_and_fit(sectors, 
             growth_rates, 
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
             growth_rates, 
             ylim = c(0,0.2), 
             ylab = "Proteome fraction",
             bquote("Growth rate"~ mu ~ (h^-1)))




