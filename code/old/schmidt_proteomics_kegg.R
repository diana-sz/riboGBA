rm(list=ls(all=TRUE))

library(RColorBrewer)
library(here)
require(dplyr)

directory <- paste0(here(), "/code")
setwd(directory)

annot_col <- "Annotated functional COG group (description)"
uniprot_col <- "Uniprot Accession"

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


annotation <- read.csv("../data/ecoli_kegg_annotation.csv")


growth_rates <-  metadata[colnames(proteomics), "Growth.rate..h.1."]
volumes <- metadata[colnames(proteomics), "Single.cell.volume..fl.1"]

for(col in c("X0", "X1", "X2")){
  new_col <- c()
  for(row in 1:nrow(annotation)){
    split <- unlist(strsplit(annotation[row, col], "', '"))
    chosen <- split[sample.int(length(split), 1)] # select random group if there are multiple
    new_col <- c(new_col, trimws(gsub("[^a-zA-Z \\s]", "", chosen)))
  }
  annotation[, col] <- new_col
}


get_sector_fraction <- function(targets, proteomics, annotation, level){
  sector_sum <- data.frame()
  for (target in targets){
    target_rows <- annotation$X[annotation[level] == target]
    target_sums <- colSums(proteomics[target_rows, ], na.rm = TRUE)
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




##### Translation & ribosomal proteins ######

met <- get_sector_fraction(c("Metabolism"), proteomics, annotation, "X0")
info <- get_sector_fraction(c("Genetic Information Processing"), proteomics, annotation, "X0")

sectors <- list("Metabolism" = met,
                "Genetic Information Processing" = info)
plot_and_fit(sectors, 
             growth_rates, 
             ylim = c(0,0.6), 
             ylab = "Proteome fraction",
             bquote("Growth rate"~ mu ~ (h^-1)))



translation <- get_sector_fraction(c("Translation", "Folding sorting and degradation"), proteomics, annotation, "X1")
ribosomes <- get_sector_fraction(c("Ribosome PATHeco"), proteomics, annotation, "X2")
trna_synthetases <- get_sector_fraction(c("AminoacyltRNA biosynthesis PATHeco"), proteomics, annotation, "X2")

sectors <- list("Translation" = translation-trna_synthetases,
                "Ribosomal proteins" = ribosomes,
                "tRNA synthetases" = trna_synthetases)
plot_and_fit(sectors, 
             growth_rates, 
             ylim = c(0,0.4), 
             ylab = "Proteome fraction",
             bquote("Growth rate"~ mu ~ (h^-1)))


carbs <- get_sector_fraction(c("Carbohydrate metabolism"), proteomics, annotation, "X1")
energy <- get_sector_fraction(c("Ribosome PATHeco"), proteomics, annotation, "X2")

sectors <- list("Carbohydrate metabolism" = carbs,
                "Energy metabolism" = energy)
plot_and_fit(sectors, 
             growth_rates, 
             ylim = c(0,0.3), 
             ylab = "Proteome fraction",
             bquote("Growth rate"~ mu ~ (h^-1)))



aas <- get_sector_fraction(c("Amino acid metabolism"), proteomics, annotation, "X1")
nts <- get_sector_fraction(c("Nucleotide metabolism"), proteomics, annotation, "X1")

sectors <- list("Amino acid metabolism" = aas,
                "Nucleotide metabolism" = nts)
plot_and_fit(sectors, 
             growth_rates, 
             ylim = c(0,0.15), 
             ylab = "Proteome fraction",
             bquote("Growth rate"~ mu ~ (h^-1)))


lip <- get_sector_fraction(c("Lipid metabolism"), proteomics, annotation, "X1")

sectors <- list("Lipid metabolism" = lip)
plot_and_fit(sectors, 
             growth_rates, 
             ylim = c(0,0.015), 
             ylab = "Proteome fraction",
             bquote("Growth rate"~ mu ~ (h^-1)))
