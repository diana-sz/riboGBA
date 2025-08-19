rm(list=ls(all=TRUE))

library(RColorBrewer)
library(here)

directory <- paste0(here(), "/code")
setwd(directory)

annot_col <- "Annotated functional COG group (description)"
uniprot_col <- "Uniprot Accession"

metadata <- read.csv("../data/Schmidt_proteomics_table_s23.csv", skip=2, nrows = 26)
metadata <- metadata[metadata$Strain == "BW25113",]
rownames(metadata) <- metadata[,1]

proteomics <- read.csv("../data/Schmidt_proteomics_table_s6.csv", skip=1)
start_index <- grep("Mass", colnames(proteomics))
end_index <- grep("Coefficient", colnames(proteomics)) - 1 
colnames(proteomics) <- proteomics[1,]
proteomics <- proteomics[2:nrow(proteomics),]
proteomics["Molecular weight (Da)"] <- lapply(proteomics["Molecular weight (Da)"], as.numeric)

masses <- proteomics[, start_index:end_index]
masses[] <- lapply(masses, as.numeric)


selected_conditions <- c("Glucose", "Acetate", "Fumarate", "Glucosamine", 
                         "Glycerol", "Pyruvate", "Chemostat 0.5", "Chemostat 0.35",
                         "Chemostat 0.20", "Chemostat 0.12", "Xylose", "Mannose",
                         "Galactose ", "Succinate", "Fructose")

masses <- masses[,selected_conditions]
growth_rates <-  metadata[colnames(masses), "Growth.rate..h.1."]
volumes <- metadata[colnames(masses), "Single.cell.volume..fl.1"]

get_sector_fraction <- function(targets, masses, proteomics, id_column){
  sector_sum <- data.frame()
  for (target in targets){
      target_rows <- proteomics[id_column] == target

      target_sums <- colSums(masses[target_rows, ], na.rm = TRUE)
      sector_sum <- rbind(sector_sum, target_sums) 
  }
  
  summed_groups <- colSums(sector_sum, na.rm=TRUE)
  fractions <- summed_groups/colSums(masses, na.rm=TRUE)
  names(fractions) <- colnames(masses)
  return(fractions)
}


get_sector_mw <- function(targets, proteomics, id_column){
  total_sum <- 0
  for (target in targets){
    target_masses <- proteomics["Molecular weight (Da)"][proteomics[id_column] == target]
    total_sum <- total_sum + sum(target_masses)
  }
  return(total_sum)
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


synthetases <- proteomics[grep("tRNA synthetase", proteomics$Description, ignore.case = TRUE), "Uniprot Accession"]
synthetase_sector <- get_sector_fraction(synthetases, masses, proteomics, "Uniprot Accession")

translation <- proteomics[grep("Translation, ribosomal structure and biogenesis",
                          proteomics[,annot_col], ignore.case = TRUE), "Uniprot Accession"]
translation_wo_synthetases <- setdiff(translation, synthetases)

trans_sector <- get_sector_fraction(translation_wo_synthetases, masses, proteomics, uniprot_col)
ptm_sector <- get_sector_fraction("Posttranslational modification, protein turnover, chaperones", masses, proteomics, annot_col)
trans_fraction <- trans_sector+ptm_sector


# get_sector_mw(ribosomal_proteins, proteomics, uniprot_col)
# get_sector_mw(translation_proteins, proteomics, annot_col)
# get_sector_mw(translation_wo_synthetases, proteomics, uniprot_col)
ribosomal_proteins <- proteomics[grep("Ribosomal protein", proteomics$Description, ignore.case = TRUE), "Uniprot Accession"]
rp_fraction <- get_sector_fraction(ribosomal_proteins, masses, proteomics, uniprot_col)

r_sector_groups <- list("Ribosomal proteins"=rp_fraction,
                        "Translation proteins (w/o tRNA synthetases)"=trans_fraction)#,
                       # "tRNA synthetases"=synthetase_sector)
plot_and_fit(r_sector_groups, growth_rates, c(0,0.4), "Proteome fraction", bquote("Growth rate"~ mu ~ (h^-1)))


mean(trans_fraction/rp_fraction) * get_sector_mw(ribosomal_proteins, proteomics, uniprot_col)


##### RNases #####
rnases <- c("P30850", "P21499", "P05055", "P21513")
rnase_fractions <- get_sector_fraction(rnases, masses, proteomics, uniprot_col)
plot_and_fit(list("RNases"=rnase_fractions), growth_rates, c(0,0.01), "Proteome fraction", bquote("Growth rate"~ mu ~ (h^-1)))



##### Transcription #####
#proteomics[proteomics[,annot_col] == "Transcription", "Description"]

polymerases <- proteomics[grep("RNA polymerase subunit", proteomics$Description, ignore.case = TRUE), uniprot_col]
rnap_mw <- get_sector_mw(polymerases, proteomics, uniprot_col)

# sigma <- proteomics[grep("sigma factor", proteomics$Description, ignore.case = TRUE), uniprot_col]
# avg_sigma <- get_sector_mw(sigma, proteomics, uniprot_col)/length(sigma)

rnap_fractions <- get_sector_fraction(polymerases, masses, proteomics, uniprot_col)

rna <- c("RNA processing and modification", "Transcription")
transcription_fractions <- get_sector_fraction(rna, masses, proteomics, annot_col)

dna <- c("Cell cycle control, cell division, chromosome partitioning",
        "Replication, recombination and repair")
dna_fractions <- get_sector_fraction(dna, masses, proteomics, annot_col)


transcription_sector <- list("RNA polymerases"=rnap_fractions,
                             "Transcription proteins"=transcription_fractions,
                             "DNA/cell_cycle"=dna_fractions)
plot_and_fit(transcription_sector, growth_rates, c(0,0.05), "Proteome fraction", bquote("Growth rate"~ mu ~ (h^-1)))



mean(transcription_fractions/rnap_fractions) * get_sector_mw(polymerases, proteomics, uniprot_col)


##### ######
#proteomics[proteomics[,annot_col] == "Amino acid transport and metabolism", "Description"]

eaa_fractions <- get_sector_fraction("Amino acid transport and metabolism", masses, proteomics, annot_col)
ent_fractions <- get_sector_fraction("Nucleotide transport and metabolism", masses, proteomics, annot_col)
anabolism_sectors <- list("AA metabolism" = eaa_fractions,
                          "NT metabolism"=ent_fractions)
plot_and_fit(anabolism_sectors, growth_rates, c(0,0.2), "Proteome fraction", bquote("Growth rate"~ mu ~ (h^-1)))


ep_fractions <- get_sector_fraction("Energy production and conversion", masses, proteomics, annot_col)
c_fractions <- get_sector_fraction("Carbohydrate transport and metabolism", masses, proteomics, annot_col)
catabolism_sectors <- list("Energy metabolism"=ep_fractions, "Carbohydrate metabolism"=c_fractions)
plot_and_fit(catabolism_sectors, growth_rates, c(0,0.4), "Proteome fraction", bquote("Growth rate"~ mu ~ (h^-1)))


l_fractions <- get_sector_fraction(c("Cell wall/membrane/envelope biogenesis",
                                     "Lipid transport and metabolism"), 
                                   masses, proteomics, annot_col)
plot_and_fit(list("Lipids/cell wall"=l_fractions), growth_rates, c(0,0.15), "Proteome fraction", bquote("Growth rate"~ mu ~ (h^-1)))




##### aggregate sectors #####

aggregated <- aggregate(masses, by = list(proteomics[,annot_col]), FUN = sum, na.rm = TRUE)

other <- c("Cell motility",
           "Defense mechanisms",
           "Function unknown",
           "General function prediction only",
           "Signal transduction mechanisms",
           "Intracellular trafficking, secretion, and vesicular transport",
           "Secondary metabolites biosynthesis, transport and catabolism",
           "Coenzyme transport and metabolism",
           "Inorganic ion transport and metabolism",
           "")
information_processing <- c("Posttranslational modification, protein turnover, chaperones",
                            "Transcription", "RNA processing and modification",
                            "Translation, ribosomal structure and biogenesis",
                            "Cell cycle control, cell division, chromosome partitioning",
                            "Replication, recombination and repair")
metabolism <- c("Carbohydrate transport and metabolism",
                "Nucleotide transport and metabolism",
                "Amino acid transport and metabolism",
                "Energy production and conversion",
                "Cell wall/membrane/envelope biogenesis",
                "Lipid transport and metabolism"
                #"Coenzyme transport and metabolism",
                #"Inorganic ion transport and metabolism"
                )

other_sector <- get_sector_fraction(other, masses, proteomics, annot_col)
info_sector <- get_sector_fraction(information_processing, masses, proteomics, annot_col)
met_sector <- get_sector_fraction(metabolism, masses, proteomics, annot_col)

sectors <- cbind(info_sector,met_sector,other_sector)

plot_and_fit(list("Other"=other_sector,
                  "Info protessing"=info_sector,
                  "met_sector"=met_sector), 
             growth_rates, c(0,0.7), "Proteome fraction", bquote("Growth rate"~ mu ~ (h^-1)))


ordered_mu <- growth_rates[order(growth_rates)]
ordered_phis <- sectors[order(growth_rates),]
colors <- brewer.pal(3, name = "Paired")
cum_comp <- t(apply(ordered_phis, 1, cumsum))
plot(NA, xlim = c(0, max(ordered_mu)), ylim=c(0,1),
     yaxs="i", xaxs="i", 
     #xaxt = "n", yaxt = "n",
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
legend('bottomright', rev(colnames(ordered_phis)), col = rev(colors), pch = 15, cex=1, border = NA)




# 
# aggregated <- aggregate(masses, by = list(proteomics$`Annotated functional COG class`), FUN = sum, na.rm = TRUE)
# rownames(aggregated) <- aggregated$Group.1
# aggregated <- aggregated[,2:ncol(aggregated)]
# aggregated["OTHER",] <- colSums(aggregated[c(1,2,5),])
# aggregated <- t(aggregated[c("INFORMATION STORAGE AND PROCESSING", "METABOLISM", "OTHER"),])
# 
# matplot(growth_rates,aggregated/volumes, pch=15:18)
# legend("topright", colnames(aggregated), pch=15:18)
# 
# # composition
# phis <- aggregated/rowSums(aggregated)
# matplot(growth_rates, phis, pch=15:18)
# legend("topright", colnames(aggregated), pch=15:18)
# 
# 
# ordered_mu <- growth_rates[order(growth_rates)]
# ordered_phis <- phis[order(growth_rates),]
# colors <- brewer.pal(3, name = "Paired")
# cum_comp <- t(apply(ordered_phis, 1, cumsum))
# plot(NA, xlim = c(0, max(ordered_mu)), ylim=c(0,1),
#      yaxs="i", xaxs="i", 
#      #xaxt = "n", yaxt = "n",
#      xlab=bquote("Growth rate"~ mu ~ (h^-1)), 
#      ylab="Relative proteome composition")
# for(i in ncol(cum_comp):1){
#   col <- colors[i]
#   polygon(c(ordered_mu, rev(ordered_mu)), c(rep(0, nrow(cum_comp)), rev(cum_comp[,i])),
#           col=col, border=col)
#   par(new=TRUE)
# }
# par(new=FALSE)
# legend('bottomright', colnames(ordered_phis), col = colors, pch = 15, cex=1, border = NA)


# plot(growth_rates, colSums(masses)/volumes)
# plot(growth_rates, volumes)
