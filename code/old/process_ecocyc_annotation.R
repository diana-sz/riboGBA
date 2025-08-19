library(here)

setwd(here())

ecocyc_annot <- read.table("data/ecoli_ecocyc_groups.txt", sep = "\t", header=TRUE)

d17_groups <- c("translation",
                "lipid biosynthetic process",          
                "amino acid biosynthetic process",
                "nucleobase biosynthetic process",
                "DNA-templated transcription",
                "RNA processing",
                "DNA biosynthetic process",
                "generation of precursor metabolites and energy",
                #"tRNA aminoacylation",
                "protein maturation")

# subgroups of translation that will be plotted separately
ribo_genes <- ecocyc_annot[ecocyc_annot$Gene.Ontology.Terms == "structural constituent of ribosome", 2]
ribo_genes <- strsplit(ribo_genes, " // ")[[1]]

trna_syn <- ecocyc_annot[ecocyc_annot$Gene.Ontology.Terms == "tRNA aminoacylation", 2]
trna_syn <- strsplit(trna_syn, " // ")[[1]]



gene_annot <- list()
for(row in 1:nrow(ecocyc_annot)){
  group <- ecocyc_annot[row,1]
  if(!group %in% d17_groups){
    next
  }
  genes <- strsplit(ecocyc_annot[row,2], " // ")[[1]]
  for(gene in genes){
    
    error_check <- try({gene_annot[[gene]] <- c(gene_annot[[gene]], group)}, silent = TRUE)
    if(class(error_check) == "try-error"){
      gene_annot[[gene]] <- c(group)
    }
  }
}


# custom 
n <- 0
for(gene in names(gene_annot)){
  if(length(gene_annot[[gene]]) > 1){
    #print(paste(gene, gene_annot[[gene]]))
    n = n+1
    }
}

print(paste(n, "/", length(gene_annot), "genes have non-unique annotation"))




