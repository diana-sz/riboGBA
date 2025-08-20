library(RColorBrewer)
library(here)

setwd(here())
#A8dekel1_dekel.csv
#A8dekel0_dekel.csv
data <- read.csv("data/A8dekel0_dekel.csv", row.names = 1)

par(mar = c(4,4,3,1), mfrow = c(1,2))

# Make a discrete palette
x_C2_vals <- sort(unique(data$x_C2))
cols <- setNames(brewer.pal(length(x_C2_vals), "Set1"), x_C2_vals)

prot_vals <- sort(unique(data$protein))
shapes <- setNames(16:(16+length(prot_vals)-1), prot_vals)

for(kcat in unique(data$kcat)){
  
  one_kcat <- data[data$kcat == kcat,]
  
  ref <- one_kcat[one_kcat$x_C2 == min(one_kcat$x_C2), "mu"]
  
  for(x_C2_val in unique(one_kcat$x_C2)){
    one_c <- one_kcat[one_kcat$x_C2 == x_C2_val, "mu"]
    one_kcat[one_kcat$x_C2 == x_C2_val, "benefit"] <- one_c - ref
  }
  
  # Assign colors based on discrete x_C2 values
  point_cols <- cols[as.character(one_kcat$x_C2)]
  point_shapes <- shapes[as.character(one_kcat$protein)]
  
  plot(mu ~ phi, data = one_kcat, ylab = "Growth rate", col = point_cols, pch = point_shapes)
  plot(benefit ~ phi, data = one_kcat, ylab = "Benefit", col = point_cols, pch = point_shapes)
  title(main = paste0("kcat ", kcat), outer = TRUE, line = -1)
  
  legend("topright", legend = x_C2_vals, col = cols, pch = 16, title = "x_C2")
}
