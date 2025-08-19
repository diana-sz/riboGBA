library(RColorBrewer)
library(here)

setwd(here())

data <- read.csv("data/A8dekel1_dekel.csv", row.names = 1)

par(mar = c(4,4,1,1), mfrow = c(1,2))

# Make a discrete palette
x_C2_vals <- sort(unique(data$x_C2))
cols <- setNames(brewer.pal(length(x_C2_vals), "Set1"), x_C2_vals)

for(kcat in unique(data$kcat)){
  
  one_kcat <- data[data$kcat == kcat,]
  
  ref <- one_kcat[one_kcat$x_C2 == min(one_kcat$x_C2), "mu"]
  
  for(x_C2_val in unique(one_kcat$x_C2)){
    one_c <- one_kcat[one_kcat$x_C2 == x_C2_val, "mu"]
    one_kcat[one_kcat$x_C2 == x_C2_val, "benefit"] <- one_c - ref
  }
  
  # Assign colors based on discrete x_C2 values
  point_cols <- cols[as.character(one_kcat$x_C2)]
  
  plot(mu ~ phi, data = one_kcat, ylab = "Growth rate", col = point_cols, pch = 16)
  plot(benefit ~ phi, data = one_kcat, ylab = "Benefit", col = point_cols, pch = 16)
  
  legend("topright", legend = x_C2_vals, col = cols, pch = 16, title = "x_C2")
}
