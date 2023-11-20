library(UpSetR)
library(MicrobiotaProcess)

##working with a phyloseq object called 'phy'

data <- get_upset(phy, factorNames = "Source")


upset(data, sets = c("Platypus", "Water", "Food"), 
      sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")

