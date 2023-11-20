#load packages
library(tidyverse)
library(phyloseq)
library(microbiomeMarker)


#assessing various methods for differential abundance testing
#LEFSe
lefse_out_platypus <- microbiomeMarker::run_lefse(platypus, taxa_rank= "Genus", wilcoxon_cutoff = 0.05, group= "Captive", kw_cutoff = 0.05, multigrp_strat = FALSE, lda_cutoff = 4)

#looking at microbiome markers - LEFSe
library(knitr)
lefse_out_platypus #6 markers identified

dat<-marker_table(lefse_out_platypus) %>% data.frame() %>% select(1:4)
head(dat)
dat %>% kable(align = "c")
write.csv(dat, "platypus_lefse_out_genus.csv")

p_abund <- plot_abundance(lefse_out_platypus, group = "Captive")
p_abund

plot_heatmap(lefse_out_platypus, transform = "log10p", group = "Captive")
plot_ef_bar(lefse_out_platypus)

plot_cladogram(lefse_out_platypus, color = c(Wild = "darkgreen", Captive = "darkblue"), clade_label_level = 4) +
  theme(plot.margin = margin(0, 0, 0, 0))

