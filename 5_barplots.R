library("ggplot2")
library("RColorBrewer")
library("microbiome")
##########################################
##                Bar plots             ##
##########################################

#Merge by groups for figure

all <- merge_samples(phy, "Category")
all <- prune_taxa((taxa_sums(all) > 0), all)# 1648 ASVs in 5 samples
Captive <- merge_samples(platypus, "Captive")
Captive <- prune_taxa((taxa_sums(Captive) > 0), Captive)# 654 ASVs in 2 samples


 ##Genus level##
# To represent at the Genus level (instead of ASV level)
Genus_all <- tax_glom(all,taxrank = "Genus")
Genus_Captive <- tax_glom(Captive,taxrank = "Genus")

# Instead of making copies of this code - I just run the same code and change the phyloseq object that is being transformed
#I also manually reorganized the colours so each group (all, Captive, stage) would match.
# Transform counts in relative abundance and select most abundant families
Genus <- transform_sample_counts(Genus_Captive, function(x) 100 * x/sum(x))
genus <- psmelt(Genus)
genus$Genus <- as.character(genus$Genus)

#rename Genera with <3% abundance
genus$Genus[genus$Abundance < 3]<- "< 3% Abundance"
write.csv(genus, "genus_plat_barplot_all.csv")
genus <- read.csv("genus_plat_barplot.csv")

#How many levels in Genus
HowMany <- length(levels(as.factor(genus$Genus)))

#reorder levels of all
genus$Sample <- factor(genus$Sample, levels = c("Platypus","Water","Earthworm","Bloodworms","Mealworm"))


ggplot(genus, aes(x = Sample, y = Abundance, fill = reorder(Genus,Abundance))) +  
  geom_bar(stat = "identity") +
  theme(legend.position="bottom", axis.title.x = element_blank()) +  geom_col(width = 0.1) +
  ylab("Relative Abundance of Bacterial Genera \n") +  
  scale_fill_manual(values = c("#9A879D","#f7fff7","#0B5351","#92140C","#7D84B2","#4B7F52","#FB9A99","#E3D8F1","#E56B70","#1B998B","#48304D",
                               "#A6CEE3","#6A3D9A","#1F78B4","#CAB2D6","#33A02C","#1F3172", "#E7298A","#B2DF8A","#FFFF99","#E31A1C",
                               "#333333","#1B9E77",  "#999999", "#1F3172","#FDBF6F")) + coord_flip()


##Family level## 
# To represent at the Family level (instead of ASV level)
Family <- tax_glom(all,taxrank = "Family")

# Transform counts in relative abundance and select most abundant families
Family <- transform_sample_counts(Family, function(x) 100 * x/sum(x))
family <- psmelt(Family)
family$Family <- as.character(family$Family)

#rename Families with <5% abundance
family$Family[family$Abundance < 5]<- "< 5% Abundance"

#How many levels in Family
HowMany <- length(levels(as.factor(family$Family)))

write.csv(family, "family_all_barplot.csv")
family <- read.csv("family_all_barplot.csv")

#reorder levels
family$Sample <- factor( family$Sample, levels = c("Platypus","Water","Earthworm","Bloodworms","Mealworm"))


ggplot(family, aes(x = Sample, y = Abundance, fill = reorder(Family,Abundance))) +  
  geom_bar(stat = "identity", width = 0.85) +
  theme(legend.position="right", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab("Relative Abundance of Bacterial Families \n")+  
  scale_fill_manual(values = c("#9A879D","#f7fff7","#0B5351","#92140C","#7D84B2","#4B7F52","#FB9A99","#E3D8F1","#E56B70","#1B998B","#48304D",
                               "#A6CEE3","#6A3D9A","#1F78B4","#FF7F00","#CAB2D6","#33A02C","#1F3172", "#E7298A","#B2DF8A","#FFFF99","#E31A1C",
                               "#333333","#1B9E77",  "#999999", "#1F3172","#FDBF6F"))

##Phylum level## 
# To represent at the Family level (instead of ASV level)
Phylum <- tax_glom(Captive,taxrank = "Phylum")

# Transform counts in relative abundance and select most abundant families
Phylum <- transform_sample_counts(Phylum, function(x) 100 * x/sum(x))
family <- psmelt(Phylum)
family$Phylum <- as.character(family$Phylum)

#rename Families with <1% abundance
family$Phylum[family$Abundance < 1]<- "< 1% Abundance"

#How many levels in Phylum
HowMany <- length(levels(as.factor(family$Phylum)))

write.csv(family, "phylum_barplot_justplaty.csv")
family <- read.csv("phylum_barplot.csv")

#reorder levels
family$Sample <- factor( family$Sample, levels = c("Platypus","Water","Earthworm","Bloodworms","Mealworm"))


ggplot(family, aes(x = Sample, y = Abundance, fill = reorder(Phylum, Abundance))) +  
  geom_bar(stat = "identity", width = 0.85) +
  theme(legend.position="right", axis.title.x = element_blank()) + 
  ylab("Relative Abundance (%) of Bacterial Phyla \n")+  
  scale_fill_manual(values = c("#9A879D","#f7fff7","#0B5351","#92140C","#7D84B2","#4B7F52","#FB9A99","#E3D8F1","#E56B70","#1B998B","#48304D"))

#for different species
family <- read.csv("phylum_barplot_other_organisms.csv")

#reorder levels
family$Sample <- factor( family$Sample, levels = c("Platypus","Tasmanian Devil","Human","Echidna"))


ggplot(family, aes(x = Sample, y = Abundance, fill = reorder(Phylum, Abundance))) +  
  geom_bar(stat = "identity", width = 0.85) +
  theme(legend.position="right", axis.title.x = element_blank()) + 
  ylab("Relative Abundance (%) of Bacterial Phyla \n")+  
  scale_fill_manual(values = c("#f7fff7","#0B5351","#92140C","#7D84B2","#1B998B","#48304D"))

ggplot(family, aes(x = Sample, y = Abundance, fill = reorder(Phylum, Abundance))) +  
  geom_bar(stat = "identity", width = 0.85) +
  theme_minimal() +  # Add this line to set a plain white background
  theme(legend.position = "right", axis.title.x = element_blank()) + 
  ylab("Relative Abundance (%) of Bacterial Phyla \n") +  
  scale_fill_manual(values = c("#48304D","#0B5351","#92140C","#7D84B2","#1B998B"))

