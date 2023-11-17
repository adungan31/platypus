library(decontam)
library(microbiome)

# identify contaminants first from PCR negatives
consList <- isContaminant(seqtab = phy, neg = "PCR_Neg", method = "prevalence")

# pull out the names of contaminants
cons <- rownames(consList)[consList$contaminant=="TRUE"]
cons <- as.character(cons) #8 contaminants identified from PCR_Negs

# - - - - - - - - - - - - - - - - - - - - - - - - - #

# to get info on the contaminants, uncomment the following code to
# run it ON THE FILE WITH THE CONTAMINANT ASVs,
# then combine the consPer.csv and taxonomy.csv file data

# subset the non PCRNeg platypus samples
#vvv <- subset_samples(phy, PCR_Neg == "FALSE")
# merge the samples
#yyy <- merge_samples(vvv, "PCR_Neg", fun = sum)
# transform counts to percentages
#yyy <- transform_sample_counts(yyy, function(x) 100 * x/sum(x))
# extract the cons percentage data
#zzz <- prune_taxa(x = yyy, taxa = cons)
# write otu table to dataframe
#xxx <- data.frame(t(zzz@otu_table))
# write xxx to csv
#write.csv(x = xxx, row.names = TRUE, file = "consPCRNeg.csv")
# subset the contaminant ASVs
#phyCons <- prune_taxa(phy, taxa = cons)
# write the contaminants to a file for reference
#contaminants <-phyCons@tax_table
#contaxa <- contaminants@.Data
#write.csv(contaxa, "contaxaPCR_Neg.csv")

# 8 contaminant ASVs in PCR negatives
# total contamination in the samples = 0.58%


# - - - - - - - - - - - - - - - - - - - - - - - - - #

# remove the contaminants from the main phy phyloseq file
phy <- remove_taxa(phy, taxa = cons)


#Remove PCR_Negs
phy <- subset_samples(phy, PCR_Neg == "FALSE") 


# identify contaminants from extract blanks
consList2 <- isContaminant(seqtab = phy, neg = "Ext_Blank", method = "prevalence")

# pull out the names of contaminants
cons2 <- rownames(consList2)[consList2$contaminant=="TRUE"]
cons2 <- as.character(cons2) #16 contaminants identified from extraction blanks

# - - - - - - - - - - - - - - - - - - - - - - - - - #

# to get info on the contaminants, uncomment the following code to
# run it ON THE FILE WITH THE CONTAMINANT ASVs,
# then combine the consPer.csv and taxonomy.csv file data

# subset the non extraction blank platypus samples
#vvv2 <- subset_samples(phy, Ext_Blank == "FALSE")
# merge the samples
#yyy2 <- merge_samples(vvv2, group = "Ext_Blank", fun = sum)
# transform counts to percentages
#yyy2 <- transform_sample_counts(yyy2, function(x) 100 * x/sum(x))
# extract the cons percentage data
#zzz2 <- prune_taxa(x = yyy2, taxa = cons2)
# write otu table to dataframe
#xxx2 <- data.frame(t(zzz2@otu_table))
# write xxx to csv
#write.csv(x = xxx2, row.names = TRUE, file = "consExtraction.csv")
# subset the contaminant ASVs
#phyCons2 <- prune_taxa(phy, taxa = cons2)
# write the contaminants to a file for reference
#contaminants2 <-phyCons2@tax_table
#contaxa2 <- contaminants2@.Data
#write.csv(contaxa2, "contaxaExtraction.csv")


# 16 contaminant ASVs in Extraction blanks
# total contamination in the samples from extraction = 1.73%
# total contamination in the samples = 2.31%

# remove the contaminants from the main phy phyloseq file
phy <- remove_taxa(phy, taxa = cons2)

# Remove extraction Blanks from phy 
phy <- subset_samples(phy, Ext_Blank=="FALSE")

#view reads by sample
sort(sample_sums(phy))
phy <- prune_samples((sample_sums(phy) > 10), phy)
phy_tree(phy) <- root(phy_tree(phy),sample(taxa_names(phy),1), resolve.root = TRUE)


platypus <- subset_samples(phy, Species == "Platypus")
platypus <- prune_taxa((taxa_sums(platypus) > 0), platypus)# 654 ASVs in 22 samples
phy_tree(platypus) <- root(phy_tree(platypus),sample(taxa_names(platypus),1), resolve.root = TRUE)
# - - - - - - - - - - - - - - - - - - - - - - - - - #

rm(phyCons)
rm(phyCons2)
rm(vvv)
rm(vvv2)
rm(xxx)
rm(xxx2)
rm(yyy)
rm(yyy2)
rm(zzz)
rm(zzz2)
rm(cons)
rm(cons2)
rm(contaminants)
rm(contaminants2)
rm(contaxa)
rm(contaxa2)
rm(consList)
rm(consList2)
