#update R in RStudio (do this only once, before running analyses)
#library(installr)
#updateR()v4.3.1

#set working directory
setwd("C:/Users/dungana/Dropbox/ARC_DP_MinimalMicrobiome/Projects/Platypus/R/input_files")

# load useful function
KillZeroRCs <- function(x) {
  x[ which( rowSums(x) != 0) , ] -> x
  x[ , which( colSums(x) != 0) ] -> x
  return(x)
}

# read in OTU table
ASV <- read.table(
  "asv.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1
)

# remove zero-sum row/columns
ASV <- KillZeroRCs(ASV)

# read in tax table
# fill logical indicates that rows have unequal lengths due to blank fields
TAXA <- read.table(
  "taxonomy.txt",
  sep = "\t",
  fill = TRUE,
  header = TRUE,
  row.names = 1
)

colnames(TAXA) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")


# read in METAdata
META <- read.table(
  "metadata.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1
)


# convert OTU and tax tables to matrices for phyloseq
ASV <- as.matrix(ASV)
TAXA <- as.matrix(TAXA)

# load libraries
library(ape)
library(phyloseq)
library(microbiome)

# convert tree file to phyloseq object
TREE <- phyloseq::read_tree("tree.nwk")

# Create the otu_table from your ASV matrix
otu_table <- otu_table(ASV, taxa_are_rows = TRUE)

# Create the tax_table from your TAXA data frame
# Assuming that the first column contains OTU names, and the rest represent taxonomic ranks
tax_table <- tax_table(TAXA)

# Create the sample_data from your META data frame (assuming you have this)
# Replace META with your actual sample metadata
sample_data <- sample_data(META)

# Create the phyloseq object
phy <- phyloseq(otu_table, tax_table, sample_data)

# merge phyloseq objects
phy <- merge_phyloseq(phy, TREE)

# Delete the matrices as they are not used in further analyses
rm(ASV)
rm(TAXA)
rm(META)
rm(TREE)
rm(sample_data)
rm(otu_table)
rm(tax_table)

# remove zero sum ASV
phy <- prune_taxa((taxa_sums(phy) > 0), phy) #2323 ASVs in 45 samples

#remove ASVs with less than 10 reads as this is not biologically feasible
phy <- prune_taxa((taxa_sums(phy) > 10), phy) #1685 ASVs in 45 samples

#adjust tree for ASVs that have been removed
phy_tree(phy) <- root(phy_tree(phy),sample(taxa_names(phy),1), resolve.root = TRUE)




