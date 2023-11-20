##############################
##          PERMANOVA       ##
##############################

library(vegan)
library(pairwiseAdonis)
library(microbiome)

###By all###

##A PERMANOVA was run on the data set with a bray distance matrix

comp_all <- microbiome::transform(phy, transform = "compositional")

#Generate bray distance matrix
bray_dist_matrix_all <- phyloseq::distance(comp_all, method = "bray")

#Dispersion test and plot; homogeneity of dispersion among groups is an assumption for adonis

#bray
dispr.bray_all <- vegan::betadisper(bray_dist_matrix_all, phyloseq::sample_data(comp_all)$Category)
dispr.bray_all
plot(dispr.bray_all)
anova(dispr.bray_all)#p=0.009848
#reject the assumption of homogeneity of dispersion by stage

Category <- phyloseq::sample_data(comp_all)$Category

#ADONIS test
vegan::adonis2(bray_dist_matrix_all ~ Category, strata = sample_sums(comp_all), permutations = 9999, permutest="dispersion") 

#         Df SumOfSqs      R2      F  Pr(>F)
#Category  4   6.1245 0.41354 5.8174  1e-04 ***
#Residual 33   8.6855 0.58646                  
#Total    37  14.8100 1.00000                                    

###Running a pairwise adonis

#pull out metadata from phyloseq object
meta <- meta(comp_all)

result <- pairwise.adonis(bray_dist_matrix_all, factors = meta$Category, perm = 9999, p.adjust.m = "holm")
write.csv(result, "pairwise-adonis-bray-all.csv")

#                     pairs Df SumsOfSqs    F.Model        R2 p.value p.adjusted sig
#1    Platypus vs Earthworm  1  1.543952   6.144988 0.2183333  0.0003     0.0027   *
#2     Platypus vs Mealworm  1  1.830067   7.692856 0.2590811  0.0006     0.0042   *
#3   Platypus vs Bloodworms  1  1.960436   8.293445 0.2737703  0.0003     0.0027   *
#4        Platypus vs Water  1  1.439520   4.671461 0.1474975  0.0001     0.0010  **
#5    Earthworm vs Mealworm  1  1.142013  12.560179 0.7584567  0.1000     0.3000    
#6  Earthworm vs Bloodworms  1  1.334127  16.145895 0.8014484  0.1000     0.3000    
#7       Earthworm vs Water  1  1.034071   2.697438 0.2306008  0.0056     0.0336   .
#8   Mealworm vs Bloodworms  1  1.479203 161.872176 0.9758850  0.1000     0.3000    
#9        Mealworm vs Water  1  1.271640   3.626107 0.2871912  0.0061     0.0336   .
#10     Bloodworms vs Water  1  1.260198   3.631652 0.2875041  0.0059     0.0336   .

Platypus - A
Earthworm - C
Water - B
Mealworm - C
Bloodworm - C

##Run again for just platypus

comp <- microbiome::transform(platypus, transform = "compositional")

#Generate bray distance matrix
bray_dist_matrix <- phyloseq::distance(comp, method = "bray")

#Dispersion test and plot; homogeneity of dispersion among groups is an assumption for adonis
#bray
dispr.bray <- vegan::betadisper(bray_dist_matrix, phyloseq::sample_data(comp)$Captive)
dispr.bray
plot(dispr.bray)
anova(dispr.bray)#p=0.4306
#fail to reject the assumption of homogeneity of dispersion by captivity status

captive <- phyloseq::sample_data(comp)$Captive

#ADONIS test
vegan::adonis2(bray_dist_matrix ~ captive, strata = sample_sums(comp), permutations = 9999) 

#         Df SumOfSqs      R2      F Pr(>F)   
#captive   1   0.4957 0.09536 2.0028 0.0333 *
#Residual 19   4.7030 0.90464                
#Total    20   5.1988 1.00000                                       

rm(comp)
rm(comp_all)
rm(dispr.bray)
rm(dispr.bray_all)
rm(meta)
rm(result)
rm(bray_dist_matrix)
rm(bray_dist_matrix_all)
rm(captive)
rm(Category)
