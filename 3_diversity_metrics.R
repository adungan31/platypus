library(plotrix)
library(microbiome)
library(gridExtra)

# check read depth & select rarefaction level
sort(sample_sums(phy))
sort(sample_sums(platypus))

# rarefy
phy_rare<- rarefy_even_depth(phy, sample.size = 7491, rngseed = 1)#minimum number of reads
platypus_rare<- rarefy_even_depth(platypus, sample.size = 30503, rngseed = 1)#minimum number of reads

#4 samples and 800 ASVs were removed (all)
#2 ASVs were removed (platypus)

### All samples ###

# extract metadata from subsetted file
divMeta <- meta(phy_rare)

# add diversity index data to metadata file
divMeta <- cbind(divMeta, estimate_richness(phy_rare))

##BoxPlot###

#reorder levels of Category
divMeta$Category <- factor(divMeta$Category, levels=c("Platypus","Water","Earthworm","Bloodworms","Mealworm"))


phy.obs1 <- ggplot(divMeta, aes(x=Category, y=Observed, fill=Category), rm.na=TRUE) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none")+ theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#B2DF8A","#1F78B4","#6A3D9A","#92140C","#333333")) + lims(y=c(0,420))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("Observed ASVs")
print(phy.obs1)

phy.sim1 <- ggplot(divMeta, aes(x=Category, y=Simpson, fill=Category)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#B2DF8A","#1F78B4","#6A3D9A","#92140C","#333333"))+ lims(y=c(0,1))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("Simpson")
print(phy.sim1)

phy.shan1 <- ggplot(divMeta, aes(x=Category, y=Shannon, fill=Category)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#B2DF8A","#1F78B4","#6A3D9A","#92140C","#333333"))+ lims(y=c(0,5))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("Shannon")
print(phy.shan1)


grid.arrange(phy.obs,phy.sim,phy.shan, ncol=3)

### platypus samples ###

# extract metadata from subsetted file
divplatypus <- meta(platypus_rare)

# add diversity index data to metadata file
divplatypus <- cbind(divplatypus, estimate_richness(platypus_rare))

##BoxPlot###

#reorder levels of Captive
divplatypus$Captive <- factor(divplatypus$Captive, levels=c("Wild","Captive"))


phy.obs2 <- ggplot(divplatypus, aes(x=Captive, y=Observed, fill=Captive), rm.na=TRUE) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none")+ theme(axis.title.x = element_blank())+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#9A879D","#D3D3D3"))+ lims(y=c(0,175))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))
  
print(phy.obs2)

phy.sim2 <- ggplot(divplatypus, aes(x=Captive, y=Simpson, fill=Captive)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") +  theme(axis.title.x = element_blank())+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#9A879D","#D3D3D3"))+ lims(y=c(0,1))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))
print(phy.sim2)

phy.shan2 <- ggplot(divplatypus, aes(x=Captive, y=Shannon, fill=Captive)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.title.x = element_blank())+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#9A879D","#D3D3D3"))+ lims(y=c(0,4))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))
print(phy.shan2)


grid.arrange(phy.obs1,phy.sim1,phy.shan1,phy.obs2,phy.sim2,phy.shan2, ncol=3, nrow=2)

rm(phy.obs2)
rm(phy.shan2)
rm(phy.sim2)
rm(platypus_rare)
rm(phy_rare)
rm(phy.obs1)
rm(phy.shan1)
rm(phy.sim1)