composition <- transform(platypus, transform = "compositional")

wunifrac <- ordinate(composition, method = "PCoA", distance = "bray")
plot_ordination(physeq = composition,
                ordination = wunifrac,
                type = "samples",
                color = "Category") + 
  theme_bw() +
  geom_point(size=3) + ggtitle(" Bray")  + scale_color_manual(values = c("#9A879D","#D3D3D3"))


composition <- transform(phy, transform = "compositional")


wunifrac <- ordinate(composition, method = "PCoA", distance = "bray")
plot_ordination(physeq = composition,
                     ordination = wunifrac,
                     type = "samples",
                     color = "Category",
                shape = "Category") + scale_shape_manual(values = c(15,15,15,16,17))+ scale_color_manual(values = c("#92140C","#6A3D9A","#333333","#B2DF8A","#1F78B4","#1B998B","#48304D",
                                                                                                                    "#A6CEE3","#CAB2D6","#1F3172", "#E7298A","#B2DF8A","#FFFF99","#E31A1C",
                                                                                                                    "#1B9E77",  "#999999", "#1F3172","#FDBF6F"))+
  theme_bw() +
  geom_point(size=3) + ggtitle("Bray") 


rm(wunifrac)
rm(composition)

