library(readxl)
AJP10_16S <- read_excel("~/Pinto/AJP10/AJP10_16S.xlsx", 
                        sheet = "taxa_full", col_types="numeric")
library(readxl)
AJP10_16S <- read_excel("~/Pinto/AJP10/16S_OTU/newest/otu_16s_seqtab_v3.xlsx", 
                        sheet = "otu_16s_seqtab_v3", col_types="numeric")
library(vegan)
library(dplyr)

nitro_meta1 <- read_excel("~/Pinto/DO Transcriptomics Project Data/Transcriptomics_Data/nitrogen_conc.xlsx", 
                          sheet = "Metadata Setup")


nitro_meta2 <-as.data.frame(nitro_meta1)
nitro_meta <-as.data.frame(nitro_meta2, row.names=nitro_meta2[,1])


taxa_rRNA <- AJP10_16S[,-1]
taxa_rRNA <-as.data.frame(taxa_rRNA)
taxa_rRNA <-t(taxa_rRNA)
taxa_rRNA <-as.data.frame(taxa_rRNA)
taxa_rRNA_per <- taxa_rRNA %>%
  mutate(across(where(is.numeric), ~ ./sum(.)))


set.seed(123)
taxa_nmds = metaMDS(taxa_rRNA_per, distance = "bray")
taxa_nmds$diss
taxa_scores = as.data.frame(scores(taxa_nmds)$sites)
taxa_scores

taxa_scores$Type = as.factor(nitro_meta$Type)
taxa_scores$Reactor = nitro_meta$Reactor
taxa_scores$Oxygen = as.factor(nitro_meta$Oxygen)
taxa_scores$NuAc = as.factor(nitro_meta$`NA`)
taxa_meta <-nitro_meta[,c(2,4,5,12)]

taxa_en = envfit(taxa_nmds, taxa_meta, permutations=999, na.rm=TRUE)
taxa_en_coord_cont = as.data.frame(scores(taxa_en, "vectors")) * ordiArrowMul(taxa_en)
taxa_en_coord_cat = as.data.frame(scores(taxa_en, "factors")) * ordiArrowMul(taxa_en)

layout(matrix(c(1,2,3), ncol=3,nrow=1))

par(mar=c(5,5,1,1))
plot(taxa_scores[,2]~taxa_scores[,1], xlim=c(-5,5), ylim=c(-2,5), cex=3, lwd=2, xlab="NMDS1 (27.6%)", ylab="NMDS2 (9.6%)",
     pch=c(24,24,24,24,21,21,21,21,22,22,22,22,24,24,24,24,21,21,21,21,22,22,22,22,24,24,24,24,21,21,21,21,22,22,22,22,17,17,19,19,15,15),
     col="black", bg=c("plum1","plum1","plum3","plum3","plum1","plum1","plum2","plum2","plum3","plum3","plum4","plum4","skyblue1","skyblue1","skyblue3","skyblue3","skyblue1","skyblue1","skyblue2","skyblue2","skyblue3","skyblue3","skyblue4","skyblue4","springgreen1","springgreen1","springgreen2","springgreen2","springgreen1","springgreen1","springgreen2","springgreen2","springgreen3","springgreen3","springgreen4","springgreen4","black","black","black","black","black","black"),
     axes="F", cex.lab=1.5 )
axis(1, lwd=3, lwd.ticks=3, cex.axis=1.5)
axis(2, lwd=3, lwd.ticks=3, cex.axis=1.5)
legend("topright", bty="n", cex=1.5, title="DNA and RNA", legend=c("SS R1 & R3","SS R4 & R5","IFAS R4 & R5","6 mg/L","4 mg/L","2 mg/L","DNA"), pch=c(21,22,24,21,21,21,21), pt.bg=c("lightgray","darkgray","lightgray","springgreen","skyblue","plum1","black") )

#legend("topright", bty="n", cex=1.5, legend=c("Anox. R1","Anox. R3","Aero. R4_1","Aero. R4_2","IFAS R4_1","IFAS R4_2","6 mg/L","4 mg/L","2 mg/L","DNA"), pch=c(21,21,22,22,24,24,21,21,21,21), pt.bg=c("lightgray","gray","darkgray","dimgray","lightgray","darkgray","springgreen","skyblue","plum1","black") )
legend("topleft", bty="n", cex=1.5, legend="a)")
#arrows(x0=0,y0=0, x1=taxa_en_coord_cont[,1]/6, y1=taxa_en_coord_cont[,2]/6, lwd=2)
#arrows(x0=0,y0=0, x1=taxa_en_coord_cat[,1]/6, y1=taxa_en_coord_cat[,2]/6, lwd=2)
#text(x=taxa_en_coord_cat[,1]/6,y=taxa_en_coord_cat[,2]/6,labels=row.names(taxa_en_coord_cat), cex=1, col="purple4")
#text(x=taxa_en_coord_cont[,1]/6,y=taxa_en_coord_cont[,2]/6,labels=row.names(taxa_en_coord_cont), cex=1, col="purple4")

ifas_rRNA <-taxa_rRNA_per[c(1:4,13:16,25:28),]
set.seed(123)
ifas_nmds = metaMDS(ifas_rRNA, distance = "bray")
ifas_nmds$diss
ifas_scores = as.data.frame(scores(ifas_nmds)$sites)


nMDS <- metaMDS(ifas_rRNA, distance = "bray", k = 2)
MDS <- cmdscale(vegdist(ifas_rRNA, method = "bray"), k = 2, eig = T, add = T )
round(MDS$eig*100/sum(MDS$eig),1)

#c(4,6,7,8,13,14,16,17)
ifas_meta <-nitro_meta[c(1:4,13:16,25:28),c(6,7,8)]
#ifas_scores$DO = as.factor(ifas_meta$`DO`)
ifas_scores$NH3 = ifas_meta$NH3avg
ifas_scores$NH3 = ifas_meta$NH3avg
ifas_scores$NO2 = ifas_meta$NO2avg
ifas_scores$NO3 = ifas_meta$NO3avg
#ifas_scores$Reactor = as.factor(ifas_meta$Reactor)


ifas_en = envfit(ifas_nmds, ifas_meta, permutations=999, na.rm=TRUE)
ifas_en_coord_cont = as.data.frame(scores(ifas_en, "vectors")) * ordiArrowMul(ifas_en)
#ifas_en_coord_cat = as.data.frame(scores(ifas_en, "factors")) * ordiArrowMul(ifas_en)


par(mar=c(5,5,1,1))
plot(ifas_scores[,2]~ifas_scores[,1], xlim=c(-0.6,0.6), ylim=c(-0.6,0.6), cex=3, lwd=2, xlab="NMDS1 (20.6%)", ylab="NMDS2 (12.3%)",
     pch=24,
     col="black", bg=c("plum1","plum1","plum3","plum3","skyblue","skyblue","skyblue3","skyblue3","springgreen","springgreen","springgreen3","springgreen3"),
     axes="F", cex.lab=1.5 )
axis(1, lwd=3, lwd.ticks=3, cex.axis=1.5)
axis(2, lwd=3, lwd.ticks=3, cex.axis=1.5)
legend("topright", bty="n", cex=1.5, title="IFAS RNA", legend=c("R4","R5","6 mg/L","4 mg/L","2 mg/L"), pch=24, pt.bg=c("lightgrey","gray","springgreen","skyblue","plum1") )
legend("topleft", bty="n", cex=1.5, legend="b)")
arrows(x0=0,y0=0, x1=ifas_en_coord_cont[,1]/7, y1=ifas_en_coord_cont[,2]/7, lwd=2)
#arrows(x0=0,y0=0, x1=ifas_en_coord_cat[,1], y1=ifas_en_coord_cat[,2], lwd=2)
#text(x=ifas_en_coord_cat[,1],y=ifas_en_coord_cat[,2],labels=row.names(ifas_en_coord_cat), cex=1, col="purple4")
text(x=ifas_en_coord_cont[,1]/8,y=ifas_en_coord_cont[,2]/8,labels=row.names(ifas_en_coord_cont), cex=1, col="purple4")




ss_rRNA <-taxa_rRNA_per[c(5:12,17:24,29:36),]
set.seed(123)
ss_nmds = metaMDS(ss_rRNA, distance = "bray")
ss_nmds$diss
ss_scores = as.data.frame(scores(ss_nmds)$sites)

ss_meta <-nitro_meta[c(5:12,17:24,29:36),c(6,7,8)]
#ss_scores$DO = as.factor(ss_meta$`DO`)
ss_scores$NH3 = ss_meta$NH3avg
ss_scores$NH3 = ss_meta$NH3avg
ss_scores$NO2 = ss_meta$NO2avg
ss_scores$NO3 = ss_meta$NO3avg
#ss_scores$Reactor = as.factor(ss_meta$Reactor)


ss_en = envfit(ss_nmds, ss_meta, permutations=999, na.rm=TRUE)
ss_en_coord_cont = as.data.frame(scores(ss_en, "vectors")) * ordiArrowMul(ss_en)

ss_ramp_col <-c("plum1","plum1","plum2","plum2","plum3","plum3","plum4","plum4","skyblue1","skyblue1","skyblue2","skyblue2","skyblue3","skyblue3","skyblue4","skyblue4","springgreen1","springgreen1","springgreen2","springgreen2","springgreen3","springgreen3","springgreen4","springgreen4")
ss_bi_col <-c("plum1","plum1","plum3","plum3","plum1","plum1","plum3","plum3","skyblue1","skyblue1","skyblue3","skyblue3","skyblue1","skyblue1","skyblue3","skyblue3","springgreen","springgreen","springgreen3","springgreen3","springgreen1","springgreen1","springgreen3","springgreen3")

par(mar=c(5,5,1,1))
plot(ss_scores[,2]~ss_scores[,1], xlim=c(-1,1), ylim=c(-0.5,0.5), cex=3, lwd=2, xlab="NMDS1 (15.6%)", ylab="NMDS2 (10.8%)",
     pch=c(21,21,21,21,22,22,22,22,21,21,21,21,22,22,22,22,21,21,21,21,22,22,22,22),
     col="black", bg=ss_bi_col,
     axes="F", cex.lab=1.5 )
axis(1, lwd=3, lwd.ticks=3, cex.axis=1.5)
axis(2, lwd=3, lwd.ticks=3, cex.axis=1.5)
legend("topright", bty="n", cex=1.5, title="Sus.Slge. RNA", legend=c("R1","R3","R4","R5","6 mg/L","4 mg/L","2 mg/L"), pch=c(22,22,21,21,21,21,21), pt.bg=c("lightgray","darkgray","lightgray","darkgray","springgreen","skyblue","plum1") )
legend("topleft", bty="n", cex=1.5, legend="c)")
arrows(x0=0,y0=0, x1=ss_en_coord_cont[,1], y1=ss_en_coord_cont[,2], lwd=2)
#arrows(x0=0,y0=0, x1=ss_en_coord_cat[,1], y1=ss_en_coord_cat[,2], lwd=2)
#text(x=ss_en_coord_cat[,1],y=ss_en_coord_cat[,2],labels=row.names(ss_en_coord_cat), cex=1, col="purple4")
text(x=ss_en_coord_cont[,1]*1.05,y=ss_en_coord_cont[,2]*1.05,labels=row.names(ss_en_coord_cont), cex=1, col="purple4")