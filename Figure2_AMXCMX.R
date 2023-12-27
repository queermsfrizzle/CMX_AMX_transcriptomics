library(readxl)
aob16s <- read_excel("~/Pinto/AJP10/16S_OTU/newest/16S_nitro.xlsx", sheet = "Monas_OTU")
layout(matrix(c(1,4,1,4,8,5,2,5,2,6,8,6,3,7,3,7), ncol=8,nrow=2))


amx16s <- read_excel("~/Pinto/AJP10/16S_OTU/newest/16S_nitro.xlsx", sheet = "AMX_OTU")
par(mar=c(5,5,2,2))
amx16savg <-as.matrix(amx16s[1:4,2:10])
amxbp <-barplot(amx16savg, col=c("chocolate1","bisque","darkgoldenrod1","red3"),
                ylab="16S rRNA Transcript Compositon [%]",cex.lab=1.5, axes="F", lwd=2, axisnames=FALSE, ylim=c(0,0.2))
legend("topleft", legend="a)", cex=1.5, bty="n", pch=26)
axis(2, at=c(0,0.05,0.1,0.15,0.2), labels=c("0%","5%","10%","15%","20%"), cex.axis=1.5)
axis(1, at=amxbp[1:3], label=c("2","4","6"), cex.axis=1.5, lwd=2)                      
axis(1, at=amxbp[4:6], label=c("2","4","6"), cex.axis=1.5, lwd=2)                      
axis(1, at=amxbp[7:9], label=c("2","4","6"), cex.axis=1.5, lwd=2)                      
axis(1, at=amxbp[2], label="SS", cex.axis=1.5, padj=1.5)
axis(1, at=amxbp[2], label="R1 & R3", cex.axis=1.5, padj=3)
axis(1, at=amxbp[5], label="SS", cex.axis=1.5, padj=1.5)
axis(1, at=amxbp[5], label="R4 & R5", cex.axis=1.5, padj=3)
axis(1, at=amxbp[8], label="IFAS", cex.axis=1.5, padj=1.5)
axis(1, at=amxbp[8], label="R4 and R5", cex.axis=1.5, padj=3)
abline(v=c(3.7,7.3), lwd=4, lty=2, col="gray")
legend("left", title="Anammox", legend=t(amx16s[1:4,1]), pt.bg=c("chocolate1","bisque","darkgoldenrod1","red3"), pch=22, cex=1.5, bty="n")

axis(4, at=0.05, label="5%", cex.axis=1.5)


nob16s <- read_excel("~/Pinto/AJP10/16S_OTU/newest/16S_nitro.xlsx", sheet = "Spira_OTU")

par(mar=c(5,5,2,2))
nob16savg <-as.matrix(nob16s[1:3,2:10])
nobbp <-barplot(nob16savg, col=c("mediumpurple","mediumpurple4","darkorchid1"),
                ylab="16S rRNA Transcript Compositon [%]",cex.lab=1.5, axes="F", lwd=2, axisnames=FALSE, ylim=c(0,0.05))
legend("topleft", legend="b)", cex=1.5, bty="n", pch=26)
axis(2, at=c(0,0.01,0.02,0.03,0.04,0.05), labels=c("0%","1%","2%","3%","4%","5%"), cex.axis=1.5)
axis(1, at=nobbp[1:3], label=c("2","4","6"), cex.axis=1.5, lwd=2)
axis(1, at=nobbp[4:6], label=c("2","4","6"), cex.axis=1.5, lwd=2)                      
axis(1, at=nobbp[7:9], label=c("2","4","6"), cex.axis=1.5, lwd=2)                      
axis(1, at=nobbp[2], label="SS", cex.axis=1.5, padj=1.5)
axis(1, at=nobbp[2], label="R1 & R3", cex.axis=1.5, padj=3)
axis(1, at=nobbp[5], label="SS", cex.axis=1.5, padj=1.5)
axis(1, at=nobbp[5], label="R4 & R5", cex.axis=1.5, padj=3)
axis(1, at=nobbp[8], label="IFAS", cex.axis=1.5, padj=1.5)
axis(1, at=nobbp[8], label="R4 and R5", cex.axis=1.5, padj=3)
abline(v=c(3.7,7.3), lwd=4, lty=2, col="gray")
legend("left", title="NOB",legend=t(nob16s[1:3,1]), pt.bg=c("mediumpurple","mediumpurple4","darkorchid"), pch=22, cex=1.5, bty="n")
axis(4, at=0.0015, label="0.15%", cex.axis=1.5)


par(mar=c(5,5,2,2))
aob16savg <-as.matrix(aob16s[1:2,2:10])
aobbp <-barplot(aob16savg, col=c("green2","forestgreen"),
                ylab="16S rRNA Transcript Compositon [%]",cex.lab=1.5, axes="F", lwd=2, axisnames=FALSE, ylim=c(0,0.0015))
legend("topleft", legend="c)", cex=1.5, bty="n", pch=26)
axis(2, at=c(0,0.0005,0.001,0.0015), labels=c("0%","0.05%","0.10%","0.15%"), cex.axis=1.5)
axis(1, at=aobbp[1:3], label=c("2","4","6"), cex.axis=1.5, lwd=2)
axis(1, at=aobbp[4:6], label=c("2","4","6"), cex.axis=1.5, lwd=2)                      
axis(1, at=aobbp[7:9], label=c("2","4","6"), cex.axis=1.5, lwd=2)                      
axis(1, at=aobbp[2], label="SS", cex.axis=1.5, padj=1.5)
axis(1, at=aobbp[2], label="R1 & R3", cex.axis=1.5, padj=3)
axis(1, at=aobbp[5], label="SS", cex.axis=1.5, padj=1.5)
axis(1, at=aobbp[5], label="R4 & R5", cex.axis=1.5, padj=3)
axis(1, at=aobbp[8], label="IFAS", cex.axis=1.5, padj=1.5)
axis(1, at=aobbp[8], label="R4 and R5", cex.axis=1.5, padj=3)
abline(v=c(3.7,7.3), lwd=4, lty=2, col="gray")
legend("left", title="AOB",legend=t(aob16s[1:2,1]), pt.bg=c("green2","forestgreen"), pch=22, cex=1.5, bty="n")






amx_rt <- read_excel("~/Pinto/AJP10/genes and otus quantified .xlsx",sheet = "AMX_Plot")

par(mar=c(5,5,2,1))
amx_rt_data <-as.matrix(amx_rt[3:8,2:10])
amxbp <-barplot(amx_rt_data, col=c("chocolate1","red3","darkgoldenrod1","bisque","coral4","pink2"),
                ylab="hzo RT-qPCR Copies [log(copies/mL)]",cex.lab=1.5, axes="F", lwd=2, axisnames=FALSE, ylim=c(0,7), xlim=c(0,13))
legend("topleft", legend="d)", cex=1.5, bty="n", pch=26)

axis(2, cex.axis=1.5, at=atc, labels=byc)

axis(1, at=amxbp[1:3], label=c("2","4","6"), cex.axis=1.5, lwd=2)                      
axis(1, at=amxbp[4:6], label=c("2","4","6"), cex.axis=1.5, lwd=2)                      
axis(1, at=amxbp[7:9], label=c("2","4","6"), cex.axis=1.5, lwd=2)                      
axis(1, at=amxbp[2], label="SS", cex.axis=1.5, padj=1.5)
axis(1, at=amxbp[2], label="R1 & R3", cex.axis=1.5, padj=3)
axis(1, at=amxbp[5], label="SS", cex.axis=1.5, padj=1.5)
axis(1, at=amxbp[5], label="R4 & R5", cex.axis=1.5, padj=3)
axis(1, at=amxbp[8], label="IFAS", cex.axis=1.5, padj=1.5)
axis(1, at=amxbp[8], label="R4 and R5", cex.axis=1.5, padj=3)
abline(v=c(3.7,7.3), lwd=4, lty=2, col="gray")

amx_rt_dna <-as.matrix(amx_rt[1,12:13])
segments(x0=0,x1=7.2, y0=amx_rt_dna[1], lwd=3, col="black", lty=2)
segments(x0=7.3,x1=12, y0=amx_rt_dna[2], lwd=3, col="black", lty=2)

legend("topright", legend=t(amx_rt[3:8,1]), pt.bg=c("chocolate1","red3","darkgoldenrod1","bisque","coral4","pink2"), pch=22, cex=1.5, bty="n")







nxr_rt <- read_excel("~/Pinto/AJP10/genes and otus quantified .xlsx",sheet = "NXR_Plot2")

par(mar=c(5,5,2,1))
nxr_rt_data <-as.matrix(nxr_rt[3:9,2:10])
nxrbp <-barplot(nxr_rt_data, col=c("violetred","orchid","lightslateblue","deeppink","magenta4","blue","mediumorchid1"),
                ylab="nxrA RT-qPCR Copies [log(copies/mL)]",cex.lab=1.5, axes="F", lwd=2, axisnames=FALSE, ylim=c(0,7), xlim=c(0,13))
legend("topleft", legend="e)", cex=1.5, bty="n", pch=26)

axis(2, cex.axis=1.5, at=atc, labels=byc)

axis(1, at=nxrbp[1:3], label=c("2","4","6"), cex.axis=1.5, lwd=2)                      
axis(1, at=nxrbp[4:6], label=c("2","4","6"), cex.axis=1.5, lwd=2)                      
axis(1, at=nxrbp[7:9], label=c("2","4","6"), cex.axis=1.5, lwd=2)                      
axis(1, at=nxrbp[2], label="SS", cex.axis=1.5, padj=1.5)
axis(1, at=nxrbp[2], label="R1 & R3", cex.axis=1.5, padj=3)
axis(1, at=nxrbp[5], label="SS", cex.axis=1.5, padj=1.5)
axis(1, at=nxrbp[5], label="R4 & R5", cex.axis=1.5, padj=3)
axis(1, at=nxrbp[8], label="IFAS", cex.axis=1.5, padj=1.5)
axis(1, at=nxrbp[8], label="R4 and R5", cex.axis=1.5, padj=3)
abline(v=c(3.7,7.3), lwd=4, lty=2, col="gray")

nxr_rt_dna <-as.matrix(nxr_rt[1,12:13])
segments(x0=0,x1=7.2, y0=nxr_rt_dna[1], lwd=3, col="black", lty=2)
segments(x0=7.3,x1=12, y0=nxr_rt_dna[2], lwd=3, col="black", lty=2)

legend("topright", legend=t(nxr_rt[3:9,1]), pt.bg=c("violetred","orchid","lightslateblue","deeppink","magenta4","blue","mediumorchid1"), pch=22, cex=1.5, bty="n")


cmx_rt <- read_excel("~/Pinto/AJP10/genes and otus quantified .xlsx",sheet = "CMX_Plot")

par(mar=c(5,5,2,1))
cmx_rt_data <-as.matrix(cmx_rt[3:8,2:10])
cmxbp <-barplot(cmx_rt_data, col=c("blue3","cornflowerblue","turquoise2","steelblue1","cyan3","lightskyblue3"),
                ylab="amoA CMX RT-qPCR Copies [log(copies/mL)]",cex.lab=1.5, axes="F", lwd=2, axisnames=FALSE, ylim=c(0,7), xlim=c(0,13))
legend("topleft", legend="f)", cex=1.5, bty="n", pch=26)

axis(2, cex.axis=1.5, at=atc, labels=byc)

axis(1, at=cmxbp[1:3], label=c("2","4","6"), cex.axis=1.5, lwd=2)                      
axis(1, at=cmxbp[4:6], label=c("2","4","6"), cex.axis=1.5, lwd=2)                      
axis(1, at=cmxbp[7:9], label=c("2","4","6"), cex.axis=1.5, lwd=2)                      
axis(1, at=cmxbp[2], label="SS", cex.axis=1.5, padj=1.5)
axis(1, at=cmxbp[2], label="R1 & R3", cex.axis=1.5, padj=3)
axis(1, at=cmxbp[5], label="SS", cex.axis=1.5, padj=1.5)
axis(1, at=cmxbp[5], label="R4 & R5", cex.axis=1.5, padj=3)
axis(1, at=cmxbp[8], label="IFAS", cex.axis=1.5, padj=1.5)
axis(1, at=cmxbp[8], label="R4 and R5", cex.axis=1.5, padj=3)
abline(v=c(3.7,7.3), lwd=4, lty=2, col="gray")

cmx_rt_dna <-as.matrix(cmx_rt[1,12:13])
segments(x0=0,x1=7.2, y0=cmx_rt_dna[1], lwd=3, col="black", lty=2)
segments(x0=7.3,x1=12, y0=cmx_rt_dna[2], lwd=3, col="black", lty=2)

legend("topright", legend=t(cmx_rt[3:8,1]), pt.bg=c("blue3","cornflowerblue","turquoise2","steelblue1","cyan3","lightskyblue3"), pch=22, cex=1.5, bty="n")








aob_rt <- read_excel("~/Pinto/AJP10/genes and otus quantified .xlsx",sheet = "AOB_Plot")

par(mar=c(5,5,2,1))
aob_rt_data <-as.matrix(aob_rt[3:8,2:10])
aobbp <-barplot(aob_rt_data, col=c("aquamarine2","aquamarine4","darkolivegreen","springgreen2","seagreen3","darkolivegreen2"),
                ylab="amoA AOB RT-qPCR Copies [log(copies/mL)]",cex.lab=1.5, axes="F", lwd=2, axisnames=FALSE, ylim=c(0,7), xlim=c(0,13))
legend("topleft", legend="g)", cex=1.5, bty="n", pch=26)

atc <-seq(from=0, to=7, by=1)
byc <-seq(from=3, to=10, by=1)
axis(2, cex.axis=1.5, at=atc, labels=byc)


axis(1, at=aobbp[1:3], label=c("2","4","6"), cex.axis=1.5, lwd=2)                      
axis(1, at=aobbp[4:6], label=c("2","4","6"), cex.axis=1.5, lwd=2)                      
axis(1, at=aobbp[7:9], label=c("2","4","6"), cex.axis=1.5, lwd=2)                      
axis(1, at=aobbp[2], label="SS", cex.axis=1.5, padj=1.5)
axis(1, at=aobbp[2], label="R1 & R3", cex.axis=1.5, padj=3)
axis(1, at=aobbp[5], label="SS", cex.axis=1.5, padj=1.5)
axis(1, at=aobbp[5], label="R4 & R5", cex.axis=1.5, padj=3)
axis(1, at=aobbp[8], label="IFAS", cex.axis=1.5, padj=1.5)
axis(1, at=aobbp[8], label="R4 and R5", cex.axis=1.5, padj=3)
abline(v=c(3.7,7.3), lwd=4, lty=2, col="gray")

aob_rt_dna <-as.matrix(aob_rt[1,12:13])
segments(x0=0,x1=7.2, y0=aob_rt_dna[1], lwd=3, col="black", lty=2)
segments(x0=7.3,x1=12, y0=aob_rt_dna[2], lwd=3, col="black", lty=2)

aob_legend <-c("DNA qPCR",t(aob_rt[3:8,1]))
legend("topright", legend=aob_legend, lwd=c(3,0,0,0,0,0,0),lty=c(2,0,0,0,0,0,0), pt.bg=c("black","aquamarine4","aquamarine2","darkolivegreen","springgreen2","seagreen3","darkolivegreen2"), pch=c(26,22,22,22,22,22,22), cex=1.5, bty="n")



