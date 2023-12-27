layout(matrix(c(1,3,2,4), ncol=2,nrow=2))

aob_rt <- read_excel("~/Pinto/AJP10/genes and otus quantified .xlsx",sheet = "AOB_Plot")

par(mar=c(5,5,2,5))
aob_rt_data <-as.matrix(aob_rt[3:8,2:10])
aobbp <-barplot(aob_rt_data, col=c("aquamarine2","aquamarine4","darkolivegreen","springgreen2","seagreen3","darkolivegreen2"),
                ylab="amoA AOB RT-qPCR Copies [log(copies/mL)]",cex.lab=1.5, axes="F", lwd=2, axisnames=FALSE, ylim=c(0,7))
legend("topleft", legend="a)", cex=1.5, bty="n", pch=26)

axis(2, cex.axis=1.5)

axis(4, at=c(0,7), labels=c("",""),cex.axis=1.5)
mtext("OTU Percent Abundance [%]", side=4, cex=1.5, padj=1)
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

legend("top", legend=t(aob_rt[3:8,1]), pt.bg=c("aquamarine4","aquamarine2","darkolivegreen","springgreen2","seagreen3","darkolivegreen2"), pch=22, cex=1.5, bty="n")
legend("topright",legend="DNA qPCR", lwd=2, cex=1.5, lty=2, bty="n")






cmx_rt <- read_excel("~/Pinto/AJP10/genes and otus quantified .xlsx",sheet = "CMX_Plot")

par(mar=c(5,5,2,5))
cmx_rt_data <-as.matrix(cmx_rt[3:8,2:10])
cmxbp <-barplot(cmx_rt_data, col=c("blue3","cornflowerblue","turquoise2","steelblue1","cyan3","lightskyblue3"),
                ylab="amoA CMX RT-qPCR Copies [log(copies/mL)]",cex.lab=1.5, axes="F", lwd=2, axisnames=FALSE, ylim=c(0,7))
legend("topleft", legend="b)", cex=1.5, bty="n", pch=26)

axis(2, cex.axis=1.5)

axis(4, at=c(0,7), labels=c("",""),cex.axis=1.5)
mtext("OTU Percent Abundance [%]", side=4, cex=1.5, padj=1)
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

legend("top", legend=t(cmx_rt[3:8,1]), pt.bg=c("blue3","cornflowerblue","turquoise2","steelblue1","cyan3","lightskyblue3"), pch=22, cex=1.5, bty="n")
legend("topright",legend="DNA qPCR", lwd=2, cex=1.5, lty=2, bty="n")







nxr_rt <- read_excel("~/Pinto/AJP10/genes and otus quantified .xlsx",sheet = "NXR_Plot")

par(mar=c(5,5,2,5))
nxr_rt_data <-as.matrix(nxr_rt[3:8,2:10])
nxrbp <-barplot(nxr_rt_data, col=c("violetred","orchid","lightslateblue","deeppink","magenta4","mediumorchid1"),
                ylab="nxrA RT-qPCR Copies [log(copies/mL)]",cex.lab=1.5, axes="F", lwd=2, axisnames=FALSE, ylim=c(0,7))
legend("topleft", legend="c)", cex=1.5, bty="n", pch=26)

axis(2, cex.axis=1.5)

axis(4, at=c(0,7), labels=c("",""),cex.axis=1.5)
mtext("OTU Percent Abundance [%]", side=4, cex=1.5, padj=1)
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

legend("top", legend=t(nxr_rt[3:8,1]), pt.bg=c("violetred","orchid","lightslateblue","deeppink","magenta4","mediumorchid1"), pch=22, cex=1.5, bty="n")
legend("topright",legend="DNA qPCR", lwd=2, cex=1.5, lty=2, bty="n")




amx_rt <- read_excel("~/Pinto/AJP10/genes and otus quantified .xlsx",sheet = "AMX_Plot")

par(mar=c(5,5,2,5))
amx_rt_data <-as.matrix(amx_rt[3:8,2:10])
amxbp <-barplot(amx_rt_data, col=c("chocolate1","red3","darkgoldenrod1","bisque","coral4","pink2"),
                ylab="hzo RT-qPCR Copies [log(copies/mL)]",cex.lab=1.5, axes="F", lwd=2, axisnames=FALSE, ylim=c(0,7))
legend("topleft", legend="d)", cex=1.5, bty="n", pch=26)

axis(2, cex.axis=1.5)

axis(4, at=c(0,7), labels=c("",""),cex.axis=1.5)
mtext("OTU Percent Abundance [%]", side=4, cex=1.5, padj=1)
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

legend("top", legend=t(amx_rt[3:8,1]), pt.bg=c("chocolate1","red3","darkgoldenrod1","bisque","coral4","pink2"), pch=22, cex=1.5, bty="n")
legend("topright",legend="DNA qPCR", lwd=2, cex=1.5, lty=2, bty="n")



