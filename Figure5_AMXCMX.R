library(readxl)
des <- read_excel("~/Pinto/DO Transcriptomics Project Data/deseq_diff_expression.xlsx", sheet = "Plotting2")



layout(matrix(c(1,2,3,4), ncol=4,nrow=1))


cmx_mags <- read_excel("~/Pinto/DO Transcriptomics Project Data/CMX_MAGS_Volc.xlsx", 
                       sheet = "MAG Plot")

cmx_mags2 <-as.matrix(t(cmx_mags[,2:3]))

par(mar=c(5,6,1,1))
bp <-barplot(cmx_mags2, beside=TRUE, col=c("plum2","springgreen2"), axes="F", ylim=c(0,1000000), lwd=2)
cmx_magssd <-as.matrix(t(cmx_mags[,4:5]))
arrows(x0=bp, y0=cmx_mags2, y1=(cmx_mags2+cmx_magssd), lwd=2, angle=90, length=0.1)
arrows(x0=bp, y0=cmx_mags2, y1=(cmx_mags2-cmx_magssd), lwd=2, angle=90, length=0.1)
axis(1, at=c(2,5), labels=c("CMX_1","CMX_2"), cex.axis=1.5)
legend("topleft", legend="a)", cex=1.5, bty="n", pch=26)
axis(2, at=c(0,200000,400000,600000,800000,1000000), labels=c("0","200k","400k","600k","800k","1,000k"), cex.axis=1.5, las=2)
axis(2,at=500000,labels="Transcripts Per Million Reads", cex.axis=1.5, padj=-4)
legend("topright", title="TPM in IFAS",legend=c("DO 2 mg/L", "DO 6 mg/L"),cex=1.5, bty="n",pch=22, pt.bg=c("plum2","springgreen2"))



nxr_cmx1 <- read_excel("~/Pinto/DO Transcriptomics Project Data/nxr expression in cmx.xlsx", sheet = "plot")



nxr_cmx2 <-as.matrix(t(nxr_cmx1[,2:3]))

par(mar=c(5,6,1,1))
bp <-barplot(nxr_cmx2, beside=TRUE, col=c("plum2","springgreen2"), axes="F", ylim=c(0,150000), lwd=2)
nxr_cmxsd <-as.matrix(t(nxr_cmx1[,4:5]))
arrows(x0=bp, y0=nxr_cmx2, y1=(nxr_cmx2+nxr_cmxsd), lwd=2, angle=90, length=0.1)
arrows(x0=bp, y0=nxr_cmx2, y1=(nxr_cmx2-nxr_cmxsd), lwd=2, angle=90, length=0.1)
axis(1, at=c(2,5), labels=c("nxrA 1","nxrA 2"), cex.axis=1.5)
legend("topleft", legend="b)", cex=1.5, bty="n", pch=26)
axis(2, at=c(0,50000,100000,150000), labels=c("0","50k","100k","150k"), cex.axis=1.5, las=2)
axis(2,at=500000,labels="Transcripts Per Million Reads", cex.axis=1.5, padj=-4)
legend("topright", title="CMX1 TPM in IFAS",legend=c("DO 2 mg/L", "DO 6 mg/L"),cex=1.5, bty="n",pch=22, pt.bg=c("plum2","springgreen2"))




par(mar=c(5,5,1,1))
plot(des$cmx1pval~des$cmx1diff, pch=21, bg="black",col="black", cex=1.5, axes=F, 
     ylab="-log(p-value)",xlab="log(Differential TPM Expression)", cex.lab=1.5,
     xlim=c(-4,4), ylim=c(0,13))

axis(1, at=c(-4,-3,-2,-1,0,1,2,3,4), labels=c("16x","8x","4x","2x","0","2x","4x","8x","16x"), cex.axis=1.5)
axis(1, at=c(-2,2), labels=c("DO 2 mg/L", "DO 6 mg/L"), cex.axis=1.5, padj=1)
axis(2, at=c(0,3,6,9,12), labels=c(0,3,6,9,12), cex.axis=1.5)

#log10(0.05) = 1.30103
#log10(0.01) = 2
volhigh <-subset(des, cmx1pval >1.30103)
points(volhigh$cmx1pval~volhigh$cmx1diff, bg="gray", pch=21, cex=1.5)

voldo2 <-subset(volhigh, cmx1diff  < (-1))
points(voldo2$cmx1pval~voldo2$cmx1diff, bg="plum2", pch=21, cex=1.5)
voldo3 <-subset(volhigh, cmx1diff  < (-2))
points(voldo3$cmx1pval~voldo3$cmx1diff, bg="plum3", pch=21, cex=1.5)
voldo4 <-subset(volhigh, cmx1diff  < (-3))
points(voldo4$cmx1pval~voldo4$cmx1diff, bg="plum4", pch=21, cex=1.5)

voldo4 <-subset(volhigh, cmx1diff  > (1))
points(voldo4$cmx1pval~voldo4$cmx1diff, bg="springgreen1", pch=21, cex=1.5)
voldo5 <-subset(volhigh, cmx1diff  > (2))
points(voldo5$cmx1pval~voldo5$cmx1diff, bg="springgreen3", pch=21, cex=1.5)
voldo6 <-subset(volhigh, cmx1diff  > (3))
points(voldo6$cmx1pval~voldo6$cmx1diff, bg="springgreen4", pch=21, cex=1.5)
abline(v=(-log10(2)), lwd=2, lty=2)
abline(v=(log10(2)), lwd=2, lty=2)
abline(h=(-log10(0.05)), lwd=2, lty=2)
legend("topleft", legend="c)", cex=1.5, bty="n", pch=26)
legend("topright", title="CMX1 in IFAS", legend=c("DO 2 mg/L","DO 6 mg/L","Diff Exp. < 2x","p > 0.05","Amm. Ox.","Nitri. Ox." ), pch=21, pt.bg=c("plum2","springgreen2","gray","black","orange","blue"), bty="n", cex=1.5)

des2 <- read_excel("~/Pinto/DO Transcriptomics Project Data/deseq_diff_expression.xlsx", sheet = "Plotting3")

points(x=des2$cmx1adiff,y=des2$cmx1apval, pch=21, cex=3, col="orange", lwd=3)
points(x=des2$cmx1ndiff,y=des2$cmx1npval, pch=21, cex=3, col="blue", lwd=3)

points(x=des2$cmx1ndiff,y=des2$cmx1npval, pch=21, cex=3, col="blue", lwd=3)

#-3.682715493	12.62032851
#-3.220770349	1.467462776
text(x=-3, y=12.2, labels="nxrA_2", cex=1.5, col="blue")
text(x=-2.5, y=1.55, labels="nxrA_1", cex=1.5, col="blue")




par(mar=c(5,5,1,1))
plot(des$cmx2pval~des$cmx2diff, pch=21, bg="black",col="black", cex=1.5, axes=F, 
     ylab="-log(p-value)",xlab="log(Differential TPM Expression)", cex.lab=1.5,
     xlim=c(-4,4), ylim=c(0,13))

axis(1, at=c(-4,-3,-2,-1,0,1,2,3,4), labels=c("16x","8x","4x","2x","0","2x","4x","8x","16x"), cex.axis=1.5)
axis(1, at=c(-2,2), labels=c("DO 2 mg/L", "DO 6 mg/L"), cex.axis=1.5, padj=1)
axis(2, at=c(0,3,6,9,12), labels=c(0,3,6,9,12), cex.axis=1.5)

#log10(0.05) = 1.30103
#log10(0.01) = 2
volhigh <-subset(des, cmx2pval >1.30103)
points(volhigh$cmx2pval~volhigh$cmx2diff, bg="gray", pch=21, cex=1.5)

voldo2 <-subset(volhigh, cmx2diff  < (-1))
points(voldo2$cmx2pval~voldo2$cmx2diff, bg="plum2", pch=21, cex=1.5)
voldo3 <-subset(volhigh, cmx2diff  < (-2))
points(voldo3$cmx2pval~voldo3$cmx2diff, bg="plum3", pch=21, cex=1.5)
voldo4 <-subset(volhigh, cmx2diff  < (-3))
points(voldo4$cmx2pval~voldo4$cmx2diff, bg="plum4", pch=21, cex=1.5)

voldo4 <-subset(volhigh, cmx2diff  > (1))
points(voldo4$cmx2pval~voldo4$cmx2diff, bg="springgreen1", pch=21, cex=1.5)
voldo5 <-subset(volhigh, cmx2diff  > (2))
points(voldo5$cmx2pval~voldo5$cmx2diff, bg="springgreen3", pch=21, cex=1.5)
voldo6 <-subset(volhigh, cmx2diff  > (3))
points(voldo6$cmx2pval~voldo6$cmx2diff, bg="springgreen4", pch=21, cex=1.5)
abline(v=(-log10(2)), lwd=2, lty=2)
abline(v=(log10(2)), lwd=2, lty=2)
abline(h=(-log10(0.05)), lwd=2, lty=2)
legend("topleft", legend="d)", cex=1.5, bty="n", pch=26)
legend("topright", title="CMX2 in IFAS", legend=c("DO 2 mg/L","DO 6 mg/L","Diff Exp. < 2x","p > 0.05","Amm. Ox.","Nitri. Ox." ), pch=21, pt.bg=c("plum2","springgreen2","gray","black","orange","blue"), bty="n", cex=1.5)

points(x=des2$cmx2adiff,y=des2$cmx2apval, pch=21, cex=3, col="orange", lwd=3)
points(x=des2$cmx2ndiff,y=des2$cmx2npval, pch=21, cex=3, col="blue", lwd=3)







