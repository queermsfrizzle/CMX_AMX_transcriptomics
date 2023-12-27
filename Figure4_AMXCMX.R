library(readxl)
des <- read_excel("~/Pinto/DO Transcriptomics Project Data/deseq_diff_expression.xlsx", sheet = "Plotting2")



layout(matrix(c(1,2,3), ncol=3,nrow=1))



amx_mags <- read_excel("~/Pinto/DO Transcriptomics Project Data/volcano plots amx.xlsx", 
                       sheet = "AMX_MAG_IFAS_TOTAL_TPM")
amx_mags2 <-as.matrix(t(amx_mags[,2:3]))

par(mar=c(5,6,1,1))
bp <-barplot(amx_mags2, beside=TRUE, col=c("plum2","springgreen2"), axes="F", ylim=c(0,1000000), lwd=2)
amx_magssd <-as.matrix(t(amx_mags[,4:5]))
arrows(x0=bp, y0=amx_mags2, y1=(amx_mags2+amx_magssd), lwd=2, angle=90, length=0.1)
arrows(x0=bp, y0=amx_mags2, y1=(amx_mags2-amx_magssd), lwd=2, angle=90, length=0.1)
axis(1, at=c(2,5,8,11), labels=c("AMX_66","AMX_281","AMX357","AMX_465"), cex.axis=1.5)
legend("topleft", legend="a)", cex=1.5, bty="n", pch=26)
axis(2, at=c(0,200000,400000,600000,800000,1000000), labels=c("0","200k","400k","600k","800k","1,000k"), cex.axis=1.5, las=2)
axis(2,at=500000,labels="Transcripts Per Million Reads", cex.axis=1.5, padj=-4)
legend("topright", title="TPM in IFAS", legend=c("DO 2 mg/L", "DO 6 mg/L"),cex=1.5, bty="n",pch=22, pt.bg=c("plum2","springgreen2"))







par(mar=c(5,5,1,1))
plot(y=des$amxipval,x=des$amxidiff, pch=21, bg="black",col="black", cex=1.5, axes=F, 
     ylab="-log(p-value)",xlab="log(Differential TPM Expression)", cex.lab=1.5,
     xlim=c(-4,4), ylim=c(0,13))

axis(1, at=c(-4,-3,-2,-1,0,1,2,3,4), labels=c("16x","8x","4x","2x","0","2x","4x","8x","16x"), cex.axis=1.5)
axis(1, at=c(-2,2), labels=c("DO 2 mg/L", "DO 6 mg/L"), cex.axis=1.5, padj=1)
axis(2, at=c(0,3,6,9,12), labels=c(0,3,6,9,12), cex.axis=1.5)

#log10(0.05) = 1.30103
#log10(0.01) = 2
volhigh <-subset(des, amxipval >1.30103)
points(volhigh$amxipval~volhigh$amxidiff, bg="gray", pch=21, cex=1.5)

voldo2 <-subset(volhigh, amxidiff  < (-1))
points(voldo2$amxipval~voldo2$amxidiff, bg="plum2", pch=21, cex=1.5)
voldo3 <-subset(volhigh, amxidiff  < (-2))
points(voldo3$amxipval~voldo3$amxidiff, bg="plum3", pch=21, cex=1.5)
voldo4 <-subset(volhigh, amxidiff  < (-3))
points(voldo4$amxipval~voldo4$amxidiff, bg="plum4", pch=21, cex=1.5)

voldo4 <-subset(volhigh, amxidiff  > (1))
points(voldo4$amxipval~voldo4$amxidiff, bg="springgreen1", pch=21, cex=1.5)
voldo5 <-subset(volhigh, amxidiff  > (2))
points(voldo5$amxipval~voldo5$amxidiff, bg="springgreen3", pch=21, cex=1.5)
voldo6 <-subset(volhigh, amxidiff  > (3))
points(voldo6$amxipval~voldo6$amxidiff, bg="springgreen4", pch=21, cex=1.5)
abline(v=(-log10(2)), lwd=2, lty=2)
abline(v=(log10(2)), lwd=2, lty=2)
abline(h=(-log10(0.05)), lwd=2, lty=2)
legend("topleft", legend="b)", cex=1.5, bty="n", pch=26)
legend("left", title="AMX_281 in IFAS", legend=c("DO 2 mg/L","DO 6 mg/L","Diff Exp. < 2x","p > 0.05","Flagella","N-Cycle" ), pch=21, pt.bg=c("plum2","springgreen2","gray","black","orange","blue"), bty="n", cex=1.5)

des2 <- read_excel("~/Pinto/DO Transcriptomics Project Data/deseq_diff_expression.xlsx", sheet = "Plotting3")

points(x=des2$amxifdiff,y=des2$amxifpval, pch=21, cex=3, col="orange", lwd=3)
points(x=des2$amxindiff,y=des2$amxinpval, pch=21, cex=3, col="blue", lwd=3)








par(mar=c(5,5,1,1))
plot(des$amxspval~des$amxsdiff, pch=21, bg="black",col="black", cex=1.5, axes=F, 
     ylab="-log(p-value)",xlab="log(Differential TPM Expression)", cex.lab=1.5,
     xlim=c(-4,4), ylim=c(0,13))

axis(1, at=c(-4,-3,-2,-1,0,1,2,3,4), labels=c("16x","8x","4x","2x","0","2x","4x","8x","16x"), cex.axis=1.5)
axis(1, at=c(-2,2), labels=c("DO 2 mg/L", "DO 6 mg/L"), cex.axis=1.5, padj=1)
axis(2, at=c(0,3,6,9,12), labels=c(0,3,6,9,12), cex.axis=1.5)

#log10(0.05) = 1.30103
#log10(0.01) = 2
volhigh <-subset(des, amxspval >1.30103)
points(volhigh$amxspval~volhigh$amxsdiff, bg="gray", pch=21, cex=1.5)

voldo2 <-subset(volhigh, amxsdiff  < (-1))
points(voldo2$amxspval~voldo2$amxsdiff, bg="plum2", pch=21, cex=1.5)
voldo3 <-subset(volhigh, amxsdiff  < (-2))
points(voldo3$amxspval~voldo3$amxsdiff, bg="plum3", pch=21, cex=1.5)
voldo4 <-subset(volhigh, amxsdiff  < (-3))
points(voldo4$amxspval~voldo4$amxsdiff, bg="plum4", pch=21, cex=1.5)

voldo4 <-subset(volhigh, amxsdiff  > (1))
points(voldo4$amxspval~voldo4$amxsdiff, bg="springgreen1", pch=21, cex=1.5)
voldo5 <-subset(volhigh, amxsdiff  > (2))
points(voldo5$amxspval~voldo5$amxsdiff, bg="springgreen3", pch=21, cex=1.5)
voldo6 <-subset(volhigh, amxsdiff  > (3))
points(voldo6$amxspval~voldo6$amxsdiff, bg="springgreen4", pch=21, cex=1.5)
abline(v=(-log10(2)), lwd=2, lty=2)
abline(v=(log10(2)), lwd=2, lty=2)
abline(h=(-log10(0.05)), lwd=2, lty=2)
legend("topleft", legend="c)", cex=1.5, bty="n", pch=26)
legend("left", title="AMX_281 in SS", legend=c("DO 2 mg/L","DO 6 mg/L","Diff Exp. < 2x","p > 0.05","Flagella","N-Cycle" ), pch=21, pt.bg=c("plum2","springgreen2","gray","black","orange","blue"), bty="n", cex=1.5)

points(x=des2$amxsfdiff,y=des2$amxsfpval, pch=21, cex=3, col="orange", lwd=3)
points(x=des2$amxsndiff,y=des2$amxsnpval, pch=21, cex=3, col="blue", lwd=3)