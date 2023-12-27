library(readxl)
library(shape)
layout(matrix(c(1,1,2,2,3,3,4), ncol=7,nrow=1))
par(mar=c(5,5,2,2), xpd=FALSE)
amohzo1 <- read_excel("~/Pinto/DO Transcriptomics Project Data/amoa cmx to hzo.xlsx", sheet = "total")
amohzo <-as.matrix(amohzo1[1:4,2:10])

plot(amohzo[2,7:9]~amohzo[1,7:9], bg="blue2", pch=21, cex=2, cex.lab=1.5, axes="F", ylim=c(2,5), xlab="", ylab="",
     xlim=c(4,6))
##plus 3 to account for uL to mL
axis(2, cex.axis=1.5, at=c(2,3,4,5), labels=c("5","6","7","8"))
axis(1, cex.axis=1.5, at=c(4,5,6), labels=c("7","8","9"))
abline(lm(amohzo[2,7:9]~amohzo[1,7:9]), lwd=2, lty=2, col="blue2")
legend("topleft", legend="a)", cex=1.5, bty="n", pch=26)
legend("top", title="hzo to amoA CMX", legend="amoA CMX", cex=1.5, bty="n", pch=21, pt.bg="blue2")

mtext("amoA CMX RT-qPCR [log(copies/mL)]", side=2, cex=1.5, padj=-2 ,col="blue4")
mtext("hzo AMX RT-qPCR [log(copies/mL)]", side=1, cex=1.5, padj=2 ,col="red4")
arrows(x0=amohzo[1,7:9], y0=amohzo[2,7:9], y1=amohzo[2,7:9]+amohzo[4,7:9], angle=90, lwd=2, length=0.1)
arrows(x0=amohzo[1,7:9], y0=amohzo[2,7:9], y1=amohzo[2,7:9]-amohzo[4,7:9], angle=90, lwd=2, length=0.1)
arrows(x0=amohzo[1,7:9], y0=amohzo[2,7:9], x1=amohzo[1,7:9]+amohzo[3,7:9], angle=90, lwd=2, length=0.1)
arrows(x0=amohzo[1,7:9], y0=amohzo[2,7:9], x1=amohzo[1,7:9]-amohzo[3,7:9], angle=90, lwd=2, length=0.1)
points(amohzo[2,7:9]~amohzo[1,7:9], bg="blue2", cex=3, pch=21)




TPM_Nitro <- read_excel("~/Pinto/DO Transcriptomics Project Data/Transcriptomics_Data/Consolidated_Nitro_MAG_Reads.xlsx", sheet = "TPM_Nitro")





x <-rep(seq(from = 1, to =8, by =1),times=7)
y <-rep(c(7,6,5,4,3,2,1), each=8)
TPM_Nitro_CMX <-subset(TPM_Nitro, MAG == "CMX")
TPM_CMX <-TPM_Nitro_CMX[,3:34]
new_cols <-c("D2_S1","D6_S1","D2_S3","D6_S3","D2_S5","D6_S5","D2_I5","D6_I5")
TPM_CMX_avg <-sapply(seq(1, ncol(TPM_CMX), 4), function(j) rowMeans(TPM_CMX[, j+(0:3)])) %>% 
  as.data.frame() %>% 
  setNames(new_cols)
y_names <-as.vector(t(TPM_Nitro_CMX[,1]))
TPM_CMX_cex <-log(TPM_CMX_avg, base=10)
m = as.matrix(TPM_CMX_cex)
q = c()
for (i in seq(1:nrow(m))){q = c(q, m[i,])}
TPM_CMX_SD <-sapply(split.default(TPM_CMX, rep(seq((ncol(TPM_CMX) / 4)), each = 4)), function(i)
  +     apply(i, 1, sd, na.rm = TRUE))
TPM_CMX_USD <-TPM_CMX_SD + TPM_CMX_avg
TPM_CMX_SD_cex <-log(TPM_CMX_USD, base=10)
m2 = as.matrix(TPM_CMX_SD_cex)
q2 = c()
for (i in seq(1:nrow(m2))){q2 = c(q2, m2[i,])}
blue50 <- rgb(0, 0, 139, max = 255, alpha = 125, names = "blue50")      
par(mar=c(5,8,1,1))
plot(y ~ x, cex=(q-1.8)^(2), axes=F, bg="blue4", pch=21, ylab="", xlab="", cex.lab=1.5, xlim=c(1,8), ylim=c(0,7.5))
points(y ~ x, cex=(q2-1.8)^(2), pch=16, col=blue50)
abline(v=2.5, lwd=2, lty=2)
abline(v=4.5, lwd=2, lty=2)
abline(v=6.5, lwd=2, lty=2)
axis(2, cex=1.5, labels=y_names, at=c(7,6,5,4,3,2,1), cex.axis=1.5, las=2)
axis(1, at=c(1,2), label=c("DO=2mg/L","DO=6mg/L"), cex.axis=1)
axis(1, at=1.5, label=c("SS R1"), cex.axis=1.5, padj = 2, lwd=0)
axis(1, at=c(3,4), label=c("DO=2mg/L","DO=6mg/L"), cex.axis=1)
axis(1, at=3.5, label=c("SS R3"), cex.axis=1.5, padj = 2, lwd=0)
axis(1, at=c(5,6), label=c("DO=2mg/L","DO=6mg/L"), cex.axis=1)
axis(1, at=5.5, label=c("SS R4"), cex.axis=1.5, padj = 2, lwd=0)
axis(1, at=c(7,8), label=c("DO=2mg/L","DO=6mg/L"), cex.axis=1)
axis(1, at=7.5, label=c("IFAS R4"), cex.axis=1.5, padj = 2, lwd=0)
mtext("Comammox Nitrospira", side=2, padj=-3.3, cex=1.5)
legend("topleft", legend="b)", cex=1.5, bty="n", pch=26)


CMX1 <-(TPM_CMX_cex[,2]-TPM_CMX_cex[,1])/2
Arrows(x0=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5),x=1.5+CMX1, y0=c(7,6,5,4,3,2,1),y1=c(7,6,5,4,3,2,1), col="skyblue2", arr.type="T", lwd=2)
axis(1, at=c(1,1.5,2), label=c("10x","=","10x"),padj=-4, cex.axis=1, tck=0.01)
axis(1, at=c(1.5), label="log diff.", cex.axis=1.5, padj=-4.5)

CMX2 <-(TPM_CMX_cex[,4]-TPM_CMX_cex[,3])/4
Arrows(x0=c(3.5,3.5,3.5,3.5,3.5,3.5,3.5),x=3.5+CMX2, y0=c(7,6,5,4,3,2,1),y1=c(7,6,5,4,3,2,1), col="skyblue2", arr.type="T", lwd=2)
axis(1, at=c(3,3.5,4), label=c("10x","=","10x"),padj=-4, cex.axis=1, tck=0.01)
axis(1, at=c(3.5), label="log diff.", cex.axis=1.5, padj=-4.5)

CMX3 <-(TPM_CMX_cex[,6]-TPM_CMX_cex[,5])/4
Arrows(x0=c(5.5,5.5,5.5,5.5,5.5,5.5,5.5),x=5.5+CMX3, y0=c(7,6,5,4,3,2,1),y1=c(7,6,5,4,3,2,1), col="skyblue2", arr.type="T", lwd=2)
axis(1, at=c(5,5.5,6), label=c("10x","=","10x"),padj=-4, cex.axis=1, tck=0.01)
axis(1, at=c(5.5), label="log diff.", cex.axis=1.5, padj=-4.5)

CMX4 <-(TPM_CMX_cex[,8]-TPM_CMX_cex[,7])/4
Arrows(x0=c(7.5,7.5,7.5,7.5,7.5,7.5,7.5),x=7.5+CMX4, y0=c(7,6,5,4,3,2,1),y1=c(7,6,5,4,3,2,1), col="skyblue2", arr.type="T", lwd=2)
axis(1, at=c(7,7.5,8), label=c("10x","=","10x"),padj=-4, cex.axis=1, tck=0.01)
axis(1, at=c(7.5), label="log diff.", cex.axis=1.5, padj=-4.5)




x <-rep(seq(from = 1, to =8, by =1),times=7)
y <-rep(c(7,6,5,4,3,2,1), each=8)
TPM_Nitro_AMX <-subset(TPM_Nitro, MAG == "AMX")
TPM_AMX <-TPM_Nitro_AMX[,3:34]
new_cols <-c("D2_S1","D6_S1","D2_S3","D6_S3","D2_S5","D6_S5","D2_I5","D6_I5")
TPM_AMX_avg <-sapply(seq(1, ncol(TPM_AMX), 4), function(j) rowMeans(TPM_AMX[, j+(0:3)])) %>% 
  as.data.frame() %>% 
  setNames(new_cols)
y_names <-as.vector(t(TPM_Nitro_AMX[,1]))
TPM_AMX_cex <-log(TPM_AMX_avg, base=10)
m = as.matrix(TPM_AMX_cex)
q = c()
for (i in seq(1:nrow(m))){q = c(q, m[i,])}
TPM_AMX_SD <-sapply(split.default(TPM_AMX, rep(seq((ncol(TPM_AMX) / 4)), each = 4)), function(i)
  +     apply(i, 1, sd, na.rm = TRUE))
TPM_AMX_USD <-TPM_AMX_SD + TPM_AMX_avg
TPM_AMX_SD_cex <-log(TPM_AMX_USD, base=10)
m2 = as.matrix(TPM_AMX_SD_cex)
q2 = c()
for (i in seq(1:nrow(m2))){q2 = c(q2, m2[i,])}
red50 <- rgb(139, 0, 0, max = 255, alpha = 125, names = "red50")      
par(mar=c(5,8,1,1))
plot(y ~ x, cex=(q-1.8)^(2), axes=F, bg="red4", pch=21, ylab="", xlab="", cex.lab=1.5, xlim=c(1,8), ylim=c(0,7.5))
points(y ~ x, cex=(q2-1.8)^(2), pch=16, col=red50)
abline(v=2.5, lwd=2, lty=2)
abline(v=4.5, lwd=2, lty=2)
abline(v=6.5, lwd=2, lty=2)
axis(2, cex=1.5, labels=y_names, at=c(7,6,5,4,3,2,1), cex.axis=1.5, las=2)
axis(1, at=c(1,2), label=c("DO=2mg/L","DO=6mg/L"), cex.axis=1)
axis(1, at=1.5, label=c("SS R1"), cex.axis=1.5, padj = 2, lwd=0)
axis(1, at=c(3,4), label=c("DO=2mg/L","DO=6mg/L"), cex.axis=1)
axis(1, at=3.5, label=c("SS R3"), cex.axis=1.5, padj = 2, lwd=0)
axis(1, at=c(5,6), label=c("DO=2mg/L","DO=6mg/L"), cex.axis=1)
axis(1, at=5.5, label=c("SS R4"), cex.axis=1.5, padj = 2, lwd=0)
axis(1, at=c(7,8), label=c("DO=2mg/L","DO=6mg/L"), cex.axis=1)
axis(1, at=7.5, label=c("IFAS R4"), cex.axis=1.5, padj = 2, lwd=0)
mtext("Anammox Bacteria", side=2, padj=-3.3, cex=1.5)
legend("topleft", legend="c)", cex=1.5, bty="n", pch=26)


AMX1 <-(TPM_AMX_cex[,2]-TPM_AMX_cex[,1])/2
Arrows(x0=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5),x=1.5+AMX1, y0=c(7,6,5,4,3,2,1),y1=c(7,6,5,4,3,2,1), col="hotpink2", arr.type="T", lwd=2)
axis(1, at=c(1,1.5,2), label=c("10x","=","10x"),padj=-4, cex.axis=1, tck=0.01)
axis(1, at=c(1.5), label="log diff.", cex.axis=1.5, padj=-4.5)

AMX2 <-(TPM_AMX_cex[,4]-TPM_AMX_cex[,3])/4
Arrows(x0=c(3.5,3.5,3.5,3.5,3.5,3.5,3.5),x=3.5+AMX2, y0=c(7,6,5,4,3,2,1),y1=c(7,6,5,4,3,2,1), col="hotpink2", arr.type="T", lwd=2)
axis(1, at=c(3,3.5,4), label=c("10x","=","10x"),padj=-4, cex.axis=1, tck=0.01)
axis(1, at=c(3.5), label="log diff.", cex.axis=1.5, padj=-4.5)

AMX3 <-(TPM_AMX_cex[,6]-TPM_AMX_cex[,5])/4
Arrows(x0=c(5.5,5.5,5.5,5.5,5.5,5.5,5.5),x=5.5+AMX3, y0=c(7,6,5,4,3,2,1),y1=c(7,6,5,4,3,2,1), col="hotpink2", arr.type="T", lwd=2)
axis(1, at=c(5,5.5,6), label=c("10x","=","10x"),padj=-4, cex.axis=1, tck=0.01)
axis(1, at=c(5.5), label="log diff.", cex.axis=1.5, padj=-4.5)

AMX4 <-(TPM_AMX_cex[,8]-TPM_AMX_cex[,7])/4
Arrows(x0=c(7.5,7.5,7.5,7.5,7.5,7.5,7.5),x=7.5+AMX4, y0=c(7,6,5,4,3,2,1),y1=c(7,6,5,4,3,2,1), col="hotpink2", arr.type="T", lwd=2)
axis(1, at=c(7,7.5,8), label=c("10x","=","10x"),padj=-4, cex.axis=1, tck=0.01)
axis(1, at=c(7.5), label="log diff.", cex.axis=1.5, padj=-4.5)


par(mar=c(0,0,0,0))
plot.new()
legend("right", legend=c("Average TPM", "St. Dev.","","100","","1,000","","10,000","", "100,000"), col=c("black","lightgray", "black","black","black","black","black","black","black","black"), pch=26, , cex=1.5, pt.cex=c(3,3,0,0.04,0,1.44,0,4.84,0,10.24), bty="n")
legend("center", legend=c("", "","","","","","","","", ""), col=c("black","lightgray", "black","black","black","black","black","black","black","black"), pch=16, , cex=1.5, pt.cex=c(3,3,0,0.04,0,1.44,0,4.84,0,10.24), bty="n")

