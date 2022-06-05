#BiocManager::install("antiProfilesData")
#install.packages("corrplot")
library(antiProfilesData)
library(preprocessCore)
library(Hmisc)
library(corrplot)
library(grDevices)
#1

class(apColonData)
pdata <- pData(apColonData)
fdata <- fData(apColonData)
edata <- exprs(apColonData)

#1.a

pLen = ncol(pdata)
eLen = ncol(edata)

sapply(pdata, class)
for (i in 1:pLen){
  print(class(pdata[,i]))
}

for (i in 1:eLen){
  # print(i)
  print(class(edata[,i]))
}
#dim(edata)
#class(pdata)
#class(edata)

#1.b
colnames(edata)
rownames(edata)

colnames(pdata)
row.names(pdata)

row.names(fdata)
#1.c
for (i in 1:pLen){
  print(summary(pdata[,i]))
}

for (i in 1:eLen){
  # print(i)
  print(summary(edata[,i]))
}

#1.d

table(pdata$Tissue, useNA = "ifany")
table(pdata$SubType, useNA = "ifany")
table(pdata$ClinicalGroup, useNA = "ifany")
table(pdata$Status, useNA = "ifany")

#1.e

cov(edata[,1:10], y = NULL, use = "everything", method = "spearman")
corr <- cor(edata[,1:10], y = NULL, use = "everything", method = "spearman")
#rcorr(edata[,1:10])
col<- colorRampPalette(c("blue", "white", "red"))(20)

corrplot(corr, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
dev.copy(device = jpeg, file = "correlationPlot.jpeg")
dev.off()

heatmap(x = corr, col = col, symm = TRUE)
dev.copy(device = jpeg, file = "correlationPlotHeat.jpeg")
dev.off()

#1.f
x <- edata[,"GSM95478"]
y <- edata[,"GSM95473"]

plot(x, y, pch = 16, col = "pink")
abline(lm(y ~ x), col = "purple", lwd = 3)
dev.copy(device = jpeg, file = "lineRelation.jpeg")
dev.off()

###############################################################################

#2

#first method of normalization
edataScaled = scale(edata)
#print(edataScaled)
pc2 = prcomp(edataScaled)
#plot(pc2$rotation[,1],svd2$v[,1])
edata_centered_2t = t(t(edataScaled) - colMeans(edataScaled))
#edata_centered_2 = edataScaled - colMeans(edataScaled)
svd_2 = svd(edata_centered_2t)
plot(pc2$rotation,svd_2$v,col=5)
dev.copy(device = jpeg, file = "svdPCAscaled.jpeg")
dev.off()

par(mfrow = c(1,2))
plot(pc2$rotation,col= 1)
dev.copy(device = jpeg, file = "svdPCAscaled_1.jpeg")
dev.off()
plot(svd_2$v, col = 2)
dev.copy(device = jpeg, file = "svdPCAscaled_2.jpeg")
dev.off()
par(mfrow = c(1,1))

#second method of normalization
edataNorm = log2(edata + 1)
edataOmit <- na.omit(edataNorm)
pc3 = prcomp(edataOmit)
edata_centered_3t = t(t(edataOmit) - colMeans(edataOmit))
#edata_centered_3 = edataOmit - colMeans(edataOmit)
svd_3 = svd(edata_centered_3)
plot(pc3$rotation,svd_3$v,col=6)
dev.copy(device = jpeg, file = "svdPCAlog.jpeg")
dev.off()
par(mfrow = c(1,2))
plot(pc3$rotation,col= 3)
dev.copy(device = jpeg, file = "svdPCAlog_1.jpeg")
dev.off()
plot(svd_3$v, col = 4)
dev.copy(device = jpeg, file = "svdPCAlog_2.jpeg")
dev.off()
par(mfrow = c(1,1))

#third method of normalization
edataQuanNorm = normalize.quantiles(as.matrix(edata))
pc4 = prcomp(edataQuanNorm)
edata_centered_4t = t(t(edataQuanNorm) - colMeans(edataQuanNorm))
#edata_centered_4 = edataQuanNorm - colMeans(edataQuanNorm)
svd_4 = svd(edata_centered_4)
plot(pc4$rotation,svd_4$v,col=3)
dev.copy(device = jpeg, file = "svdPCAquantile.jpeg")
dev.off()
par(mfrow = c(1,2))
plot(pc4$rotation,col= 5)
dev.copy(device = jpeg, file = "svdPCAquantile_1.jpeg")
dev.off()
plot(svd_4$v, col = 6)
dev.copy(device = jpeg, file = "svdPCAquantile_2.jpeg")
dev.off()

print("Similar Plots and Values using Scale normalisation")
print(pc2$rotation)
print(svd_2$v)
      
plot(pc2$rotation,col= 1)
plot(svd_2$v, col = 2)

print(length(pc2$rotation))
print(sum(pc2$rotation == svd_2$v))
print(sum(pc2$rotation != svd_2$v))
par(mfrow = c(1,1))
###############################################################################

#3


phenotypes <- as.factor(c(rep("Aries",29), rep("Taurus", 24), rep("Gemini", 22), 
                          rep("Cancer", 19), rep("Leo", 21), rep("Virgo",18), 
                          rep("Libra", 19), rep("Scorpio", 20), rep("Sagittarius", 23), 
                          rep("Capricorn", 18), rep("Aquarius", 20), rep("Pisces",23)))
p <- c(1,1,1,1,1,1,1,1,1,1,1,1)
p <-  p/sum(p)
tPheno <- table(phenotypes)

chisq.test(tPheno, p = p)


cat("H0 Hypothesis: zodiac signs are evenly distributed across visual artists.\n",
    "Since the p-value of the chisq-test turned out to be higher than 0.05,\n",
    "therefore, the H0","(","null hypothesis",")", "is accepted")
cat("H1 Hypothesis: zodiac signs are unevenly distributed across visual artists.\n",
    "Since the p-value of the chisq-test turned out to be higher than 0.05,",
    "therefore, \nthe H0 is accepted, and the H1","(","alternative hypothesis",
    ")", "is rejected")

###############################################################################

#4

ecuDistO = dist(t(edataOmit[,1:10]))
#ecuDistS = dist(t(edataScaled[,1:10]))
#ecuDistQ = dist(t(edataQuanNorm[,1:10]))
hierClust = hclust(ecuDistO)
plot(hierClust,hang = -1)
dev.copy(device = jpeg, file = "hierClusters.jpeg")
dev.off()

kCluster = kmeans(edataOmit,centers=4)
dim(kCluster$centers)
clust1Cent = kCluster$centers[1,]
clust2Cent = kCluster$centers[2,]
clust3Cent = kCluster$centers[3,]
clust4Cent = kCluster$centers[4,]
cat("Centroid of first cluster:",clust1Cent,
    "\nCentroid of second cluster:",clust2Cent,
    "\nCentroid of third cluster:",clust3Cent,
    "\nCentroid of forth cluster:",clust4Cent)
table(kCluster$cluster)
