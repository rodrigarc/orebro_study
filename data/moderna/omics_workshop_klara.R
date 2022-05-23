#Exercise 3: Data analysis using R
#Start by creating a folder on your desktop, name it "Data analysis Course R" and put the datafiles for the exercise in that folder 

#Setting the working directory
setwd("~/Box Sync/PhD courses/3102 Omics data analysis/Omics workshop")
library(dplyr)
library(matrixStats)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(wesanderson)
library(limma)
library(gridExtra)
# ~ equals klaralenart folder in this case
MA <- read.csv(file = "MA data.csv", header = T, row.names = 1, check.names = FALSE)
MAinfo <- read.csv(file = "MAinfo.csv", header = T, row.names = 1, check.names = FALSE)
ProbeBG <- read.csv(file = "ProbeOverBackground.csv", header = T, row.names = 1, check.names = FALSE)
ProbeInfo <- read.csv(file = "ProbeInfo.csv", header = T, row.names = 1, check.names = FALSE)

identical(rownames(MA), rownames(ProbeInfo)) #T
identical(rownames(MAinfo), colnames(MA)) #T
identical(rownames(ProbeInfo), rownames(ProbeBG)) #T
identical(colnames(MA), colnames(ProbeBG)) #T

# Log2 transform the MA data
MA <- log2(MA)

# Use only probes with associated gene ID
MA.clean <- subset(MA[ProbeInfo$GeneID != 0,])        #43603 genes to 26776 genes
ProbeBG <- subset(ProbeBG[ProbeInfo$GeneID != 0,])
# Use only probes with expression over background in all samples of at least one group
ProbeBG$ProbeID <- row.names(ProbeBG)
ProbeBG <- ProbeBG[,c(36, 1:35)]
ProbeBG <- data.frame(lapply(ProbeBG, function(x) {gsub("TRUE", as.numeric(1), x)}))
ProbeBG <- data.frame(lapply(ProbeBG, function(x) {gsub("FALSE", as.numeric(0), x)}))
PBG2 <- ProbeBG[,-1]
rownames(PBG2) <- ProbeBG[,1]
ProbeBG <- PBG2
ProbeBG$sum <- rowSums(data.matrix(ProbeBG))
ProbeBG$sum <- ProbeBG$sum - 35
ProbeBG <- ProbeBG %>% rename("I01_W4-24H" = "I01_W4.24H", 
                              "I04_W4-24H" = "I04_W4.24H", 
                              "I06_W4-24H" = "I06_W4.24H",
                              "I14_W4-24H" = "I14_W4.24H",
                              "I16_W4-24H" = "I16_W4.24H",
                              "I02_W4-24H" = "I02_W4.24H",
                              "I13_W4-24H" = "I13_W4.24H")


  # data.matrix() replaces factors with integers. TRUE = 2, FALSE = 1. 
  # I subtracted 35 from all cells for easier visualization, and then TRUE = 1.
  # I can remove all that are = 0, because they have FALSE (not expressed above background) for all samples
hist(ProbeBG$sum, breaks = 36)  # histogram
count(ProbeBG$sum > 0 & ProbeBG$sum < 35) # 10490 genes are expressed above background in some but not all samples
count(ProbeBG$sum > 0 & ProbeBG$sum < 10) # 4081 genes
count(ProbeBG$sum > 9 & ProbeBG$sum < 20) # 1587 genes
count(ProbeBG$sum > 19 & ProbeBG$sum < 30) # 2617 genes
count(ProbeBG$sum > 29 & ProbeBG$sum < 35) # 2205 genes

ProbeBG <- subset(ProbeBG[ProbeBG$sum != 0,])

A <- rownames(MAinfo[(MAinfo$Vaccine =="mRNA"),])
B <- rownames(MAinfo[(MAinfo$Vaccine =="Vaxigrip"),])
C <- rownames(MAinfo[(MAinfo$Vaccine =="Fluad"),])

Asum <- ProbeBG %>% select(., matches(A))
Asum$sum <- (rowSums(data.matrix(Asum)) - 15)
A.good <- subset(Asum[Asum$sum == 15,])

Bsum <- ProbeBG %>% select(., matches(B))
Bsum$sum <- (rowSums(data.matrix(Bsum)) - 9)
B.good <- subset(Bsum[Bsum$sum == 9,])

Csum <- ProbeBG %>% select(., matches(C))
Csum$sum <- (rowSums(data.matrix(Csum)) - 11)
C.good <- subset(Csum[Csum$sum == 11,])

X <- unique(rownames(A.good), rownames(B.good), rownames(C.good), fromLast = FALSE)
Good.all <- as.data.frame(unique(X))
Good.all <- Good.all %>% rename(Probe = "unique(X)")
rownames(Good.all) <- Good.all$Probe

# Use only Probes that are expressed above background in all samples from 1 group
MA.clean <- subset(MA.clean[rownames(Good.all),])    # 26776 genes to 14204 genes

# PCA plots - group
MA.pca <- prcomp(t(MA.clean), scale. = TRUE)
score.df <- as.data.frame(MA.pca$x)
score.df$Vaccine <- MAinfo[rownames(score.df),"Vaccine"]
score.df$Timepoint <- MAinfo[rownames(score.df),"Timepoint"]
score.df$NHP <- MAinfo[rownames(score.df),"AnimalID"]

#Skree plot
Contributions <- round(summary(MA.pca)$importance[2,] * 100, 2)
barplot(Contributions[1:20], las=3, main="Screeplot MA PCA", xlab="component", ylab="Contribution(%)", col="lightblue")
PC1contr <- Contributions[1]
PC2contr <- Contributions[2]

# PCA colored by vaccine
PCA1 <- ggplot(data = score.df, aes(x = PC1, y = PC2, label = rownames(score.df))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(aes(color=Vaccine), size = 2) +
  scale_color_manual(values=wes_palette(n=3, "GrandBudapest1")) +
  xlab(sprintf("PC1 (%s%%)", PC1contr)) + ylab(sprintf("PC2 (%s%%)", PC2contr))+ 
  ggtitle("PC1 vs PC2 - Vaccine group")

PCA2 <- ggplot(data = score.df, aes(x = PC1, y = PC2, label = rownames(score.df))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(aes(color=Timepoint), size = 2) +
  scale_color_manual(values=wes_palette(n=3, "BottleRocket2")) +
  xlab(sprintf("PC1 (%s%%)", PC1contr)) + ylab(sprintf("PC2 (%s%%)", PC2contr))+ 
  ggtitle("PC1 vs PC2 - Timepoint")

my_col <- colorRampPalette(brewer.pal(n = 8, name = "Dark2")) (15)
PCA3 <- ggplot(data = score.df, aes(x = PC1, y = PC2, label = rownames(score.df))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(aes(color=NHP), size = 2) +
  scale_color_manual(values=my_col) +
  xlab(sprintf("PC1 (%s%%)", PC1contr)) + ylab(sprintf("PC2 (%s%%)", PC2contr))+ 
  ggtitle("PC1 vs PC2 - NHP")

ggsave("PCA_new.pdf", arrangeGrob(PCA1, PCA2, PCA3, ncol = 3), height = 5, width = 15)

# Heatmap by Euclidian distance
my_palette <- colorRampPalette(c("deepskyblue3", "white", "violetred2"))(n = 299)
gr.col <- as.numeric(MAinfo$Vaccine)
col.group <- wes_palette(n=3, "GrandBudapest1")
tp.col <- as.numeric(MAinfo$Timepoint)
col.tp <- wes_palette(n=3, "BottleRocket2")

par(mfrow = c(1,2))
pdf("Heatmap_new.pdf", height = 5, width = 10)
heatmap.2(data.matrix(MA.clean), main="Heatmap by vaccine group", sub="Euclidian, scale by row", trace="none", symbreaks=TRUE, margins = c(10, 20), 
          key=FALSE, col=my_palette, na.color="gray80", cexCol = 0.3, cexRow = 0.5, density.info="none", 
          dendrogram="column", ColSideColors=col.group[gr.col], scale="row", labRow = F, labCol = F)
legend("topright", legend=unique(MAinfo$Vaccine),fill=col.group[unique(gr.col)])

heatmap.2(data.matrix(MA.clean), main="Heatmap by timepoint", sub="Euclidian, scale by row", trace="none", symbreaks=TRUE, margins = c(10, 20), 
                key=FALSE, col=my_palette, na.color="gray80", cexCol = 0.3, cexRow = 0.5, density.info="none", 
                dendrogram="column", ColSideColors=col.tp[tp.col], scale="row", labRow = F, labCol = F)
legend("topright", legend=unique(MAinfo$Timepoint),fill=col.tp[unique(tp.col)])
dev.off()

# Heatmap by Pearson-Ward
rowv <- as.dendrogram(hclust(as.dist(1-cor(t(MA.clean), method ="pearson", use = "pairwise.complete.obs")), method="ward.D2"))
colv <- as.dendrogram(hclust(as.dist(1-cor(MA.clean , method ="pearson", use = "pairwise.complete.obs")), method="ward.D2"))

par(mfrow = c(1,2))
pdf("Heatmap_PW.pdf", height = 5, width = 10)
heatmap.2(data.matrix(MA.clean), main="Heatmap by vaccine group", sub="Pearson,WardD2, scale by row", trace="none", symbreaks=TRUE, margins = c(10, 20), 
          key=FALSE, col=my_palette, na.color="gray80", cexCol = 0.3, cexRow = 0.5, density.info="none", 
          Rowv = rowv, Colv = colv, dendrogram="column", ColSideColors=col.group[gr.col], scale="row", labRow = F, labCol = F)
legend("topright", legend=unique(MAinfo$Vaccine),fill=col.group[unique(gr.col)])

heatmap.2(data.matrix(MA.clean), main="Heatmap by timepoint", sub="Pearson,WardD2, scale by row", trace="none", symbreaks=TRUE, margins = c(10, 20), 
          key=FALSE, col=my_palette, na.color="gray80", cexCol = 0.3, cexRow = 0.5, density.info="none", 
          Rowv = rowv, Colv = colv, dendrogram="column", ColSideColors=col.tp[tp.col], scale="row", labRow = F, labCol = F)
legend("topright", legend=unique(MAinfo$Timepoint),fill=col.tp[unique(tp.col)])
dev.off()

# Limma analysis for DEGs - all groups together
# Setting up contrasts (first extract relevant sample IDs)
Comp1 <- rownames(MAinfo[(MAinfo$Timepoint == "0H" | MAinfo$Timepoint == "24H"),])
Comp2 <- rownames(MAinfo[(MAinfo$Timepoint == "0H" | MAinfo$Timepoint == "W4-24H"),])
Comp2 <- Comp2[-c(5, 8, 11, 12, 13, 14, 21)]

# Prepare files including only relevant sample IDs
MAinfo24H <- MAinfo[Comp1,]
MAinfoW4 <- MAinfo[Comp2,]
MA.clean24H <- MA.clean[,Comp1]
MA.cleanW4 <- MA.clean[,Comp2]

# Calcs for 24H vs 0H
NHP24H <- factor(MAinfo24H$AnimalID)
TP24H <- factor(MAinfo24H$Timepoint, levels=c("0H","24H"))
design24H <- model.matrix(~NHP24H+TP24H)
fit24H <- lmFit(MA.clean24H, design24H)
fit24H <- eBayes(fit24H)
topTable(fit24H, coef="TP24H24H")  # Top10 hits
DEG24H <- topTable(fit24H, coef="TP24H24H", sort="none", n=Inf)
DEG24H$neg.adj.p <- -log10(DEG24H$adj.P.Val)

plot(DEG24H$logFC, DEG24H$neg.adj.p, xlim=c(-6,6), ylim = c(0, 6), main="24H DEG analysis", xlab="fold difference (log2)",ylab="adj.p-value (-log10)",pch=20, cex=0.3, col="black")
p.threshold = -log10(0.01)
abline(h=p.threshold,v=c(-1,1),lty=2,col="red")

# Calcs for W4-24H vs 0H
NHPW4 <- factor(MAinfoW4$AnimalID)
TPW4 <- factor(MAinfoW4$Timepoint, levels=c("0H","W4-24H"))
designW4 <- model.matrix(~NHPW4+TPW4)
fitW4 <- lmFit(MA.cleanW4, designW4)
fitW4 <- eBayes(fitW4)
topTable(fitW4, coef="TPW4W4-24H")  # Top10 hits
DEGW4 <- topTable(fitW4, coef="TPW4W4-24H", sort="none", n=Inf)
DEGW4$neg.adj.p <- -log10(DEGW4$adj.P.Val)

plot(DEGW4$logFC, DEGW4$neg.adj.p, xlim=c(-6,6), ylim = c(0, 6), main="W4-24H DEG analysis", xlab="fold difference (log2)",ylab="adj.p-value (-log10)",pch=20, cex=0.3, col="black")
p.threshold = -log10(0.01)
abline(h=p.threshold,v=c(-1,1),lty=2,col="red")

# Make analysis per group
mRNA <- MAinfo[(MAinfo$Vaccine =="mRNA"),]
Vaxigrip <- MAinfo[(MAinfo$Vaccine =="Vaxigrip"),]
Fluad <- MAinfo[(MAinfo$Vaccine =="Fluad"),]

mRNA24H <- mRNA[mRNA$Timepoint == "0H" | mRNA$Timepoint == "24H",]
Vaxigrip24H <- Vaxigrip[Vaxigrip$Timepoint == "0H" | Vaxigrip$Timepoint == "24H",]
Fluad24H <- Fluad[Fluad$Timepoint == "0H" | Fluad$Timepoint == "24H",]

Aa <- rownames(mRNA24H)
Bb <- rownames(Vaxigrip24H)
Cc <- rownames(Fluad24H)

mRNA.MA.24H <- MA.clean[,Aa]
Vaxi.MA.24H <- MA.clean[,Bb]
Flua.MA.24H <- MA.clean[,Cc]

par(mfrow = c(1,1))
p.threshold = -log10(0.05)   # p.adj < 0.05 is considered significant

# mRNA, 24H vs 0H
NHP.A <- factor(mRNA24H$AnimalID)
TP.A <- factor(mRNA24H$Timepoint, levels=c("0H","24H"))
design.A <- model.matrix(~NHP.A+TP.A)
fit.A <- lmFit(mRNA.MA.24H, design.A)
fit.A <- eBayes(fit.A)
topTable(fit.A, coef="TP.A24H")  # Top10 hits
DEG.mRNA.24H <- topTable(fit.A, coef="TP.A24H", sort="none", n=Inf)
DEG.mRNA.24H$neg.adj.p <- -log10(DEG.mRNA.24H$adj.P.Val)
DEG.mRNA.24H$neg.p <- -log10(DEG.mRNA.24H$P.Value)

#Vaxigrip, 24H vs 0H
NHP.B <- factor(Vaxigrip24H$AnimalID)
TP.B <- factor(Vaxigrip24H$Timepoint, levels=c("0H","24H"))
design.B <- model.matrix(~NHP.B+TP.B)
fit.B <- lmFit(Vaxi.MA.24H, design.B)
fit.B <- eBayes(fit.B)
topTable(fit.B, coef="TP.B24H")  # Top10 hits
DEG.Vaxi.24H <- topTable(fit.B, coef="TP.B24H", sort="none", n=Inf)
DEG.Vaxi.24H$neg.adj.p <- -log10(DEG.Vaxi.24H$adj.P.Val)
DEG.Vaxi.24H$neg.p <- -log10(DEG.Vaxi.24H$P.Value)

#Fluad, 24H vs 0H
NHP.C <- factor(Fluad24H$AnimalID)
TP.C <- factor(Fluad24H$Timepoint, levels=c("0H","24H"))
design.C <- model.matrix(~NHP.C+TP.C)
fit.C<- lmFit(Flua.MA.24H, design.C)
fit.C <- eBayes(fit.C)
topTable(fit.C, coef="TP.C24H")  # Top10 hits
DEG.Flua.24H <- topTable(fit.C, coef="TP.C24H", sort="none", n=Inf)
DEG.Flua.24H$neg.adj.p <- -log10(DEG.Flua.24H$adj.P.Val)
DEG.Flua.24H$neg.p <- -log10(DEG.Flua.24H$P.Value)

pdf("Volcano plots DEG - raw p.pdf", height = 5, width = 15)       # Volcano plot function included in the limma package
par(mfrow = c(1,3))
volcanoplot(fit.A, coef = "TP.A24H", style = "p-value", highlight = 0, col = "grey", main = "DEG 24H vs 0H - mRNA",
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35, xlim = c(-5, 9), ylim = c(0, 8.5))
abline(h=p.threshold,v=c(-1,1),lty=2,col="red")

volcanoplot(fit.B, coef = "TP.B24H", style = "p-value", highlight = 0, col = "grey", main = "DEG 24H vs 0H - Vaxigrip",
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35, xlim = c(-5, 9), ylim = c(0, 8.5))
abline(h=p.threshold,v=c(-1,1),lty=2,col="red")

volcanoplot(fit.C, coef = "TP.C24H", style = "p-value", highlight = 0, col = "grey", main = "DEG 24H vs 0H - Fluad",
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35, xlim = c(-5, 9), ylim = c(0, 8.5))
abline(h=p.threshold,v=c(-1,1),lty=2,col="red")
dev.off()

pdf("Volcano plots DEG - FDR or p.pdf", height = 15, width = 10)  #Comparison of plots if used FDR or raw p-values (different thresholds used)
par(mfrow = c(3, 2))
plot(DEG.mRNA.24H$logFC, DEG.mRNA.24H$neg.adj.p, xlim=c(-5,9), ylim = c(0, 8.5), main="mRNA - adj.p", xlab="fold difference (log2)",ylab="adj.p-value (-log10)",pch=20, cex=0.3, col="black")
p.threshold = -log10(0.05)
abline(h=p.threshold,v=c(-1,1),lty=2,col="red")

plot(DEG.mRNA.24H$logFC, DEG.mRNA.24H$neg.p, xlim=c(-5,9), ylim = c(0, 8.5), main="mRNA - raw p", xlab="fold difference (log2)",ylab="p-value (-log10)",pch=20, cex=0.3, col="black")
p.threshold = -log10(0.01)
abline(h=p.threshold,v=c(-1,1),lty=2,col="red")

plot(DEG.Vaxi.24H$logFC, DEG.Vaxi.24H$neg.adj.p, xlim=c(-5,9), ylim = c(0, 8.5), main="Vaxigrip - adj.p", xlab="fold difference (log2)",ylab="adj.p-value (-log10)",pch=20, cex=0.3, col="black")
p.threshold = -log10(0.05)
abline(h=p.threshold,v=c(-1,1),lty=2,col="red")

plot(DEG.Vaxi.24H$logFC, DEG.Vaxi.24H$neg.p, xlim=c(-5,9), ylim = c(0, 8.5), main="Vaxigrip - raw p", xlab="fold difference (log2)",ylab="p-value (-log10)",pch=20, cex=0.3, col="black")
p.threshold = -log10(0.01)
abline(h=p.threshold,v=c(-1,1),lty=2,col="red")

plot(DEG.Flua.24H$logFC, DEG.Flua.24H$neg.adj.p, xlim=c(-5,9), ylim = c(0, 8.5), main="Fluad - adj.p", xlab="fold difference (log2)",ylab="adj.p-value (-log10)",pch=20, cex=0.3, col="black")
p.threshold = -log10(0.05)
abline(h=p.threshold,v=c(-1,1),lty=2,col="red")

plot(DEG.Flua.24H$logFC, DEG.Flua.24H$neg.p, xlim=c(-5,9), ylim = c(0, 8.5), main="Fluad - raw p", xlab="fold difference (log2)",ylab="p-value (-log10)",pch=20, cex=0.3, col="black")
p.threshold = -log10(0.01)
abline(h=p.threshold,v=c(-1,1),lty=2,col="red")
dev.off()

# Extract lists of DEGs for functional analysis - raw p-value < 0.01 and LogFC > 1 | LogFC < -1
Func.A <- DEG.mRNA.24H[abs(DEG.mRNA.24H$logFC) > 1 & DEG.mRNA.24H$P.Value < 0.01,] # 3916 genes
Func.B <- DEG.Vaxi.24H[abs(DEG.Vaxi.24H$logFC) > 1 & DEG.Vaxi.24H$P.Value < 0.01,] # 2122 genes
Func.C <- DEG.Flua.24H[abs(DEG.Flua.24H$logFC) > 1 & DEG.Flua.24H$P.Value < 0.01,] # 4809 genes

ProbeInfo$ProbeID <- row.names(ProbeInfo)   # names need to be in columns, not as row.names to be able to merge tables
Func.A$ProbeID <- row.names(Func.A)
Func.B$ProbeID <- row.names(Func.B)
Func.C$ProbeID <- row.names(Func.C)

Func.A <- left_join(Func.A, ProbeInfo, by = "ProbeID")
Func.B <- left_join(Func.B, ProbeInfo, by = "ProbeID")
Func.C <- left_join(Func.C, ProbeInfo, by = "ProbeID")

Func.A.RNK <- Func.A[,c(10,1:9)]
Func.A.RNK <- Func.A.RNK[order(-Func.A.RNK$logFC),]
Func.B.RNK <- Func.B[,c(10,1:9)]
Func.B.RNK <- Func.B.RNK[order(-Func.B.RNK$logFC),]
Func.C.RNK <- Func.C[,c(10,1:9)]
Func.C.RNK <- Func.C.RNK[order(-Func.C.RNK$logFC),]

write.table(Func.A.RNK, file = "Func.A.rnk", row.names = F, col.names = T, sep = "\t")
write.table(Func.B.RNK, file = "Func.B.rnk", row.names = F, col.names = T, sep = "\t")
write.table(Func.C.RNK, file = "Func.C.rnk", row.names = F, col.names = T, sep = "\t")

Func.A.RNK <- Func.A.RNK[,c(2, 10)] %>% rename(A_LogFC = logFC)
Func.B.RNK <- Func.B.RNK[,c(2, 10)] %>% rename(B_LogFC = logFC)
Func.C.RNK <- Func.C.RNK[,c(2, 10)] %>% rename(C_LogFC = logFC)

Func.RNK <- left_join(Func.A.RNK, Func.B.RNK, by = "ProbeID")
Func.RNK <- left_join(Func.RNK, Func.C.RNK, by = "ProbeID")
Func.RNK[is.na(Func.RNK)] <- 0

DEG.pca <- prcomp(t(Func.RNK[,c(2:4)]), scale. = TRUE)
score.df <- as.data.frame(DEG.pca$x)
score.df$Vaccine <- c("mRNA", "Vaxigrip", "Fluad")

Contributions <- round(summary(DEG.pca)$importance[2,] * 100, 2)
barplot(Contributions[1:20], las=3, main="Screeplot DEG PCA", xlab="component", ylab="Contribution(%)", col="lightblue")
PC1contr <- Contributions[1]
PC2contr <- Contributions[2]

# PCA colored by vaccine
ggplot(data = score.df, aes(x = PC1, y = PC2, label = rownames(score.df))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(aes(color=Vaccine), size = 2) +
  scale_color_manual(values=wes_palette(n=3, "GrandBudapest1")) +
  xlab(sprintf("PC1 (%s%%)", PC1contr)) + ylab(sprintf("PC2 (%s%%)", PC2contr))+ 
  ggtitle("PC1 vs PC2 - Vaccine group")
# This gives one point per vaccine group. I want one point per NHP.
# Need to calculate DEGs for each NHP separately, at 24H

I01.FC <- MA.clean24H$I01_24H/MA.clean24H$I01_D0
I02.FC <- MA.clean24H$I02_24H/MA.clean24H$I02_D0
I03.FC <- MA.clean24H$I03_24H/MA.clean24H$I03_D0
I04.FC <- MA.clean24H$I04_24H/MA.clean24H$I04_D0
I05.FC <- MA.clean24H$I05_24H/MA.clean24H$I05_D0
I06.FC <- MA.clean24H$I06_24H/MA.clean24H$I06_D0
I07.FC <- MA.clean24H$I07_24H/MA.clean24H$I07_D0
I09.FC <- MA.clean24H$I09_24H/MA.clean24H$I09_D0
I10.FC <- MA.clean24H$I10_24H/MA.clean24H$I10_D0
I12.FC <- MA.clean24H$I12_24H/MA.clean24H$I12_D0
I13.FC <- MA.clean24H$I13_24H/MA.clean24H$I13_D0
I14.FC <- MA.clean24H$I14_24H/MA.clean24H$I14_D0
I16.FC <- MA.clean24H$I16_24H/MA.clean24H$I16_D0
I17.FC <- MA.clean24H$I17_24H/MA.clean24H$I17_D0

LogFC <- cbind(I01.FC, I02.FC, I03.FC, I04.FC, I05.FC, I06.FC, I07.FC, 
               I09.FC, I10.FC, I12.FC, I13.FC, I14.FC, I16.FC, I17.FC)
rownames(LogFC) <- row.names(MA.clean24H)
LogFC <- as.data.frame(LogFC)
LogFC$ProbeID <- rownames(LogFC)
LogFC <- left_join(LogFC, ProbeInfo, by = "ProbeID")

DEG.RNK <- filter(LogFC, ProbeID %in% Func.RNK$ProbeID)
DEG.RNK <- DEG.RNK[,c(1:14)]

NHP.DEG <- colnames(DEG.RNK)
NHP.DEG <- NHP.DEG[c(1:14)]
Vaccine <- c("mRNA", "Vaxigrip", "Fluad", "mRNA", "Fluad", "mRNA", "Fluad", "Vaxigrip",
             "Vaxigrip", "Fluad", "Fluad", "mRNA", "mRNA", "Vaxigrip")

DEG.RNK.info <- as.data.frame(NHP.DEG, Vaccine)
DEG.RNK.info$Vaccine <- row.names(DEG.RNK.info)
row.names(DEG.RNK.info) <- DEG.RNK.info$NHP.DEG
DEG.RNK.info <- DEG.RNK.info %>% rename(NHP = NHP.DEG)
DEG.RNK.info$NHP <- DEG.RNK.info$NHP %>% gsub("\\.FC", "", .)

DEG.pca2 <- prcomp(t(DEG.RNK), scale. = TRUE)
score.df2 <- as.data.frame(DEG.pca2$x)
score.df2$Vaccine <- DEG.RNK.info[rownames(score.df2),"Vaccine"]
score.df2$NHP <- DEG.RNK.info[rownames(score.df2),"NHP"]

Contributions <- round(summary(DEG.pca2)$importance[2,] * 100, 2)
barplot(Contributions[1:20], las=3, main="Screeplot DEG PCA", xlab="component", ylab="Contribution(%)", col="lightblue")
PC1contr <- Contributions[1]
PC2contr <- Contributions[2]

q <- ggplot(data = score.df2, aes(x = PC1, y = PC2, label = rownames(score.df2)), margin = c(5,5)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(aes(color=Vaccine), size = 2) +
  scale_color_manual(values=wes_palette(n=3, "GrandBudapest1")) +
  xlab(sprintf("PC1 (%s%%)", PC1contr)) + ylab(sprintf("PC2 (%s%%)", PC2contr))+ 
  ggtitle("PC1 vs PC2 - Vaccine group")

ggsave("PCA DEG.pdf", plot = q, width = 6, height = 5,)

# Extract lists of DEGs for GSEA - all genes, ordered by descending logFC values
GSEA.A <- DEG.mRNA.24H
GSEA.B <- DEG.Vaxi.24H
GSEA.C <- DEG.Flua.24H

GSEA.A$ProbeID <- row.names(GSEA.A)
GSEA.B$ProbeID <- row.names(GSEA.B)
GSEA.C$ProbeID <- row.names(GSEA.C)

GSEA.A <- left_join(GSEA.A, ProbeInfo, by = "ProbeID")
GSEA.B <- left_join(GSEA.B, ProbeInfo, by = "ProbeID")
GSEA.C <- left_join(GSEA.C, ProbeInfo, by = "ProbeID")

GSEA.A.RNK <- GSEA.A[,c(10,1)]
GSEA.A.RNK <- GSEA.A.RNK[order(-GSEA.A.RNK$logFC),]
GSEA.A.RNK$Rank <- 1:14204    # Third column needed for clusterProfiler:GSEA
GSEA.B.RNK <- GSEA.B[,c(10,1)]
GSEA.B.RNK <- GSEA.B.RNK[order(-GSEA.B.RNK$logFC),]
GSEA.B.RNK$Rank <- 1:14204
GSEA.C.RNK <- GSEA.C[,c(10,1)]
GSEA.C.RNK <- GSEA.C.RNK[order(-GSEA.C.RNK$logFC),]
GSEA.B.RNK$Rank <- 1:14204

write.table(GSEA.A.RNK, file = "GSEA.A.rnk", row.names = F, col.names = F, sep = "\t")
write.table(GSEA.B.RNK, file = "GSEA.B.rnk", row.names = F, col.names = F, sep = "\t")
write.table(GSEA.C.RNK, file = "GSEA.C.rnk", row.names = F, col.names = F, sep = "\t")

library(GSA)
library(qusage)
library(tidyr)
BTM <- read.gmt("BTM_for_GSEA_20131008.gmt")
BTM$ont <- BTM$ont %>% gsub("\\(", "", .) %>% gsub("\\)", "", .)
BTM$ID <- BTM$ont
BTM <- BTM[,c(1,3,2)]
BTM$ont <- BTM$ont %>% gsub("M\\d*.\\d*$", "", .)
BTM$ID <- BTM$ID %>% gsub(".*M", "M", .) %>% gsub(".*S", "S", .)
BTM2name <- BTM %>% select(ID, gene)

BTMs <- BTM[,c(1,2)] %>% unique(.)

C7 <- GSA.read.gmt(file = "c7.all.v7.2.symbols.gmt")


Intersect <- intersect(GSEA.A.RNK$GeneID, BTMgenes) #1532
Intersect <- intersect(GSEA.A.RNK$GeneID, C7genes) #8983
?intersect

library(clusterProfiler)
browseVignettes("clusterProfiler")

# Perform GSEA within clusterProfiler
GSEA.list.1 <- GSEA.A.RNK$logFC
names(GSEA.list.1) <- as.character(GSEA.A.RNK$GeneID)
GSEA.list.1 <- sort(GSEA.list.1, decreasing = TRUE)
GSEA.1 <- GSEA(GSEA.list.1, TERM2NAME = BTM2name, TERM2GENE = BTM2name, pAdjustMethod = "fdr", pvalueCutoff = 0.1)
GSEA.1 <- as.data.frame(GSEA.1)

GSEA.list.2 <- GSEA.B.RNK$logFC
names(GSEA.list.2) <- as.character(GSEA.B.RNK$GeneID)
GSEA.list.2 <- sort(GSEA.list.2, decreasing = TRUE)
GSEA.2 <- GSEA(GSEA.list.2, TERM2NAME = BTM2name, TERM2GENE = BTM2name, pAdjustMethod = "fdr", pvalueCutoff = 0.1)
GSEA.2 <- as.data.frame(GSEA.2)

GSEA.list.3 <- GSEA.C.RNK$logFC
names(GSEA.list.3) <- as.character(GSEA.C.RNK$GeneID)
GSEA.list.3 <- sort(GSEA.list.3, decreasing = TRUE)
GSEA.3 <- GSEA(GSEA.list.3, TERM2NAME = BTM2name, TERM2GENE = BTM2name, pAdjustMethod = "fdr", pvalueCutoff = 0.1)
GSEA.3 <- as.data.frame(GSEA.3)

GSEA.1 <- remove_rownames(GSEA.1)
GSEA.2 <- remove_rownames(GSEA.2)
GSEA.3 <- remove_rownames(GSEA.3)

GSEA_Comb <- GSEA.1[,c(1, 5, 7)] %>% 
  full_join(., GSEA.2[,c(1,5,7)], by = "ID", suffix = c(".mRNA", ".Vaxi")) %>%
  full_join(., GSEA.3[,c(1,5,7)], by = "ID", suffix = c("", ".Flua")) %>%
  rename(NES.Flua = NES, p.adjust.Flua = p.adjust, BTMx = ID) %>%
  left_join(., BTMs, by = "BTMx", suffix = c("", ""))
row.names(GSEA_Comb) <- GSEA_Comb[,c(1)]
GSEA_Comb <- GSEA_Comb[,-c(1)]
GSEA_Comb <- subset(GSEA_Comb, ont != "TBA ")

my_palette <- colorRampPalette(c("navyblue", "white", "red4"))(n = 299)
heatmap.2(data.matrix(GSEA_Comb[,c(1, 3, 5)]), main="Heatmap of DEGs", dendrogram = "none",     
          trace="none", na.color="white", scale="none", Rowv = NA, Colv = NA, margin = c(5,15), 
          col=my_palette, key=FALSE, cexCol = 0.3, cexRow = 0.5, density.info="none", 
          labRow = F, labCol = F, na.rm = FALSE, 
          sepwidth=c(0.0005, 0.0005), sepcolor="black", colsep = 1:4, rowsep = 1:88)

# Error msg: NA/NaN/Inf in foreign function call (arg 10)
is.infinite(data.matrix(GSEA_Comb[,c(1, 3, 5)]))  # No infinite values in the matrix
is.na(data.matrix(GSEA_Comb[,c(1, 3, 5)])) # Has NA values
is.nan(data.matrix(GSEA_Comb[,c(1, 3, 5)])) # No NaN values in the matrix
# Resolved using " Rowv = NA, Colv = NA "

# Try making a bubble/baloon plot
library(reshape2)  # Need to reformat data first
GSEA_Comb$BTMx <- row.names(GSEA_Comb)
GSEA_Comb <- GSEA_Comb[,c(8, 1:6)]
Reform1 <- melt(GSEA_Comb[,c(1, 2, 4, 6)]) %>%
  rename(Group = variable, NES = value)
Reform1$Group <- gsub("NES\\.", "", Reform1$Group)

Reform2 <- melt(GSEA_Comb[,c(1, 3, 5, 7)]) %>%
  rename(Group = variable, FDR = value)
Reform2$Group <- gsub("p\\.adjust\\.", "", Reform2$Group)
Baloon <- full_join(Reform1, Reform2, by = c("BTMx", "Group"))
Baloon$Group <- Baloon$Group %>% gsub("Flua", "Fluad", .) %>%
  gsub("Vaxi", "Vaxigrip", .)

# Import higher annotation file
BTM_HA <- read.csv(file = "btm_higher_annotation.csv") %>% rename(BTMx = ID)
BTM_HA <- BTM_HA[, -c(4:24)]

BTMs <- BTMs %>% rename(BTMx = ID)
BTMs$ont <- BTMs$ont %>% gsub("S\\d\\d$", "", .) %>% gsub("S\\d$", "", .)
BTMs <- BTMs %>% mutate(BTM = paste(ont, "(", BTMx, ")", sep = ""))
Baloon <- left_join(Baloon, BTMs[,c(2,3)], by = "BTMx")
Baloon <- left_join(Baloon, BTM_HA[,c(1,3)], by = "BTMx")
Baloon <- Baloon[order(Baloon$Higher.annotation),]
BTMorder <- Reform1 %>% 
  left_join(., BTM_HA[,c(1,3)], by = "BTMx")
BTMorder <- BTMorder %>% inner_join(., Baloon, by = "BTMx") %>% .[,c(1,4,8)]
BTMorder <- unique(BTMorder)
BTMorder <- BTMorder[order(BTMorder$Higher.annotation),]
BTMorder <- BTMorder$BTM

pdf(file = "GSEA.pdf", height = 9, width = 5.5)
ggballoonplot(Baloon, x = "Group", y = "BTM", size = "FDR",
              fill = "NES", size.range = c(1, 10)) +
  scale_fill_gradientn(colors = my_palette, limits = c(-3.1, 3.1)) +
  scale_size_continuous(trans = 'reverse') +
  scale_x_discrete(limits = c("mRNA", "Vaxigrip", "Fluad")) +
  scale_y_discrete(limits = BTMorder)
dev.off()

min <- min(Baloon$FDR, na.rm = T)
max <- max(Baloon$FDR, na.rm = T)


library(ggpubr)
?ggballoonplot

#EXERCISE 1: Basic plotting of quant data

#Read Proteomics data from timecourse experiment of drug treatment.
#The data contains relative quantities of 12161 proteins (gene centric) across 10 samples (ctrlx3, 2hx2, 6hx2, 24hx3).
#All quantities are relative to the average of Ctrls (log2)
QCmatrixNN <- as.matrix(read.table("proteomicsQCmatrix.not.normalised.txt",header=T, row.names=1, check.names = FALSE, na.strings=c("", "NA", "Inf"), sep="\t"))

#Get the dimensions of the data matrix. Number of rows will be indicated first and then the number of columns
dim(QCmatrixNN)
  # rows: 12161
  # columns: 10

#View the first(head) or last(tail) six rows of the matrix
head(QCmatrixNN)
tail(QCmatrixNN)
  # genes (rows) are in the alphabetical order

#View specified columns and rows (matrix[rows,columns])
QCmatrixNN[1:10,1:5]

#Boxplot data to view distribution of quant values
boxplot(QCmatrixNN)

#Define range for y axis
boxplot(QCmatrixNN, ylim=c(-1,1))

#Annotate and color
boxplot(QCmatrixNN, ylim=c(-1,1), main="Drug timecourse boxplots", ylab="relative quantity (log2)", col="red")

#Color by sample group
#Read meta data with sample information for annotation
QCInfo <- read.table("QCInfo.txt",header=T, row.names=1, check.names = FALSE, na.strings=c("", "NA", "Inf"), sep="\t")

#Define sample groups for annotation
gr.col <- as.numeric(QCInfo$SampleType)

#Define colors for annotation. Number of colors should be equal to the number of sample groups
col.QC <- c("greenyellow", "yellow", "orange", "red")

#Create boxplot with coloring based on defined annotation
boxplot(QCmatrixNN, ylim=c(-1,1), main="Drug timecourse boxplots", ylab="relative quantity (log2)", col=col.QC[gr.col])

#Add annotation legend
legend("topright", legend=unique(QCInfo$SampleType),fill=col.QC[unique(gr.col)], cex=0.8)

#Normalise data based on sample medians
library(matrixStats)
#Create vector of sample medians
Norm.Factors <- colMedians(QCmatrixNN)
#Use sample medians to center the data
QCmatrix.Norm <- scale(QCmatrixNN, center = Norm.Factors, scale = FALSE)
?scale  # info on scale function. We use Norm.Factors as it represensts column medians.
        # If you ise center = TRUE, it will scale based on column means

#Boxplot normalised data
boxplot(QCmatrix.Norm, ylim=c(-1,1), main="Drug timecourse boxplots", ylab="relative quantity (log2)", col=col.QC[gr.col])
legend("topright", legend=unique(QCInfo$SampleType),fill=col.QC[unique(gr.col)], cex=0.3)

#Save normalised data as txt file
write.table(QCmatrix.Norm, file="QCmatrix.Normalised.txt",row.names=T,quote=F,sep="\t")

#Histogram of ctrlA quantifications
hist(QCmatrix.Norm[,"CtrlA"], breaks=100, col="gold", xlim=c(-1,1), main="Histogram CtrlA")

#95% interval of ctrlA quantifications
spread95ctrlA <- quantile(QCmatrix.Norm[,"CtrlA"], probs=c(.025,.975))

#Add lines for 95%spread
hist(QCmatrix.Norm[,"CtrlA"], breaks=100, col="gold", xlim=c(-1,1), main="Histogram CtrlA")
abline(v=spread95ctrlA, col="red")

#Add text to describe lines
text(0.5,2000,labels="95%interval", col="red")

#Histograms for CtrlA and 24hA for comparison of quant distributions
#Define number of plots (rows,columns) that are displayed together
par(mfrow=c(1,2))
#Create histograms
hist(QCmatrix.Norm[,"CtrlA"], breaks=100, col="gold", xlim=c(-1,1), main="Histogram CtrlA")
abline(v=spread95ctrlA, col="red")
hist(QCmatrix.Norm[,"24hA"], breaks=100, col="blue", xlim=c(-1,1), main="Histogram 24hA")
spread95.24hA <- quantile(QCmatrix.Norm[,"24hA"], probs=c(.025,.975))
abline(v=spread95.24hA, col="red")

#Scatter plot of quant values from CtrlA vs CtrlB compared to 24hA vs 24hB 
par(mfrow=c(1,2))
plot(QCmatrix.Norm[,"CtrlA"],QCmatrix.Norm[,"CtrlB"], xlim=c(-2,2), ylim=c(-2,2), pch=20, cex=0.2, col="purple", main="CtrlA vs CtrlB")
abline(v=c(-1,1), h=c(-1,1), col="red")
plot(QCmatrix.Norm[,"24hA"],QCmatrix.Norm[,"24hB"], xlim=c(-2,2), ylim=c(-2,2), pch=20, cex=0.2, col="purple", main="24hA vs 24hB")
abline(v=c(-1,1), h=c(-1,1), col="red")

#Select genes that are 2 fold up or downregulated in 24hA and 24hB samples
Up24hAandB <- intersect(rownames(QCmatrix.Norm[QCmatrix.Norm[,"24hA"]> 1,]),rownames(QCmatrix.Norm[QCmatrix.Norm[,"24hB"]> 1,]))
Down24hAndB <- intersect(rownames(QCmatrix.Norm[QCmatrix.Norm[,"24hA"]< -1,]),rownames(QCmatrix.Norm[QCmatrix.Norm[,"24hB"]< -1,]))
?intersect  # intersect gives you the values that are common to sets A and B (i.e. only rownames (= genes) that are upregulated (>1) both in 24hA and 24hB

#Create matrix of up and downregulated genes
QCmatrixReg24 <- QCmatrix.Norm[c(Up24hAandB,Down24hAndB),]

#Barplot of specific genes quantification over the samples
par(mfrow=c(1,2))
barplot(QCmatrixReg24["CDKN1B",], main="CDKN1B", col="green", ylab = "relative quantity(log2)", cex.names=0.5)
barplot(QCmatrixReg24["CCND1",], main="CCND1", col="red", ylab = "relative quantity(log2)", cex.names=0.5)

#Create pdf-file of barplot stored in working directory
pdf(file="CDKN1BandCCND1.pdf", width=10, height=5)
par(mfrow=c(1,2))
barplot(QCmatrix.Norm["CDKN1B",], main="CDKN1B", col="green", ylab = "relative quantity(log2)", cex.names=0.5)
barplot(QCmatrix.Norm["CCND1",], main="CCND1", col="red", ylab = "relative quantity(log2)", cex.names=0.5)
dev.off()




#EXERCISE 2: Breast Cancer (BC) PAM 50 data set

#Read gene expression data (PAM50) for a cohort of BC samples
BCmatrix <- as.matrix(read.table("BCmatrixR.txt",header=T, row.names=1, check.names = FALSE, na.strings=c("", "NA", "Inf"), sep="\t"))
BCSampleInfo <- read.table("BreastInfoR.txt",header=T, row.names=1, check.names = FALSE, na.strings=c("", "NA", "Inf"), sep="\t")
#control that data matrix and sample info is matched
identical(colnames(BCmatrix), rownames(BCSampleInfo))
?identical # tests two objects (colnames vs rownames is this case) if they are exactly equal

#BC Heatmap
BiocManager::install("gplots")
library("gplots")
library("RColorBrewer")
my_palette <- colorRampPalette(c("deepskyblue3", "white", "violetred2"))(n = 299)
gr.col <- as.numeric(BCSampleInfo$subtype)
col.SubTypes <- c("pink","aquamarine3","greenyellow", "yellow", "orange", "red", "darkred")

#Euclidian
pdf(file="BC_Heatmap_Euclidian.pdf", width=10, height=7)
heatmap.2(BCmatrix, main="BC_Heatmap", sub="Euclidian, scale by row", trace="none", symbreaks=TRUE, margins = c(10, 20), 
          lhei = c(2, 8), lwid = c(8, 20), symkey=TRUE, col=my_palette, na.color="gray80", cexCol = 0.3, cexRow = 0.5, density.info="none", 
          dendrogram="both", ColSideColors=col.SubTypes[gr.col], scale="row")
legend("topright", legend=unique(BCSampleInfo$subtype),fill=col.SubTypes[unique(gr.col)])
dev.off()

#Pearson Ward
rowv <- as.dendrogram(hclust (as.dist(1-cor(t(BCmatrix), method ="pearson", use = "pairwise.complete.obs")), method="ward.D2"))
colv <- as.dendrogram(hclust (as.dist(1-cor(BCmatrix , method ="pearson", use = "pairwise.complete.obs")), method="ward.D2"))
pdf(file="BC_Heatmap_Pearson.pdf", width=10, height=7)
heatmap.2(BCmatrix, main="BC_Heatmap", sub="Pearson,WardD2, scale by row", trace="none", symbreaks=TRUE, margins = c(10, 20), 
          lhei = c(2, 8), lwid = c(8, 20), symkey=TRUE, col=my_palette, na.color="gray80", cexCol = 0.1, cexRow = 0.1, Rowv=rowv, density.info="none", 
          dendrogram="both", Colv=colv,ColSideColors=col.SubTypes[gr.col], scale="row")
legend("topright", legend=unique(BCSampleInfo$subtype),fill=col.SubTypes[unique(gr.col)])
dev.off()

#Define LuminalA samples and non-LuminalA samples
LumA <- rownames(BCSampleInfo[BCSampleInfo$subtype=="LumA",])
NonLumA <- setdiff(rownames(BCSampleInfo),LumA)
library(matrixStats)
?setdiff # returns a list of elements/vectors that are present in A but not in B. The order matters!

#Calculate fold difference between LuminalA and non-LuminalA samples
LumAFold <- rowMeans(BCmatrix[,LumA],na.rm = TRUE)-rowMeans(BCmatrix[,NonLumA],na.rm = TRUE)

#Calculate p-values and adjusted p-values (Benjamini-Hochberg) for LuminalA vs non-LuminalA comparison
t.result <- apply(BCmatrix, 1, function(x) t.test(x[LumA],x[NonLumA]))
 ?apply # Effectively returns a vector or array or list (PAM50 GENE NAMES) of values obtained by applying a function to 
        # margins of an array or matrix (BCmatrix). For a matrix, a margin=1 indicates rows. https://www.guru99.com/r-apply-sapply-tapply.html
        # Basically, it takes the BC matrix, performs the calculation by rows and applies a t-test on LumA vs NonLumA columns.
?t.test # Performs one and two sample t-tests on vectors of data.
t.test(BCmatrix[,LumA],BCmatrix[,NonLumA])  #doesn't work

p_values <- unlist(lapply(t.result, function(x) x$p.value))
 ?lapply # lapply(X, FUN, ...): Returns a list of the same length as X (=t.result), each element of which is the result of applying FUN to the corresponding element of X.
 ?unlist # changes a list into a vector?
         # Effectively, extracts p-values from t.result (using lapply) and writes them down as a vector (unlist)
adj.p_values <- p.adjust(p_values, method = "fdr")
  ?p.adjust # Requires a numeric vector of p-values, returns adjusted p-values using the chosen method.
            # Available methods: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
            # "fdr" = "BH"
 
#Create results table of analysis
LumA_results <- data.frame(fc=LumAFold, p.value=p_values, adj.p.value=adj.p_values,row.names=rownames(BCmatrix))
#Sort results table by fc
LumA_results_sorted <- LumA_results[order(LumA_results$fc,decreasing=TRUE),]

#View Luminal A top5 genes
LumA_results_sorted[1:5,]

#Volcano Plot
fc <- LumA_results$fc
p.values <- LumA_results$adj.p.value
neg.log10.pvalue <- -log10(p.values)
max.fc <- round(max(abs(fc),na.rm=TRUE))
max.p <- round(max(neg.log10.pvalue,na.rm=TRUE))
pdf(file="LumA_Volcano.pdf", width=10, height=7)
plot(fc,neg.log10.pvalue,xlim=c(-max.fc,max.fc),main="LuminalA BC DE analysis", xlab="fold difference (log2)",ylab="adj.p-value (-log10)",pch=20, cex=1, col="green4")
p.threshold = -log10(0.01)
abline(h=p.threshold,v=c(-1,1),lty=2,col="red")
dev.off()

#HER2 samples diff analysis
Her2 <- rownames(BCSampleInfo[BCSampleInfo$subtype=="Her2",])
NonHer2 <- setdiff(rownames(BCSampleInfo),Her2)
library(matrixStats)

#Calculate fold difference between Her2 and non-Her2 samples
Her2Fold <- rowMeans(BCmatrix[,Her2],na.rm = TRUE)-rowMeans(BCmatrix[,NonHer2],na.rm = TRUE)

#Calculate p-values and adjusted p-values (Benjamini-Hochberg) for Her2 vs non-Her2 comparison
t.result <- apply(BCmatrix, 1, function (x) t.test(x[Her2],x[NonHer2]))
p_values <- unlist(lapply(t.result, function(x) x$p.value))
adj.p_values <- p.adjust(p_values, method = "fdr")

t.result.df <- data.frame("pvalue" = apply(BCmatrix, 1, function(x) t.test(x[LumA],x[NonLumA])$p.value))

#Results Table
Her2_results <- data.frame(fc=Her2Fold, p.value=p_values, adj.p.value=adj.p_values,row.names=rownames(BCmatrix))
#Sort results table by fc
Her2_results_sorted <- Her2_results[order(Her2_results$fc,decreasing=TRUE),]

#View Her2 top5 genes
Her2_results_sorted[1:5,]

#Volcano Plot
fc <- Her2_results$fc
p.values <- Her2_results$adj.p.value
neg.log10.pvalue <- -log10(p.values)
max.fc <- round(max(abs(fc),na.rm=TRUE))
max.p <- round(max(neg.log10.pvalue,na.rm=TRUE))
pdf(file="Her2_Volcano.pdf", width=10, height=7)
plot(fc,neg.log10.pvalue,xlim=c(-max.fc,max.fc),main="Her2 BC DE analysis", xlab="fold difference (log2)",ylab="adj.p-value (-log10)",pch=20, cex=1, col="deepskyblue3")
p.threshold = -log10(0.01)
abline(h=p.threshold,v=c(-1,1),lty=2,col="red")
dev.off()






#EXERCISE 3: Cancer Cell Line Encyclopedia (CCLE) gene expression analysis  

#read CCLE gene expression data 
CCLEmatrix <- as.matrix(read.table("CCLEmatrix.txt",header=T, row.names=1, check.names = FALSE, na.strings=c("", "NA", "Inf"), sep="\t"))
CCLEinfo <- read.table("CCLEinfo.txt",header=T, row.names=1, check.names = FALSE, na.strings=c("", "NA", "Inf"), sep="\t")

#Clustering and heatmap of cell lines based on gene expression
my_palette <- colorRampPalette(c("deepskyblue3", "white", "violetred2"))(n = 299)
gr.col <- as.numeric(CCLEinfo$Type)
col.CLType <- c("deepskyblue3","darkorchid3","greenyellow", "yellow", "orange", "red")

#Pearson Ward
rowv <- as.dendrogram(hclust (as.dist(1-cor(t(CCLEmatrix), method ="pearson")), method="ward.D2"))
colv <- as.dendrogram(hclust (as.dist(1-cor(CCLEmatrix , method ="pearson")), method="ward.D2"))
pdf(file="CCLE_Heatmap_Pearson.pdf", width=10, height=7)
heatmap.2(CCLEmatrix, main="CCLE_Heatmap", sub="Pearson,WardD2", trace="none", symbreaks=TRUE, margins = c(10, 20), 
          lhei = c(2, 8), lwid = c(8, 20), symkey=TRUE, col=my_palette, cexCol = 0.1, cexRow = 0.1, Rowv=rowv, density.info="none", 
          dendrogram="both", Colv=colv,ColSideColors=col.CLType[gr.col], scale="row")
legend("topright", legend=unique(CCLEinfo$Type),fill=col.CLType[unique(gr.col)], cex=0.5)
dev.off()
?hclust # hierarchical clustering. 
        # Different methods available: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid"

#Principal component analysis (PCA)
CCLE.pca <- prcomp(t(CCLEmatrix), scale. = TRUE)
score.df <- as.data.frame(CCLE.pca$x)
score.df$type <- CCLEinfo[rownames(score.df),"Type"]
?prcomp # Performs a principal components analysis on the given data matrix and returns the results as an object of class prcomp

#Skree plot
Contributions <- round(summary(CCLE.pca)$importance[2,] * 100, 2)
barplot(Contributions[1:20], las=3, main="Screeplot CCLE PCA", xlab="component", ylab="Contribution(%)", col="lightblue")
PC1contr <- Contributions[1]
PC2contr <- Contributions[2]

#Create PCA plot for first two components using ggplot2 package (needs to be installed)
library(ggplot2)
ggplot(data = score.df, aes(x = PC1, y = PC2, label = rownames(score.df))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(aes(color=type), size = 2) +
  xlab(sprintf("PC1 (%s%%)", PC1contr)) + ylab(sprintf("PC2 (%s%%)", PC2contr))+ 
  ggtitle("CCLE PC1 vs PC2")
?sprintf # print function. % is an escape character. s indicates a string (in this case, PC1contr - numeric). 
         # %% is an escape character+% that gives you "%"

biplot(CCLE.pca, cex=0.3)
? biplot # part of stats package. Slow and unclear. Plots observations on x and variables on y. 




#Extra task: What is the difference between the two main Lung_NSCLC clusters (see heatmap)?
#Extract NSCLC subset of the data
NSCLCInfo <- CCLEinfo[CCLEinfo$Type=="LUNG_NSCLC",]
NSCLCmatrix <- CCLEmatrix[,rownames(NSCLCInfo)]

#Heatmap to cluster NSCLC cellines based on gene expression
my_palette <- colorRampPalette(c("deepskyblue3", "white", "violetred2"))(n = 299)
rowv <- as.dendrogram(hclust (as.dist(1-cor(t(NSCLCmatrix), method ="pearson")), method="ward.D2"))
colv <- as.dendrogram(hclust (as.dist(1-cor(NSCLCmatrix, method ="pearson")), method="ward.D2"))
pdf(file="NSCLC_Heatmap_Pearson.pdf", width=10, height=7)
heatmap.2(NSCLCmatrix, main="NSCLC_Heatmap", sub="Pearson,WardD2", trace="none", symbreaks=TRUE, margins = c(10, 20), 
          lhei = c(2, 8), lwid = c(8, 20), symkey=TRUE, col=my_palette, cexCol = 0.1, cexRow = 0.1, Rowv=rowv, density.info="none", 
          dendrogram="both", Colv=colv, scale="row")
dev.off()

#divide NSCLC cell lines into two subgroups based on clustering
tree <- hclust (as.dist(1-cor(NSCLCmatrix , method ="pearson")), method="ward.D2")
treeCuts <- cutree(tree, 2)                  # cut the hierarchical clustering tree in 2 parts
Cluster1 <- names(treeCuts)[treeCuts==1]
Cluster2 <- names(treeCuts)[treeCuts==2]

#Differential expression analysis between the two NSCLC clusters
library(matrixStats)
C1vsC2Fold <- rowMeans(NSCLCmatrix[,Cluster1])-rowMeans(NSCLCmatrix[,Cluster2])

t.result <- apply(NSCLCmatrix, 1, function (x) t.test(x[Cluster1],x[Cluster2]))
p_values <- unlist(lapply(t.result, function(x) x$p.value))
adj.p_values <- p.adjust(p_values, method = "fdr")

#Results Table
NSCLC_DE_results <- data.frame(fc=C1vsC2Fold, p.value=p_values, adj.p.value=adj.p_values,row.names=rownames(NSCLCmatrix))
#results table sorted by adjusted p-value
NSCLC_DE_results_sorted <- NSCLC_DE_results[order(NSCLC_DE_results$adj.p.value,decreasing=FALSE),]

#Volcano Plot
fc <- NSCLC_DE_results$fc
p.values <- NSCLC_DE_results$adj.p.value
neg.log10.pvalue <- -log10(p.values)
max.fc <- round(max(abs(fc),na.rm=TRUE))
max.p <- round(max(neg.log10.pvalue,na.rm=TRUE))
pdf(file="NSCLC_DE_Clusters_Volcano.pdf", width=10, height=7)
plot(fc,neg.log10.pvalue,xlim=c(-max.fc,max.fc),main="NSCLC clusters DE analysis", xlab="fold difference (log2)",ylab="adj.p-value (-log10)",pch=20, cex=1, col="orange1")
p.threshold = -log10(0.01)
abline(h=p.threshold,v=c(-1,1),lty=2,col="red")
dev.off()

#View top10 genes
NSCLC_DE_results_sorted[1:10,]

#Save results table as txt file
write.table(NSCLC_DE_results_sorted, file="NSCLC_DE_analysis.txt",row.names=T,quote=F,sep="\t")


