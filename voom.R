library("limma")
library("dplyr")
library("edgeR")
library("ggplot2")


##For ChAT

counts <- read.table("ChATcounts.txt", header = TRUE)

newCounts <- counts[apply(counts[,-1], 1, function(x) !all(x==0)),]
head(newCounts)
str(newCounts)
rownames(newCounts) <- newCounts[,1]
newCounts[,1] <- NULL
head(newCounts)

newCountsMatrix <- data.matrix(newCounts)
str(newCountsMatrix)
head(newCountsMatrix)

group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=newCountsMatrix,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)

dge <- DGEList(counts=newCountsMatrix)
dge <- calcNormFactors(dge)
isexpr <- rowSums(cpm(dge) >10) >= 2
dge <- dge[isexpr,]

vwts <- voomWithQualityWeights(dge, design, normalization="none", plot=TRUE)
vfit2 <- lmFit(vwts)
vfit2 <- eBayes(vfit2)
topTable(vfit2,coef=2, sort.by="P")


##v <- voom(dge,design,plot=TRUE)

##v <- voom(dge,design,plot=TRUE,normalize="quantile")

##fit <- lmFit(v,design)
##fit <- eBayes(fit)

## topTable(fit,coef=ncol(design))
top <- topTable(vfit2,coef=2,number=Inf,sort.by="P")

head(top)
sum(top$adj.P.Val<0.1)
plotMDS(vwts, xlim=c(-2,2), ylim=c(-1,1))
write.csv(top, "ChAT_voom.csv")

##Construct the plot object
g = ggplot(data=top, aes(x=logFC, y=-log10(P.Value), color=top$P.Value<=0.05)) +
      geom_point(alpha=0.4, size=5) +
      theme(legend.position = "none") +
      xlim(c(-1.8, 1.8)) + ylim(c(0, 3.5)) +
      xlab("log2 fold change") + ylab("-log10 p-value")
g
ggsave("ChAT_voom.png")


## for Drd1

counts <- read.table("Drd1counts.txt", header = TRUE)

newCounts <- counts[apply(counts[,-1], 1, function(x) !all(x==0)),]
head(newCounts)
str(newCounts)
rownames(newCounts) <- newCounts[,1]
newCounts[,1] <- NULL
head(newCounts)

newCountsMatrix <- data.matrix(newCounts)
str(newCountsMatrix)
head(newCountsMatrix)

group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=newCountsMatrix,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)

dge <- DGEList(counts=newCountsMatrix)
dge <- calcNormFactors(dge)
isexpr <- rowSums(cpm(dge) >10) >= 2
dge <- dge[isexpr,]

vwts <- voomWithQualityWeights(dge, design, normalization="none", plot=TRUE)
vfit2 <- lmFit(vwts)
vfit2 <- eBayes(vfit2)
topTable(vfit2,coef=2, sort.by="P")


##v <- voom(dge,design,plot=TRUE)

##v <- voom(dge,design,plot=TRUE,normalize="quantile")

##fit <- lmFit(v,design)
##fit <- eBayes(fit)

## topTable(fit,coef=ncol(design))
top <- topTable(vfit2,coef=2,number=Inf,sort.by="P")

head(top)
sum(top$adj.P.Val<0.1)
plotMDS(vwts, xlim=c(-2,2), ylim=c(-1,1))
write.csv(top, "Drd1_voom.csv")

##Construct the plot object
g = ggplot(data=top, aes(x=logFC, y=-log10(P.Value), color=top$P.Value<=0.05)) +
      geom_point(alpha=0.4, size=5) +
      theme(legend.position = "none") +
      xlim(c(-1.8, 1.8)) + ylim(c(0, 3.5)) +
      xlab("log2 fold change") + ylab("-log10 p-value")
g
ggsave("Drd1_voom.png")



## For Drd2

counts <- read.table("Drd2counts.txt", header = TRUE)

newCounts <- counts[apply(counts[,-1], 1, function(x) !all(x==0)),]
head(newCounts)
str(newCounts)
rownames(newCounts) <- newCounts[,1]
newCounts[,1] <- NULL
head(newCounts)

newCountsMatrix <- data.matrix(newCounts)
str(newCountsMatrix)
head(newCountsMatrix)

group <- factor(c(1,1,1,2,2,2,2,2))
y <- DGEList(counts=newCountsMatrix,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)

dge <- DGEList(counts=newCountsMatrix)
dge <- calcNormFactors(dge)
isexpr <- rowSums(cpm(dge) >10) >= 2
dge <- dge[isexpr,]

vwts <- voomWithQualityWeights(dge, design, normalization="none", plot=TRUE)
vfit2 <- lmFit(vwts)
vfit2 <- eBayes(vfit2)
topTable(vfit2,coef=2, sort.by="P")


##v <- voom(dge,design,plot=TRUE)

##v <- voom(dge,design,plot=TRUE,normalize="quantile")

##fit <- lmFit(v,design)
##fit <- eBayes(fit)

## topTable(fit,coef=ncol(design))
top <- topTable(vfit2,coef=2,number=Inf,sort.by="P")

head(top)
sum(top$adj.P.Val<0.1)
plotMDS(vwts, xlim=c(-2,2), ylim=c(-1,1))
write.csv(top, "Drd2_voom.csv")

##Construct the plot object
g = ggplot(data=top, aes(x=logFC, y=-log10(P.Value), color=top$P.Value<=0.05)) +
      geom_point(alpha=0.4, size=5) +
      theme(legend.position = "none") +
      xlim(c(-1.8, 1.8)) + ylim(c(0, 3.5)) +
      xlab("log2 fold change") + ylab("-log10 p-value")
g

ggsave("Drd2_voom.png")



## For Cogalt2
counts <- read.table("Cogalt2counts.txt", header = TRUE)

newCounts <- counts[apply(counts[,-1], 1, function(x) !all(x==0)),]
head(newCounts)
str(newCounts)
rownames(newCounts) <- newCounts[,1]
newCounts[,1] <- NULL
head(newCounts)

newCountsMatrix <- data.matrix(newCounts)
str(newCountsMatrix)
head(newCountsMatrix)

group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=newCountsMatrix,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)

dge <- DGEList(counts=newCountsMatrix)
dge <- calcNormFactors(dge)
isexpr <- rowSums(cpm(dge) >10) >= 2
dge <- dge[isexpr,]

vwts <- voomWithQualityWeights(dge, design, normalization="none", plot=TRUE)
vfit2 <- lmFit(vwts)
vfit2 <- eBayes(vfit2)
topTable(vfit2,coef=2, sort.by="P")


##v <- voom(dge,design,plot=TRUE)

##v <- voom(dge,design,plot=TRUE,normalize="quantile")

##fit <- lmFit(v,design)
##fit <- eBayes(fit)

## topTable(fit,coef=ncol(design))
top <- topTable(vfit2,coef=2,number=Inf,sort.by="P")

head(top)
sum(top$adj.P.Val<0.1)
plotMDS(vwts, xlim=c(-2,2), ylim=c(-1,1))
write.csv(top, "Cogalt2_voom.csv")

##Construct the plot object
g = ggplot(data=top, aes(x=logFC, y=-log10(P.Value), color=top$P.Value<=0.05)) +
      geom_point(alpha=0.4, size=5) +
      theme(legend.position = "none") +
      xlim(c(-1.8, 1.8)) + ylim(c(0, 3.5)) +
      xlab("log2 fold change") + ylab("-log10 p-value")
g
ggsave("Cogalt2_voom.png")

## For S100a10
counts <- read.table("S100a10counts.txt", header = TRUE)

newCounts <- counts[apply(counts[,-1], 1, function(x) !all(x==0)),]
head(newCounts)
str(newCounts)
rownames(newCounts) <- newCounts[,1]
newCounts[,1] <- NULL
head(newCounts)

newCountsMatrix <- data.matrix(newCounts)
str(newCountsMatrix)
head(newCountsMatrix)

group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=newCountsMatrix,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)

dge <- DGEList(counts=newCountsMatrix)
dge <- calcNormFactors(dge)
isexpr <- rowSums(cpm(dge) >10) >= 2
dge <- dge[isexpr,]

vwts <- voomWithQualityWeights(dge, design, normalization="none", plot=TRUE)
vfit2 <- lmFit(vwts)
vfit2 <- eBayes(vfit2)
topTable(vfit2,coef=2, sort.by="P")


##v <- voom(dge,design,plot=TRUE)

##v <- voom(dge,design,plot=TRUE,normalize="quantile")

##fit <- lmFit(v,design)
##fit <- eBayes(fit)

## topTable(fit,coef=ncol(design))
top <- topTable(vfit2,coef=2,number=Inf,sort.by="P")

head(top)
sum(top$adj.P.Val<0.1)
plotMDS(vwts, xlim=c(-2,2), ylim=c(-1,1))
write.csv(top, "S100a10_voom.csv")

##Construct the plot object
g = ggplot(data=top, aes(x=logFC, y=-log10(P.Value), color=top$P.Value<=0.05)) +
      geom_point(alpha=0.4, size=5) +
      theme(legend.position = "none") +
      xlim(c(-1.8, 1.8)) + ylim(c(0, 3.5)) +
      xlab("log2 fold change") + ylab("-log10 p-value")
g
ggsave("S100a10_voom.png")
