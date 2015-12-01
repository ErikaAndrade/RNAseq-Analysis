library("limma")
library("dplyr")
library("edgeR")
library("ggplot2")
## For ChAT

counts <- read.table("ChATcounts.txt", header = TRUE)

newCounts <- counts[apply(counts[,-1], 1, function(x) !all(x<=1)),]
head(newCounts)
str(newCounts)
rownames(newCounts) <- newCounts[,1]
newCounts[,1] <- NULL


newCountsMatrix <- data.matrix(newCounts)
str(newCountsMatrix)
head(newCountsMatrix)

dge <- DGEList(counts=newCountsMatrix)
dge <- calcNormFactors(dge)
isexpr <- rowSums(cpm(dge) >10) >= 2
dge <- dge[isexpr,]


group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=newCountsMatrix,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(dge,design)

fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)

##fit <- glmQLFit(y,design)
##qlf <- glmQLFTest(fit,coef=2)
##topTags(qlf)

plotMDS(y, xlim=c(-2,2), ylim=c(-1,1))
?topTags
top <- topTags(lrt, 300000)
TopTable <- top$table

write.csv(top$table, "ChAT_edgeR.csv")

##Construct the plot object
g = ggplot(data=TopTable, aes(x=logFC, y=-log10(PValue), color=TopTable$PValue<=0.05)) +
      geom_point(alpha=0.4, size=5) +
      theme(legend.position = "none") +
      xlim(c(-2, 2)) + ylim(c(0, 4)) +
      xlab("log2 fold change") + ylab("-log10 p-value")
g
ggsave("ChAT_edgeR.png")

## For Drd1

counts <- read.table("Drd1counts.txt", header = TRUE)

newCounts <- counts[apply(counts[,-1], 1, function(x) !all(x<=1)),]
head(newCounts)
str(newCounts)
rownames(newCounts) <- newCounts[,1]
newCounts[,1] <- NULL

newCountsMatrix <- data.matrix(newCounts)
str(newCountsMatrix)
head(newCountsMatrix)

dge <- DGEList(counts=newCountsMatrix)
dge <- calcNormFactors(dge)
isexpr <- rowSums(cpm(dge) >10) >= 2
dge <- dge[isexpr,]

group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=newCountsMatrix,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(dge,design)

fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)

##fit <- glmQLFit(y,design)
##qlf <- glmQLFTest(fit,coef=2)
##topTags(qlf)

plotMDS(y, xlim=c(-2,2), ylim=c(-1,1))
?topTags
top <- topTags(lrt, 300000)
TopTable <- top$table

write.csv(top$table, "Drd1_edgeR.csv")

##Construct the plot object
g = ggplot(data=TopTable, aes(x=logFC, y=-log10(PValue), color=TopTable$PValue<=0.05)) +
      geom_point(alpha=0.4, size=5) +
      theme(legend.position = "none") +
      xlim(c(-2, 2)) + ylim(c(0, 4)) +
      xlab("log2 fold change") + ylab("-log10 p-value")
g
ggsave("Drd1_edgeR.png")

## For Drd2

counts <- read.table("Drd2counts.txt", header = TRUE)

newCounts <- counts[apply(counts[,-1], 1, function(x) !all(x<=1)),]
head(newCounts)
str(newCounts)
rownames(newCounts) <- newCounts[,1]
newCounts[,1] <- NULL


newCountsMatrix <- data.matrix(newCounts)
str(newCountsMatrix)
head(newCountsMatrix)

dge <- DGEList(counts=newCountsMatrix)
dge <- calcNormFactors(dge)
isexpr <- rowSums(cpm(dge) >10) >= 2
dge <- dge[isexpr,]


group <- factor(c(1,1,1,2,2,2,2,2))
y <- DGEList(counts=newCountsMatrix,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(dge,design)

fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)

##fit <- glmQLFit(y,design)
##qlf <- glmQLFTest(fit,coef=2)
##topTags(qlf)

plotMDS(y, xlim=c(-2,2), ylim=c(-1,1))
?topTags
top <- topTags(lrt, 300000)
TopTable <- top$table

write.csv(top$table, "Drd2_edgeR.csv")

##Construct the plot object
g = ggplot(data=TopTable, aes(x=logFC, y=-log10(PValue), color=TopTable$PValue<=0.05)) +
      geom_point(alpha=0.4, size=5) +
      theme(legend.position = "none") +
      xlim(c(-2, 2)) + ylim(c(0, 4)) +
      xlab("log2 fold change") + ylab("-log10 p-value")
g
ggsave("Drd2_edgeR.png")

## Cogalt2
counts <- read.table("Cogalt2counts.txt", header = TRUE)

newCounts <- counts[apply(counts[,-1], 1, function(x) !all(x<=1)),]
head(newCounts)
str(newCounts)
rownames(newCounts) <- newCounts[,1]
newCounts[,1] <- NULL


newCountsMatrix <- data.matrix(newCounts)
str(newCountsMatrix)
head(newCountsMatrix)

dge <- DGEList(counts=newCountsMatrix)
dge <- calcNormFactors(dge)
isexpr <- rowSums(cpm(dge) >10) >= 2
dge <- dge[isexpr,]


group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=newCountsMatrix,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(dge,design)

fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)

##fit <- glmQLFit(y,design)
##qlf <- glmQLFTest(fit,coef=2)
##topTags(qlf)

plotMDS(y, xlim=c(-2,2), ylim=c(-1,1))
?topTags
top <- topTags(lrt, 300000)
TopTable <- top$table

write.csv(top$table, "Cogalt2_edgeR.csv")

##Construct the plot object
g = ggplot(data=TopTable, aes(x=logFC, y=-log10(PValue), color=TopTable$PValue<=0.05)) +
      geom_point(alpha=0.4, size=5) +
      theme(legend.position = "none") +
      xlim(c(-2, 2)) + ylim(c(0, 4)) +
      xlab("log2 fold change") + ylab("-log10 p-value")
g
ggsave("Cogalt2_edgeR.png")

## S100a10

counts <- read.table("S100a10counts.txt", header = TRUE)

newCounts <- counts[apply(counts[,-1], 1, function(x) !all(x<=1)),]
head(newCounts)
str(newCounts)
rownames(newCounts) <- newCounts[,1]
newCounts[,1] <- NULL

newCountsMatrix <- data.matrix(newCounts)
str(newCountsMatrix)
head(newCountsMatrix)

dge <- DGEList(counts=newCountsMatrix)
dge <- calcNormFactors(dge)
isexpr <- rowSums(cpm(dge) >10) >= 2
dge <- dge[isexpr,]

group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=newCountsMatrix,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(dge,design)

fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)

##fit <- glmQLFit(y,design)
##qlf <- glmQLFTest(fit,coef=2)
##topTags(qlf)

plotMDS(y, xlim=c(-2,2), ylim=c(-1,1))
?topTags
top <- topTags(lrt, 300000)
TopTable <- top$table

write.csv(top$table, "S100a10_edgeR.csv")

##Construct the plot object
g = ggplot(data=TopTable, aes(x=logFC, y=-log10(PValue), color=TopTable$PValue<=0.05)) +
      geom_point(alpha=0.4, size=5) +
      theme(legend.position = "none") +
      xlim(c(-2, 2)) + ylim(c(0, 4)) +
      xlab("log2 fold change") + ylab("-log10 p-value")
g
ggsave("S100a10_edgeR.png")
