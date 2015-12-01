library("limma")
library("dplyr")
library("edgeR")
library("ggplot2")
library("xlsx")
?read.csv
setwd("~/Documents/Cell type R/")
ChAT_CuffDiff2 <- read.table("gene_exp.diff", header = TRUE)
head(ChAT_CuffDiff2)
ChAT_CuffDiffNOzeros <- subset(ChAT_CuffDiff2, ChAT_CuffDiff2$value_2 != "0")
ChAT_CuffDiffSIG <- subset(ChAT_CuffDiff2, ChAT_CuffDiff2$q_value <=0.05)
ChAT_CuffDiffSIGnoZEROs <- subset(ChAT_CuffDiffNOzeros, ChAT_CuffDiffNOzeros$q_value <=0.05)

##Construct the plot object
g = ggplot(data=ChAT_CuffDiffNOzeros, aes(x=ChAT_CuffDiffNOzeros$log2.fold_change., y=-log10(ChAT_CuffDiffNOzeros$p_value), color=ChAT_CuffDiffNOzeros$p_value<=0.05)) +
      geom_point(alpha=0.3, size=5) +
      theme(legend.position = "none") +
      xlim(c(-5, 5)) + ylim(c(0.1, 5)) +
      xlab("log2 fold change") + ylab("-log10 p-value")
g
ggsave("ChAT_CuffDiff2.png")

ChAT_CuffDiff <- read.table("Galaxy17-[Chat-Cuffdiff-default-Refseq-_on_data_2,_data_7,_and_others__transcript_differential_expression_testing].tabular", header = TRUE)
head(ChAT_CuffDiff)


##Construct the plot object
g = ggplot(data=ChAT_CuffDiff, aes(x=log2.fold_change., y=-log10(p_value), color=p_value<=0.05)) +
      geom_point(alpha=0.3, size=5) +
      theme(legend.position = "none") +
      xlim(c(-5, 5)) + ylim(c(0.1, 5)) +
      xlab("log2 fold change") + ylab("-log10 p-value")
g
ggsave("ChAT_CuffDiff.png")

Drd1_CuffDiff <- read.table("Galaxy31-Drd1[Cuffdiff_on_data_4,_data_2,_and_others__gene_differential_expression_testing].tabular", header = TRUE)
head(Drd1_CuffDiff)


##Construct the plot object
g = ggplot(data=Drd1_CuffDiff, aes(x=log2.fold_change., y=-log10(p_value), color=p_value<=0.05)) +
      geom_point(alpha=0.3, size=5) +
      theme(legend.position = "none") +
      xlim(c(-5, 5)) + ylim(c(0.1, 5)) +
      xlab("log2 fold change") + ylab("-log10 p-value")
g
ggsave("Drd1_CuffDiff.png")


Drd2_CuffDiff <- read.table("Galaxy42-Drd2[Cuffdiff_on_data_8,_data_10,_and_others__gene_differential_expression_testing].tabular", header = TRUE)
head(Drd2_CuffDiff)


##Construct the plot object
g = ggplot(data=Drd2_CuffDiff, aes(x=log2.fold_change., y=-log10(p_value), color=p_value<=0.05)) +
      geom_point(alpha=0.3, size=5) +
      theme(legend.position = "none") +
      xlim(c(-5, 5)) + ylim(c(0.1, 5)) +
      xlab("log2 fold change") + ylab("-log10 p-value")
g
ggsave("Drd2_CuffDiff.png")

S100a10_CuffDiff <- read.table("Galaxy131-[S100a10-Cuffdiff_-default-Refseq-on_data_7,_data_123,_and_others__gene_differential_expression_testing].tabular", header = TRUE)
head(S100a10_CuffDiff)


##Construct the plot object
g = ggplot(data=S100a10_CuffDiff, aes(x=log2.fold_change., y=-log10(p_value), color=p_value<=0.05)) +
      geom_point(alpha=0.3, size=5) +
      theme(legend.position = "none") +
      xlim(c(-5, 5)) + ylim(c(0.1, 5)) +
      xlab("log2 fold change") + ylab("-log10 p-value")
g
ggsave("S100a10_CuffDiff.png")


Cogalt2_CuffDiff <- read.table("Galaxy37-[Glt25d2-Cuffdiff-default-refseq-_on_data_6,_data_7,_and_others__gene_differential_expression_testing].tabular", header = TRUE)
head(Cogalt2_CuffDiff)


##Construct the plot object
g = ggplot(data=Cogalt2_CuffDiff, aes(x=log2.fold_change., y=-log10(p_value), color=p_value<=0.05)) +
      geom_point(alpha=0.3, size=5) +
      theme(legend.position = "none") +
      xlim(c(-5, 5)) + ylim(c(0.1, 5)) +
      xlab("log2 fold change") + ylab("-log10 p-value")
g
ggsave("Cogalt2_CuffDiff.png")
