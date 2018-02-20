# Josue Vega
# B. cinerea transcripts on Arabidopsis
# Spearman Correlations, test statistic, p-value, + bonferroni corrected p values
# Then only p<.05 for bonferroni

rm(list=ls())
library(dplyr)
library(readr)
setwd("~/../Desktop/B. cinera/Hwaviness/")
HwaviBC <- read.csv("WavyGWAS_lsmeans.fxmod1_R_output.csv")

rawTranscripts <- read.csv("../Hwaviness/data/BcAtTranscripts/lsmeans_zscale_allreads.csv")

BcAt.corr <- NULL

for (i in 1:nrow(AntOverrep)){
  myrow <- AntOverrep[i,]
  myA <- myrow$AllWavyGenFreq
  myB <- myrow$WGGenFreq
  myCmA <- sum(AntOverrep[,3])-myA
  myDmB <- sum(AntOverrep$WGGenFreq)-myB
  
  PvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$anac055, method="spearman", exact=FALSE)$p.value
  PvalueSpea <- format(round(PvalueSpea, 5), nsmall = 4)
  rhovalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$anac055, method="spearman", exact=FALSE)$estimate
  rhovalueSpea <- format(round(rhovalueSpea, 5), nsmall = 4)
  
  fisher.p.up <- fisher.test(matrix(c(myA, myB, myCmA, myDmB), nrow=2, ncol=2), alternative="greater")$p.value
  fisher.p.down <- fisher.test(matrix(c(myA, myB, myCmA, myDmB), nrow=2, ncol=2), alternative="less")$p.value
  BcAt.corr <- append(BcAt.corr, PvalueSpea)
  fisher.p.under <- append(fisher.p.under, fisher.p.down)
}
assign(paste("fisher.p.over", names(AntOverrep)[3], sep="."),fisher.p.over)
assign(paste("fisher.p.under", names(AntOverrep)[3], sep="."),fisher.p.under)