# Josue Vega
# B. cinerea transcripts on Arabidopsis
# Spearman Correlations, test statistic, p-value, + bonferroni corrected p values
# Then only p<.05 for bonferroni

rm(list=ls())
library(dplyr)
library(readr)
setwd("~/../Desktop/B. cinera/Hwaviness/")
HwaviBC <- read.csv("../WavyGWAS_lsmeans.fxmod1_R_output.csv") #contains lsmeans instead of average
HwaviBC <- HwaviBC [,c("Isolate", "HwaviEstimate")] #elimating excess columns for merging

rawTranscripts <- read.csv("../Hwaviness/data/BcAtTranscripts/lsmeans_zscale_allreads.csv")#Phentype data (over9k entries)

mergeTranscripts <- merge(HwaviBC, rawTranscripts, by = "Isolate")

BcAt.corr <- NULL
BcAt.corr <- as.data.frame(BcAt.corr)
BcAt.corr <-data.frame(mergeTranscripts$Isolate)
names(BcAt.corr)<-c("p.value","rho.value")

for (i in 4:ncol(mergeTranscripts)){ #start at 4 since merge scooted it over
  mycol <- mergeTranscripts[i,]
  
  PvalueSpea <- cor.test(mergeTranscripts$HwaviEstimate,mergeTranscripts[[i]], method="spearman", exact=FALSE)$p.value
  PvalueSpea <- format(round(PvalueSpea, 5), nsmall = 4)
  rhovalueSpea <- cor.test(mergeTranscripts$HwaviEstimate,mergeTranscripts[[i]], method="spearman", exact=FALSE)$estimate
  rhovalueSpea <- format(round(rhovalueSpea, 5), nsmall = 4)

  BcAt.corr <- paste(BcAt.corr, PvalueSpea)
}
