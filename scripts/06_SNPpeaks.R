#Josue Vega
#plotting from bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
## set correct working directory
#setwd("~/Documents/GitRepos/BcSolGWAS/")
#setwd("~/../Desktop/B. cinera/Hwaviness/")

## update these file names!
#Input File: data/04_bigRRoutput/trueMAF20_20NA/HWavi_trueMAF20_20NA.HEM.PlotFormat.csv
#Output File: data/05_genes/Hwavi_TopSNPs_long.csv, data/05_genes/Hwavi_TopSNPs_wide.csv
#this goes into gene annotation and then venn diagrams
#Plots: NONE
############################################################################

#Load packages
## install any that are missing
library(plyr); library(ggplot2); library(grid)

#Import data
## check file names
HEM.plotdata <- read.csv("data/04_bigRRoutput/trueMAF20_20NA/HWavi_trueMAF20_20NA.HEM.PlotFormat.csv")

#get threshhold values 
HEM.thresh <- read.csv("data/04_bigRRoutput/trueMAF20_20NA/HWavi_trueMAF20_20NA.HEM.Thresh.csv")

#take the top 50 over the threshold for each phenotype
TH95pos <- HEM.thresh[1,]
TH95neg <- HEM.thresh[5,]
TH99pos <- HEM.thresh[3,]
TH99neg <- HEM.thresh[7,]
TH999pos <- HEM.thresh[4,]
TH999neg <- HEM.thresh[8,]


##check if you need to remove first column
#remove empty first column
HEM.plotdata <- HEM.plotdata[,-c(1)]

names(HEM.plotdata)

#conditionally replace values < threshold with zero
#for 99% Threshold
HEM.topSNPs_99Thr <- HEM.plotdata
HEM.topSNPs_99Thr$Estimate[HEM.topSNPs_99Thr$Estimate < TH99pos$Estimate & HEM.topSNPs_99Thr$Estimate > 0] <- 0
HEM.topSNPs_999Thr$Estimate[HEM.topSNPs_99Thr$Estimate > TH99neg$Estimate & HEM.topSNPs_99Thr$Estimate < 0] <- 0
#remove rows if estimate = 0 
#this keeps only the N rows with SNPs > 99% threshold
HEM.topSNPs_99Thr <- HEM.topSNPs_99Thr[!(HEM.topSNPs_99Thr$Estimate==0),]

#for 99.9% Threshold
HEM.topSNPs_999Thr <- HEM.plotdata
HEM.topSNPs_999Thr$Estimate[HEM.topSNPs_999Thr$Estimate < TH999pos$Estimate & HEM.topSNPs_999Thr$Estimate > 0] <- 0
HEM.topSNPs_999Thr$Estimate[HEM.topSNPs_999Thr$Estimate > TH999neg$Estimate & HEM.topSNPs_999Thr$Estimate < 0] <- 0
#remove rows if estimate = 0 
#this keeps only the N rows with SNPs > 99% threshold
HEM.topSNPs_999Thr <- HEM.topSNPs_999Thr[!(HEM.topSNPs_999Thr$Estimate==0),]

library(ggplot2)
plot1 <- ggplot(HEM.topSNPs_999Thr, aes(x=Index, y=Estimate))
plot1 + geom_point()+
  theme_bw()
  
#save it

##check these file names
write.csv(HEM.topSNPs_999Thr, "data/05_genes/Hwavi_TopSNPs999.csv")
write.csv(HEM.topSNPs_99Thr, "data/05_genes/Hwavi_TopSNPs99.csv")