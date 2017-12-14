#Josue Vega
#plotting from bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
## customize to your working directory
setwd("~/../Desktop/B. cinera/Hwaviness/")
#setwd("~/Projects/BcSolGWAS/")

## customize these notes to match your file names
#Input file: Sl_DomesticationLS_MAF20.HEM.PlotFormat.csv AND Sl_DomesticationLS_MAF20.HEM.Thresh.csv
#Output file: NONE
#Plots: Basic and greyscale domestication Manhattan plots
#go to script 09_bigRR_meta for domestication color plot (fig R8)
###########################################################################
#Plotting the HEM results

#Load plotting package
## install these packages if you do not have them
library(ggplot2); library(grid); library(plyr)

#Import data (reorganized from script ReformatBigRRouts.R)
## customize this to the output files from script 04
HEM.plotdata <- read.csv("data/04_bigRRoutput/trueMAF_10NA/SlBc_domest_trueMAF20_10NA.HEM.PlotFormat.csv")
## if the first column of HEM.plotdata is a nameless column of numbers, remove it here:
HEM.plotdata <- HEM.plotdata[,-c(1)]

#get threshhold values 
## customize this to the threshold files from script 04
HEM.thresh <- read.csv("data/04_bigRRoutput/trueMAF_10NA/SlBc_domest_trueMAF20_10NA.HEM.Thresh.csv")
## if the first column of HEM.thresh is a nameless column of numbers, remove it here:
HEM.thresh <- HEM.thresh[,-c(1)]

## visually check your threshold file: you want the 99% positive threshold row saved here:
TH99pos <- HEM.thresh[3,]
for (i in 2:ncol(TH99pos)){
  assign(paste("TH99pos_", names(TH99pos[i]), sep=""),as.numeric(TH99pos[i]))
}

## visually check your threshold file: you want the 99.9% positive threshold row saved here:
TH999pos <- HEM.thresh[4,]
for (i in 2:ncol(TH999pos)){
  assign(paste("TH999pos_", names(TH999pos[i]), sep=""),as.numeric(TH999pos[i]))
}

## visually check your threshold file: you want the 95% positive threshold row saved here:
TH95pos <- HEM.thresh[1,]
for (i in 2:ncol(TH95pos)){
  assign(paste("TH95pos_", names(TH95pos[i]), sep=""),as.numeric(TH95pos[i]))
}

## visually check your threshold file: you want the 99% negative threshold row saved here:
TH99neg <- HEM.thresh[7,]
for (i in 2:ncol(TH99neg)){
  assign(paste("TH99neg_", names(TH99neg[i]), sep=""),as.numeric(TH99neg[i]))
}

## visually check your threshold file: you want the 99.9% negative threshold row saved here:
TH999neg <- HEM.thresh[8,]
for (i in 2:ncol(TH999neg)){
  assign(paste("TH999neg_", names(TH999neg[i]), sep=""),as.numeric(TH999neg[i]))
}

## visually check your threshold file: you want the 95% negative threshold row saved here:
TH95neg <- HEM.thresh[5,]
for (i in 2:ncol(TH95neg)){
  assign(paste("TH95neg_", names(TH95neg[i]), sep=""),as.numeric(TH95neg[i]))
}

#create a custom color scale
## this will make each chromosome alternating colors of light and dark grey
myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60")
names(myColors) <- levels(HEM.plotdata$Chrom)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

## here you draw your first plot, using the 95% and 99% thresholds
#lowTR
#jpeg(paste("plots/paper/SlBc_trueMAF20_10NA_lowTR_manhattan",names(HEM.plotdata[i]),".jpg", sep=""), width=8, height=5, units='in', res=600)

##replace [,i] and [i] with the column number for "estimate" aka the effect size estimate for each SNP
print(ggplot(HEM.plotdata, aes(x=Index, y=HEM.plotdata[,i]))+
  theme_bw()+
  colScale+
  geom_point(aes(color = factor(Chrom)))+
  labs(list(y="SNP Effect Estimate", title=("YOUR TITLE HERE")))+
  guides(col = guide_legend(nrow = 8, title="Chromosome"))+
  geom_hline(yintercept=get(paste("TH95neg_", names(HEM.plotdata[i]), sep="")), colour = "black", lty=2) +
  geom_hline(yintercept=get(paste("TH95pos_", names(HEM.plotdata[i]), sep="")), colour = "black", lty=2) +
  geom_text(aes(0,get(paste("TH95neg_", names(HEM.plotdata[i]), sep="")), label = "95% Threshold", vjust = 1.2, hjust = .05), col = "black")+
  geom_hline(yintercept=get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), lty=2) +
  geom_hline(yintercept=get(paste("TH99pos_", names(HEM.plotdata[i]), sep="")), lty=2) +
  geom_text(aes(0,get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), label = "99% Threshold", vjust = 1.5, hjust = .05), col = "black")+
  theme(legend.position="none")+
  #scale_x_continuous(name="Chromosome", breaks = c(1674920, 5242762, 8999082, 11057880, 13580364, 17181767, 20009659, 22371413, 24388631, 26765466, 28559407, 30104939, 31864287, 33980281, 35775028, 38876579), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
  expand_limits(y=0))
#dev.off()

## here you draw your second plot, using the 99.9% threshold
#highTR

##replace [,i] and [i] with the column number for "estimate" aka the effect size estimate for each SNP
print(ggplot(HEM.plotdata, aes(x=Index, y=HEM.plotdata[,i]))+
        theme_bw()+
        colScale+
        geom_point(aes(color = factor(Chrom)))+
        labs(list(y="SNP Effect Estimate", title=("YOUR TITLE HERE")))+
        guides(col = guide_legend(nrow = 8, title="Chromosome"))+
        geom_hline(yintercept=get(paste("TH999neg_", names(HEM.plotdata[i]), sep="")), colour = "black", lty=2) +
        geom_hline(yintercept=get(paste("TH999pos_", names(HEM.plotdata[i]), sep="")), colour = "black", lty=2) +
        geom_text(aes(0,get(paste("TH999neg_", names(HEM.plotdata[i]), sep="")), label = "99.9% Threshold", vjust = 1.2, hjust = .05), col = "black")+
        theme(legend.position="none")+
        #scale_x_continuous(name="Chromosome", breaks = c(1674920, 5242762, 8999082, 11057880, 13580364, 17181767, 20009659, 22371413, 24388631, 26765466, 28559407, 30104939, 31864287, 33980281, 35775028, 38876579), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
        expand_limits(y=0))
#dev.off()