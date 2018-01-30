#Josue Vega
# Trying new graphical packages
# focus on regression (pearson/spearman)
#################################
#################################

# Building Table of Top 30 SNP Effect Size
## Row 1 = Positive or Negative SNP Effects
## Row 2 = SNP Effect Sizes (absolute Value)
## Column 1 = Isolate Name
### Chart filled with 1 or 0 if the isolate contains the trait 
rm(list=ls())
library(ggplot2)
library(dplyr)

setwd("~/../Desktop/B. cinera/")
#Column Descending Values of SNP effect estimate
EffectEstimate <- read.csv("Hwaviness/data/04_bigRRoutput/trueMAF20_20NA/HWavi_trueMAF20_20NA.HEM.PlotFormat.csv")

EffectEstimate <- EffectEstimate[c(2,3,5,9)]
#EffectEstimateAbs <- abs(EffectEstimate)
EffectEstimateAbsOrg <- arrange(abs(EffectEstimate), desc(Estimate))
EffectEstimateAbsOrg <- EffectEstimateAbsOrg[(1:30),]
NewVar <- merge(EffectEstimate, EffectEstimateAbsOrg, by = "Index", all.y=TRUE)
#Adding gene IDs
geneIDs <- read.csv("Hwaviness/geneID_fxn_estimate.csv")
geneIDs <- geneIDs[c(2, 12)]
NewVar <- merge(geneIDs, NewVar, by = "Index", all.y=TRUE)
EffectEstimateCom <- NewVar[,-c(6:7)]
EffectEstimate.list <- as.list(as.data.frame(EffectEstimateCom))

#row of Isolates
IsolateNames <- read.csv("Hwaviness/data/03_bigRRinput/Domestication/hpbinSNP_bigRR_trueMAF20_50NA.csv")
IsolateNames <- IsolateNames[(5:97)]#

testest <- merge(IsolateNames, EffectEstimateCom, by="Index", all.y = TRUE)
Top30 <- t(testest)
write.csv(Top30, "Hwaviness/Top30Est_Bc_JV.csv")

View(testest)
View(Top30)


######
  









#################################
#################################
rm(list=ls())
library(ggplot2)
library(dplyr)

setwd("~/../Desktop/B. cinera/")
HwaviBC <- read.csv("WavyGWAS_lsmeans.fxmod1_R_output.csv")
Hwavi <- HwaviBC
AvgNorm <- read.csv("BR_JV_ManualHyphalDat_032817(edit).csv")

Isolate <- AvgNorm %>% 
  group_by(Isolate) %>%
  summarise(avg_pheno = mean(Phenotype, na.rm = TRUE), 
            min_pheno = min(Phenotype, na.rm = TRUE), 
            max_pheno = max(Phenotype, na.rm = TRUE),
            sd_pheno = sd(Phenotype, na.rm = TRUE),
            total = n())

#Dot Plot of Isolates (increasing distr)
p <- ggplot(Hwavi, aes(x = reorder(Isolate, HwaviEstimate), y=HwaviEstimate))  +
  theme_bw() + 
  geom_point(colour="black", stat="identity") +
  xlab("Isolate") + ylab("Hyphal Waviness") +
  theme(axis.text.x=element_text(face = "bold", color = "black", size = 10.5, angle = -60, hjust = 0))
print(p)

#Violin Plot of LsMeans
p <- ggplot(Hwavi, aes(x=NA, y=HwaviEstimate))  +
  theme_bw() + 
  theme(axis.text.x=element_blank()) +
  geom_violin(colour="black", fill = "grey", scale="area") +
  geom_boxplot(width=.1) +
  xlab("Isolate") + ylab("Hyphal Waviness") 
print(p)

#Density Plot of LsMeans
plot(density(Hwavi$HwaviEstimate), main="Frequency Distribution of LsMeans of Hyphal Waviness")
#plot(density(Isolate$avg_pheno))

