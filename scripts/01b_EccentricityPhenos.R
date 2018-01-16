#01b_EccentricityPhenos.R
#NES
#01/12/18

#prepare lesion eccentricity phenotypes of B. cinerea on A. thaliana for correlation to B. cinerea hyphal waviness
#-------------------------------------------------------------------------
rm(list=ls())
setwd("~/../Desktop/B. cinera/Hwaviness/data/BcAtPhenos/")

isolist <- read.csv("IsolateKey.csv")
myphenos <- read.csv("BcAt_LesionImageData.csv")
names(myphenos)[6] <- "Sample"

allphenos <- merge(myphenos,isolist, by="Sample")

allphenos <- allphenos[,c("Experiment","GrowingFlat","AgarFlat","Plant","Isolate","Lesion.0.m.eccentricity","Lesion.red.m.eccentricity","Lesion.grn.m.eccentricity", "Lesion.blu.m.eccentricity", "Lesion.redgrn.m.eccentricity","Lesion.redblu.m.eccentricity", "Lesion.grnblu.m.eccentricity", "Lesion.Bred.m.eccentricity", "Lesion.Bgrn.m.eccentricity", "Lesion.Bblu.m.eccentricity", "Lesion.Bredgrn.m.eccentricity", "Lesion.Bredblu.m.eccentricity", "Lesion.Bgrnblu.m.eccentricity" )]

write.csv(allphenos, "BcAt_LesionEccentricity.csv")
View(allphenos)

###### Finding Linear Model 

library(lme4)
library(lmer)
#Experiment -> block of experiment in which it was done 

fullmod <- aov(Lesion.0.m.eccentricity ~ Isolate + Plant + Experiment + (Isolate*GrowingFlat) + (Isolate*AgarFlat), allphenos) # all terms included with basic interaction 
fullmod


##Current Most Accurate model with dropped terms
Cmodel <- lm(Lesion.0.m.eccentricity ~ Isolate + Plant + Experiment + (Isolate*Plant) + (Plant*GrowingFlat) + (Isolate*AgarFlat), allphenos) ## Dropped (Plant*AgarFlat) term b/c change was insignificant 

currmod <- Cmodel # Change for whatever model
summary(currmod) # provide summary stat4 over the parameters of the lm
anova(currmod)# get the anova table for the linear model
shapiro.test(residuals(currmod))


###### Finding LsMeans
f=NULL
library(data.table)
library(lsmeans)
library(lme4)
library(lmerTest)
BcAtPhenos.lm <- lmer(Lesion.0.m.eccentricity ~ Isolate + Plant + Experiment + Isolate*Plant + (1|GrowingFlat) + (1|AgarFlat), allphenos)
# BcAtPhenos.lm <- lm(Lesion.0.m.eccentricity ~ Isolate + Plant + Experiment + (Isolate*Plant) + (Plant*GrowingFlat) + (Isolate*AgarFlat), allphenos)

anova(BcAtPhenos.lm) #check to make sure model is working
#Wavy.lsm <- lsmeans(Wavy.lm, "Isolate")
BcAtPhenos.lsm <- lsmeansLT(BcAtPhenos.lm) #lsmeans is deprecaed -> lsmeansLT works (recommended by R)
df <- as.data.frame(print(BcAtPhenos.lsm))
setDT(df, keep.rownames = T)[]

write.csv(df, "C:/Users/vegaj/Desktop/B. cinera/Hwaviness/data/BcAtPhenos/BcAtPhenosGWAS_lsmeans.fxmod1.csv")
