#01b_EccentricityPhenos.R
#NES
#01/12/18

#prepare lesion eccentricity phenotypes of B. cinerea on A. thaliana for correlation to B. cinerea hyphal waviness
#-------------------------------------------------------------------------
rm(list=ls())
setwd("~/../Desktop/B. cinera/Hwaviness/data/BcAtPhenos/")

#read in the files
isolist <- read.csv("IsolateKey.csv")
myphenos <- read.csv("BcAt_LesionImageData.csv")
#rename Sample to match isolist$Sample
names(myphenos)[6] <- "Sample"

#add in correct isolate names
allphenos <- merge(myphenos,isolist, by="Sample")

#keep only the columns we need: independent variables and all eccentricity phenotypes
allphenos <- allphenos[,c("Experiment","GrowingFlat","AgarFlat","Plant","Isolate","Lesion.0.m.eccentricity","Lesion.red.m.eccentricity","Lesion.grn.m.eccentricity", "Lesion.blu.m.eccentricity", "Lesion.redgrn.m.eccentricity","Lesion.redblu.m.eccentricity", "Lesion.grnblu.m.eccentricity", "Lesion.Bred.m.eccentricity", "Lesion.Bgrn.m.eccentricity", "Lesion.Bblu.m.eccentricity", "Lesion.Bredgrn.m.eccentricity", "Lesion.Bredblu.m.eccentricity", "Lesion.Bgrnblu.m.eccentricity" )]

#Here's the independent variable info for designing a linear model:
# Experiment = 1 or 2, which of 2 replicate experimental BLOCKS the data are from
# GrowingFlat = which BLOCK the plants came from (this is pretty similar to PlateBlock in your Hwaviness data)
# AgarFlat = which BLOCK the plant * isolate interaction was measured in (this is most similar to PlateBlock in Hwaviness)
# Plant = which plant GENOTYPE was used
# Isolate = which Botrytis ISOLATE was used (this is the same as Isolate in Hwaviness)

#and to start, only use Lesion.0.m.eccentricity as your dependent variable. Ignore the rest.

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
summary(currmod) # provide summary stat over the parameters of the lm
anova(currmod)# get the anova table for the linear model
shapiro.test(residuals(currmod))

