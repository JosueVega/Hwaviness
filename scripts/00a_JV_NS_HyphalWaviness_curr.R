JV_HyphalWaviness_NS.R
#Created on R by Josue Vega and Nicole Soltis Jul 2017 
#Current Purpose
#Individual Scripts for Data Analysis of Hyphal Waviness

library(readr)

mydata <- read_csv("C:/Users/vegaj/Desktop/B. cinera/BR_JV_ManualHyphalDat_032817(edit).csv") #finds and opens csv file
#mydata <- read.csv("~/Projects/WavyHyphae/data/BR_JV_ManualHyphalDat_032817.csv")
#recodedate <- read.csv("C:/Users/vegaj/Desktop/scripts_forJV/")

library(dplyr)
library(ggplot2)
library(beanplot)
library(RColorBrewer) 
library(lme4)
library(tidyverse)
#check library open first

#Phenotype (1-10) ->
#10 = shortest wavelength
#1 = no waviness (longest wavelength)
#PlateBlock - 5 Blocks per plate
#DotRep - Numbering order of PlateBlock spore
#Rep - PlateBlock + DotRep (Concatinated)


ls(mydata) #lists variables

#mydata$Date <- paste(mydata$Date, "-16", sep="") 
#mydata <- merge(mydata, recodedate, by="Date")
#unique(mydata$NewDate)
mydata$NewDate <- as.Date(mydata$NewDate) # Orders dates properly

mydata$Rep <- paste(mydata$PlateBlock, mydata$DotRep) #finishes Rep by concatinating PlateBlock+DotRep

#------------------------------------ 
#Summary Data Tables -> Graphs


#Creates table of averages under PlateBlock (mean, min, max, sd, n)
plateBlock<- mydata %>% 
  group_by(PlateBlock) %>%
  summarise(avg_pheno = mean(Phenotype, na.rm = TRUE), 
            min_pheno = min(Phenotype, na.rm = TRUE), 
            max_pheno = max(Phenotype, na.rm = TRUE),
            sd_pheno = sd(Phenotype, na.rm = TRUE),
            total = n())

#Creates table of averages under Isolate (mean, min, max, n)
isolate <- mydata %>% 
  group_by(Isolate) %>%
  summarise(avg_pheno = mean(Phenotype, na.rm = TRUE), 
            min_pheno = min(Phenotype, na.rm = TRUE), 
            max_pheno = max(Phenotype, na.rm = TRUE),
            sd_pheno = sd(Phenotype, na.rm = TRUE),
            total = n())

#Creates table of averages under Date (mean, min, max, n)
Date <- mydata %>% 
  group_by(NewDate) %>%
  summarise(avg_pheno = mean(Phenotype, na.rm = TRUE), 
            min_pheno = min(Phenotype, na.rm = TRUE), 
            max_pheno = max(Phenotype, na.rm = TRUE),
            sd_pheno = sd(Phenotype, na.rm = TRUE),
            total = n())


#--------------------------------- 
#Trend Graphs for Visualization 

# GI added this- it looks cleaner than the default grey background
p <- ggplot(date, aes(x=(Date), y=avg_pheno, fill=Date)) + 
  geom_bar(colour="black", stat="identity") + 
  xlab("Date") + ylab("Average Phenotype (1 - 10)") +
  ggtitle("Average Hyphal Phenotype based on Date") + 
  geom_errorbar(aes(ymin = avg_pheno-sd_pheno, ymax = avg_pheno+sd_pheno), width=.2, position=position_dodge(.9)) +        guides(fill=FALSE) + #for removal of the legend
  theme(axis.text.x=element_text(angle = -90, hjust = 0))
print(p)

#Graph of Isolate p
p <- ggplot(isolate, aes(x=(Isolate), y=avg_pheno, fill=Isolate)) + 
  geom_bar(colour="black", stat="identity") + 
  xlab("Isolate") + ylab("Average Phenotype (1 - 10)") +
  ggtitle("Average Hyphal Phenotype based on Isolate") + 
  geom_errorbar(aes(ymin = avg_pheno-sd_pheno, ymax = avg_pheno+sd_pheno), width=.2, position=position_dodge(.9)) +        guides(fill=FALSE) + #for removal of the legend
  theme(axis.text.x=element_text(angle = -90, hjust = 0))
print(p)

###look at the function transform {base} and "reorder" within that to work on sorting by isolate and by date
  #Unsure what you meant here, marked my reorder comment that I thought you meant


# Graph of Phenotype based on Isolate
p <- ggplot(isolate, aes(x=(Isolate), y=avg_pheno, fill=Isolate)) + 
  geom_bar(colour="black", stat="identity") +
  xlab("Isolate") + ylab("Average Phenotype (1 - 10)") +
  #reorder(isolate, avg_pheno) + #Unsure what you meant
  ggtitle("Average Hyphal Phenotype based on Isolate") + 
  geom_errorbar(aes(ymin = avg_pheno-sd_pheno, ymax = avg_pheno+sd_pheno), width=.2, position=position_dodge(.9)) +        guides(fill=FALSE) + #for removal of the legend
  theme(axis.text.x=element_text(angle = -90, hjust = 0))
print(p)

#Graph of average Phenptype based on PlateBlock
p <- ggplot(plateBlock, aes(x=reorder(PlateBlock, avg_pheno), y=avg_pheno, fill=PlateBlock)) + 
  geom_bar(colour="black", stat="identity") +
  xlab("Isolate") + ylab("Average Phenotype (1 - 10)") +
  ggtitle("Average Hyphal Phenotype based on PlateBlock") + 
  geom_errorbar(aes(ymin = avg_pheno-sd_pheno, ymax = avg_pheno+sd_pheno), width=.2, position=position_dodge(.9)) +        guides(fill=FALSE) + #for removal of the legend
  theme(axis.text.x=element_text(angle = -90, hjust = 0))
print(p)

#----------------------------------------------- 
#Assumption Checks

#Reminder
#IndVar: Date, Isolate, Ordering, PDAConc, Rep
#DepVar: Phenotype

#Checking for normality
mydataCheck <- mydata
View(mydataCheck)

xtabs(~ Isolate + Date, mydataCheck) # Creates a contingency Table


#Checking normality of dependent variable
attach(mydataCheck) #set path
hist(Phenotype) #gives histogram of distribution of Phenotype

mydataCheck$Phenotype.t <- mydataCheck$Phenotype + 1 #? check for alterations if shifted?

shapiro.test(Phenotype) #if p-value < .05 --> Normally distributed


# I am missing the package? Or I do not have the function
qqp(mydataCheck$Phenotype.t, "norm") 
qqp(mydataCheck$Phenotype.t, "norm") 

# Should I continue with transformations to check for normalization?
  #Such as square and log transformations

#Double Checking for normality using log transformation
mydatatransf <- (log((Phenotype)))
hist(mydatatransf) #still normal?
shapiro.test(mydatatransf) #still p<.05


#check for homoscedasticity
boxplot(Phenotype~Isolate + Ordering, ylab="YTITLE", main="PLOTTITLE", las=3)
#Doable? The size is too big and stops R (not enough memory)


#statistic
bartlett.test(Phenotype~Isolate, mydataCheck) #only works with 1 indepVar?
leveneTest(Phenotype~Isolate) # Not a function I have (am I missing a library?)
var.test(Phenotype~Isolate) #Missing 2nd level?


#----------------------------------------------- 
#ANOVA + Power Tests

#Calculating One Way ANOVA (dep=Pheno, indep=Isolate
hyphalAOV = aov(Phenotype ~ Isolate, mydata)
summary(hyphalAOV)
TukeyHSD(hyphalAOV)
#Should this be represented in a table?


#Power analysis: Sufficient Data Collection?
library(pwr)
#one-way anova
pwr.anova.test(k=5,f=0.25,sig.level=0.05,power=0.8)
#Used Nicole's Notes for Template 
#two-way anova
#cohen's f2, related to R-squared (model fit)
#u is df of numerator
#= ((levels of factor 1)-1) * ((levels of factor 2)-1)
#v is total number of subjects across all categories
#MINUS (levels of factor 1)*(levels of factor 2)
#these values are given in anova table
pwr.f2.test(u=2, v=294, sig.level=0.05, power=0.8)
#f2 is Cohen's f-squared value, ~R-squared, so effect size low.
#by convention, f2 0.02 is small, 0.15 is medium, 0.35 is large

#----------------------------------------------- 
#Linear Model Explanation 

lm(Phenotype ~ Isolate + Date, mydata) # Isolate and Date as seperate factors
lm(Phenotype ~ Isolate * Date, mydata) #date and isolate interacting
summary(aov(Phenotype ~ Isolate * Date, mydata))


#----------------------------------------------- 
# Linear 2-Way ANOVA -> Parametric Test 
mydataANOVA <- anova(lm(Phenotype ~ Isolate + Isolate*PlateBlock + Date, mydata)) # 2Way ANOVA of a linear model
summary(mydataANOVA) #tabular form of the data, normal distr of information 
mydataANOVA # Shows Response (Phenotype) and ANOVA Table

# Mixed Model 2-Way ANOVA -> Parametric Test
mydataANOVA2 <- anova(lmer(Phenotype ~ Isolate + (1|PlateBlock) + Date, mydata))
summary(mydataANOVA2) #tabular form of the data, normal distr of information 
mydataANOVA2 # Shows Response (Phenotype) and ANOVA Table


# Analysis of Variance --> Graphical 
mydata.AOV <- aov(Phenotype~ Isolate + (1|PlateBlock) + Date, mydata)#invalid list type for Date (is it because Date is a different type of information)
plot(mydata.AOV) # Ordering produces QQ plot + Residual Plot
summary(mydata.AOV)
 
# Residual Normally Distributed 
MYResid <-lm(mydataCheck~Isolate*Date)
#run shapiro-wilk goodness of fit test on the residuals
shapiro.test(residuals(MYResid))


#----------------------------------------------------
##Here's where we develop the model to extract least-squared means of the data and feed into bigRR

library(lme4)

# DotRep Excluded
fullmod <- aov(Phenotype ~ Isolate + PlateBlock + NewDate + PDAConc, mydata) # all
fullmod

mymod1 <- lm(Phenotype ~ Isolate + Date + PlateBlock, mydata) 
mymod2 <- lm(Phenotype ~ Isolate + Date + Isolate*PlateBlock, mydata) 
mymod3 <- lm(Phenotype ~ Isolate + Date + Rep + PlateBlock, mydata)
mymod4 <- lm(Phenotype ~ Isolate + Date, mydata)
mymod5 <- lmer(Phenotype ~ Isolate + (1|Date) + Isolate*PlateBlock, mydata) 
mymod6 <- lmer(Phenotype ~ Isolate + Date + (1|PlateBlock), mydata)
mymod7 <- lm(Phenotype ~ Isolate + Date + Isolate*PlateBlock , mydata)


mymodnull <- lm(Phenotype ~ Date + PlateBlock + Rep, mydata)

# # # # # Test Script
currmod <- mymod6 # Change for whatever model
summary(currmod) # provide summary stat over the parameters of the lm
anova(currmod)# get the anova table for the linear model
 

shapiro.test(residuals(currmod))

anova(mymod2, mymod5)
rand(currmod)
# extract data for the random effect


#----------------------------------------------------------------
##here's how to run lsmeans and make a nice data table. This only works for a mixed model (fixed effects AND random effects)

##include random effects as (1|term)
##include fixed effects as term

df=NULL
library(data.table)
library(lsmeans)
library(lme4)
library(lmerTest)
Wavy.lm <- lmer(Phenotype ~ Isolate + (Isolate/PlateBlock) + (1|Date) + (1|PlateBlock) + (Date/PlateBlock), data=mydata)
#Wavy.lm <- lmer(Phenotype ~ Isolate + (Isolate*PlateBlock) + Date + (Date/PlateBlock), data=mydata)#date made random for lsmeans to work, NewDate is treated as a random effect
 
anova(Wavy.lm) #check to make sure model is working
#Wavy.lsm <- lsmeans(Wavy.lm, "Isolate")
Wavy.lsm <- lsmeansLT(Wavy.lm, "Isolate") #lsmeans is deprecaed -> lsmeansLT works (recommended by R)
df <- as.data.frame(print(Wavy.lsm))
setDT(df, keep.rownames = T)[]

write.csv(df, "C:/Users/vegaj/Desktop/B. cinera/HWavinessVisuals/WavyGWAS_2ndTry_lsmeans.fxmod1.csv")

#------------------------------------------------------------------
# bigRR
