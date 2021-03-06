#reformat bigRR output data
#Josue Vega

#--------------------------------------------------------
rm(list=ls())
library(tidyr)
#setwd("~/Documents/GitRepos/BcSolGWAS/data/SNP_files")

setwd("~/../Desktop/B. cinera/Hwaviness/data/04_bigRRoutput/trueMAF20_20NA/")
#Import data
#reorganize file Sl_LesionSize.HEM.csv
HEMdat <- read.csv("HWavi_trueMAF20_20NA.HEM.csv")
View(HEMdat)
#first remove first 4 rows (threshold data)
## check that you are removing just the threshold data
HEMdat <- HEMdat[,-c(1)]
HEMthresh <- HEMdat[1:8,]
HEMdat <- HEMdat[-c(1:8),]
HEMdat2 <- HEMdat

#split chromosome and Contig into 2 variables
names(HEMdat)
unique(HEMdat$X.1)

HEMdat$X.1 <- gsub(pattern = "Chromosome1\\.", replacement = "Chromosome1.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome2\\.", replacement = "Chromosome2.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome3\\.", replacement = "Chromosome3.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome4\\.", replacement = "Chromosome4.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome5\\.", replacement = "Chromosome5.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome6\\.", replacement = "Chromosome6.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome7\\.", replacement = "Chromosome7.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome8\\.", replacement = "Chromosome8.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome9\\.", replacement = "Chromosome9.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome10\\.", replacement = "Chromosome10.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome11\\.", replacement = "Chromosome11.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome12\\.", replacement = "Chromosome12.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome13\\.", replacement = "Chromosome13.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome14\\.", replacement = "Chromosome14.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome15\\.", replacement = "Chromosome15.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.", replacement = "Chromosome16.0.", HEMdat$X.1)
unique(HEMdat$X.1)

#Then everything with 1.0.(1.)number I'll change to 1.(1.)number
HEMdat$X.1 <- gsub(pattern = "Chromosome1\\.0\\.1\\.", replacement = "Chromosome1.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome2\\.0\\.1\\.", replacement = "Chromosome2.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome2\\.0\\.2\\.", replacement = "Chromosome2.2.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome2\\.0\\.3\\.", replacement = "Chromosome2.3.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome2\\.0\\.4\\.", replacement = "Chromosome2.4.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome2\\.0\\.5\\.", replacement = "Chromosome2.5.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome3\\.0\\.1\\.", replacement = "Chromosome3.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome3\\.0\\.2\\.", replacement = "Chromosome3.2.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome4\\.0\\.1\\.", replacement = "Chromosome4.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome5\\.0\\.1\\.", replacement = "Chromosome5.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome6\\.0\\.1\\.", replacement = "Chromosome6.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome6\\.0\\.2\\.", replacement = "Chromosome6.2.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome6\\.0\\.3\\.", replacement = "Chromosome6.3.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome7\\.0\\.1\\.", replacement = "Chromosome7.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome7\\.0\\.2\\.", replacement = "Chromosome7.2.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome8\\.0\\.1\\.", replacement = "Chromosome8.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome8\\.0\\.2\\.", replacement = "Chromosome8.2.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome9\\.0\\.1\\.", replacement = "Chromosome9.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome10\\.0\\.1\\.", replacement = "Chromosome10.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome11\\.0\\.1\\.", replacement = "Chromosome11.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome12\\.0\\.1\\.", replacement = "Chromosome12.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome13\\.0\\.1\\.", replacement = "Chromosome13.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome13\\.0\\.2\\.", replacement = "Chromosome13.2.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome14\\.0\\.1\\.", replacement = "Chromosome14.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome14\\.0\\.2\\.", replacement = "Chromosome14.2.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome15\\.0\\.1\\.", replacement = "Chromosome15.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome15\\.0\\.2\\.", replacement = "Chromosome15.2.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome15\\.0\\.3\\.", replacement = "Chromosome15.3.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome15\\.0\\.4\\.", replacement = "Chromosome15.4.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.1\\.", replacement = "Chromosome16.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.2\\.", replacement = "Chromosome16.2.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.3\\.", replacement = "Chromosome16.3.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.4\\.", replacement = "Chromosome16.4.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.5\\.", replacement = "Chromosome16.5.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.6\\.", replacement = "Chromosome16.6.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.7\\.", replacement = "Chromosome16.7.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.8\\.", replacement = "Chromosome16.8.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.9\\.", replacement = "Chromosome16.9.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.10\\.", replacement = "Chromosome16.10.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.11\\.", replacement = "Chromosome16.11.", HEMdat$X.1)

HEMdat2 <- separate (HEMdat, X.1, into = c("Chrom", "Contig", "Pos") )
#double check
## make sure that chromosome and contig are now 2 separate terms
unique(HEMdat2$Chrom)
unique(HEMdat2$Contig)
HEMdat3 <- HEMdat2
HEMdat3$Chr.Cont <- paste(HEMdat2$Chrom, HEMdat2$Contig, sep=".")

HEM.plotdata <- HEMdat2

HEM.plotdata$Pos <- as.character(HEM.plotdata$Pos)#ensure that position data is not in scientific notation

#Reformat Chromosomes and Positions
HEM.plotdata$Chrom <- gsub("Chromosome", "", HEM.plotdata$Chrom)
HEM.plotdata$Chrom <- as.numeric(as.character(HEM.plotdata$Chrom))
HEM.plotdata$Contig <- as.numeric(as.character(HEM.plotdata$Contig))
HEM.plotdata$Pos <- as.numeric(as.character(HEM.plotdata$Pos))

#sort dataframe rows in order of Chrom, then Pos
HEM.plotdata <- HEM.plotdata[with(HEM.plotdata, order(Chrom, Contig, Pos)), ]

#now make Contigs line up consecutively
HEM.plotdata$Chrom.Cont <- paste(HEM.plotdata$Chrom, HEM.plotdata$Contig, sep=".")
HEM.plotdata$Chrom.Cont <- as.numeric(HEM.plotdata$Chrom.Cont)

#let's try making the chrom.Cont integers so that R isn't confused
unique(HEM.plotdata$Chrom.Cont)
HEM.plotdata$Chrom.Cont.F <- as.factor(HEM.plotdata$Chrom.Cont)
unique(HEM.plotdata$Chrom.Cont.F)
recode.vars <- data.frame(OGvals = c(1,1.1,2,2.1,2.2,2.3,2.4,2.5, 3, 3.1, 3.2, 4,  4.1, 5,  5.1, 6, 6.2, 6.3, 7,  7.1, 7.2, 8,  8.1, 8.2, 9,  9.1, 10.1, 11, 11.1, 12, 12.1, 13, 13.1, 14, 14.1, 14.2, 15, 15.1, 15.2, 15.3, 15.4, 16, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, 16.11), newvals = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52))

##make sure that these are paired numbers that make sense
print(recode.vars)

HEM.plotdata$Chrom.Cont.Int <- recode.vars$newvals[match(HEM.plotdata$Chrom.Cont.F, recode.vars$OGvals)]
## check, do these values make sense?
unique(HEM.plotdata$Chrom.Cont.Int)

#Make plotting variables -- a continuous count from Chromosome 1, Contig 1, Position 1 to the end of the last Contig of Chromosome 16.
HEM.plotdata$Index = NA
lastbase = 0

#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(HEM.plotdata$Chrom.Cont.Int)) {
  print(i)
  #for chromosome 1
  if (i==1) {
    #for the subset of HEM.plotdata rows with Chromosome 1, set Index variable for each row to equal Pos.
    HEM.plotdata[HEM.plotdata$Chrom.Cont.Int==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom.Cont.Int==i, ]$Pos
    #for all other chromosomes: 
  }	else {
    #lastbase for chromosome i is the greater of:
    #current lastbase counter plus the maxiumum position of chromosome i-1
    #OR 1
    #changed lastbase+tail to lastbase+max
    lastbase=lastbase+max(subset(HEM.plotdata,HEM.plotdata$Chrom.Cont.Int==i-1)$Pos, 1)
    #and then for the subset of HEM.plotdata rows with Chromosome i, set Index variable for each row to equal Pos + lastbase
    HEM.plotdata[HEM.plotdata$Chrom.Cont.Int==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom.Cont.Int==i, ]$Pos+lastbase
  }
}

## write this out to a file name that makes sense
write.csv(HEM.plotdata, "HWavi_trueMAF20_20NA.HEM.PlotFormat.csv") 
write.csv(HEMthresh, "Hwavi_trueMAF20_20NA.HEM.Thresh.csv")

#read in to 06_plots