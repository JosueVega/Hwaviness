#Josue Vega
#02242017
#10_GeneAnnotations.R

##update input file names
#Goal: summarize gene annotation files and draw venn diagrams
#Input Files: TopSNPs_Long_annot.csv

#and draw SNP to gene histogram
#------------------------------------------------------
rm(list=ls())
## install these if needed
library(dplyr); library(ggplot2)
## set correct directory
setwd("~/../Desktop/B. cinera/Hwaviness/")
#read in genes file: see 07_SNPdat_annot.R for which
##AllGenes <- read.csv("data/05_genes/TopSNPs_Long_trueMAF20_10NA.csv")
#this is just the list of genes, from SNPdat. SNP info etc. must be reattached.
##myAnnot <- read.csv("data/06_snpdat/name.FORPERL.output.csv")
#names(myAnnot)

#SNP list
## read in correct file
## can run this separately for 99% threshold and 99.9% threshold

mySNPs <- read.csv("data/05_genes/Hwavi_TopSNPs99.csv")
myGenes <- read.csv("data/06_snpdat/snpdat99/Hwavi_TopSNPs99.FORPERL.csv")

# mySNPs <- read.csv("data/05_genes/Hwavi_TopSNPs999.csv")
# myGenes <- read.csv("data/06_snpdat/snpdat999/Hwavi_TopSNPs999.FORPERL.csv")

#merge SNPs with genes
mySNPs$Chrom2 <- paste("CHROMOSOME",mySNPs$Chrom, sep='')
mySNPs$Chrom.Pos <- paste(mySNPs$Chrom2, mySNPs$Pos, sep='.')
myGenes$Chrom.Pos <- paste(myGenes$Chromosome.Number, myGenes$SNP.Position, sep='.')
#now, need to narrow down myGenes list to only include genes WITHIN 2kb WINDOW of SNP
#window options we can try:
#1kb (500 bp each side), 2kb (1kb each side), 4kb (2kb each side)
myGenes$Distance.to.nearest.feature <- as.numeric(myGenes$Distance.to.nearest.feature)
myGenes$Distance.to.nearest.feature[is.na(myGenes$Distance.to.nearest.feature)] <-0
myGenes2kb <- myGenes[myGenes$Distance.to.nearest.feature<1000,]

#sub in IPgen1, IPgen2, IPgen4
mygen <- myGenes2kb[,c("Chrom.Pos","gene.ID.containing.the.current.feature")]
#this includes multiple SNPs per gene
mygen <- unique(mygen)
#remove X and Chrom2
## check which columns to remove
mySNPs <- mySNPs[,-c(1,10)]
mySNPwGen <- merge(mySNPs, mygen, by="Chrom.Pos")

#then count number of phenotypes 
names(mySNPwGen)
colnames(mySNPwGen)[10] <- "geneID"

write.csv(mySNPwGen, "data/06_snpdat/99Thr_snpGenes_2kbwin.csv")
# write.csv(mySNPwGen, "data/06_snpdat/999Thr_snpGenes_2kbwin.csv")

#if any SNP is > 0 for a given gene, will get sum > 0
#else, sum = 0
myGeneSummary <- mySNPwGen %>%
  group_by(geneID) %>%
  summarize(pheno_count = sum((abs(Estimate))>0, na.rm = TRUE))

table(myGeneSummary$pheno_count)

##name according to threshold
write.csv(myGeneSummary, "data/06_snpdat/99Thr_genes_2kbwin.csv")
# write.csv(myGeneSummary, "data/06_snpdat/999Thr_genes_2kbwin.csv")
