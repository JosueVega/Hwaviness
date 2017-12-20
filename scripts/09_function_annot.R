#11_GeneAnnot_function
#Josue Vega
#052217

##customize this part
#Input: SNP data/GWAS_files/05_annotation/Domest_topSNPs_genestoAnnot.csv 
#and Gene data/GWAS_files/05_annotation/Domest_TopSNPs_10NA_intoAnt.csv
#from script 10_GeneAnnot_10NA_venns_figR8.R
#Output:
#Figures: None
#---------------------------------------------------------------------------------
rm(list=ls())

##set directory
setwd("~/../Desktop/B. cinera/Hwaviness")
library(dplyr)
myGenes <- read.csv("data/06_snpdat/99Thr_genes_2kbwin.csv")
allFuncs <- read.csv("data/07_SNPdat_annot/botrytis_cinerea__t4__1_pfam_to_genes_mycleaned.csv")

head(allFuncs)

#keep only the columns we need
myFuncs <- allFuncs[,c("GENE","PFAM_NAME","PFAM_DESCRIPTION")]
#myFuncs <- myFuncs[!duplicated(myFuncs[,"GENE"]),]

#now within each gene, merge back onto AnnotGenes
colnames(myFuncs)[1] <- "geneID"
myGenes <- myGenes[,-c(1)]
AnnotGenes<- merge(myGenes, myFuncs, by="geneID")
write.csv(AnnotGenes, "data/07_SNPdat_annot/AllAnnots_99Thr_byGene.csv")

#summary DF for WHOLE GENOME
#remove empty levels of the variable (in case function appears 0 times in whole genome)
myFuncs$PFAM_DESCRIPTION <- droplevels(myFuncs$PFAM_DESCRIPTION)
#how many times does a given function occur in the whole genome?
WGAntCats <- as.data.frame(table(myFuncs$PFAM_DESCRIPTION))
colnames(WGAntCats)[2] <-"WGGenFreq"
colnames(WGAntCats)[1] <- "Function"
WGAntCats <- WGAntCats[-c(1,2),]
#need to add droplevels() to correctly count the ocurrences of each Function, but this doesn't work for underrepresentation
AnnotGenes$PFAM_DESCRIPTION <- droplevels(AnnotGenes$PFAM_DESCRIPTION)

#-----------------------------------------------------------------
WavyAntCats <- as.data.frame(table(AnnotGenes$PFAM_DESCRIPTION))
colnames(WavyAntCats)[2] <-"AllWavyGenFreq"
colnames(WavyAntCats)[1]<- "Function"
AntOverrep <- merge(WGAntCats, WavyAntCats, by="Function", all=T)

#with all=T, can make underrepresentation work by filling in zeroes
AntOverrep[is.na(AntOverrep)] <- 0

#remove blank rows
AntOverrep <- AntOverrep[-c(1),]

#test for overrepresentation of a function
#https://stats.stackexchange.com/questions/72553/which-statistical-test-should-be-used-to-test-for-enrichment-of-gene-lists
#fisher's exact test
#http://www.biostathandbook.com/fishers.html
#function my list, function whole list, not-function my list, not-function whole list

#next: modify this to run for each phenotype

  fisher.p.over <- NULL
  fisher.p.under <- NULL
  for (i in 1:nrow(AntOverrep)){
    myrow <- AntOverrep[i,]
    myA <- myrow$AllWavyGenFreq
    myB <- myrow$WGGenFreq
    myCmA <- sum(AntOverrep[,3])-myA
    myDmB <- sum(AntOverrep$WGGenFreq)-myB
    fisher.p.up <- fisher.test(matrix(c(myA, myB, myCmA, myDmB), nrow=2, ncol=2), alternative="greater")$p.value
    fisher.p.down <- fisher.test(matrix(c(myA, myB, myCmA, myDmB), nrow=2, ncol=2), alternative="less")$p.value
    fisher.p.over <- append(fisher.p.over, fisher.p.up)
    fisher.p.under <- append(fisher.p.under, fisher.p.down)
  }
  assign(paste("fisher.p.over", names(AntOverrep)[3], sep="."),fisher.p.over)
  assign(paste("fisher.p.under", names(AntOverrep)[3], sep="."),fisher.p.under)

#write fisher values into df
names(AntOverrep)
AntOverrep$fisher.up.All <- fisher.p.over.AllWavyGenFreq
AntOverrep$fisher.dn.All <- fisher.p.under.AllWavyGenFreq

write.csv(AntOverrep, "data/07_SNPdat_annot/Thr99_FuncAnnotatedGenes.csv")
