#Josue Vega
#format SNP peaks for SNPdat annotation
#100217
#-----------------------------------------------------------
## run this for both TopSNPs999 and TopSNPs99
rm(list=ls())
##correct the working directory
#setwd("~/path/to/dir")

#load files
setwd("~/../Desktop/B. cinera/Hwaviness/")
myTopSNPs <- read.csv("data/05_genes/Hwavi_TopSNPs999.csv")
#myTopSNPs <- read.csv("data/05_genes/Hwavi_TopSNPs99.csv")

#format files for SNPdat 
names(myTopSNPs)
myTopSNPs$chromosome.id <- paste("Chromosome",myTopSNPs$Chrom.Cont, sep="")
myTopSNPs$position <- myTopSNPs$Pos
myTopSNPs$mutation <- "A"
myTopSNPs.snpdat <- myTopSNPs[,c("chromosome.id", "position", "mutation")]
#get rid of ".0"s
myTopSNPs.snpdat$chromosome.id <- gsub("\\.0$", "", myTopSNPs.snpdat$chromosome.id)
unique(myTopSNPs.snpdat$chromosome.id)
write.table(myTopSNPs.snpdat, file="data/05_genes/snpdat999/Hwavi_TopSNPs999.FORPERL.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
#write.table(myTopSNPs.snpdat, file="data/05_genes/snpdat99/Hwavi_TopSNPs99.FORPERL.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
