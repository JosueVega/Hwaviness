#Josue Vega
#plot of SNPs along gene of interest

#12_singleGeneManhattan.R
#---------------------------------------------
rm(list=ls())
setwd("~/../Desktop/B. cinera/Hwaviness/")
#setwd("~/Projects/BcSolGWAS/")
library(plyr); library(ggplot2); library(grid)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("data/04_bigRRoutput/trueMAF20_20NA/HWavi_trueMAF20_20NA.HEM.PlotFormat.csv")

#get threshhold values 
HEM.thresh <- read.csv("data/04_bigRRoutput/trueMAF20_20NA/HWavi_trueMAF20_20NA.HEM.Thresh.csv")

#take the SNPs over the threshold for each phenotype

TH99pos <- HEM.thresh[3,]
for (i in 2:ncol(TH99pos)){
  assign(paste("TH99pos_", names(TH99pos[i]), sep=""),as.numeric(TH99pos[i]))
}
TH99neg <- HEM.thresh[7,]
for (i in 2:ncol(TH99neg)){
  assign(paste("TH99neg_", names(TH99neg[i]), sep=""),as.numeric(TH99neg[i]))
}


names(HEM.plotdata)
HEM.plotdata <- HEM.plotdata[,-c(1)]
#only look at chromosome 16
#this includes all contigs (chrom and seg are separated)
HEM.plotdata <- HEM.plotdata[which(HEM.plotdata$Chrom.Seg=='2.2'),]

#get the start position of chromosome 16
min(HEM.plotdata$Index)
max(HEM.plotdata$Index)
max(HEM.plotdata$Index) - min(HEM.plotdata$Index)

#narrow window: +- 1 kb
HEM.plotdata$Chr2.2Index <- HEM.plotdata$Index - min(HEM.plotdata$Index) + 1
min(HEM.plotdata$Chr2.2Index)
#now get target region within chromosome 2.2
#my gene/feature: about 1kb
#825306 to 826345

#and I'll add 2kb on each side
HEM.plotdataSM <- HEM.plotdata[which(HEM.plotdata$Pos > 823306),]
HEM.plotdataSM <- HEM.plotdataSM[which(HEM.plotdataSM$Pos < 828345),]

#trying a bigger window (8kb) to find missing phenos
#HEM.plotdataSM <- HEM.plotdata[which(HEM.plotdata$Pos < 355042),]
#HEM.plotdataSM <- HEM.plotdataSM[which(HEM.plotdataSM$Pos > 337285),]

#and remove any duplicated POS here
SNPlist_8a <- HEM.plotdataSM[!duplicated(HEM.plotdataSM$Pos), ]

#SNPlist <- as.data.frame(SNPlist_8a$Pos)
#write.csv(SNPlist,"data/genome/chr16_analysis/SNPlistFig8a.csv")

HEM.plotdata.OG <- HEM.plotdata
HEM.plotdata <- HEM.plotdataSM

#add a new variable so we can plot in order of seg.pos
HEM.plotdata$Seg.Pos <- paste(HEM.plotdata$Segment, HEM.plotdata$Pos, sep='.')
HEM.plotdata <- HEM.plotdata[order(HEM.plotdata$Segment, HEM.plotdata$Pos),]

#All groups (4:6)
#keep only: SNPs over 99% Threshold
#now very few over 99.9% Thr
for (i in c(4:15)){
  assign(paste("HEMpos.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] > get(paste("TH99pos_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,i)))
  assign(paste("HEMneg.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] < get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,i)))
}

#combine pos and neg by group
for (i in c(4:15)){
  assign(paste("HEM.", names(HEM.plotdata[i]), sep=""), rbind(get(paste("HEMpos.", names(HEM.plotdata[i]), sep="")),get(paste("HEMneg.", names(HEM.plotdata[i]), sep=""))))
}

#also keep non-sig SNPs as grey dots
for (i in c(4:15)){
  assign(paste("ns.HEMpos.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] < get(paste("TH99pos_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,i)))
  assign(paste("ns.HEMneg.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] > get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,i)))
}
for (i in c(4:15)){
  assign(paste("ns.HEM.", names(HEM.plotdata[i]), sep=""), rbind(get(paste("ns.HEMpos.", names(HEM.plotdata[i]), sep="")),get(paste("ns.HEMneg.", names(HEM.plotdata[i]), sep=""))))
}
for (i in c(4:15)){
  mydf <- paste("ns.HEM.", names(HEM.plotdata[i]), sep="")
  renamedf <- get(mydf)
  colnames(renamedf)[4] <- "Effect"
  assign(mydf, renamedf)
  myblob <- rep("NS", nrow(get(mydf)))
  assign(mydf, cbind(get(mydf), Trait = myblob))
}
HEM.NS <- rbind(ns.HEM.LA410, ns.HEM.LA480, ns.HEM.LA1547, ns.HEM.LA1589, ns.HEM.LA1684, ns.HEM.LA2093, ns.HEM.LA2176, ns.HEM.LA2706, ns.HEM.LA3008, ns.HEM.LA3475, ns.HEM.LA4345, ns.HEM.LA4355)

#then combine
#4:15
for (i in c(4:15)){
  mydf <- paste("HEM.", names(HEM.plotdata[i]), sep="")
  renamedf <- get(mydf)
  colnames(renamedf)[4] <- "Effect"
  assign(mydf, renamedf)
  myblob <- rep(names(HEM.plotdata[i]), nrow(get(mydf)))
  assign(mydf, cbind(get(mydf), Trait = myblob))
}
HEM.topSNPsB <- rbind(HEM.LA410, HEM.LA480, HEM.LA1547, HEM.LA1589, HEM.LA1684, HEM.LA2093, HEM.LA2176, HEM.LA2706, HEM.LA3008, HEM.LA3475, HEM.LA4345, HEM.LA4355)


HEM.topSNPs <- rbind(HEM.NS, HEM.topSNPsB)

#create a custom color scale
#myColors <- c("grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60")
myColors <- c("gray85", "#999999", "#292929","#684800" ,"#CBA22A", "#63B2D3", "#1FA69D", "#57B761", "#DAD94C","#2B869D", "#EE82EE", "#D2652D", "#CC79A7")
names(myColors) <- levels(HEM.topSNPs$Trait)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

#modify colors so that wild are oranges and domesticated are blues
levels(HEM.topSNPs$Trait)
#"NS"     "LA410"  "LA480"  "LA1547" "LA1589" "LA1684" "LA2093"
#"LA2176" "LA2706" "LA3008" "LA3475" "LA4345" "LA4355"
#NS, D, W, W, W, W, W, W, D, D, D, D, D
#oranges: coral1, deeppink3, yellow1, darkorange, red4, pink1
#blues: darkblue, dodgerblue1, blueviolet, lawngreen, seagreen4, mediumorchid1
myColors <- c("gray85", "darkblue", "coral1", "deeppink3", "goldenrod2", "darkorange", "red4", "pink1", "dodgerblue1", "blueviolet", "lawngreen", "seagreen4", "mediumorchid1")

names(myColors) <- levels(HEM.topSNPs$Trait)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

HEM.topSNPs.P <- HEM.topSNPs

#add lines for plot
myblocks <- read.csv("data/genome/chr16_analysis/BlockBoundaries.csv")
block2snp <- read.csv("data/genome/chr16_analysis/plink/fig8aMatch/MatchDrawLines.csv")
names(myblocks)[1] <- "SNPnum"
mylines <- merge(myblocks, block2snp, by="SNPnum")
#add indexing here!
mylines <- mylines[,c(1,3)]
mylines <- merge(mylines, getmyindex, by="Pos")

names(HEM.topSNPs.P)

SNPlist <- HEM.plotdataSM$Pos
write.csv(SNPlist, "data/genome/chr2_analysis/SNPlistFig8a.csv")

#how many genotypes per SNP?
HEM.SNP.list <- HEM.topSNPs[,c("Pos","Trait")]
HEM.SNP.w <-  as.data.frame(table(HEM.SNP.list))

#get SNPs for breaks to match with 8b haplotype plot
unique(HEM.topSNPs.P$Pos)

#SNP blocks from SNPlist32Fig8b.csv

#jpeg("paper/plots/FigR8/Sl_LesionSize_trueMAF20_NA10_lowTR.gene01Chr16.ManhattanPlot.jpg", width=7, height=5, units='in', res=600)
jpeg("paper/plots/FigR8/Sl_LesionSize_trueMAF20_NA10_lowTR.gene01Chr2.2.ManhattanPlot.jpg", width=7, height=5, units='in', res=600)
ggplot(HEM.topSNPs.P, aes(x=Pos, y=100*Effect))+
  theme_bw()+
  colScale+
  geom_point(aes(color = factor(Trait)))+
  labs(list(y=expression(paste("Estimated SNP Effect Size (",mm^{2},")"))))+
  guides(col = guide_legend(nrow = 8, title="SNP position"))+
  theme(legend.position="none")+
  scale_y_continuous(breaks=c(0.015, 0.01, 0.005, 0, -0.005, -0.01, -0.015), limits=c(-0.015, 0.015))+
  scale_x_continuous(name="SNP position on Chr 2.2 (kb)", 
                     #this is by Pos
                     breaks=c(823000, 824000, 825000, 826000, 827000, 828000, 829000), 
                     
                     labels=c("823", "824", "825", "826", "827", "828", "829"), limits=c(823000, 829000))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  #floated above chart: ymin=0.009, ymax=0.011
  #exon 1 825306	826178
  geom_rect(mapping=aes(ymin=-0.001, ymax=0.001, xmin=825306, xmax=826178), alpha=0.01, fill="darkturquoise")+
  #exon 2 826235	826345
  geom_rect(mapping=aes(ymin=-0.001, ymax=0.001, xmin=826235, xmax=826345), alpha=0.01, fill="darkturquoise")+
  #block 1
  geom_rect(mapping=aes(ymin=0.010, ymax=0.012, xmin=823323, xmax=823506), alpha=0.01, fill="red")+
  #block 2
  geom_rect(mapping=aes(ymin=0.010, ymax=0.012, xmin=823848, xmax=824176), alpha=0.01, fill="red")+
  #block 3
  geom_rect(mapping=aes(ymin=0.010, ymax=0.012, xmin=824908, xmax=827148), alpha=0.01, fill="red")+
  #block 4
  geom_rect(mapping=aes(ymin=0.010, ymax=0.012, xmin=827641, xmax=828305), alpha=0.01, fill="red")+
  #geom_vline(xintercept=c(x),linetype="dotted")+
  expand_limits(y=0)
dev.off()