#Josue Vega
# Separate visual projects for HWaviness

#################################
# Table 1: ANOVA table for hyphal waviness

rm(list=ls())
library(lme4)
library(readr)
setwd("~/../Desktop/B. cinera/Hwaviness")
mydata <- read.csv("../BR_JV_ManualHyphalDat_032817(edit).csv") 

# Mixed Model 2-Way ANOVA -> Parametric Test
mydataANOVA1 <- anova(lmer(Phenotype ~ Isolate + (1|PlateBlock) + Date, mydata))
mydataANOVA2 <- aov(Phenotype ~ Isolate + (Isolate*PlateBlock) + Date, mydata)
mydataANOVA2 <- summary(mydataANOVA2) #tabular form of the data, normal distr of information 
mydataANOVA2
write.csv(mydataANOVA2, "../../HWaviness Parts/Completed/ANOVA_Hwavi.csv") # Shows Response (Phenotype) and ANOVA Table



#################################
## Histogram of phenotype Counts

rm(list=ls())
library(ggplot2)
setwd("~/../Desktop/B. cinera/Hwaviness")

#999 Tresh
pheno <- read.csv("data/07_SNPdat_annot/AllAnnots_999Thr_byGene.csv")
hist(pheno$pheno_count, axes=TRUE, breaks = 5, plot=TRUE, labels=TRUE, xlab="Number of SNPs per Gene", main="99.9% Treshold: SNP Frequency Distribution")

#99 Tresh
pheno2 <- read.csv("data/07_SNPdat_annot/AllAnnots_99Thr_byGene.csv")
hist(pheno2$pheno_count, axes=TRUE, breaks = 5, plot=TRUE, labels=TRUE, xlab="Number of SNPs per Gene", main="99% Threshold: SNP Frequency Distribution")

#################################
#Merging GeneID and Function and Estimate

rm(list=ls())
library(ggplot2)
setwd("~/../Desktop/B. cinera/Hwaviness")

IDgene_99_Estimate <- read.csv("data/06_snpdat/99Thr_snpGenes_2kbwin.csv")
IDgene_99_Estimate <- IDgene_99_Estimate [,c("Estimate","geneID", "Chrom")]
IDgene_99_Estimate <- unique(IDgene_99_Estimate, by = "geneID")

# View(IDgene_99_Estimate)
annots_99_Func <- read.csv("data/07_SNPdat_annot/AllAnnots_99Thr_byGene.csv")
annots_99_Func <- annots_99_Func [,c("PFAM_NAME","PFAM_DESCRIPTION", "geneID")]
annots_99_Func <- unique(annots_99_Func, by = "geneID")
View(annots_99_Func)
# View(annots_99_Func)
overall <- merge(IDgene_99_Estimate,annots_99_Func,by="geneID")
View(overall)
# annots_99_Func <- annots_99_Func [,c("PFAM_NAME","PFAM_DESCRIPTION")]
# count(unique(annots_99_Func))
# IDgene_99_Estimate <- IDgene_99_Estimate [,c("Estimate")]
# count(unique(IDgene_99_Estimate))

write.csv(overall, "../HWavinessVisuals/HWaviGeneID_fxn_estimateThresh99.csv")

#################################

#  scatter plot of hyphal waviness vs. pectin growth on sugar agar

rm(list=ls())
library(ggplot2)
library(dplyr)

setwd("~/../Desktop/B. cinera/")

library(readr)
SugarPectin <- read_delim("C:/Users/vegaj/Desktop/B. cinera/HWavinessVisuals/SugarPectin_Lsmeans.csv", ";", escape_double = FALSE, trim_ws = TRUE)
SugarPectin <- SugarPectin [,c("GenoRename", "P72", "P48", "S72", "S48")]
colnames(SugarPectin)<-c("Isolate", "P72", "P48", "S72", "S48")
HwaviBC <- read.csv("BR_JV_ManualHyphalDat_032817(edit).csv")
Hwavi <- HwaviBC


Isolate <- Hwavi %>% 
  group_by(Isolate) %>%
  summarise(avg_pheno = mean(Phenotype, na.rm = TRUE), 
            min_pheno = min(Phenotype, na.rm = TRUE), 
            max_pheno = max(Phenotype, na.rm = TRUE),
            sd_pheno = sd(Phenotype, na.rm = TRUE),
            total = n())

SugarPectinWaviness <- merge(SugarPectin, Isolate, by = "Isolate")
#View(SugarPectinWaviness)

#Scatter Plot of Pectin/Sugar(72+48hrs) v Waviness
##P72

pvaluePear <- summary(lm(avg_pheno~P72, SugarPectinWaviness))$coefficients["P72","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(SugarPectinWaviness$avg_pheno,SugarPectinWaviness$P72, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(SugarPectinWaviness$avg_pheno,SugarPectinWaviness$P72, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/Pectin72hScatter_Waviness.pdf")
plot1 <- ggplot(SugarPectinWaviness, aes(P72,avg_pheno)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  xlab("Isolate Growth on Pectin after 72hr") + ylab("Average Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Comparison of Isolate growth on Pectin after 72hr to Isolate Hyphal Waviness") +
  geom_smooth(method='lm')
plot1
dev.off()

##P48

pvaluePear <- summary(lm(avg_pheno~P48, SugarPectinWaviness))$coefficients["P48","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(SugarPectinWaviness$avg_pheno,SugarPectinWaviness$P48, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(SugarPectinWaviness$avg_pheno,SugarPectinWaviness$P48, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/Pectin48hScatter_Waviness.pdf")
plot2 <- ggplot(SugarPectinWaviness, aes(P48,avg_pheno)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  xlab("Isolate Growth on Pectin after 48hr") + ylab("Average Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Comparison of Isolate growth on Pectin after 48hr to Isolate Hyphal Waviness") +
  geom_smooth(method='lm')
plot2
dev.off()
##S72

pvaluePear <- summary(lm(avg_pheno~S72, SugarPectinWaviness))$coefficients["S72","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(SugarPectinWaviness$avg_pheno,SugarPectinWaviness$S72, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(SugarPectinWaviness$avg_pheno,SugarPectinWaviness$S72, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/Sugar72hScatter_Waviness.pdf")
plot3 <- ggplot(SugarPectinWaviness, aes(S72,avg_pheno)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  xlab("Isolate Growth on Sugar after 72hr") + ylab("Average Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Comparison of Isolate growth on Sugar after 72hr to Isolate Hyphal Waviness") +
  geom_smooth(method='lm')
plot3
dev.off()
##S48

pvaluePear <- summary(lm(avg_pheno~S48, SugarPectinWaviness))$coefficients["S48","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(SugarPectinWaviness$avg_pheno,SugarPectinWaviness$S48, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(SugarPectinWaviness$avg_pheno,SugarPectinWaviness$S48, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/Sugar48hScatter_Waviness.pdf")
plot4 <- ggplot(SugarPectinWaviness, aes(S48,avg_pheno)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  xlab("Isolate Growth on Sugar after 48hr") + ylab("Average Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Comparison of Isolate growth on Sugar after 48hr to Isolate Hyphal Waviness") +
  geom_smooth(method='lm')
plot4
dev.off()

##Multiple Plots in one
library(cowplot)
gridPlot <- plot_grid( plot2, plot1,plot4, plot3, labels = c("A", "B", "C", "D"), scale = .75)
gridPlot

#################################
#################################
# scatter plot of hyphal waviness vs. lesion size on tomato, by isolate

rm(list=ls())
library(ggplot2)
library(dplyr)

setwd("~/../Desktop/B. cinera/")

BC_tomato <- read.csv("HWavinessVisuals/Bc_tomatoLesion.csv")
HwaviBC <- read.csv("BR_JV_ManualHyphalDat_032817(edit).csv")
Hwavi <- HwaviBC


Isolate <- Hwavi %>% 
  group_by(Isolate) %>%
  summarise(avg_pheno = mean(Phenotype, na.rm = TRUE), 
            min_pheno = min(Phenotype, na.rm = TRUE), 
            max_pheno = max(Phenotype, na.rm = TRUE),
            sd_pheno = sd(Phenotype, na.rm = TRUE),
            total = n())

BC_TomatoWaviness <- merge(BC_tomato, Isolate, by = "Isolate")

#Scatter Plot mean Lesion v Waviness

pvaluePear <- summary(lm(avg_pheno~meanLesion, BC_TomatoWaviness))$coefficients["meanLesion","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BC_TomatoWaviness$avg_pheno,BC_TomatoWaviness$meanLesion, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BC_TomatoWaviness$avg_pheno,BC_TomatoWaviness$meanLesion, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/TomatoScatter_Lesion_Waviness.pdf")
plotTom <- ggplot(BC_TomatoWaviness, aes(meanLesion,avg_pheno)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  xlab("Mean Lesion Size on Tomato") + ylab("Average Hyphal Waviness") +
  ggtitle(NULL) +
  #ggtitle("Mean Lesion Size on Tomato against Hyphal Waviness per Isolate") +
  geom_smooth(method='lm', se = FALSE)
plotTom
dev.off()

# 2nd scatter plot of hyphal waviness vs. lesion size on different interaction types, by isolate

library(ggplot2)
library(dplyr)

setwd("~/../Desktop/B. cinera/")

BC_Interact <- read.csv("HWavinessVisuals/DataFor100IsoCollection.csv")
BC_Interact <- BC_Interact [,c("Isolate","Cendivia","Brapa","Cintybus","Glycine",	"Helianthus",	"Solanum")]
HwaviBC <- read.csv("BR_JV_ManualHyphalDat_032817(edit).csv")
Hwavi <- HwaviBC


Isolate <- Hwavi %>% 
  group_by(Isolate) %>%
  summarise(avg_pheno = mean(Phenotype, na.rm = TRUE), 
            min_pheno = min(Phenotype, na.rm = TRUE), 
            max_pheno = max(Phenotype, na.rm = TRUE),
            sd_pheno = sd(Phenotype, na.rm = TRUE),
            total = n())

BC_InteractWaviness <- merge(BC_Interact, Isolate, by = "Isolate")

#Scatter Plot mean Interactions Lesion v Waviness
##C. endivia - dicot

pvaluePear <- summary(lm(avg_pheno~Cendivia, BC_InteractWaviness))$coefficients["Cendivia","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$Cendivia, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$Cendivia, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/CendiviaScatter_Lesion_Waviness.pdf")
plot1 <- ggplot(BC_InteractWaviness, aes(Cendivia,avg_pheno)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) + #geom_text(label = (rvalueSpea), parse = TRUE) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  xlab("Mean Lesion Size on C. endivia") + ylab("Average Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Mean Lesion Size on Cichorium endivia against Hyphal Waviness per Isolate") +
  geom_smooth(method='lm')
plot1
dev.off()

##B. rapa - dicot

pvaluePear <- summary(lm(avg_pheno~Brapa, BC_InteractWaviness))$coefficients["Brapa","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$Brapa, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$Brapa, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/BrapaScatter_Lesion_Waviness.pdf")
plot2 <- ggplot(BC_InteractWaviness, aes(Brapa,avg_pheno)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  geom_text(x = 2, y = 6, label = (rvalueSpea), parse = TRUE) +
  xlab("Mean Lesion Size on B. rapa") + ylab("Average Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Mean Lesion Size on Brassica rapa against Hyphal Waviness per Isolate") +
  geom_smooth(method='lm')
plot2
dev.off()
##C. intybus - dicot

pvaluePear <- summary(lm(avg_pheno~Cintybus, BC_InteractWaviness))$coefficients["Cintybus","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$Cintybus, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$Cintybus, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/CintybusScatter_Lesion_Waviness.pdf")
plot3 <- ggplot(BC_InteractWaviness, aes(Cintybus,avg_pheno)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  geom_text(x = 2, y = 6, label = (rvalueSpea), parse = TRUE) +
  xlab("Mean Lesion Size on C. intybus") + ylab("Average Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Mean Lesion Size on Cichorium intybus against Hyphal Waviness per Isolate") +
  geom_smooth(method='lm')
plot3
dev.off()
##Glycine max

pvaluePear <- summary(lm(avg_pheno~Glycine, BC_InteractWaviness))$coefficients["Glycine","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$Glycine, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$Glycine, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/GlycineScatter_Lesion_Waviness.pdf")
plot4 <- ggplot(BC_InteractWaviness, aes(Glycine,avg_pheno)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  geom_text(x = 2, y = 6, label = (rvalueSpea), parse = TRUE) +
  xlab("Mean Lesion Size on G. max") + ylab("Average Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Mean Lesion Size on Glycine max against Hyphal Waviness per Isolate") +
  geom_smooth(method='lm')
plot4
dev.off()
##Helianthus - dicot

pvaluePear <- summary(lm(avg_pheno~Helianthus, BC_InteractWaviness))$coefficients["Helianthus","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$Helianthus, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$Helianthus, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/HelianthusScatter_Lesion_Waviness.pdf")
plot5 <- ggplot(BC_InteractWaviness, aes(Helianthus,avg_pheno)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  xlab("Mean Lesion Size on H. annuus") + ylab("Average Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Mean Lesion Size on Helianthus annuus against Hyphal Waviness per Isolate") +
  geom_smooth(method='lm')
plot5
dev.off()
##Solanum - dicot

pvaluePear <- summary(lm(avg_pheno~Cendivia, BC_InteractWaviness))$coefficients["Cendivia","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$Cendivia, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$Cendivia, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/SolanumScatter_Lesion_Waviness.pdf")
plot6 <- ggplot(BC_InteractWaviness, aes(Solanum,avg_pheno)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  geom_text(x = 2, y = 6, label = (rvalueSpea), parse = TRUE) +
  xlab("Mean Lesion Size on Solanum") + ylab("Average Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Mean Lesion Size on Solanum against Hyphal Waviness per Isolate") +
  geom_smooth(method='lm')
plot6
dev.off()

##Average Eudicot Data




##Multiple Plots in one
library(cowplot)
gridPlot <- plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, plotTom, labels = c("A", "B", "C", "D", "E", "F", "G"), scale = .75)
gridPlot

#################################
#BOXPLOT: comparison of hyphal waviness across isolates (try violin or box and whiskers plot)

rm(list=ls())
library(ggplot2)
library(dplyr)

setwd("~/../Desktop/B. cinera/")
HwaviBC <- read.csv("BR_JV_ManualHyphalDat_032817(edit).csv")
Hwavi <- HwaviBC
pdf("../HWaviness Parts/Completed/Boxplot_Hwavi_Isolate.pdf", width = 15, height = 5)
ggplot(Hwavi, aes(reorder(Isolate, Phenotype, mean),Phenotype)) +
  theme_bw() + 
  geom_boxplot() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  xlab("Isolate") + ylab("Phenotype Distribution") 

dev.off()

#################################
# VIOLIN PLOT: comparison of hyphal waviness across isolates

rm(list=ls())
library(ggplot2)
library(dplyr)

setwd("~/../Desktop/B. cinera/Hwaviness/")
HwaviBC <- read.csv("../BR_JV_ManualHyphalDat_032817(edit).csv")
Hwavi <- HwaviBC
pdf("../../HWaviness Parts/Completed/Violin_Hwavi_Isolate.pdf", width = 15, height = 5)
ggplot(Hwavi, aes(reorder(Isolate, Phenotype, mean),Phenotype)) + 
  theme_bw() + 
  geom_violin(fill='#56B4E9', trim = FALSE) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  xlab("Isolate") + ylab("Hyphal Waviness")
dev.off()


#################################

#Scatter of Eccentricity of Isolates

rm(list=ls())
library(ggplot2)
library(dplyr)

setwd("~/../Desktop/B. cinera/")

Ecc <- read.csv("Hwaviness/data/BcAtPhenos/BcAtPhenosGWAS_lsmeans.fxmod1_clean.csv")

HwaviBC <- read.csv("BR_JV_ManualHyphalDat_032817(edit).csv")
Hwavi <- HwaviBC


Isolate <- Hwavi %>% 
  group_by(Isolate) %>%
  summarise(avg_pheno = mean(Phenotype, na.rm = TRUE), 
            min_pheno = min(Phenotype, na.rm = TRUE), 
            max_pheno = max(Phenotype, na.rm = TRUE),
            sd_pheno = sd(Phenotype, na.rm = TRUE),
            total = n())

BC_EccWaviness <- merge(Ecc, Isolate, by = "Isolate")

pvaluePear <- summary(lm(avg_pheno~Estimate, BC_EccWaviness))$coefficients["Estimate","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BC_EccWaviness$avg_pheno,BC_EccWaviness$Estimate, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BC_EccWaviness$avg_pheno,BC_EccWaviness$Estimate, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/EccenScatter_Lesion_Waviness.pdf")
plot6 <- ggplot(BC_EccWaviness, aes(Estimate,avg_pheno)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  geom_text(x = 2, y = 6, label = (rvalueSpea), parse = TRUE) +
  xlab("Mean Eccentricity") + ylab("Average Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Mean Lesion Size on Solanum against Hyphal Waviness per Isolate") +
  geom_smooth(method='lm')
plot6
dev.off()
##################################################################
#LsMean comparisons of the plots instead of the averages

##Pectin+Sugar x LsMeans
rm(list=ls())
library(ggplot2)
library(dplyr)

setwd("~/../Desktop/B. cinera/")

library(readr)
SugarPectin <- read_delim("C:/Users/vegaj/Desktop/B. cinera/HWavinessVisuals/SugarPectin_Lsmeans.csv", ";", escape_double = FALSE, trim_ws = TRUE)
SugarPectin <- SugarPectin [,c("GenoRename", "P72", "P48", "S72", "S48")]
colnames(SugarPectin)<-c("Isolate", "P72", "P48", "S72", "S48")
HwaviBC <- read.csv("WavyGWAS_lsmeans.fxmod1_R_output.csv")
Isolate <- HwaviBC

SugarPectinWaviness <- merge(SugarPectin, Isolate, by = "Isolate")
View(SugarPectinWaviness)

#Scatter Plot of Pectin/Sugar(72+48hrs) v Waviness
##P72

rvaluePear <- summary(lm(HwaviEstimate~P72, SugarPectinWaviness))$r.squared 
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(HwaviEstimate~P72, SugarPectinWaviness))$coefficients["P72","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(SugarPectinWaviness$HwaviEstimate,SugarPectinWaviness$P72, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(SugarPectinWaviness$HwaviEstimate,SugarPectinWaviness$P72, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/LsMeanScatters/Pectin72hScatter_LsWaviness.pdf")
plot1 <- ggplot(SugarPectinWaviness, aes(P72,HwaviEstimate)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  xlab("Isolate Growth on Pectin after 72hr") + ylab("LsMeans Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Comparison of Isolate growth on Pectin after 72hr to Isolate Hyphal Waviness") +
  geom_smooth(method='lm')
plot1
dev.off()

##P48

rvaluePear <- summary(lm(HwaviEstimate~P48, SugarPectinWaviness))$r.squared 
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(HwaviEstimate~P48, SugarPectinWaviness))$coefficients["P48","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(SugarPectinWaviness$HwaviEstimate,SugarPectinWaviness$P48, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(SugarPectinWaviness$HwaviEstimate,SugarPectinWaviness$P48, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/LsMeanScatters/Pectin48hScatter_LsWaviness.pdf")
plot2 <- ggplot(SugarPectinWaviness, aes(P48,HwaviEstimate)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  xlab("Isolate Growth on Pectin after 48hr") + ylab("LsMeans Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Comparison of Isolate growth on Pectin after 48hr to Isolate Hyphal Waviness") +
  geom_smooth(method='lm')
plot2
dev.off()
##S72

rvaluePear <- summary(lm(HwaviEstimate~S72, SugarPectinWaviness))$r.squared 
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(HwaviEstimate~S72, SugarPectinWaviness))$coefficients["S72","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(SugarPectinWaviness$HwaviEstimate,SugarPectinWaviness$S72, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(SugarPectinWaviness$HwaviEstimate,SugarPectinWaviness$S72, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/LsMeanScatters/Sugar72hScatter_LsWaviness.pdf")
plot3 <- ggplot(SugarPectinWaviness, aes(S72,HwaviEstimate)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  xlab("Isolate Growth on Sugar after 72hr") + ylab("LsMeans Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Comparison of Isolate growth on Sugar after 72hr to Isolate Hyphal Waviness") +
  geom_smooth(method='lm')
plot3
dev.off()
##S48

rvaluePear <- summary(lm(HwaviEstimate~S48, SugarPectinWaviness))$r.squared 
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(HwaviEstimate~S48, SugarPectinWaviness))$coefficients["S48","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(SugarPectinWaviness$HwaviEstimate,SugarPectinWaviness$S48, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(SugarPectinWaviness$HwaviEstimate,SugarPectinWaviness$S48, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/LsMeanScatters/Sugar48hScatter_LsWaviness.pdf")
plot4 <- ggplot(SugarPectinWaviness, aes(S48,HwaviEstimate)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  xlab("Isolate Growth on Sugar after 48hr") + ylab("LsMeans Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Comparison of Isolate growth on Sugar after 48hr to Isolate Hyphal Waviness") +
  geom_smooth(method='lm')
plot4
dev.off()

##Multiple Plots in one
library(cowplot)
gridPlot <- plot_grid( plot2, plot1,plot4, plot3, labels = c("A", "B", "C", "D"), scale = .75)
gridPlot


#Scatter of Eccentricity of Isolates x LsMeans

rm(list=ls())
library(ggplot2)
library(dplyr)

setwd("~/../Desktop/B. cinera/")

Ecc <- read.csv("Hwaviness/data/BcAtPhenos/BcAtPhenosGWAS_lsmeans.fxmod1_clean.csv")

HwaviBC <- read.csv("WavyGWAS_lsmeans.fxmod1_R_output.csv")
Isolate <- HwaviBC

BC_EccWaviness <- merge(Ecc, Isolate, by = "Isolate")

rvaluePear <- summary(lm(HwaviEstimate~Estimate, BC_EccWaviness))$r.squared 
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(HwaviEstimate~Estimate, BC_EccWaviness))$coefficients["Estimate","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BC_EccWaviness$HwaviEstimate,BC_EccWaviness$Estimate, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BC_EccWaviness$HwaviEstimate,BC_EccWaviness$Estimate, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/LsMeanScatters/EccenScatter_Lesion_LsWaviness.pdf")
plot6 <- ggplot(BC_EccWaviness, aes(Estimate,HwaviEstimate)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  geom_text(x = 2, y = 6, label = (rvalueSpea), parse = TRUE) +
  xlab("Mean Eccentricity") + ylab("LsMeans Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Mean Lesion Size on Solanum against Hyphal Waviness per Isolate") +
  geom_smooth(method='lm')
plot6
dev.off()

# scatter plot of Lsmeanshyphal waviness vs. lesion size on tomato, by isolate

rm(list=ls())
library(ggplot2)
library(dplyr)

setwd("~/../Desktop/B. cinera/")

BC_tomato <- read.csv("HWavinessVisuals/Bc_tomatoLesion.csv")
HwaviBC <- read.csv("WavyGWAS_lsmeans.fxmod1_R_output.csv")
Isolate <- HwaviBC

BC_TomatoWaviness <- merge(BC_tomato, Isolate, by = "Isolate")

#Scatter Plot mean Lesion v Waviness
rvaluePear <- summary(lm(HwaviEstimate~meanLesion, BC_TomatoWaviness))$r.squared
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(HwaviEstimate~meanLesion, BC_TomatoWaviness))$coefficients["meanLesion","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BC_TomatoWaviness$HwaviEstimate,BC_TomatoWaviness$meanLesion, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BC_TomatoWaviness$HwaviEstimate,BC_TomatoWaviness$meanLesion, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/LsMeanScatters/TomatoScatter_Lesion_Waviness.pdf")
plotTom <- ggplot(BC_TomatoWaviness, aes(meanLesion,HwaviEstimate)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  xlab("Mean Lesion Size on Tomato") + ylab("LsMeans Hyphal Waviness") +
  ggtitle(NULL) +
  #ggtitle("Mean Lesion Size on Tomato against Hyphal Waviness per Isolate") +
  geom_smooth(method='lm', se = FALSE)
plotTom
dev.off()

# 2nd scatter plot of hyphal waviness vs. lesion size on different interaction types, by isolate

setwd("~/../Desktop/B. cinera/")

BC_Interact <- read.csv("HWavinessVisuals/DataFor100IsoCollection.csv")
BC_Interact <- BC_Interact [,c("Isolate","Cendivia","Brapa","Cintybus","Glycine",	"Helianthus",	"Solanum")]

BC_InteractWaviness <- merge(BC_Interact, Isolate, by = "Isolate")

#Scatter Plot mean Interactions Lesion v Waviness
##C. endivia - dicot

rvaluePear <- summary(lm(HwaviEstimate~Cendivia, BC_InteractWaviness))$r.squared 
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(HwaviEstimate~Cendivia, BC_InteractWaviness))$coefficients["Cendivia","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$Cendivia, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$Cendivia, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/LsMeanScatters/CendiviaScatter_Lesion_Waviness.pdf")
plot1 <- ggplot(BC_InteractWaviness, aes(Cendivia,HwaviEstimate)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) + #geom_text(label = (rvalueSpea), parse = TRUE) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  xlab("Mean Lesion Size on C. endivia") + ylab("LsMeans Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Mean Lesion Size on Cichorium endivia against Hyphal Waviness per Isolate") +
  geom_smooth(method='lm')
plot1
dev.off()

##B. rapa - dicot

rvaluePear <- summary(lm(HwaviEstimate~Brapa, BC_InteractWaviness))$r.squared 
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(HwaviEstimate~Brapa, BC_InteractWaviness))$coefficients["Brapa","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$Brapa, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$Brapa, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/LsMeanScatters/BrapaScatter_Lesion_Waviness.pdf")
plot2 <- ggplot(BC_InteractWaviness, aes(Brapa,HwaviEstimate)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  geom_text(x = 2, y = 6, label = (rvalueSpea), parse = TRUE) +
  xlab("Mean Lesion Size on B. rapa") + ylab("LsMeans Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Mean Lesion Size on Brassica rapa against Hyphal Waviness per Isolate") +
  geom_smooth(method='lm')
plot2
dev.off()
##C. intybus - dicot

rvaluePear <- summary(lm(HwaviEstimate~Cintybus, BC_InteractWaviness))$r.squared 
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(HwaviEstimate~Cintybus, BC_InteractWaviness))$coefficients["Cintybus","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$Cintybus, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$Cintybus, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/LsMeanScatters/CintybusScatter_Lesion_Waviness.pdf")
plot3 <- ggplot(BC_InteractWaviness, aes(Cintybus,HwaviEstimate)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  geom_text(x = 2, y = 6, label = (rvalueSpea), parse = TRUE) +
  xlab("Mean Lesion Size on C. intybus") + ylab("LsMeans Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Mean Lesion Size on Cichorium intybus against Hyphal Waviness per Isolate") +
  geom_smooth(method='lm')
plot3
dev.off()
##Glycine max

rvaluePear <- summary(lm(HwaviEstimate~Glycine, BC_InteractWaviness))$r.squared
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(HwaviEstimate~Glycine, BC_InteractWaviness))$coefficients["Glycine","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$Glycine, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$Glycine, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/LsMeanScatters/GlycineScatter_Lesion_Waviness.pdf")
plot4 <- ggplot(BC_InteractWaviness, aes(Glycine,HwaviEstimate)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  geom_text(x = 2, y = 6, label = (rvalueSpea), parse = TRUE) +
  xlab("Mean Lesion Size on G. max") + ylab("LsMeans Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Mean Lesion Size on Glycine max against Hyphal Waviness per Isolate") +
  geom_smooth(method='lm')
plot4
dev.off()
##Helianthus - dicot

rvaluePear <- summary(lm(HwaviEstimate~Helianthus, BC_InteractWaviness))$r.squared 
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(HwaviEstimate~Helianthus, BC_InteractWaviness))$coefficients["Helianthus","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$Helianthus, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$Helianthus, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/LsMeanScatters/HelianthusScatter_Lesion_Waviness.pdf")
plot5 <- ggplot(BC_InteractWaviness, aes(Helianthus,HwaviEstimate)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  xlab("Mean Lesion Size on H. annuus") + ylab("LsMeans Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Mean Lesion Size on Helianthus annuus against Hyphal Waviness per Isolate") +
  geom_smooth(method='lm')
plot5
dev.off()
##Solanum - dicot

rvaluePear <- summary(lm(HwaviEstimate~Solanum, BC_InteractWaviness))$r.squared 
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(HwaviEstimate~Solanum, BC_InteractWaviness))$coefficients["Solanum","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$Solanum, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$Solanum, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/LsMeanScatters/SolanumScatter_Lesion_Waviness.pdf")
plot6 <- ggplot(BC_InteractWaviness, aes(Solanum,HwaviEstimate)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  geom_text(x = 2, y = 6, label = (rvalueSpea), parse = TRUE) +
  xlab("Mean Lesion Size on Solanum") + ylab("LsMeans Hyphal Waviness") + 
  ggtitle(NULL) +
  #ggtitle("Mean Lesion Size on Solanum against Hyphal Waviness per Isolate") +
  geom_smooth(method='lm')
plot6
dev.off()

##Multiple Plots in one
library(cowplot)
gridPlot <- plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, plotTom, labels = c("A", "B", "C", "D", "E", "F", "G"), scale = .75)
gridPlot
##################################################################

# scatter plot of Lsmeanshyphal waviness vs. Ecc per plant, by isolate

rm(list=ls())
library(ggplot2)
library(dplyr)

setwd("~/../Desktop/B. cinera/")

BcAt <- read.csv("Hwaviness/data/BcAtPhenos/BcAtPhenosGWAS_lsmeans.fxmod1_organiz.csv")
HwaviBC <- read.csv("WavyGWAS_lsmeans.fxmod1_R_output.csv")
Isolate <- HwaviBC

BcAt_LsWaviness <- merge(BcAt, Isolate, by = "Isolate")

#Scatter Plot mean Interactions Lesion v Waviness
## anac088 - dicot

rvaluePear <- summary(lm(HwaviEstimate~anac055, BcAt_LsWaviness))$r.squared 
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(HwaviEstimate~anac055, BcAt_LsWaviness))$coefficients["anac055","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$anac055, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$anac055, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/LsMeanScatters/anac088Scatter_Lesion_LsWaviness.pdf")
plot1 <- ggplot(BcAt_LsWaviness, aes(anac055,HwaviEstimate)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) + #geom_text(label = (rvalueSpea), parse = TRUE) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  xlab("Lesion Eccentricity on anac088") + ylab("LsMeans Hyphal Waviness") + 
  ggtitle(NULL) +
  geom_smooth(method='lm')
plot1
dev.off()

##coi1 - dicot

rvaluePear <- summary(lm(HwaviEstimate~coi1, BcAt_LsWaviness))$r.squared 
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(HwaviEstimate~coi1, BcAt_LsWaviness))$coefficients["coi1","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$coi1, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$coi1, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/LsMeanScatters/coi1Scatter_Lesion_LsWaviness.pdf")
plot2 <- ggplot(BcAt_LsWaviness, aes(coi1,HwaviEstimate)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  geom_text(x = 2, y = 6, label = (rvalueSpea), parse = TRUE) +
  xlab("Lesion Eccentricity on coi1") + ylab("LsMeans Hyphal Waviness") + 
  ggtitle(NULL) +
  geom_smooth(method='lm')
plot2
dev.off()
##col0 - dicot

rvaluePear <- summary(lm(HwaviEstimate~col0, BcAt_LsWaviness))$r.squared 
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(HwaviEstimate~col0, BcAt_LsWaviness))$coefficients["col0","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$col0, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$col0, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/LsMeanScatters/col0Scatter_Lesion_LsWaviness.pdf")
plot3 <- ggplot(BcAt_LsWaviness, aes(col0,HwaviEstimate)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  geom_text(x = 2, y = 6, label = (rvalueSpea), parse = TRUE) +
  xlab("Lesion Eccentricity on col0") + ylab("LsMeans Hyphal Waviness") + 
  ggtitle(NULL) +
  geom_smooth(method='lm')
plot3
dev.off()
##npr1

rvaluePear <- summary(lm(HwaviEstimate~npr1, BcAt_LsWaviness))$r.squared
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(HwaviEstimate~npr1, BcAt_LsWaviness))$coefficients["npr1","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$npr1, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$npr1, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/LsMeanScatters/npr1Scatter_Lesion_LsWaviness.pdf")
plot4 <- ggplot(BcAt_LsWaviness, aes(npr1,HwaviEstimate)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  geom_text(x = 2, y = 6, label = (rvalueSpea), parse = TRUE) +
  xlab("Lesion Eccentricity on npr1") + ylab("LsMeans Hyphal Waviness") + 
  ggtitle(NULL) +
  geom_smooth(method='lm')
plot4
dev.off()
##pad3 - dicot

rvaluePear <- summary(lm(HwaviEstimate~pad3, BcAt_LsWaviness))$r.squared 
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(HwaviEstimate~pad3, BcAt_LsWaviness))$coefficients["pad3","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$pad3, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$pad3, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/LsMeanScatters/pad3Scatter_Lesion_LsWaviness.pdf")
plot5 <- ggplot(BcAt_LsWaviness, aes(pad3,HwaviEstimate)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  xlab("Lesion Eccentricity on pad3") + ylab("LsMeans Hyphal Waviness") + 
  ggtitle(NULL) +
  geom_smooth(method='lm')
plot5
dev.off()

##tga3 - dicot

rvaluePear <- summary(lm(HwaviEstimate~tga3, BcAt_LsWaviness))$r.squared 
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(HwaviEstimate~tga3, BcAt_LsWaviness))$coefficients["tga3","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$tga3, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$HwaviEstimate,BcAt_LsWaviness$tga3, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/LsMeanScatters/tga3Scatter_Lesion_LsWaviness.pdf")
plot6 <- ggplot(BcAt_LsWaviness, aes(tga3,HwaviEstimate)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  xlab("Lesion Eccentricity on tga3") + ylab("LsMeans Hyphal Waviness") + 
  ggtitle(NULL) +
  geom_smooth(method='lm')
plot5
dev.off()

##Multiple Plots in one
library(cowplot)
gridPlot <- plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, labels = c("A", "B", "C", "D", "E", "F"), scale = .75)
gridPlot

##################################################################

# scatter plot of Lsmeanshyphal waviness vs. Ecc per plant, by isolate

rm(list=ls())
library(ggplot2)
library(dplyr)

setwd("~/../Desktop/B. cinera/")

BcAt <- read.csv("Hwaviness/data/BcAtPhenos/BcAtPhenosGWAS_lsmeans.fxmod1_organiz.csv")
HwaviBC <- read.csv("BR_JV_ManualHyphalDat_032817(edit).csv")
Hwavi <- HwaviBC

Isolate <- Hwavi %>% 
  group_by(Isolate) %>%
  summarise(avg_pheno = mean(Phenotype, na.rm = TRUE), 
            min_pheno = min(Phenotype, na.rm = TRUE), 
            max_pheno = max(Phenotype, na.rm = TRUE),
            sd_pheno = sd(Phenotype, na.rm = TRUE),
            total = n())

BcAt_LsWaviness <- merge(BcAt, Isolate, by = "Isolate")

#Scatter Plot mean Interactions Lesion v Waviness
## anac088 - dicot

rvaluePear <- summary(lm(avg_pheno~anac055, BcAt_LsWaviness))$r.squared 
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(avg_pheno~anac055, BcAt_LsWaviness))$coefficients["anac055","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$anac055, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$anac055, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)


pdf("../HWaviness Parts/Completed/anac088Scatter_Lesion_LsWaviness.pdf")
plot1 <- ggplot(BcAt_LsWaviness, aes(anac055,avg_pheno)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) + #geom_text(label = (rvalueSpea), parse = TRUE) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  xlab("Lesion Eccentricity on anac088") + ylab("Average Hyphal Waviness") + 
  ggtitle(NULL) +
  geom_smooth(method='lm')
plot1
dev.off()

##coi1 - dicot

rvaluePear <- summary(lm(avg_pheno~coi1, BcAt_LsWaviness))$r.squared 
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(avg_pheno~coi1, BcAt_LsWaviness))$coefficients["coi1","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$coi1, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$coi1, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/coi1Scatter_Lesion_LsWaviness.pdf")
plot2 <- ggplot(BcAt_LsWaviness, aes(coi1,avg_pheno)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  geom_text(x = 2, y = 6, label = (rvalueSpea), parse = TRUE) +
  xlab("Lesion Eccentricity on coi1") + ylab("Average Hyphal Waviness") + 
  ggtitle(NULL) +
  geom_smooth(method='lm')
plot2
dev.off()
##col0 - dicot

rvaluePear <- summary(lm(avg_pheno~col0, BcAt_LsWaviness))$r.squared 
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(avg_pheno~col0, BcAt_LsWaviness))$coefficients["col0","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$col0, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$col0, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/col0Scatter_Lesion_LsWaviness.pdf")
plot3 <- ggplot(BcAt_LsWaviness, aes(col0,avg_pheno)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  geom_text(x = 2, y = 6, label = (rvalueSpea), parse = TRUE) +
  xlab("Lesion Eccentricity on col0") + ylab("Average Hyphal Waviness") + 
  ggtitle(NULL) +
  geom_smooth(method='lm')
plot3
dev.off()
##npr1

rvaluePear <- summary(lm(avg_pheno~npr1, BcAt_LsWaviness))$r.squared
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(avg_pheno~npr1, BcAt_LsWaviness))$coefficients["npr1","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$npr1, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$npr1, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/npr1Scatter_Lesion_LsWaviness.pdf")
plot4 <- ggplot(BcAt_LsWaviness, aes(npr1,avg_pheno)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  geom_text(x = 2, y = 6, label = (rvalueSpea), parse = TRUE) +
  xlab("Lesion Eccentricity on npr1") + ylab("Average Hyphal Waviness") + 
  ggtitle(NULL) +
  geom_smooth(method='lm')
plot4
dev.off()
##pad3 - dicot

rvaluePear <- summary(lm(avg_pheno~pad3, BcAt_LsWaviness))$r.squared 
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(avg_pheno~pad3, BcAt_LsWaviness))$coefficients["pad3","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$pad3, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$pad3, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/pad3Scatter_Lesion_LsWaviness.pdf")
plot5 <- ggplot(BcAt_LsWaviness, aes(pad3,avg_pheno)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  xlab("Lesion Eccentricity on pad3") + ylab("Average Hyphal Waviness") + 
  ggtitle(NULL) +
  geom_smooth(method='lm')
plot5
dev.off()

##tga3 - dicot

rvaluePear <- summary(lm(avg_pheno~tga3, BcAt_LsWaviness))$r.squared 
rvaluePear <- format(round(rvaluePear, 5), nsmall = 4)
pvaluePear <- summary(lm(avg_pheno~tga3, BcAt_LsWaviness))$coefficients["tga3","Pr(>|t|)"]
pvaluePear <- format(round(pvaluePear, 5), nsmall = 4)
rvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$tga3, method="spearman", exact=FALSE)$p.value
rvalueSpea <- format(round(rvalueSpea, 5), nsmall = 4)
pvalueSpea <- cor.test(BcAt_LsWaviness$avg_pheno,BcAt_LsWaviness$tga3, method="spearman", exact=FALSE)$estimate
pvalueSpea <- format(round(pvalueSpea, 5), nsmall = 4)

pdf("../HWaviness Parts/Completed/tga3Scatter_Lesion_Waviness.pdf")
plot6 <- ggplot(BcAt_LsWaviness, aes(tga3,avg_pheno)) + 
  geom_point(color='red') + #geom_text(aes(label=Isolate), position = position_nudge(y = -0.1),  size=3) +
  annotate(geom = 'text', label = paste('r2 = ', rvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -1) +
  annotate(geom = 'text', label = paste('p = ', pvalueSpea), x = Inf, y = -Inf, hjust = 1, vjust = -2.5) +
  xlab("Lesion Eccentricity on tga3") + ylab("Average Hyphal Waviness") + 
  ggtitle(NULL) +
  geom_smooth(method='lm')
plot5
dev.off()

##Multiple Plots in one
library(cowplot)
gridPlot <- plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, labels = c("A", "B", "C", "D", "E", "F"), scale = .75)
gridPlot