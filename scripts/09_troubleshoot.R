genes <- read.csv("geneID_99thr_fxn_estimate.csv")
genes <- read.csv("HWaviGeneID_fxn_estimateThresh99.csv")
names(genes)

genes2 <- genes[genes$Chrom.Pos=="CHROMOSOME7.363510",]

genes3 <- genes[genes$geneID=="BcT4_6156",]
unique(genes3$Pos)
unique(genes3$PFAM_NAME)
genes3$Pos.PFAM <- paste(genes3$Pos, genes3$PFAM_NAME, sep="")
unique(genes3$Pos.PFAM) #720 unique SNP-PFAM combos
genes3$Pos.PFAM.Desc <- paste(genes3$Pos, genes3$PFAM_NAME, genes3$PFAM_DESCRIPTION, sep="")
hist(table(genes3$Pos.PFAM.Desc))


genes3 <- genes3[,-c(1)]
genes3$Est.PFAM.Desc <- paste(genes3$Estimate, genes3$PFAM_NAME,genes3$PFAM_DESCRIPTION, sep="")
unique(genes3$Est.PFAM.Desc)
genes3 <- unique(genes3)
