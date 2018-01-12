#01b_EccentricityPhenos.R
#NES
#01/12/18

#prepare lesion eccentricity phenotypes of B. cinerea on A. thaliana for correlation to B. cinerea hyphal waviness
#-------------------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/Hwaviness/data/BcAtPhenos")

isolist <- read.csv("IsolateKey.csv")
myphenos <- read.csv("BcAt_LesionImageData.csv")
names(myphenos)[6] <- "Sample"

allphenos <- merge(myphenos,isolist, by="Sample")

allphenos <- allphenos[,c("Experiment","GrowingFlat","AgarFlat","Plant","Isolate","Lesion.0.m.eccentricity","Lesion.red.m.eccentricity","Lesion.grn.m.eccentricity", "Lesion.blu.m.eccentricity", "Lesion.redgrn.m.eccentricity","Lesion.redblu.m.eccentricity", "Lesion.grnblu.m.eccentricity", "Lesion.Bred.m.eccentricity", "Lesion.Bgrn.m.eccentricity", "Lesion.Bblu.m.eccentricity", "Lesion.Bredgrn.m.eccentricity", "Lesion.Bredblu.m.eccentricity", "Lesion.Bgrnblu.m.eccentricity" )]

write.csv(allphenos, "BcAt_LesionEccentricity.csv")
