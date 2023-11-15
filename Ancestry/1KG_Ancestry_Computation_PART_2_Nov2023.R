############################################
#### Author: Dan Gustavson & Aysu Okbay ####
#### Date:   November 13 2023           ####
############################################


###################################################################################
#### PART 2  (COMMANDS IN R )                                                   ###
### Plot CATSLife against 1KG ancestry groups, then plot computed ancestries    ###
###     against 1KG. Finally, combine all IDs together for ancestry group file  ###
###################################################################################
### NOTE: Some AMR IDs were also labelled EUR. This script labels EUR first, then # 
###   adds AMR groups for those not already labelled EUR.                       ###
### NOTE: There are no South Asian's identified in CATSLife                     ###
###   so some changes to this script may be needed if your sample has SAS       ###
###################################################################################

library(data.table)
library(ggplot2)
library(tidyverse)
#args=commandArgs(trailingOnly=TRUE)
#PCs=args[1]
#out=args[2]

#####################################################
### PLOTS AGAINST 1KG BEFORE ANCESTRY COMPUTATION ###
#####################################################

data <-fread("CATS_1kG_HM3_PCs.eigenvec.annotated",header=T)
pdf(file="output_PCplots.pdf")  

PC12<-ggplot(data, aes(x=PC1, y=PC2, color=SUPERPOP)) + geom_point()
PC13<-ggplot(data, aes(x=PC1, y=PC3, color=SUPERPOP)) + geom_point()
PC14<-ggplot(data, aes(x=PC1, y=PC4, color=SUPERPOP)) + geom_point()
PC23<-ggplot(data, aes(x=PC2, y=PC3, color=SUPERPOP)) + geom_point()
PC24<-ggplot(data, aes(x=PC2, y=PC4, color=SUPERPOP)) + geom_point()
PC34<-ggplot(data, aes(x=PC3, y=PC4, color=SUPERPOP)) + geom_point()

print(PC12)
print(PC13)
print(PC14)
print(PC23)
print(PC24)
print(PC34)

dev.off()


##########################################
### PLOT ANCESTRY GROUPS TO CHECK THEM ###
##########################################


#### EUROPEAN ####

data <-fread("CATS_EUR_1kG_HM3_PCs.eigenvec.annotated",header=T)

pdf(file="output_PCplots_EUR.pdf")  

PC12<-ggplot(data, aes(x=PC1, y=PC2, color=SUPERPOP)) + geom_point()
PC13<-ggplot(data, aes(x=PC1, y=PC3, color=SUPERPOP)) + geom_point()
PC14<-ggplot(data, aes(x=PC1, y=PC4, color=SUPERPOP)) + geom_point()
PC23<-ggplot(data, aes(x=PC2, y=PC3, color=SUPERPOP)) + geom_point()
PC24<-ggplot(data, aes(x=PC2, y=PC4, color=SUPERPOP)) + geom_point()
PC34<-ggplot(data, aes(x=PC3, y=PC4, color=SUPERPOP)) + geom_point()

print(PC12)
print(PC13)
print(PC14)
print(PC23)
print(PC24)
print(PC34)

dev.off()



#### SOUTH ASIAN ####
### NOTE: None identified
data <-fread("CATS_SAS_1kG_HM3_PCs.eigenvec.annotated",header=T)
pdf(file="output_PCplots_SAS.pdf")  

PC12<-ggplot(data, aes(x=PC1, y=PC2, color=SUPERPOP)) + geom_point()
PC13<-ggplot(data, aes(x=PC1, y=PC3, color=SUPERPOP)) + geom_point()
PC14<-ggplot(data, aes(x=PC1, y=PC4, color=SUPERPOP)) + geom_point()
PC23<-ggplot(data, aes(x=PC2, y=PC3, color=SUPERPOP)) + geom_point()
PC24<-ggplot(data, aes(x=PC2, y=PC4, color=SUPERPOP)) + geom_point()
PC34<-ggplot(data, aes(x=PC3, y=PC4, color=SUPERPOP)) + geom_point()

print(PC12)
print(PC13)
print(PC14)
print(PC23)
print(PC24)
print(PC34)

dev.off()

#### EAST ASIAN ####
data <-fread("CATS_EAS_1kG_HM3_PCs.eigenvec.annotated",header=T)
pdf(file="output_PCplots_EAS.pdf")  

PC12<-ggplot(data, aes(x=PC1, y=PC2, color=SUPERPOP)) + geom_point()
PC13<-ggplot(data, aes(x=PC1, y=PC3, color=SUPERPOP)) + geom_point()
PC14<-ggplot(data, aes(x=PC1, y=PC4, color=SUPERPOP)) + geom_point()
PC23<-ggplot(data, aes(x=PC2, y=PC3, color=SUPERPOP)) + geom_point()
PC24<-ggplot(data, aes(x=PC2, y=PC4, color=SUPERPOP)) + geom_point()
PC34<-ggplot(data, aes(x=PC3, y=PC4, color=SUPERPOP)) + geom_point()

print(PC12)
print(PC13)
print(PC14)
print(PC23)
print(PC24)
print(PC34)

dev.off()


#### AMERICAN ####
data <-fread("CATS_AMR_1kG_HM3_PCs.eigenvec.annotated",header=T)
pdf(file="output_PCplots_AMR.pdf")  

PC12<-ggplot(data, aes(x=PC1, y=PC2, color=SUPERPOP)) + geom_point()
PC13<-ggplot(data, aes(x=PC1, y=PC3, color=SUPERPOP)) + geom_point()
PC14<-ggplot(data, aes(x=PC1, y=PC4, color=SUPERPOP)) + geom_point()
PC23<-ggplot(data, aes(x=PC2, y=PC3, color=SUPERPOP)) + geom_point()
PC24<-ggplot(data, aes(x=PC2, y=PC4, color=SUPERPOP)) + geom_point()
PC34<-ggplot(data, aes(x=PC3, y=PC4, color=SUPERPOP)) + geom_point()

print(PC12)
print(PC13)
print(PC14)
print(PC23)
print(PC24)
print(PC34)

dev.off()

#### AFRICAN ####
data <-fread("CATS_AFR_1kG_HM3_PCs.eigenvec.annotated",header=T)
pdf(file="output_PCplots_AFR.pdf")  

PC12<-ggplot(data, aes(x=PC1, y=PC2, color=SUPERPOP)) + geom_point()
PC13<-ggplot(data, aes(x=PC1, y=PC3, color=SUPERPOP)) + geom_point()
PC14<-ggplot(data, aes(x=PC1, y=PC4, color=SUPERPOP)) + geom_point()
PC23<-ggplot(data, aes(x=PC2, y=PC3, color=SUPERPOP)) + geom_point()
PC24<-ggplot(data, aes(x=PC2, y=PC4, color=SUPERPOP)) + geom_point()
PC34<-ggplot(data, aes(x=PC3, y=PC4, color=SUPERPOP)) + geom_point()

print(PC12)
print(PC13)
print(PC14)
print(PC23)
print(PC24)
print(PC34)

dev.off()


########################################
### Pull in ancestry files and merge ###
########################################
# Note: AMR seems to be a pretty wide group so we only classify as AMR
#       when they aren't in other groups

EUR <- read.table("CATS_EUR_FID_IID.txt")
colnames(EUR) <- c("FID","IID")
EUR$Ancestry <- "EUR"
EUR1 <- select(EUR, IID, Ancestry)

AFR <- read.table("CATS_AFR_FID_IID.txt")
colnames(AFR) <- c("FID","IID")
AFR$Ancestry <- "AFR"
AFR1 <- select(AFR, IID, Ancestry)

# NO IDs IN DATAFILE - CHECK AGAIN IN FUTURE
#SAS <- read.table("CATS_SAS_FID_IID.txt")
#colnames(SAS) <- c("FID","IID")
#SAS$Ancestry <- "SAS"
#SAS1 <- select(SAS, IID, Ancestry)

EAS <- read.table("CATS_EAS_FID_IID.txt")
colnames(EAS) <- c("FID","IID")
EAS$Ancestry <- "EAS"
EAS1 <- select(EAS, IID, Ancestry)

AMR <- read.table("CATS_AMR_FID_IID.txt")
colnames(AMR) <- c("FID","IID")
AMR$Ancestry <- "AMR"
AMR1 <- select(AMR, IID, Ancestry)

EUR_AFR_EAS <- rbind(EUR1, AFR1, EAS1)
table(EUR_AFR_EAS$Ancestry)


length(EUR_AFR_EAS$IID)
length(unique(EUR_AFR_EAS$IID))

# American has overlap with other groups so exclude those first
AMR2 <- filter(AMR1, !IID %in% unique(EUR_AFR_EAS$IID))
length(AMR1$IID)
length(AMR2$IID)

ALL <- rbind(EUR1, AFR1, EAS1, AMR2)
length(ALL$IID)
length(unique(ALL$IID))

write.table(ALL, "CATSLife_Ancestry_Nov2023.txt", quote=F, sep="\t")

