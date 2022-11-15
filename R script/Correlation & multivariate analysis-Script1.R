#Author: Numan Ibne Asad
#Project: Microbial based predictive modeling


#Correlation between microbiome data vs quality

library(magrittr)
library(car)
library(broom)
library(ggplot2)
library(tidyverse)
library(MASS)
library(car)
library(readxl)
library(tidyr)
library(dplyr)
library(data.table)
library(datasets)
library(ggpubr)
library(rstatix)
library (BiodiversityR)
library(ggplot2)
library(here)
#install.packages("grid", dependencies = TRUE, repos = "http://cran.us.r-project.org")
##ORDINATIONS

##Import and arrange objects

#Data pre processing
otu.16S=read.table(file=here("data","16s_OTU_field.txt"), row.names = 1,  header=T, comment.char = "", sep="\t")
taxa.16S=read.table(file=here("data","16s_tax.txt"), row.names = 1, header=T, comment.char = "", sep="\t")
treat.16S=read.table(file=here("data","mapping_file.field.txt"),row.names=1, header=T, sep="\t", comment.char="")

otu.16S=otu.16S[order(row.names(otu.16S)),] 
taxa.16S.S= taxa.16S[order(row.names(taxa.16S)),] 
table(row.names(otu.16S)==row.names(taxa.16S.S))

otu.16S.S =cbind(otu.16S, taxa.16S.S, by=0)


#Fungal community

otu.ITS=read.table(file=here("data","ITS.OTU.field.txt"), row.names=1, header=T, comment.char = "", sep="\t")
taxa.ITS=read.table(file=here("data","ITS_taxa.txt"),  header=T, comment.char = "", sep="\t")
taxa.ITS=taxa.ITS[-c(4869,4870),]

otu.ITS = otu.ITS[order(row.names(otu.ITS)),] 
rownames(taxa.ITS)=taxa.ITS[,1]
taxa.ITS= taxa.ITS[,-1]

taxa.ITS = taxa.ITS[order(row.names(taxa.ITS)),] 
row.names(otu.ITS)==row.names(taxa.ITS)
otu.ITS.s =cbind(otu.ITS, taxa.ITS, by = 0)

treat.ITS=read.table(file=here("data","mapping_file.field.txt"),row.names=1, header=T, sep="\t", comment.char="")
#treat.ITS=treat.ITS[-c(26,86,155,162,193,194),]


#Trimming
tax.16s=data.frame(otu.16S.S[,337:(dim(otu.16S.S)[2]-1)], row.names = paste("X",rownames(otu.16S),sep="")) #Last column, contain taxonomy string
otu.16S.t=otu.16S.S[,1:(dim(otu.16S.S)[2]-8)]

#Normalization
otu.16s.norm=t(otu.16S.t)
otu.16S.rel=(otu.16s.norm/rowSums(otu.16s.norm)) #Normalization
rowSums(otu.16S.rel) #Sanity check

#Trimming
tax.ITS=data.frame(otu.ITS.s[,338:dim(otu.ITS.s)[2]-1], row.names = paste("X",rownames(otu.ITS),sep="")) #Last column, contain taxonomy string
otu.ITS.t=otu.ITS.s[,1:(dim(otu.ITS.s)[2]-8)]

#Normalization
otu.ITS.norm=t(otu.ITS.t)
otu.ITS.rel=(otu.ITS.norm/rowSums(otu.ITS.norm)) #Normalization
rowSums(otu.ITS.rel)

#Sort by row names
treat.bact=treat.16S[order(row.names(treat.16S)),]
otu.16S.rels=otu.16S.rel[order(row.names(otu.16S.rel)),]
row.names(otu.16S.rels)==row.names(treat.bact) #Sanity check


treat.fungi=treat.ITS[order(row.names(treat.ITS)),]
otu.ITS.rels=otu.ITS.rel[order(row.names(otu.ITS.rel)),]
row.names(otu.ITS.rels)==row.names(treat.fungi) #Sanity check


## Analysis of beta diversity (Bray Curtis) and PCoA

#Calculate PCoA for whole data
library(vegan)
library(tidyverse)
bray.ITS.exp.field=vegdist(otu.ITS.rels,method = "bray")
pcoa.ITS.exp.field=cmdscale(sqrt(bray.ITS.exp.field),eig=T)

bray.16S.exp.field=vegdist(otu.16S.rels,method = "bray")
pcoa.16S.exp.field=cmdscale(sqrt(bray.16S.exp.field),eig=T)


#PERMANOVA analysis

perm.ITS.exp.field=adonis2(bray.ITS.exp.field~Cultivar*date*Treatment*block, group=block, data=treat.fungi)
perm.ITS.exp.field


write.table(perm.ITS.exp.field,file = here("output/tables", "perm.ITS.exp.treatment.txt"),eol = "\n",sep="\t")

perm.16S.exp.field=adonis2(bray.16S.exp.field~Cultivar*date*Treatment*block, group=block, data=treat.bact)
perm.16S.exp.field
write.table(perm.16S.exp.field,file = here("output/tables","perm.16S.exp.treatment.txt"),eol = "\n", sep="\t")


# Beta diversity analysis of Biolog data

Biolog.exp=read.table(file=here("data", "Biolog.Exp.txt"), row.names=1, header=T, comment.char = "", sep="\t")
Biolog.exp=na.omit(Biolog.exp)
treat.biolog=read.table(file=here("data","mapping_file.field.txt"),row.names=1, header=T, sep="\t", comment.char="")
treat.biolog=treat.biolog[-c(320),]

#Sort by row names
treat.b=treat.biolog[order(row.names(treat.biolog)),]
biolog.exp=Biolog.exp[order(row.names(Biolog.exp)),]
row.names(biolog.exp)==row.names(treat.b) #Sanity check

# Beta diversity analysis
bray.biolog.exp.field=vegdist(Biolog.exp,method = "bray")
pcoa.biolog.exp.field=cmdscale(sqrt(bray.biolog.exp.field),eig=T)

#PERMANOVA analysis
perm.biolog.exp.field=adonis2(bray.biolog.exp.field~Cultivar*date*Treatment*block, group=block, data=treat.biolog)
perm.biolog.exp.field
write.table(perm.biolog.exp.field,file = here("output/tables","perm.biolog.exp.txt"), eol = "\n", sep="\t")



# Preparation for data sorting and downstream analysis

#Sub-setting for grain qualities
metafile=read.table(file=here("data","metafile.txt"),row.names=1, header=T, sep="\t", comment.char="")
qual.data=read.table(file=here("data","model.quality.txt"),row.names = 1, header=T, sep="\t")

#Diversity index of bacteria

bact.div.field=read.table(file=here("data","bact.diversity.txt"), row.names=1,header=T, sep="\t")
fun.div.field=read.table(file=here("data","fungal.diversity.txt"), row.names=1,header=T, sep="\t")

#Carbon Utilization pattern
Biolog.exp=read.table(file=here("data","Biolog.Exp.txt"), row.names=1, header=T, comment.char = "", sep="\t")
Biolog.exp=(Biolog.exp/rowSums(Biolog.exp))#normalization


#Functional genes
bact.AOB=read.table(file=here("data","AOB.txt"),row.names = 1, header=T, sep="\t")
arch.AOA=read.table(file=here("data","AOA.field.txt"),row.names = 1, header=T, sep="\t")
nirk.bact=read.table(file=here("data","nirk.txt"),row.names = 1, header=T, sep="\t")
bact.nosz=read.table(file=here("data","nosZ.field.txt"),row.names = 1, header=T, sep="\t")
AOA.AOB=read.table(file=here("data","AOA.AOB.txt"),row.names=1, header=T, sep="\t", comment.char="")
FB.ratio=read.table(file=here("data","F.B ratio.txt"),row.names=1, header=T, sep="\t", comment.char="")
FB.ratio=log(FB.ratio)
qpcr.field=cbind(arch.AOA,bact.AOB,nirk.bact, bact.nosz,FB.ratio)

#Transformation of data
#qpcr.field$AOA =log10(qpcr.field$AOA)
#qpcr.field$AOB =log10(qpcr.field$AOB)
#qpcr.field$nirk=log10(qpcr.field$nirk)
#qpcr.field$nosZ=log10(qpcr.field$nosZ)
#qpcr.field$AOA.AOB=log(qpcr.field$AOA.AOB)


# Rearranging order
metafile=metafile[order(row.names(metafile)),]
qual.data=qual.data[order(row.names(qual.data)),]
OTU.16S.field=otu.16S.rels[order(row.names(otu.16S.rels)),]
OTU.ITS.field=otu.ITS.rels[order(row.names(otu.ITS.rels)),]
div.16S=bact.div.field[order(row.names(bact.div.field)),]
div.ITS=fun.div.field[order(row.names(bact.div.field)),]
biolog.field= Biolog.exp[order(row.names(Biolog.exp)),]
qpcr.field=qpcr.field[order(row.names(qpcr.field)),]

#Sorting out the OTUs with non-normalize data

otu.16S.non.normalize<-otu.16s.norm[order(row.names(otu.16s.norm)),]
otu.ITS.non.normalize<-otu.ITS.norm[order(row.names(otu.ITS.norm)),]



row.names(metafile)==row.names(qual.data)
row.names(metafile)==row.names(OTU.16S.field)
row.names(qual.data)==row.names(OTU.16S.field)
row.names(metafile)==row.names(OTU.ITS.field)
row.names(qual.data)==row.names(OTU.ITS.field)
row.names(div.16S)==row.names(qual.data)
row.names(div.ITS)==row.names(qual.data)
row.names(metafile)==row.names(div.16S)
row.names(metafile)==row.names(div.ITS)
row.names(metafile)==row.names(biolog.field)
row.names(metafile)==row.names(qpcr.field)
row.names(qual.data)==row.names(qpcr.field)
row.names(qual.data)==row.names(otu.16S.non.normalize)
row.names(qual.data)==row.names(otu.ITS.non.normalize)
row.names(metafile)==row.names(otu.16S.non.normalize)
row.names(metafile)==row.names(otu.ITS.non.normalize)

#Sub-setting 16S data according to cultivar
otu.16S.DT=OTU.16S.field[which(metafile$Cultivar=="Strongfield"),]
otu.16S.DS=OTU.16S.field[which(metafile$Cultivar=="AC Nass"),]


#Sub-setting ITS data according to cultivar

otu.ITS.DT=OTU.ITS.field[which(metafile$Cultivar=="Strongfield"),]
otu.ITS.DS=OTU.ITS.field[which(metafile$Cultivar=="AC Nass"),]


#Sub-setting biolog data according to cultivar

biolog.DT=biolog.field[which(metafile$Cultivar=="Strongfield"),]
biolog.DS=biolog.field[which(metafile$Cultivar=="AC Nass"),]

# Sub-setting qPCR data according to cultivar
qPCR.DT=qpcr.field[which(metafile$Cultivar=="Strongfield"),]
qPCR.DS=qpcr.field[which(metafile$Cultivar=="AC Nass"),]


# Sub-setting diversity data according to cultivar

div.16S.DT=div.16S[which(metafile$Cultivar=="Strongfield"),]
div.16S.DS=div.16S[which(metafile$Cultivar=="AC Nass"),]

div.ITS.DT=div.ITS[which(metafile$Cultivar=="Strongfield"),]
div.ITS.DS=div.ITS[which(metafile$Cultivar=="AC Nass"),]

#Sub-setting quality data according to the cultivar
quality.timepoint=cbind(metafile,qual.data)

strongfield.qual=quality.timepoint[which(metafile$Cultivar=="Strongfield"),]
AC.Nass.qual=quality.timepoint[which(metafile$Cultivar=="AC Nass"),]

# Sub-setting meta-file according to the Cultivar
metafile.DT=metafile[which(metafile$Cultivar=="Strongfield"),]
metafile.DS=metafile[which(metafile$Cultivar=="AC Nass"),]

metafile.DT=metafile.DT[order(row.names(metafile.DT)),]
metafile.DS=metafile.DS[order(row.names(metafile.DS)),]

#Sub-setting meta file data according to the dates
metafile.DT.T1=metafile.DT[which(metafile.DT$date=="10_May"),]
metafile.DT.T2=metafile.DT[which(metafile.DT$date=="24_May"),]
metafile.DT.T3=metafile.DT[which(metafile.DT$date=="7_June"),]
metafile.DT.T4=metafile.DT[which(metafile.DT$date=="21_June"),]
metafile.DT.T5=metafile.DT[which(metafile.DT$date=="5_July"),]
metafile.DT.T6=metafile.DT[which(metafile.DT$date=="19_July"),]
metafile.DT.T7=metafile.DT[which(metafile.DT$date=="1_August"),]

##Sub-setting meta file data according to the dates
metafile.DS.T1=metafile.DS[which(metafile.DS$date=="10_May"),]
metafile.DS.T2=metafile.DS[which(metafile.DS$date=="24_May"),]
metafile.DS.T3=metafile.DS[which(metafile.DS$date=="7_June"),]
metafile.DS.T4=metafile.DS[which(metafile.DS$date=="21_June"),]
metafile.DS.T5=metafile.DS[which(metafile.DS$date=="5_July"),]
metafile.DS.T6=metafile.DS[which(metafile.DS$date=="19_July"),]
metafile.DS.T7=metafile.DS[which(metafile.DS$date=="1_August"),]

## Sub-setting quality according to cultivar
qual.DT.T1=strongfield.qual[which(strongfield.qual$date=="10_May"),]
qual.DT.T2=strongfield.qual[which(strongfield.qual$date=="24_May"),]
qual.DT.T3=strongfield.qual[which(strongfield.qual$date=="7_June"),]
qual.DT.T4=strongfield.qual[which(strongfield.qual$date=="21_June"),]
qual.DT.T5=strongfield.qual[which(strongfield.qual$date=="5_July"),]
qual.DT.T6=strongfield.qual[which(strongfield.qual$date=="19_July"),]
qual.DT.T7=strongfield.qual[which(strongfield.qual$date=="1_August"),]


qual.DS.T1=AC.Nass.qual[which(AC.Nass.qual$date=="10_May"),]
qual.DS.T2=AC.Nass.qual[which(AC.Nass.qual$date=="24_May"),]
qual.DS.T3=AC.Nass.qual[which(AC.Nass.qual$date=="7_June"),]
qual.DS.T4=AC.Nass.qual[which(AC.Nass.qual$date=="21_June"),]
qual.DS.T5=AC.Nass.qual[which(AC.Nass.qual$date=="5_July"),]
qual.DS.T6=AC.Nass.qual[which(AC.Nass.qual$date=="19_July"),]
qual.DS.T7=AC.Nass.qual[which(AC.Nass.qual$date=="1_August"),]

# Sorting quality data from the meta file
qual.DT.S1=qual.DT.T1[,which(names(qual.DT.T1) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
qual.DT.S2=qual.DT.T2[,which(names(qual.DT.T2) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
qual.DT.S3=qual.DT.T3[,which(names(qual.DT.T3) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
qual.DT.S4=qual.DT.T4[,which(names(qual.DT.T4) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
qual.DT.S5=qual.DT.T5[,which(names(qual.DT.T5) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
qual.DT.S6=qual.DT.T6[,which(names(qual.DT.T6) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
qual.DT.S7=qual.DT.T7[,which(names(qual.DT.T7) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]

## Sorting quality data from the meta file
qual.DS.S1=qual.DS.T1[,which(names(qual.DS.T1) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
qual.DS.S2=qual.DS.T2[,which(names(qual.DS.T2) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
qual.DS.S3=qual.DS.T3[,which(names(qual.DS.T3) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
qual.DS.S4=qual.DS.T4[,which(names(qual.DS.T4) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
qual.DS.S5=qual.DS.T5[,which(names(qual.DS.T5) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
qual.DS.S6=qual.DS.T6[,which(names(qual.DS.T6) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
qual.DS.S7=qual.DS.T7[,which(names(qual.DS.T7) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]


## Sorting 16S ASVs according TO DT cultivar
metafile.DT=metafile.DT[order(row.names(metafile.DT)),]
row.names(metafile.DT)==row.names(otu.16S.DT)
row.names(metafile.DT)==row.names(otu.ITS.DT)


otu.16S.DT.T1=otu.16S.DT[which(metafile.DT$date=="10_May"),]
otu.16S.DT.T1=otu.16S.DT.T1[,!colSums(otu.16S.DT.T1)==0]
otu.16S.DT.T2=otu.16S.DT[which(metafile.DT$date=="24_May"),]
otu.16S.DT.T2=otu.16S.DT.T2[,!colSums(otu.16S.DT.T2)==0]
otu.16S.DT.T3=otu.16S.DT[which(metafile.DT$date=="7_June"),]
otu.16S.DT.T3=otu.16S.DT.T3[,!colSums(otu.16S.DT.T3)==0]
otu.16S.DT.T4=otu.16S.DT[which(metafile.DT$date=="21_June"),]
otu.16S.DT.T4=otu.16S.DT.T4[,!colSums(otu.16S.DT.T4)==0]
otu.16S.DT.T5=otu.16S.DT[which(metafile.DT$date=="5_July"),]
otu.16S.DT.T5=otu.16S.DT.T5[,!colSums(otu.16S.DT.T5)==0]
otu.16S.DT.T6=otu.16S.DT[which(metafile.DT$date=="19_July"),]
otu.16S.DT.T6=otu.16S.DT.T6[,!colSums(otu.16S.DT.T6)==0]
otu.16S.DT.T7=otu.16S.DT[which(metafile.DT$date=="1_August"),]
otu.16S.DT.T7=otu.16S.DT.T7[,!colSums(otu.16S.DT.T7)==0]

otu.ITS.DT.T1=otu.ITS.DT[which(metafile.DT$date=="10_May"),]
otu.ITS.DT.T1=otu.ITS.DT.T1[,!colSums(otu.ITS.DT.T1)==0]
otu.ITS.DT.T2=otu.ITS.DT[which(metafile.DT$date=="24_May"),]
otu.ITS.DT.T2=otu.ITS.DT.T2[,!colSums(otu.ITS.DT.T2)==0]
otu.ITS.DT.T3=otu.ITS.DT[which(metafile.DT$date=="7_June"),]
otu.ITS.DT.T3=otu.ITS.DT.T3[,!colSums(otu.ITS.DT.T3)==0]
otu.ITS.DT.T4=otu.ITS.DT[which(metafile.DT$date=="21_June"),]
otu.ITS.DT.T4=otu.ITS.DT.T4[,!colSums(otu.ITS.DT.T4)==0]
otu.ITS.DT.T5=otu.ITS.DT[which(metafile.DT$date=="5_July"),]
otu.ITS.DT.T5=otu.ITS.DT.T5[,!colSums(otu.ITS.DT.T5)==0]
otu.ITS.DT.T6=otu.ITS.DT[which(metafile.DT$date=="19_July"),]
otu.ITS.DT.T6=otu.ITS.DT.T6[,!colSums(otu.ITS.DT.T6)==0]
otu.ITS.DT.T7=otu.ITS.DT[which(metafile.DT$date=="1_August"),]
otu.ITS.DT.T7=otu.ITS.DT.T7[,!colSums(otu.ITS.DT.T7)==0]


## Sorting 16S diversity according TO DT cultivar
row.names(metafile.DT)==row.names(div.16S.DT)
row.names(metafile.DS)==row.names(div.16S.DS)

div.16S.DT.T1=div.16S.DT[which(metafile.DT$date=="10_May"),]
div.16S.DT.T2=div.16S.DT[which(metafile.DT$date=="24_May"),]
div.16S.DT.T3=div.16S.DT[which(metafile.DT$date=="7_June"),]
div.16S.DT.T4=div.16S.DT[which(metafile.DT$date=="21_June"),]
div.16S.DT.T5=div.16S.DT[which(metafile.DT$date=="5_July"),]
div.16S.DT.T6=div.16S.DT[which(metafile.DT$date=="19_July"),]
div.16S.DT.T7=div.16S.DT[which(metafile.DT$date=="1_August"),]


div.ITS.DT.T1=div.ITS.DT[which(metafile.DT$date=="10_May"),]
div.ITS.DT.T2=div.ITS.DT[which(metafile.DT$date=="24_May"),]
div.ITS.DT.T3=div.ITS.DT[which(metafile.DT$date=="7_June"),]
div.ITS.DT.T4=div.ITS.DT[which(metafile.DT$date=="21_June"),]
div.ITS.DT.T5=div.ITS.DT[which(metafile.DT$date=="5_July"),]
div.ITS.DT.T6=div.ITS.DT[which(metafile.DT$date=="19_July"),]
div.ITS.DT.T7=div.ITS.DT[which(metafile.DT$date=="1_August"),]

colnames(div.ITS.DT.T1)[which(colnames(div.ITS.DT.T1) %in% c("Shannon", "Simpson","OTU.richness","Chao1","ACE","PD"))] <- c("Shannon.fun", "Simpson.fun","OTU.richness.fun","Chao1.fun","ACE.fun","PD.fun")
colnames(div.ITS.DT.T2)[which(colnames(div.ITS.DT.T2) %in% c("Shannon", "Simpson","OTU.richness","Chao1","ACE","PD"))] <- c("Shannon.fun", "Simpson.fun","OTU.richness.fun","Chao1.fun","ACE.fun","PD.fun")
colnames(div.ITS.DT.T3)[which(colnames(div.ITS.DT.T3) %in% c("Shannon", "Simpson","OTU.richness","Chao1","ACE","PD"))] <- c("Shannon.fun", "Simpson.fun","OTU.richness.fun","Chao1.fun","ACE.fun","PD.fun")
colnames(div.ITS.DT.T4)[which(colnames(div.ITS.DT.T4) %in% c("Shannon", "Simpson","OTU.richness","Chao1","ACE","PD"))] <- c("Shannon.fun", "Simpson.fun","OTU.richness.fun","Chao1.fun","ACE.fun","PD.fun")
colnames(div.ITS.DT.T5)[which(colnames(div.ITS.DT.T5) %in% c("Shannon", "Simpson","OTU.richness","Chao1","ACE","PD"))] <- c("Shannon.fun", "Simpson.fun","OTU.richness.fun","Chao1.fun","ACE.fun","PD.fun")
colnames(div.ITS.DT.T6)[which(colnames(div.ITS.DT.T6) %in% c("Shannon", "Simpson","OTU.richness","Chao1","ACE","PD"))] <- c("Shannon.fun", "Simpson.fun","OTU.richness.fun","Chao1.fun","ACE.fun","PD.fun")
colnames(div.ITS.DT.T7)[which(colnames(div.ITS.DT.T7) %in% c("Shannon", "Simpson","OTU.richness","Chao1","ACE","PD"))] <- c("Shannon.fun", "Simpson.fun","OTU.richness.fun","Chao1.fun","ACE.fun","PD.fun")

## Sorting microbial function according TO DT cultivar
row.names(metafile.DT)==row.names(biolog.DT)
row.names(metafile.DS)==row.names(biolog.DS)

biolog.DT.T1=biolog.DT[which(metafile.DT$date=="10_May"),]
biolog.DT.T2=biolog.DT[which(metafile.DT$date=="24_May"),]
biolog.DT.T3=biolog.DT[which(metafile.DT$date=="7_June"),]
biolog.DT.T4=biolog.DT[which(metafile.DT$date=="21_June"),]
biolog.DT.T5=biolog.DT[which(metafile.DT$date=="5_July"),]
biolog.DT.T6=biolog.DT[which(metafile.DT$date=="19_July"),]
biolog.DT.T7=biolog.DT[which(metafile.DT$date=="1_August"),]

biolog.DS.T1=biolog.DS[which(metafile.DS$date=="10_May"),]
biolog.DS.T2=biolog.DS[which(metafile.DS$date=="24_May"),]
biolog.DS.T3=biolog.DS[which(metafile.DS$date=="7_June"),]
biolog.DS.T4=biolog.DS[which(metafile.DS$date=="21_June"),]
biolog.DS.T5=biolog.DS[which(metafile.DS$date=="5_July"),]
biolog.DS.T6=biolog.DS[which(metafile.DS$date=="19_July"),]
biolog.DS.T7=biolog.DS[which(metafile.DS$date=="1_August"),]

## Sorting microbial function according TO DT cultivar
row.names(metafile.DT)==row.names(qPCR.DT)
row.names(metafile.DS)==row.names(qPCR.DS)

##Sorting microbial function genes according TO DT cultivar
qpcr.T1.DT=qPCR.DT[metafile.DT$date=="10_May",]
qpcr.T2.DT=qPCR.DT[metafile.DT$date=="24_May",]
qpcr.T3.DT=qPCR.DT[metafile.DT$date=="7_June",]
qpcr.T4.DT=qPCR.DT[metafile.DT$date=="21_June",]
qpcr.T5.DT=qPCR.DT[metafile.DT$date=="5_July",]
qpcr.T6.DT=qPCR.DT[metafile.DT$date=="19_July",]
qpcr.T7.DT=qPCR.DT[metafile.DT$date=="1_August",]

##Sorting microbial function genes according TO DT cultivar
qpcr.T1.DS=qPCR.DS[metafile.DS$date=="10_May",]
qpcr.T2.DS=qPCR.DS[metafile.DS$date=="24_May",]
qpcr.T3.DS=qPCR.DS[metafile.DS$date=="7_June",]
qpcr.T4.DS=qPCR.DS[metafile.DS$date=="21_June",]
qpcr.T5.DS=qPCR.DS[metafile.DS$date=="5_July",]
qpcr.T6.DS=qPCR.DS[metafile.DS$date=="19_July",]
qpcr.T7.DS=qPCR.DS[metafile.DS$date=="1_August",]


## Sorting 16S ASVs according TO DS cultivar
row.names(metafile.DS)==row.names(otu.16S.DS)
row.names(metafile.DS)==row.names(otu.ITS.DS)


otu.16S.DS.T1=otu.16S.DS[which(metafile.DS$date=="10_May"),]
otu.16S.DS.T1=otu.16S.DS.T1[,!colSums(otu.16S.DS.T1)==0]
otu.16S.DS.T2=otu.16S.DS[which(metafile.DS$date=="24_May"),]
otu.16S.DS.T2=otu.16S.DS.T2[,!colSums(otu.16S.DS.T2)==0]
otu.16S.DS.T3=otu.16S.DS[which(metafile.DS$date=="7_June"),]
otu.16S.DS.T3=otu.16S.DS.T3[,!colSums(otu.16S.DS.T3)==0]
otu.16S.DS.T4=otu.16S.DS[which(metafile.DS$date=="21_June"),]
otu.16S.DS.T4=otu.16S.DS.T4[,!colSums(otu.16S.DS.T4)==0]
otu.16S.DS.T5=otu.16S.DS[which(metafile.DS$date=="5_July"),]
otu.16S.DS.T5=otu.16S.DS.T5[,!colSums(otu.16S.DS.T5)==0]
otu.16S.DS.T6=otu.16S.DS[which(metafile.DS$date=="19_July"),]
otu.16S.DS.T6=otu.16S.DS.T6[,!colSums(otu.16S.DS.T6)==0]
otu.16S.DS.T7=otu.16S.DS[which(metafile.DS$date=="1_August"),]
otu.16S.DS.T7=otu.16S.DS.T7[,!colSums(otu.16S.DS.T7)==0]


otu.ITS.DS.T1=otu.ITS.DS[which(metafile.DS$date=="10_May"),]
otu.ITS.DS.T1=otu.ITS.DS.T1[,!colSums(otu.ITS.DS.T1)==0]
otu.ITS.DS.T2=otu.ITS.DS[which(metafile.DS$date=="24_May"),]
otu.ITS.DS.T2=otu.ITS.DS.T2[,!colSums(otu.ITS.DS.T2)==0]
otu.ITS.DS.T3=otu.ITS.DS[which(metafile.DS$date=="7_June"),]
otu.ITS.DS.T3=otu.ITS.DS.T3[,!colSums(otu.ITS.DS.T3)==0]
otu.ITS.DS.T4=otu.ITS.DS[which(metafile.DS$date=="21_June"),]
otu.ITS.DS.T4=otu.ITS.DS.T4[,!colSums(otu.ITS.DS.T4)==0]
otu.ITS.DS.T5=otu.ITS.DS[which(metafile.DS$date=="5_July"),]
otu.ITS.DS.T5=otu.ITS.DS.T5[,!colSums(otu.ITS.DS.T5)==0]
otu.ITS.DS.T6=otu.ITS.DS[which(metafile.DS$date=="19_July"),]
otu.ITS.DS.T6=otu.ITS.DS.T6[,!colSums(otu.ITS.DS.T6)==0]
otu.ITS.DS.T7=otu.ITS.DS[which(metafile.DS$date=="1_August"),]
otu.ITS.DS.T7=otu.ITS.DS.T7[,!colSums(otu.ITS.DS.T7)==0]


## Sorting 16S diversity according TO DT cultivar

div.16S.DS.T1=div.16S.DS[which(metafile.DS$date=="10_May"),]
div.16S.DS.T2=div.16S.DS[which(metafile.DS$date=="24_May"),]
div.16S.DS.T3=div.16S.DS[which(metafile.DS$date=="7_June"),]
div.16S.DS.T4=div.16S.DS[which(metafile.DS$date=="21_June"),]
div.16S.DS.T5=div.16S.DS[which(metafile.DS$date=="5_July"),]
div.16S.DS.T6=div.16S.DS[which(metafile.DS$date=="19_July"),]
div.16S.DS.T7=div.16S.DS[which(metafile.DS$date=="1_August"),]


div.ITS.DS.T1=div.ITS.DS[which(metafile.DS$date=="10_May"),]
div.ITS.DS.T2=div.ITS.DS[which(metafile.DS$date=="24_May"),]
div.ITS.DS.T3=div.ITS.DS[which(metafile.DS$date=="7_June"),]
div.ITS.DS.T4=div.ITS.DS[which(metafile.DS$date=="21_June"),]
div.ITS.DS.T5=div.ITS.DS[which(metafile.DS$date=="5_July"),]
div.ITS.DS.T6=div.ITS.DS[which(metafile.DS$date=="19_July"),]
div.ITS.DS.T7=div.ITS.DS[which(metafile.DS$date=="1_August"),]

colnames(div.ITS.DS.T1)[which(colnames(div.ITS.DS.T1) %in% c("Shannon", "Simpson","OTU.richness","Chao1","ACE","PD"))] <- c("Shannon.fun", "Simpson.fun","OTU.richness.fun","Chao1.fun","ACE.fun","PD.fun")
colnames(div.ITS.DS.T2)[which(colnames(div.ITS.DS.T2) %in% c("Shannon", "Simpson","OTU.richness","Chao1","ACE","PD"))] <- c("Shannon.fun", "Simpson.fun","OTU.richness.fun","Chao1.fun","ACE.fun","PD.fun")
colnames(div.ITS.DS.T3)[which(colnames(div.ITS.DS.T3) %in% c("Shannon", "Simpson","OTU.richness","Chao1","ACE","PD"))] <- c("Shannon.fun", "Simpson.fun","OTU.richness.fun","Chao1.fun","ACE.fun","PD.fun")
colnames(div.ITS.DS.T4)[which(colnames(div.ITS.DS.T4) %in% c("Shannon", "Simpson","OTU.richness","Chao1","ACE","PD"))] <- c("Shannon.fun", "Simpson.fun","OTU.richness.fun","Chao1.fun","ACE.fun","PD.fun")
colnames(div.ITS.DS.T5)[which(colnames(div.ITS.DS.T5) %in% c("Shannon", "Simpson","OTU.richness","Chao1","ACE","PD"))] <- c("Shannon.fun", "Simpson.fun","OTU.richness.fun","Chao1.fun","ACE.fun","PD.fun")
colnames(div.ITS.DS.T6)[which(colnames(div.ITS.DS.T6) %in% c("Shannon", "Simpson","OTU.richness","Chao1","ACE","PD"))] <- c("Shannon.fun", "Simpson.fun","OTU.richness.fun","Chao1.fun","ACE.fun","PD.fun")
colnames(div.ITS.DS.T7)[which(colnames(div.ITS.DS.T7) %in% c("Shannon", "Simpson","OTU.richness","Chao1","ACE","PD"))] <- c("Shannon.fun", "Simpson.fun","OTU.richness.fun","Chao1.fun","ACE.fun","PD.fun")


# Now we will go for the  correlation between individual time point and the cultivar

# for make sure we will again check the sanity with quality data

row.names(qual.DT.S1)==row.names(otu.16S.DT.T1)
row.names(qual.DT.S2)==row.names(otu.16S.DT.T2)
row.names(qual.DT.S3)==row.names(otu.16S.DT.T3)
row.names(qual.DT.S4)==row.names(otu.16S.DT.T4)
row.names(qual.DT.S5)==row.names(otu.16S.DT.T5)
row.names(qual.DT.S6)==row.names(otu.16S.DT.T6)
row.names(qual.DT.S7)==row.names(otu.16S.DT.T7)

row.names(qual.DT.S1)==row.names(otu.ITS.DT.T1)
row.names(qual.DT.S2)==row.names(otu.ITS.DT.T2)
row.names(qual.DT.S3)==row.names(otu.ITS.DT.T3)
row.names(qual.DT.S4)==row.names(otu.ITS.DT.T4)
row.names(qual.DT.S5)==row.names(otu.ITS.DT.T5)
row.names(qual.DT.S6)==row.names(otu.ITS.DT.T6)
row.names(qual.DT.S7)==row.names(otu.ITS.DT.T7)

#Sanity for diversity for DT

row.names(qual.DT.S1)==row.names(div.16S.DT.T1)
row.names(qual.DT.S2)==row.names(div.16S.DT.T2)
row.names(qual.DT.S3)==row.names(div.16S.DT.T3)
row.names(qual.DT.S4)==row.names(div.16S.DT.T4)
row.names(qual.DT.S5)==row.names(div.16S.DT.T5)
row.names(qual.DT.S6)==row.names(div.16S.DT.T6)
row.names(qual.DT.S7)==row.names(div.16S.DT.T7)

row.names(qual.DT.S1)==row.names(div.ITS.DT.T1)
row.names(qual.DT.S2)==row.names(div.ITS.DT.T2)
row.names(qual.DT.S3)==row.names(div.ITS.DT.T3)
row.names(qual.DT.S4)==row.names(div.ITS.DT.T4)
row.names(qual.DT.S5)==row.names(div.ITS.DT.T5)
row.names(qual.DT.S6)==row.names(div.ITS.DT.T6)
row.names(qual.DT.S7)==row.names(div.ITS.DT.T7)



row.names(qual.DS.S1)==row.names(otu.16S.DS.T1)
row.names(qual.DS.S2)==row.names(otu.16S.DS.T2)
row.names(qual.DS.S3)==row.names(otu.16S.DS.T3)
row.names(qual.DS.S4)==row.names(otu.16S.DS.T4)
row.names(qual.DS.S5)==row.names(otu.16S.DS.T5)
row.names(qual.DS.S6)==row.names(otu.16S.DS.T6)
row.names(qual.DS.S7)==row.names(otu.16S.DS.T7)

row.names(qual.DS.S1)==row.names(otu.ITS.DS.T1)
row.names(qual.DS.S2)==row.names(otu.ITS.DS.T2)
row.names(qual.DS.S3)==row.names(otu.ITS.DS.T3)
row.names(qual.DS.S4)==row.names(otu.ITS.DS.T4)
row.names(qual.DS.S5)==row.names(otu.ITS.DS.T5)
row.names(qual.DS.S6)==row.names(otu.ITS.DS.T6)
row.names(qual.DS.S7)==row.names(otu.ITS.DS.T7)


#Sanity for diversity for DS

row.names(qual.DS.S1)==row.names(div.16S.DS.T1)
row.names(qual.DS.S2)==row.names(div.16S.DS.T2)
row.names(qual.DS.S3)==row.names(div.16S.DS.T3)
row.names(qual.DS.S4)==row.names(div.16S.DS.T4)
row.names(qual.DS.S5)==row.names(div.16S.DS.T5)
row.names(qual.DS.S6)==row.names(div.16S.DS.T6)
row.names(qual.DS.S7)==row.names(div.16S.DS.T7)

row.names(qual.DS.S1)==row.names(div.ITS.DS.T1)
row.names(qual.DS.S2)==row.names(div.ITS.DS.T2)
row.names(qual.DS.S3)==row.names(div.ITS.DS.T3)
row.names(qual.DS.S4)==row.names(div.ITS.DS.T4)
row.names(qual.DS.S5)==row.names(div.ITS.DS.T5)
row.names(qual.DS.S6)==row.names(div.ITS.DS.T6)
row.names(qual.DS.S7)==row.names(div.ITS.DS.T7)


# Cor.test between 16S ASVs vs quality=> Cultivar DT

i=1; unlink ("cor.DT.T1-16S-qual.pearson.txt");
while (i<=ncol(otu.16S.DT.T1)) {j=1; while(j<=ncol(qual.DT.S1)) 
{cor=cor.test(otu.16S.DT.T1[,i], qual.DT.S1[,j], method="kendall"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.16S.DT.T1)[i], " \t", 
colnames(qual.DT.S1)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(taxa.16S.S)[i,],"\n", file= here("output/tables","cor.DT.T1-16S-qual-pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DT.T2-16S-qual.pearson.txt");
while (i<=ncol(otu.16S.DT.T2)) {j=1; while(j<=ncol(qual.DT.S2)) 
{cor=cor.test(otu.16S.DT.T2[,i], qual.DT.S2[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.16S.DT.T2)[i], " \t", 
colnames(qual.DT.S2)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(tax.16s)[i,],"\n", file=here("output/tables","cor.DT.T2-16S-qual-pearson.txt"), append=T)};j=j+1} ;i=i+1}



i=1; unlink ("cor.DT.T3-16S-qual.pearson.txt");
while (i<=ncol(otu.16S.DT.T3)) {j=1; while(j<=ncol(qual.DT.S3)) 
{cor=cor.test(otu.16S.DT.T3[,i], qual.DT.S3[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.16S.DT.T3)[i], " \t", 
colnames(qual.DT.S3)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(tax.16s)[i,],"\n", file=here("output/tables","cor.DT.T3-16S-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}



i=1; unlink ("cor.DT.T4-16S-qual.pearson.txt");
while (i<=ncol(otu.16S.DT.T4)) {j=1; while(j<=ncol(qual.DT.S4)) 
{cor=cor.test(otu.16S.DT.T4[,i], qual.DT.S4[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.16S.DT.T4)[i], " \t", 
colnames(qual.DT.S4)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(tax.16s)[i,],"\n", file=here("output/tables","cor.DT.T4-16S-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DT.T5-16S-qual.pearson.txt");
while (i<=ncol(otu.16S.DT.T5)) {j=1; while(j<=ncol(qual.DT.S5)) 
{cor=cor.test(otu.16S.DT.T5[,i], qual.DT.S5[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.16S.DT.T5)[i], " \t", 
colnames(qual.DT.S5)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(tax.16s)[i,],"\n", file=here("output/tables","cor.DT.T5-16S-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DT.T6-16S-qual.pearson.txt");
while (i<=ncol(otu.16S.DT.T6)) {j=1; while(j<=ncol(qual.DT.S6)) 
{cor=cor.test(otu.16S.DT.T6[,i], qual.DT.S6[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.16S.DT.T6)[i], " \t", 
colnames(qual.DT.S6)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(tax.16s)[i,],"\n", file=here("output/tables","cor.DT.T6-16S-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DT.T7-16S-qual.pearson.txt");
while (i<=ncol(otu.16S.DT.T7)) {j=1; while(j<=ncol(qual.DT.S7)) 
{cor=cor.test(otu.16S.DT.T7[,i], qual.DT.S7[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.16S.DT.T7)[i], " \t", 
colnames(qual.DT.S7)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(tax.16s)[i,],"\n", file=here("output/tables","cor.DT.T7-16S-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}




# Cor.test between 16S OTUs vs quality=> Cultivar Ds
i=1; unlink ("cor.DS.T1-16S-qual.pearson.txt");
while (i<=ncol(otu.16S.DS.T1)) {j=1; while(j<=ncol(qual.DS.S1)) 
{cor=cor.test(otu.16S.DS.T1[,i], qual.DS.S1[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.16S.DS.T1)[i], " \t", 
colnames(qual.DT.S1)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(taxa.16S.S)[i,],"\n", file=here("output/tables","cor.DS.T1-16S-qual-pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DS.T2-16S-qual.pearson.txt");
while (i<=ncol(otu.16S.DS.T2)) {j=1; while(j<=ncol(qual.DS.S2)) 
{cor=cor.test(otu.16S.DS.T2[,i], qual.DS.S2[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.16S.DS.T2)[i], " \t", 
colnames(qual.DS.S2)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(taxa.16S.S)[i,],"\n", file=here("output/tables","cor.DS.T2-16S-qual-pearson.txt"), append=T)};j=j+1} ;i=i+1}



i=1; unlink ("cor.DS.T3-16S-qual.pearson.txt");
while (i<=ncol(otu.16S.DS.T3)) {j=1; while(j<=ncol(qual.DS.S3)) 
{cor=cor.test(otu.16S.DS.T3[,i], qual.DS.S3[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.16S.DS.T3)[i], " \t", 
colnames(qual.DS.S3)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(taxa.16S.S)[i,],"\n", file=here("output/tables","cor.DS.T3-16S-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}



i=1; unlink ("cor.DS.T4-16S-qual.pearson.txt");
while (i<=ncol(otu.16S.DS.T4)) {j=1; while(j<=ncol(qual.DS.S4)) 
{cor=cor.test(otu.16S.DS.T4[,i], qual.DS.S4[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.16S.DS.T4)[i], " \t", 
colnames(qual.DS.S4)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(taxa.16S.S)[i,],"\n", file=here("output/tables","cor.DS.T4-16S-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DS.T5-16S-qual.pearson.txt");
while (i<=ncol(otu.16S.DS.T5)) {j=1; while(j<=ncol(qual.DS.S5)) 
{cor=cor.test(otu.16S.DS.T5[,i], qual.DS.S5[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.16S.DS.T5)[i], " \t", 
colnames(qual.DS.S5)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(taxa.16S.S)[i,],"\n", file=here("output/tables","cor.DS.T5-16S-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DS.T6-16S-qual.pearson.txt"); 
while (i<=ncol(otu.16S.DS.T6)) {j=1; while(j<=ncol(qual.DS.S6)) 
{cor=cor.test(otu.16S.DS.T6[,i], qual.DS.S6[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.16S.DS.T6)[i], " \t", 
colnames(qual.DS.S6)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(taxa.16S.S)[i,],"\n", file=here("output/tables","cor.DS.T6-16S-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DS.T7-16S-qual.pearson.txt");
while (i<=ncol(otu.16S.DS.T7)) {j=1; while(j<=ncol(qual.DS.S7)) 
{cor=cor.test(otu.16S.DS.T7[,i], qual.DS.S7[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.16S.DS.T7)[i], " \t", 
colnames(qual.DS.S7)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(taxa.16S.S)[i,],"\n", file=here("output/tables","cor.DS.T7-16S-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}



# Cor.test between diversity vs quality=> Cultivar DT
i=1; unlink ("cor.DT.T1-16S-div-qual.pearson.txt");
while (i<=ncol(div.16S.DT.T1)) {j=1; while(j<=ncol(qual.DT.S1)) 
{cor=cor.test(div.16S.DT.T1[,i], qual.DT.S1[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.16S.DT.T1)[i], " \t", 
colnames(qual.DT.S1)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DT.T1-16S-div-qual-pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DT.T2-16S-div-qual.pearson.txt");
while (i<=ncol(div.16S.DT.T2)) {j=1; while(j<=ncol(qual.DT.S2)) 
{cor=cor.test(div.16S.DT.T2[,i], qual.DT.S2[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.1){cat (colnames(div.16S.DT.T2)[i], " \t", 
colnames(qual.DT.S2)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DT.T2-16S-div-qual-pearson.txt"), append=T)};j=j+1} ;i=i+1}



i=1; unlink ("cor.DT.T3-16S-div-qual.pearson.txt");
while (i<=ncol(div.16S.DT.T3)) {j=1; while(j<=ncol(qual.DT.S3)) 
{cor=cor.test(div.16S.DT.T3[,i], qual.DT.S3[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.1){cat (colnames(div.16S.DT.T3)[i], " \t", 
colnames(qual.DT.S3)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DT.T3-16S-div-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}



i=1; unlink ("cor.DT.T4-16S-div-qual.pearson.txt");
while (i<=ncol(div.16S.DT.T4)) {j=1; while(j<=ncol(qual.DT.S4)) 
{cor=cor.test(div.16S.DT.T4[,i], qual.DT.S4[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.16S.DT.T4)[i], " \t", 
colnames(qual.DT.S4)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DT.T4-16S-div-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DT.T5-16S-div-qual.pearson.txt");
while (i<=ncol(div.16S.DT.T5)) {j=1; while(j<=ncol(qual.DT.S5)) 
{cor=cor.test(div.16S.DT.T5[,i], qual.DT.S5[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.16S.DT.T5)[i], " \t", 
colnames(qual.DT.S5)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DT.T5-16S-div-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DT.T6-16S-qual.pearson.txt");
while (i<=ncol(div.16S.DT.T6)) {j=1; while(j<=ncol(qual.DT.S6)) 
{cor=cor.test(div.16S.DT.T6[,i], qual.DT.S6[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.16S.DT.T6)[i], " \t", 
colnames(qual.DT.S6)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file="cor.DT.T6-16S-div-qual.pearson.txt", append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DT.T7-16S-div-qual.pearson.txt");
while (i<=ncol(div.16S.DT.T7)) {j=1; while(j<=ncol(qual.DT.S7)) 
{cor=cor.test(div.16S.DT.T7[,i], qual.DT.S7[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.16S.DT.T7)[i], " \t", 
colnames(qual.DT.S7)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DT.T7-16S-div-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}


# Cor.test between diversity vs quality=> Cultivar DS
i=1; unlink ("cor.DS.T1-16S-div-qual.pearson.txt");
while (i<=ncol(div.16S.DS.T1)) {j=1; while(j<=ncol(qual.DS.S1)) 
{cor=cor.test(div.16S.DS.T1[,i], qual.DS.S1[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.16S.DS.T1)[i], " \t", 
colnames(qual.DS.S1)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DS.T1-16S-div-qual-pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DS.T2-16S-div-qual.pearson.txt");
while (i<=ncol(div.16S.DS.T2)) {j=1; while(j<=ncol(qual.DS.S2)) 
{cor=cor.test(div.16S.DS.T2[,i], qual.DS.S2[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.16S.DS.T2)[i], " \t", 
colnames(qual.DS.S2)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DS.T2-16S-div-qual-pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DS.T3-16S-div-qual.pearson.txt");
while (i<=ncol(div.16S.DS.T3)) {j=1; while(j<=ncol(qual.DS.S3)) 
{cor=cor.test(div.16S.DS.T3[,i], qual.DS.S3[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.16S.DS.T3)[i], " \t", 
colnames(qual.DS.S3)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DS.T3-16S-div-qual-pearson.txt"), append=T)};j=j+1} ;i=i+1}

i=1; unlink ("cor.DS.T4-16S-div-qual.pearson.txt");
while (i<=ncol(div.16S.DS.T4)) {j=1; while(j<=ncol(qual.DS.S4)) 
{cor=cor.test(div.16S.DS.T4[,i], qual.DS.S4[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.16S.DS.T4)[i], " \t", 
colnames(qual.DS.S4)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DS.T4-16S-div-qual-pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DS.T5-16S-div-qual.pearson.txt");
while (i<=ncol(div.16S.DS.T5)) {j=1; while(j<=ncol(qual.DS.S5)) 
{cor=cor.test(div.16S.DS.T5[,i], qual.DS.S5[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.16S.DS.T5)[i], " \t", 
colnames(qual.DS.S5)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DS.T5-16S-div-qual-pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DS.T6-16S-div-qual.pearson.txt");
while (i<=ncol(div.16S.DS.T6)) {j=1; while(j<=ncol(qual.DS.S6)) 
{cor=cor.test(div.16S.DS.T6[,i], qual.DS.S6[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.16S.DS.T6)[i], " \t", 
colnames(qual.DS.S6)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DS.T6-16S-div-qual-pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DS.T7-16S-div-qual.pearson.txt");
while (i<=ncol(div.16S.DS.T7)) {j=1; while(j<=ncol(qual.DS.S7)) 
{cor=cor.test(div.16S.DS.T7[,i], qual.DS.S7[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.16S.DS.T7)[i], " \t", 
colnames(qual.DS.S7)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DS.T7-16S-div-qual-pearson.txt"), append=T)};j=j+1} ;i=i+1}



# Cor.test between ITS ASVs vs quality=> Cultivar DT

i=1; unlink ("cor.DT.T1-ITS-div-qual.txt");
while (i<=ncol(div.ITS.DT.T1)) {j=1; while(j<=ncol(qual.DT.S1)) 
{cor=cor.test(div.ITS.DT.T1[,i], qual.DT.S1[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.16S.DT.T1)[i], " \t", 
colnames(qual.DT.S1)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DT.T1-ITS-div-qual.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DT.T2-ITS-div-qual.txt");
while (i<=ncol(div.ITS.DT.T2)) {j=1; while(j<=ncol(qual.DT.S2)) 
{cor=cor.test(div.ITS.DT.T2[,i], qual.DT.S2[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.ITS.DT.T2)[i], " \t", 
colnames(qual.DT.S2)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DT.T2-ITS-div-qual.txt"), append=T)};j=j+1} ;i=i+1}

i=1; unlink ("cor.DT.T3-ITS-div-qual.txt");
while (i<=ncol(div.ITS.DT.T3)) {j=1; while(j<=ncol(qual.DT.S3)) 
{cor=cor.test(div.ITS.DT.T3[,i], qual.DT.S3[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.ITS.DT.T3)[i], " \t", 
colnames(qual.DT.S3)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DT.T3-ITS-div-qual.txt"), append=T)};j=j+1} ;i=i+1}

i=1; unlink ("cor.DT.T4-ITS-div-qual.txt");
while (i<=ncol(div.ITS.DT.T4)) {j=1; while(j<=ncol(qual.DT.S4)) 
{cor=cor.test(div.ITS.DT.T4[,i], qual.DT.S4[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.ITS.DT.T4)[i], " \t", 
colnames(qual.DT.S4)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DT.T4-ITS-div-qual.txt"), append=T)};j=j+1} ;i=i+1}

i=1; unlink ("cor.DT.T5-ITS-div-qual.txt");
while (i<=ncol(div.ITS.DT.T5)) {j=1; while(j<=ncol(qual.DT.S5)) 
{cor=cor.test(div.ITS.DT.T5[,i], qual.DT.S5[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.ITS.DT.T5)[i], " \t", 
colnames(qual.DT.S5)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DT.T5-ITS-div-qual.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DT.T6.ITS-div-qual.txt");
while (i<=ncol(div.ITS.DT.T6)) {j=1; while(j<=ncol(qual.DT.S6)) 
{cor=cor.test(div.ITS.DT.T6[,i], qual.DT.S6[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.ITS.DT.T6)[i], " \t", 
colnames(qual.DT.S6)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DT.T6-ITS.div-qual.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DT.T7.ITS-div-qual.txt");
while (i<=ncol(div.ITS.DT.T7)) {j=1; while(j<=ncol(qual.DT.S7)) 
{cor=cor.test(div.ITS.DT.T7[,i], qual.DT.S7[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.ITS.DT.T7)[i], " \t", 
colnames(qual.DT.S7)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DT.T7.ITS-div-qual.txt"), append=T)};j=j+1} ;i=i+1}

# Cor.test between ITS Diversity vs quality=> Cultivar DS

i=1; unlink ("cor.DS.T1-ITS-div-qual.txt");
while (i<=ncol(div.ITS.DS.T1)) {j=1; while(j<=ncol(qual.DS.S1)) 
{cor=cor.test(div.ITS.DS.T1[,i], qual.DS.S1[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.16S.DS.T1)[i], " \t", 
colnames(qual.DS.S1)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DS.T1-ITS-div-qual.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DS.T2-ITS-div-qual.txt");
while (i<=ncol(div.ITS.DS.T2)) {j=1; while(j<=ncol(qual.DS.S2)) 
{cor=cor.test(div.ITS.DS.T2[,i], qual.DS.S2[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.ITS.DS.T2)[i], " \t", 
colnames(qual.DS.S2)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DS.T2.ITS-div-qual.txt"), append=T)};j=j+1} ;i=i+1}

i=1; unlink ("cor.DS.T3-ITS-div-qual.txt");
while (i<=ncol(div.ITS.DS.T3)) {j=1; while(j<=ncol(qual.DS.S3)) 
{cor=cor.test(div.ITS.DS.T3[,i], qual.DS.S3[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.ITS.DS.T3)[i], " \t", 
colnames(qual.DS.S3)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DS.T3-ITS-div-qual.txt"), append=T)};j=j+1} ;i=i+1}

i=1; unlink ("cor.DS.T4-ITS-div-qual.txt");
while (i<=ncol(div.ITS.DS.T4)) {j=1; while(j<=ncol(qual.DS.S4)) 
{cor=cor.test(div.ITS.DS.T4[,i], qual.DT.S4[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.ITS.DS.T4)[i], " \t", 
colnames(qual.DS.S4)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DS.T4-ITS-div-qual.txt"), append=T)};j=j+1} ;i=i+1}

i=1; unlink ("cor.DS.T5-ITS-div-qual.txt");
while (i<=ncol(div.ITS.DS.T5)) {j=1; while(j<=ncol(qual.DS.S5)) 
{cor=cor.test(div.ITS.DS.T5[,i], qual.DS.S5[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.ITS.DS.T5)[i], " \t", 
colnames(qual.DS.S5)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DS.T5-ITS-div-qual.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DS.T6-ITS-div-qual.txt");
while (i<=ncol(div.ITS.DS.T6)) {j=1; while(j<=ncol(qual.DS.S6)) 
{cor=cor.test(div.ITS.DS.T6[,i], qual.DS.S6[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.ITS.DS.T6)[i], " \t", 
colnames(qual.DS.S6)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DS.T6-ITS-div-qual.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DS.T7-ITS-div-qual.txt");
while (i<=ncol(div.ITS.DS.T7)) {j=1; while(j<=ncol(qual.DS.S7)) 
{cor=cor.test(div.ITS.DS.T7[,i], qual.DS.S7[,j], method="spearman"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(div.ITS.DS.T7)[i], " \t", 
colnames(qual.DS.S7)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DS.T7-ITS-div-qual.txt"), append=T)};j=j+1} ;i=i+1}




# Cor.test between ITS ASVs vs quality=> Cultivar DT

i=1; unlink ("cor.DT.T1-ITS-qual.pearson.txt");
while (i<=ncol(otu.ITS.DT.T1)) {j=1; while(j<=ncol(qual.DT.S1)) 
{cor=cor.test(otu.ITS.DT.T1[,i], qual.DT.S1[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.ITS.DT.T1)[i], " \t", 
colnames(qual.DT.S1)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(tax.ITS)[i,],"\n", file=here("output/tables","cor.DT.T1-ITS-qual-pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DT.T2-ITS-qual.pearson.txt");
while (i<=ncol(otu.ITS.DT.T2)) {j=1; while(j<=ncol(qual.DT.S2)) 
{cor=cor.test(otu.ITS.DT.T2[,i], qual.DT.S2[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.ITS.DT.T2)[i], " \t", 
colnames(qual.DT.S2)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(tax.ITS)[i,],"\n", file=here("output/tables","cor.DT.T2-ITS-qual-pearson.txt"), append=T)};j=j+1} ;i=i+1}



i=1; unlink ("cor.DT.T3-ITS-qual.pearson.txt");
while (i<=ncol(otu.ITS.DT.T3)) {j=1; while(j<=ncol(qual.DT.S3)) 
{cor=cor.test(otu.ITS.DT.T3[,i], qual.DT.S3[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.ITS.DT.T3)[i], " \t", 
colnames(qual.DT.S3)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(tax.ITS)[i,],"\n", file=here("output/tables","cor.DT.T3-ITS-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}



i=1; unlink ("cor.DT.T4-ITS-qual.pearson.txt");
while (i<=ncol(otu.ITS.DT.T4)) {j=1; while(j<=ncol(qual.DT.S4)) 
{cor=cor.test(otu.ITS.DT.T4[,i], qual.DT.S4[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.ITS.DT.T4)[i], " \t", 
colnames(qual.DT.S4)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(tax.ITS)[i,],"\n", file=here("output/tables","cor.DT.T4-ITS-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}




i=1; unlink ("cor.DT.T5-ITS-qual.pearson.txt");
while (i<=ncol(otu.ITS.DT.T5)) {j=1; while(j<=ncol(qual.DT.S5)) 
{cor=cor.test(otu.ITS.DT.T5[,i], qual.DT.S5[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.ITS.DT.T5)[i], " \t", 
colnames(qual.DT.S5)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(tax.ITS)[i,],"\n", file=here("output/tables","cor.DT.T5-ITS-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DT.T6-ITS-qual.pearson.txt");
while (i<=ncol(otu.ITS.DT.T6)) {j=1; while(j<=ncol(qual.DT.S6)) 
{cor=cor.test(otu.ITS.DT.T6[,i], qual.DT.S6[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.ITS.DT.T6)[i], " \t", 
colnames(qual.DT.S6)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(tax.ITS)[i,],"\n", file=here("output/tables","cor.DT.T6-ITS-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DT.T7-ITS-qual.pearson.txt");
while (i<=ncol(otu.ITS.DT.T7)) {j=1; while(j<=ncol(qual.DT.S7)) 
{cor=cor.test(otu.ITS.DT.T7[,i], qual.DT.S7[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.ITS.DT.T7)[i], " \t", 
colnames(qual.DT.S7)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(tax.ITS)[i,],"\n", file=here("output/tables","cor.DT.T7-ITS-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}




# Cor.test between 16S ASVs vs quality=> Cultivar DS

i=1; unlink ("cor.DS.T1-ITS-qual.pearson.txt");
while (i<=ncol(otu.ITS.DS.T1)) {j=1; while(j<=ncol(qual.DS.S1)) 
{cor=cor.test(otu.ITS.DS.T1[,i], qual.DS.S1[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.ITS.DS.T1)[i], " \t", 
colnames(qual.DS.S1)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(tax.ITS)[i,],"\n", file=here("output/tables","cor.DS.T1-ITS-qual-pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DS.T2-ITS-qual.pearson.txt");
while (i<=ncol(otu.ITS.DS.T2)) {j=1; while(j<=ncol(qual.DS.S2)) 
{cor=cor.test(otu.ITS.DS.T2[,i], qual.DS.S2[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.ITS.DS.T2)[i], " \t", 
colnames(qual.DS.S2)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(tax.ITS)[i,],"\n", file=here("output/tables","cor.DS.T2-ITS-qual-pearson.txt"), append=T)};j=j+1} ;i=i+1}



i=1; unlink ("cor.DS.T3-ITS-qual.pearson.txt");
while (i<=ncol(otu.ITS.DS.T3)) {j=1; while(j<=ncol(qual.DS.S3)) 
{cor=cor.test(otu.ITS.DS.T3[,i], qual.DS.S3[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.ITS.DS.T3)[i], " \t", 
colnames(qual.DS.S3)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(tax.ITS)[i,],"\n", file=here("output/tables","cor.DS.T3-ITS-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DS.T4-ITS-qual.pearson.txt");
while (i<=ncol(otu.ITS.DS.T4)) {j=1; while(j<=ncol(qual.DS.S4)) 
{cor=cor.test(otu.ITS.DS.T4[,i], qual.DS.S4[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.ITS.DS.T4)[i], " \t", 
colnames(qual.DS.S4)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(tax.ITS)[i,],"\n", file=here("output/tables","cor.DS.T4-ITS-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}

i=1; unlink ("cor.DS.T5-ITS-qual.pearson.txt");
while (i<=ncol(otu.ITS.DS.T5)) {j=1; while(j<=ncol(qual.DS.S5)) 
{cor=cor.test(otu.ITS.DS.T5[,i], qual.DS.S5[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.ITS.DS.T5)[i], " \t", 
colnames(qual.DS.S5)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(tax.ITS)[i,],"\n", file=here("output/tables","cor.DS.T5-ITS-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}

i=1; unlink ("cor.DS.T6-ITS-qual.pearson.txt");
while (i<=ncol(otu.ITS.DS.T6)) {j=1; while(j<=ncol(qual.DS.S6)) 
{cor=cor.test(otu.ITS.DS.T6[,i], qual.DS.S6[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.ITS.DS.T6)[i], " \t", 
colnames(qual.DS.S6)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(tax.ITS)[i,],"\n", file=here("output/tables","cor.DS.T6-ITS-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DS.T7-ITS-qual.pearson.txt");
while (i<=ncol(otu.ITS.DS.T7)) {j=1; while(j<=ncol(qual.DS.S7)) 
{cor=cor.test(otu.ITS.DS.T7[,i], qual.DS.S7[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.001){cat (colnames(otu.ITS.DS.T7)[i], " \t", 
colnames(qual.DS.S7)[j], "\t", cor$p.value, "\t", cor$estimate,"\t",as.matrix(tax.ITS)[i,],"\n", file=here("output/tables","cor.DS.T7-ITS-qual.pearson.txt"), append=T)};j=j+1} ;i=i+1}


#Sanity for bio-log data set=>DT

row.names(qual.DT.S1)==row.names(biolog.DT.T1)
row.names(qual.DT.S2)==row.names(biolog.DT.T2)
row.names(qual.DT.S3)==row.names(biolog.DT.T3)
row.names(qual.DT.S4)==row.names(biolog.DT.T4)
row.names(qual.DT.S5)==row.names(biolog.DT.T5)
row.names(qual.DT.S6)==row.names(biolog.DT.T6)
row.names(qual.DT.S7)==row.names(biolog.DT.T7)


#Sanity for bio-log data set=>DS

row.names(qual.DS.S1)==row.names(biolog.DS.T1)
row.names(qual.DS.S2)==row.names(biolog.DS.T2)
row.names(qual.DS.S3)==row.names(biolog.DS.T3)
row.names(qual.DS.S4)==row.names(biolog.DS.T4)
row.names(qual.DS.S5)==row.names(biolog.DS.T5)
row.names(qual.DS.S6)==row.names(biolog.DS.T6)
row.names(qual.DS.S7)==row.names(biolog.DS.T7)

# Cor.test between biolog vs quality=> Cultivar DT

i=1; unlink ("cor.DT.T1-biolog-qual.txt");
while (i<=ncol(biolog.DT.T1)) {j=1; while(j<=ncol(qual.DT.S1)) 
{cor=cor.test(biolog.DT.T1[,i], qual.DT.S2[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(biolog.DT.T1)[i], " \t", 
colnames(qual.DT.S1)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DT.T1-biolog-qual.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DT.T2-biolog-qual.txt");
while (i<=ncol(biolog.DT.T2)) {j=1; while(j<=ncol(qual.DT.S2)) 
{cor=cor.test(biolog.DT.T2[,i], qual.DT.S2[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(biolog.DT.T2)[i], " \t", 
colnames(qual.DT.S2)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DT.T2-biolog-qual.txt"), append=T)};j=j+1} ;i=i+1}



i=1; unlink ("cor.DT.T3-biolog-qual.txt");
while (i<=ncol(biolog.DT.T3)) {j=1; while(j<=ncol(qual.DT.S3)) 
{cor=cor.test(biolog.DT.T3[,i], qual.DT.S3[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(biolog.DT.T3)[i], "\t", 
colnames(qual.DT.S3)[j],"\t", cor$p.value,"\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DT.T3-biolog-qual.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DT.T4-biolog-qual.txt");
while (i<=ncol(biolog.DT.T4)) {j=1; while(j<=ncol(qual.DT.S4)) 
{cor=cor.test(biolog.DT.T4[,i], qual.DT.S4[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(biolog.DT.T4)[i], "\t", 
colnames(qual.DT.S4)[j],"\t", cor$p.value,"\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DT.T4-biolog-qual.txt"), append=T)};j=j+1} ;i=i+1}

i=1; unlink ("cor.DT.T5-biolog-qual.txt");
while (i<=ncol(biolog.DT.T5)) {j=1; while(j<=ncol(qual.DT.S5)) 
{cor=cor.test(biolog.DT.T5[,i],qual.DT.S5[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(biolog.DT.T5)[i], "\t", 
colnames(qual.DT.S5)[j],"\t", cor$p.value,"\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DT.T5-biolog-qual.txt"), append=T)};j=j+1} ;i=i+1}

i=1; unlink ("cor.DT.T6-biolog-qual.txt");
while (i<=ncol(biolog.DT.T6)) {j=1; while(j<=ncol(qual.DT.S6)) 
{cor=cor.test(biolog.DT.T6[,i], qual.DT.S6[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(biolog.DT.T6)[i], "\t", 
colnames(qual.DT.S6)[j],"\t", cor$p.value,"\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DT.T6-biolog-qual.txt"), append=T)};j=j+1} ;i=i+1}

i=1; unlink ("cor.DT.T7-biolog-qual.txt");
while (i<=ncol(biolog.DT.T7)) {j=1; while(j<=ncol(qual.DT.S7)) 
{cor=cor.test(biolog.DT.T7[,i], qual.DT.S7[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(biolog.DT.T7)[i], "\t", 
colnames(qual.DT.S7)[j],"\t", cor$p.value,"\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DT.T7-biolog-qual.txt"), append=T)};j=j+1} ;i=i+1}

# Cor.test between biolog vs quality=> Cultivar DS

i=1; unlink ("cor.DS.T1-biolog-qual.txt");
while (i<=ncol(biolog.DS.T1)) {j=1; while(j<=ncol(qual.DS.S1)) 
{cor=cor.test(biolog.DS.T1[,i], qual.DS.S2[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(biolog.DS.T1)[i], " \t", 
colnames(qual.DS.S1)[j],"\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DS.T1-biolog-qual.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DS.T2-biolog-qual.txt");
while (i<=ncol(biolog.DS.T2)) {j=1; while(j<=ncol(qual.DS.S2)) 
{cor=cor.test(biolog.DS.T2[,i], qual.DS.S2[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(biolog.DS.T2)[i], " \t", 
colnames(qual.DS.S2)[j], "\t", cor$p.value, "\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DS.T2-biolog-qual.txt"), append=T)};j=j+1} ;i=i+1}



i=1; unlink ("cor.DS.T3-biolog-qual.txt");
while (i<=ncol(biolog.DS.T3)) {j=1; while(j<=ncol(qual.DS.S3)) 
{cor=cor.test(biolog.DS.T3[,i], qual.DS.S3[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(biolog.DS.T3)[i], "\t", 
colnames(qual.DS.S3)[j],"\t", cor$p.value,"\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DS.T3-biolog-qual.txt"), append=T)};j=j+1} ;i=i+1}


i=1; unlink ("cor.DS.T4-biolog-qual.txt");
while (i<=ncol(biolog.DS.T4)) {j=1; while(j<=ncol(qual.DS.S4)) 
{cor=cor.test(biolog.DS.T4[,i], qual.DS.S4[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(biolog.DS.T4)[i], "\t", 
colnames(qual.DS.S4)[j],"\t", cor$p.value,"\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DS.T4-biolog-qual.txt"), append=T)};j=j+1} ;i=i+1}

i=1; unlink ("cor.DS.T5-biolog-qual.txt");
while (i<=ncol(biolog.DS.T5)) {j=1; while(j<=ncol(qual.DS.S5)) 
{cor=cor.test(biolog.DS.T5[,i], qual.DS.S5[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(biolog.DS.T5)[i], "\t", 
colnames(qual.DS.S5)[j],"\t", cor$p.value,"\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DS.T5-biolog-qual.txt"), append=T)};j=j+1} ;i=i+1}

i=1; unlink ("cor.DS.T6-biolog-qual.txt");
while (i<=ncol(biolog.DS.T6)) {j=1; while(j<=ncol(qual.DS.S6)) 
{cor=cor.test(biolog.DS.T6[,i], qual.DT.S6[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(biolog.DS.T6)[i], "\t", 
colnames(qual.DS.S6)[j],"\t", cor$p.value,"\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DS.T6-biolog-qual.txt"), append=T)};j=j+1} ;i=i+1}

i=1; unlink ("cor.DS.T7-biolog-qual.txt");
while (i<=ncol(biolog.DS.T7)) {j=1; while(j<=ncol(qual.DS.S7)) 
{cor=cor.test(biolog.DS.T7[,i], qual.DS.S7[,j], method="pearson"); p.value=cor$p.value; 
if(is.na(p.value)){p.value=1}; if(p.value<0.05){cat (colnames(biolog.DS.T7)[i], "\t", 
colnames(qual.DS.S7)[j],"\t", cor$p.value,"\t", cor$estimate,"\t","\n", file=here("output/tables","cor.DS.T7-biolog-qual.txt"), append=T)};j=j+1} ;i=i+1}


# Bray-Curtis dissimilarity index for drought tolerant

library(vegan)
bray.16S.DT.T1=vegdist(otu.16S.DT.T1,method = "bray")
pcoa.16S.DT.T1=cmdscale(sqrt(bray.16S.DT.T1),eig=T)
bray.16S.DT.T2=vegdist(otu.16S.DT.T2,method = "bray")
pcoa.16S.DT.T2=cmdscale(sqrt(bray.16S.DT.T2),eig=T)
bray.16S.DT.T3=vegdist(otu.16S.DT.T3,method = "bray")
pcoa.16S.DT.T3=cmdscale(sqrt(bray.16S.DT.T3),eig=T)
bray.16S.DT.T4=vegdist(otu.16S.DT.T4,method = "bray")
pcoa.16S.DT.T4=cmdscale(sqrt(bray.16S.DT.T4),eig=T)
bray.16S.DT.T5=vegdist(otu.16S.DT.T5,method = "bray")
pcoa.16S.DT.T5=cmdscale(sqrt(bray.16S.DT.T5),eig=T)
bray.16S.DT.T6=vegdist(otu.16S.DT.T6,method = "bray")
pcoa.16S.DT.T6=cmdscale(sqrt(bray.16S.DT.T6),eig=T)
bray.16S.DT.T7=vegdist(otu.16S.DT.T7,method = "bray")
pcoa.16S.DT.T7=cmdscale(sqrt(bray.16S.DT.T7),eig=T)


# PCoA point as data.frame for DT
pcoa.point.DT.T1=data.frame(pcoa.16S.DT.T1$points)#Axis1 8.23% #Axis2 5.44%
pcoa.point.DT.T2=data.frame(pcoa.16S.DT.T2$points)#Axis1 9.10% #Axis2 5.95%
pcoa.point.DT.T3=data.frame(pcoa.16S.DT.T3$points)#Axis1 8.19% #Axis2 5.81%
pcoa.point.DT.T4=data.frame(pcoa.16S.DT.T4$points)#Axis1 7.49% #Axis2 5.52%
pcoa.point.DT.T5=data.frame(pcoa.16S.DT.T5$points)#Axis1 7.68% #Axis2 5.83%
pcoa.point.DT.T6=data.frame(pcoa.16S.DT.T6$points)#Axis1 16.51% #Axis2 8.67%
pcoa.point.DT.T7=data.frame(pcoa.16S.DT.T7$points)#Axis1 16.81% #Axis2 7.18%

colnames(pcoa.point.DT.T1)[which(colnames(pcoa.point.DT.T1) %in% c("X1", "X2"))] <- c("bact.axis1", "bact.axis2")
colnames(pcoa.point.DT.T2)[which(colnames(pcoa.point.DT.T2) %in% c("X1", "X2"))] <- c("bact.axis1", "bact.axis2")
colnames(pcoa.point.DT.T3)[which(colnames(pcoa.point.DT.T3) %in% c("X1", "X2"))] <- c("bact.axis1", "bact.axis2")
colnames(pcoa.point.DT.T4)[which(colnames(pcoa.point.DT.T4) %in% c("X1", "X2"))] <- c("bact.axis1", "bact.axis2")
colnames(pcoa.point.DT.T5)[which(colnames(pcoa.point.DT.T5) %in% c("X1", "X2"))] <- c("bact.axis1", "bact.axis2")
colnames(pcoa.point.DT.T6)[which(colnames(pcoa.point.DT.T6) %in% c("X1", "X2"))] <- c("bact.axis1", "bact.axis2")
colnames(pcoa.point.DT.T7)[which(colnames(pcoa.point.DT.T7) %in% c("X1", "X2"))] <- c("bact.axis1", "bact.axis2")



#Bray-Curtis dissimilarity index for drought tolerant

library(vegan)
bray.16S.DS.T1=vegdist(otu.16S.DS.T1,method = "bray")
pcoa.16S.DS.T1=cmdscale(sqrt(bray.16S.DT.T1),eig=T)
bray.16S.DS.T2=vegdist(otu.16S.DS.T2,method = "bray")
pcoa.16S.DS.T2=cmdscale(sqrt(bray.16S.DS.T2),eig=T)
bray.16S.DS.T3=vegdist(otu.16S.DS.T3,method = "bray")
pcoa.16S.DS.T3=cmdscale(sqrt(bray.16S.DS.T3),eig=T)
bray.16S.DS.T4=vegdist(otu.16S.DS.T4,method = "bray")
pcoa.16S.DS.T4=cmdscale(sqrt(bray.16S.DS.T4),eig=T)
bray.16S.DS.T5=vegdist(otu.16S.DS.T5,method = "bray")
pcoa.16S.DS.T5=cmdscale(sqrt(bray.16S.DS.T5),eig=T)
bray.16S.DS.T6=vegdist(otu.16S.DS.T6,method = "bray")
pcoa.16S.DS.T6=cmdscale(sqrt(bray.16S.DS.T6),eig=T)
bray.16S.DS.T7=vegdist(otu.16S.DS.T7,method = "bray")
pcoa.16S.DS.T7=cmdscale(sqrt(bray.16S.DS.T7),eig=T)

# PCoA point as data.frame for DS
pcoa.point.DS.T1=data.frame(pcoa.16S.DS.T1$points)#Axis1 8.23% #Axis2 5.44%
pcoa.point.DS.T2=data.frame(pcoa.16S.DS.T2$points)#Axis1 6.92% #Axis2 6.64%
pcoa.point.DS.T3=data.frame(pcoa.16S.DS.T3$points)#Axis1 7.36% #Axis2 5.36%
pcoa.point.DS.T4=data.frame(pcoa.16S.DS.T4$points)#Axis1 7.66% #Axis2 5.36%
pcoa.point.DS.T5=data.frame(pcoa.16S.DS.T5$points)#Axis1 6.81% #Axis2 5.95%
pcoa.point.DS.T6=data.frame(pcoa.16S.DS.T6$points)#Axis1 15.61% #Axis2 10.28%
pcoa.point.DS.T7=data.frame(pcoa.16S.DS.T7$points)#Axis1 18.07% #Axis2 11.15%

colnames(pcoa.point.DS.T1)[which(colnames(pcoa.point.DS.T1) %in% c("X1", "X2"))] <- c("bact.axis1", "bact.axis2")
colnames(pcoa.point.DS.T2)[which(colnames(pcoa.point.DS.T2) %in% c("X1", "X2"))] <- c("bact.axis1", "bact.axis2")
colnames(pcoa.point.DS.T3)[which(colnames(pcoa.point.DS.T3) %in% c("X1", "X2"))] <- c("bact.axis1", "bact.axis2")
colnames(pcoa.point.DS.T4)[which(colnames(pcoa.point.DS.T4) %in% c("X1", "X2"))] <- c("bact.axis1", "bact.axis2")
colnames(pcoa.point.DS.T5)[which(colnames(pcoa.point.DS.T5) %in% c("X1", "X2"))] <- c("bact.axis1", "bact.axis2")
colnames(pcoa.point.DS.T6)[which(colnames(pcoa.point.DS.T6) %in% c("X1", "X2"))] <- c("bact.axis1", "bact.axis2")
colnames(pcoa.point.DS.T7)[which(colnames(pcoa.point.DS.T7) %in% c("X1", "X2"))] <- c("bact.axis1", "bact.axis2")




# Bray-Curtis dissimilarity index for drought tolerant

library(vegan)
bray.ITS.DT.T1=vegdist(otu.ITS.DT.T1,method = "bray")
pcoa.ITS.DT.T1=cmdscale(sqrt(bray.ITS.DT.T1),eig=T)
bray.ITS.DT.T2=vegdist(otu.ITS.DT.T2,method = "bray")
pcoa.ITS.DT.T2=cmdscale(sqrt(bray.ITS.DT.T2),eig=T)
bray.ITS.DT.T3=vegdist(otu.ITS.DT.T3,method = "bray")
pcoa.ITS.DT.T3=cmdscale(sqrt(bray.ITS.DT.T3),eig=T)
bray.ITS.DT.T4=vegdist(otu.ITS.DT.T4,method = "bray")
pcoa.ITS.DT.T4=cmdscale(sqrt(bray.ITS.DT.T4),eig=T)
bray.ITS.DT.T5=vegdist(otu.ITS.DT.T5,method = "bray")
pcoa.ITS.DT.T5=cmdscale(sqrt(bray.ITS.DT.T5),eig=T)
bray.ITS.DT.T6=vegdist(otu.ITS.DT.T6,method = "bray")
pcoa.ITS.DT.T6=cmdscale(sqrt(bray.ITS.DT.T6),eig=T)
bray.ITS.DT.T7=vegdist(otu.ITS.DT.T7,method = "bray")
pcoa.ITS.DT.T7=cmdscale(sqrt(bray.ITS.DT.T7),eig=T)


# PCoA point as data.frame for DT
pcoa.point.ITS.DT.T1=data.frame(pcoa.ITS.DT.T1$points)#Axis1 11.72 % #Axis2 10.78%
pcoa.point.ITS.DT.T2=data.frame(pcoa.ITS.DT.T2$points)#Axis1 14.50% #Axis2 10.48%
pcoa.point.ITS.DT.T3=data.frame(pcoa.ITS.DT.T3$points)#Axis1 12.78% #Axis2 9.54%
pcoa.point.ITS.DT.T4=data.frame(pcoa.ITS.DT.T4$points)#Axis1 18.62% #Axis2 9.20%
pcoa.point.ITS.DT.T5=data.frame(pcoa.ITS.DT.T5$points)#Axis1 14.91% #Axis2 9.39%
pcoa.point.ITS.DT.T6=data.frame(pcoa.ITS.DT.T6$points)#Axis1 11.19% #Axis2 9.42%
pcoa.point.ITS.DT.T7=data.frame(pcoa.ITS.DT.T7$points)#Axis1 12.61% #Axis2 %9.10%

colnames(pcoa.point.ITS.DT.T1)[which(colnames(pcoa.point.ITS.DT.T1) %in% c("X1", "X2"))] <- c("fun.axis1", "fun.axis2")
colnames(pcoa.point.ITS.DT.T2)[which(colnames(pcoa.point.ITS.DT.T2) %in% c("X1", "X2"))] <- c("fun.axis1", "fun.axis2")
colnames(pcoa.point.ITS.DT.T3)[which(colnames(pcoa.point.ITS.DT.T3) %in% c("X1", "X2"))] <- c("fun.axis1", "fun.axis2")
colnames(pcoa.point.ITS.DT.T4)[which(colnames(pcoa.point.ITS.DT.T4) %in% c("X1", "X2"))] <- c("fun.axis1", "fun.axis2")
colnames(pcoa.point.ITS.DT.T5)[which(colnames(pcoa.point.ITS.DT.T5) %in% c("X1", "X2"))] <- c("fun.axis1", "fun.axis2")
colnames(pcoa.point.ITS.DT.T6)[which(colnames(pcoa.point.ITS.DT.T6) %in% c("X1", "X2"))] <- c("fun.axis1", "fun.axis2")
colnames(pcoa.point.ITS.DT.T7)[which(colnames(pcoa.point.ITS.DT.T7) %in% c("X1", "X2"))] <- c("fun.axis1", "fun.axis2")



# Bray-Curtis dissimilarity index for drought sensitive

library(vegan)
bray.ITS.DS.T1=vegdist(otu.ITS.DS.T1,method = "bray")
pcoa.ITS.DS.T1=cmdscale(sqrt(bray.ITS.DS.T1),eig=T)
bray.ITS.DS.T2=vegdist(otu.ITS.DS.T2,method = "bray")
pcoa.ITS.DS.T2=cmdscale(sqrt(bray.ITS.DS.T2),eig=T)
bray.ITS.DS.T3=vegdist(otu.ITS.DS.T3,method = "bray")
pcoa.ITS.DS.T3=cmdscale(sqrt(bray.ITS.DS.T3),eig=T)
bray.ITS.DS.T4=vegdist(otu.ITS.DS.T4,method = "bray")
pcoa.ITS.DS.T4=cmdscale(sqrt(bray.ITS.DS.T4),eig=T)
bray.ITS.DS.T5=vegdist(otu.ITS.DS.T5,method = "bray")
pcoa.ITS.DS.T5=cmdscale(sqrt(bray.ITS.DS.T5),eig=T)
bray.ITS.DS.T6=vegdist(otu.ITS.DS.T6,method = "bray")
pcoa.ITS.DS.T6=cmdscale(sqrt(bray.ITS.DS.T6),eig=T)
bray.ITS.DS.T7=vegdist(otu.ITS.DS.T7,method = "bray")
pcoa.ITS.DS.T7=cmdscale(sqrt(bray.ITS.DS.T7),eig=T)


# PCoA point as data.frame for DT
pcoa.point.ITS.DS.T1=data.frame(pcoa.ITS.DS.T1$points)#Axis1 13.15% #Axis2 9.78%
pcoa.point.ITS.DS.T2=data.frame(pcoa.ITS.DS.T2$points)#Axis1 12.37% #Axis2 10.21%
pcoa.point.ITS.DS.T3=data.frame(pcoa.ITS.DS.T3$points)#Axis1 13.06% #Axis2 10.59%
pcoa.point.ITS.DS.T4=data.frame(pcoa.ITS.DS.T4$points)#Axis1 13.08% #Axis2 10.69%
pcoa.point.ITS.DS.T5=data.frame(pcoa.ITS.DS.T5$points)#Axis1 13.91% #Axis2 10.38%
pcoa.point.ITS.DS.T6=data.frame(pcoa.ITS.DS.T6$points)#Axis1 13.40% #Axis2 8.95%
pcoa.point.ITS.DS.T7=data.frame(pcoa.ITS.DS.T7$points)#Axis1 11.71% #Axis2 8.66%

colnames(pcoa.point.ITS.DS.T1)[which(colnames(pcoa.point.ITS.DS.T1) %in% c("X1", "X2"))] <- c("fun.axis1", "fun.axis2")
colnames(pcoa.point.ITS.DS.T2)[which(colnames(pcoa.point.ITS.DS.T2) %in% c("X1", "X2"))] <- c("fun.axis1", "fun.axis2")
colnames(pcoa.point.ITS.DS.T3)[which(colnames(pcoa.point.ITS.DS.T3) %in% c("X1", "X2"))] <- c("fun.axis1", "fun.axis2")
colnames(pcoa.point.ITS.DS.T4)[which(colnames(pcoa.point.ITS.DS.T4) %in% c("X1", "X2"))] <- c("fun.axis1", "fun.axis2")
colnames(pcoa.point.ITS.DS.T5)[which(colnames(pcoa.point.ITS.DS.T5) %in% c("X1", "X2"))] <- c("fun.axis1", "fun.axis2")
colnames(pcoa.point.ITS.DS.T6)[which(colnames(pcoa.point.ITS.DS.T6) %in% c("X1", "X2"))] <- c("fun.axis1", "fun.axis2")
colnames(pcoa.point.ITS.DS.T7)[which(colnames(pcoa.point.ITS.DS.T7) %in% c("X1", "X2"))] <- c("fun.axis1", "fun.axis2")



#Ordination axis for ITS

pcoa.ITS.DT.T7$eig[1]/sum(pcoa.ITS.DT.T7$eig)*100#Axis1
pcoa.ITS.DT.T7$eig[2]/sum(pcoa.ITS.DT.T7$eig)*100#Axis2

pcoa.ITS.DS.T1$eig[1]/sum(pcoa.ITS.DS.T1$eig)*100#Axis1
pcoa.ITS.DS.T1$eig[2]/sum(pcoa.ITS.DS.T1$eig)*100#Axis2

#Ordination testing
library(ggplot2)
library(viridis)
pcoa.16S.DT.T1$eig[1]/sum(pcoa.16S.DT.T1$eig)*100#Axis1
pcoa.16S.DT.T1$eig[2]/sum(pcoa.16S.DT.T1$eig)*100#Axis2

pcoa.16S.DS.T7$eig[1]/sum(pcoa.16S.DS.T7$eig)*100#Axis1
pcoa.16S.DS.T7$eig[2]/sum(pcoa.16S.DS.T7$eig)*100#Axis2



#Keep only top 10 OTU with strongest correlations for Drought  Sensitive (Output correlation file)
taxa.DS.T1=c("OTU_10214","OTU_10278", "OTU_10378","OTU_11143","OTU_11397","OTU_11475","OTU_1183","OTU_12150","OTU_12183","OTU_12189")
bact.DS.T1=data.frame(otu.16S.DS.T1[,colnames(otu.16S.DS.T1) %in% taxa.DS.T1])

taxa.DS.T2=c("OTU_10002","OTU_10233","OTU_10425","OTU_10948","OTU_1270","OTU_1325","OTU_1343","OTU_14681","OTU_16054", "OTU_516")
bact.DS.T2=data.frame(otu.16S.DS.T2[,colnames(otu.16S.DS.T2) %in% taxa.DS.T2])

taxa.DS.T3=c("OTU_10331", "OTU_10632", "OTU_10742","OTU_13085","OTU_9449","OTU_11673","OTU_1174",  "OTU_1187" , "OTU_11878","OTU_12307")
bact.DS.T3=data.frame(otu.16S.DS.T3[,colnames(otu.16S.DS.T3) %in% taxa.DS.T3])

taxa.DS.T4=c("OTU_10024","OTU_10214", "OTU_10258",  "OTU_10298","OTU_10470","OTU_10497","OTU_10644","OTU_10676","OTU_10740","OTU_10981")
bact.DS.T4=data.frame(otu.16S.DS.T4[,colnames(otu.16S.DS.T4) %in% taxa.DS.T4])

taxa.DS.T5=c("OTU_10038","OTU_10282", "OTU_10553","OTU_11058","OTU_11242","OTU_1143","OTU_11844", "OTU_1217" , "OTU_1242",  "OTU_13629")
bact.DS.T5=data.frame(otu.16S.DS.T5[,colnames(otu.16S.DS.T5) %in% taxa.DS.T5])

taxa.DS.T6=c("OTU_10050","OTU_1007","OTU_10416","OTU_1053", "OTU_1067", "OTU_1079","OTU_1082","OTU_1091", "OTU_11025","OTU_1188" )
bact.DS.T6=data.frame(otu.16S.DS.T6[,colnames(otu.16S.DS.T6) %in% taxa.DS.T6])

taxa.DS.T7=c("OTU_10371","OTU_10391","OTU_10409","OTU_12375","OTU_12465","OTU_13582","OTU_14604","OTU_15282","OTU_18314","OTU_18504")
bact.DS.T7=data.frame(otu.16S.DS.T7[,colnames(otu.16S.DS.T7) %in% taxa.DS.T7])


#Keep only 10 ITS OTU with the strongest correlations for Drought  Sensitive cultivar (Output correlation file)
taxa.ITS.DS.T1=c("OTU_18","OTU_1077", "OTU_1198","OTU_133","OTU_151","OTU_2199","OTU_330","OTU_3622","OTU_3898","OTU_638")
ITS.DS.T1=data.frame(otu.ITS.DS.T1[,colnames(otu.ITS.DS.T1) %in% taxa.ITS.DS.T1])

taxa.ITS.DS.T2=c("OTU_1184","OTU_1277","OTU_1331","OTU_1474","OTU_161","OTU_2760","OTU_2956","OTU_3382","OTU_3499", "OTU_379")
ITS.DS.T2=data.frame(otu.ITS.DS.T2[,colnames(otu.ITS.DS.T2) %in% taxa.ITS.DS.T2])

taxa.ITS.DS.T3=c("OTU_1088","OTU_1142", "OTU_1330","OTU_1541","OTU_156","OTU_1703","OTU_2361","OTU_2557","OTU_3008","OTU_3172")
ITS.DS.T3=data.frame(otu.ITS.DS.T3[,colnames(otu.ITS.DS.T3) %in% taxa.ITS.DS.T3])

taxa.ITS.DS.T4=c("OTU_105","OTU_1090","OTU_132","OTU_2597","OTU_4728","OTU_5266","OTU_589","OTU_6120","OTU_7030", "OTU_813")
ITS.DS.T4=data.frame(otu.ITS.DS.T4[,colnames(otu.ITS.DS.T4) %in% taxa.ITS.DS.T4])

taxa.ITS.DS.T5=c("OTU_1032","OTU_948","OTU_1416","OTU_1585","OTU_1617","OTU_1659","OTU_1965","OTU_2062","OTU_2133", "OTU_2788")
ITS.DS.T5=data.frame(otu.ITS.DS.T5[,colnames(otu.ITS.DS.T5) %in% taxa.ITS.DS.T5])

taxa.ITS.DS.T6=c("OTU_1106","OTU_1236","OTU_1247","OTU_1250","OTU_303","OTU_94","OTU_1667","OTU_281","OTU_5084", "OTU_5873")
ITS.DS.T6=data.frame(otu.ITS.DS.T6[,colnames(otu.ITS.DS.T6) %in% taxa.ITS.DS.T6])

taxa.ITS.DS.T7=c("OTU_1694","OTU_2165","OTU_330","OTU_393","OTU_4025","OTU_4059","OTU_4372","OTU_5879","OTU_77", "OTU_791")
ITS.DS.T7=data.frame(otu.ITS.DS.T7[,colnames(otu.ITS.DS.T7) %in% taxa.ITS.DS.T7])



#Keep only 10 OTU with strongest correlations for drought tolerance cultivar (Output correlation file)
taxa.DT.T1=c("OTU_111", "OTU_1882","OTU_10146", "OTU_33589", "OTU_34132","OTU_34227", "OTU_34381","OTU_34480","OTU_34516","OTU_41461")
bact.DT.T1=data.frame(otu.16S.DT.T1[,colnames(otu.16S.DT.T1) %in% taxa.DT.T1])

taxa.DT.T2=c("OTU_111","OTU_6321","OTU_1000","OTU_10060","OTU_380","OTU_10454","OTU_10488","OTU_1651","OTU_17999","OTU_11759")
bact.DT.T2=data.frame(otu.16S.DT.T2[,colnames(otu.16S.DT.T2) %in% taxa.DT.T2])
  
taxa.DT.T3=c("OTU_10712","OTU_10782","OTU_10866","OTU_10957","OTU_2","OTU_2317","OTU_2341","OTU_2348","OTU_23745","OTU_43433")
bact.DT.T3=data.frame(otu.16S.DT.T3[,colnames(otu.16S.DT.T3) %in% taxa.DT.T3])

taxa.DT.T4=c("OTU_10000","OTU_25804","OTU_192","OTU_430","OTU_122","OTU_14774","OTU_14816","OTU_14926","OTU_1554","OTU_1605")
bact.DT.T4=data.frame(otu.16S.DT.T4[,colnames(otu.16S.DT.T4) %in% taxa.DT.T4])

taxa.DT.T5=c("OTU_1043","OTU_10631","OTU_12550","OTU_12560","OTU_11670","OTU_1304","OTU_1680","OTU_17011","OTU_18380","OTU_30617")
bact.DT.T5=data.frame(otu.16S.DT.T5[,colnames(otu.16S.DT.T5) %in% taxa.DT.T5])

taxa.DT.T6=c("OTU_10365","OTU_1037","OTU_1044","OTU_43710","OTU_10687","OTU_11003","OTU_11321","OTU_1150","OTU_1078","OTU_1148" ) 
bact.DT.T6=data.frame(otu.16S.DT.T6[,colnames(otu.16S.DT.T6) %in% taxa.DT.T6])

taxa.DT.T7=c("OTU_21928","OTU_1038","OTU_10633","OTU_23744","OTU_40119","OTU_620","OTU_11394","OTU_13293","OTU_14325","OTU_14336" ) 
bact.DT.T7=data.frame(otu.16S.DT.T7[,colnames(otu.16S.DT.T7) %in% taxa.DT.T7])

  
#Keep only 10 fungal OTU with strongest correlations for Drought (Output correlation file)  
 
taxa.ITS.DT.T1=c("OTU_6761","OTU_1032", "OTU_1498","OTU_1592","OTU_1936","OTU_2077","OTU_2238","OTU_2717","OTU_46","OTU_466")
ITS.DT.T1=data.frame(otu.ITS.DT.T1[,colnames(otu.ITS.DT.T1) %in% taxa.ITS.DT.T1])

taxa.ITS.DT.T2=c("OTU_17","OTU_87","OTU_6894","OTU_1831","OTU_2184","OTU_2339", "OTU_2399","OTU_257","OTU_2641","OTU_455")
ITS.DT.T2=data.frame(otu.ITS.DT.T2[,colnames(otu.ITS.DT.T2) %in% taxa.ITS.DT.T2])

taxa.ITS.DT.T3=c("OTU_226","OTU_2278","OTU_1029","OTU_1032","OTU_1040","OTU_1123","OTU_1148","OTU_1246","OTU_1278","OTU_1313")
ITS.DT.T3=data.frame(otu.ITS.DT.T3[,colnames(otu.ITS.DT.T3) %in% taxa.ITS.DT.T3])

taxa.ITS.DT.T4=c("OTU_1061","OTU_107","OTU_1071","OTU_1074","OTU_1076","OTU_1307","OTU_1472","OTU_1656","OTU_2042", "OTU_281")
ITS.DT.T4=data.frame(otu.ITS.DT.T4[,colnames(otu.ITS.DT.T4) %in% taxa.ITS.DT.T4])
 
  
taxa.ITS.DT.T5=c("OTU_16","OTU_2","OTU_27","OTU_338","OTU_4","OTU_5","OTU_52","OTU_62","OTU_6408","OTU_6576")
ITS.DT.T5=data.frame(otu.ITS.DT.T5[,colnames(otu.ITS.DT.T5) %in% taxa.ITS.DT.T5])

 
taxa.ITS.DT.T6=c("OTU_136","OTU_5749","OTU_1375","OTU_142","OTU_150","OTU_1568","OTU_181","OTU_276","OTU_2995", "OTU_3537")
ITS.DT.T6=data.frame(otu.ITS.DT.T6[,colnames(otu.ITS.DT.T6) %in% taxa.ITS.DT.T6])

taxa.ITS.DT.T7=c("OTU_6811","OTU_590","OTU_4","OTU_276","OTU_140","OTU_1491","OTU_1515","OTU_161","OTU_1640","OTU_1924")
ITS.DT.T7=data.frame(otu.ITS.DT.T7[,colnames(otu.ITS.DT.T7) %in% taxa.ITS.DT.T7])

# Double Sanity check 
row.names(qual.DT.S1)==row.names(bact.DT.T1)
row.names(metafile.DT.T1)==row.names(bact.DT.T1)
row.names(qual.DT.S1)==row.names(metafile.DT.T1)
row.names(qual.DT.S2)==row.names(bact.DT.T2)
row.names(metafile.DT.T2)==row.names(bact.DT.T2)
row.names(qual.DT.S2)==row.names(metafile.DT.T2)
row.names(qual.DT.S3)==row.names(bact.DT.T3)
row.names(metafile.DT.T3)==row.names(bact.DT.T3)
row.names(qual.DT.S3)==row.names(metafile.DT.T3)
row.names(qual.DT.S4)==row.names(bact.DT.T4)
row.names(metafile.DT.T4)==row.names(bact.DT.T4)
row.names(qual.DT.S4)==row.names(metafile.DT.T4) 
row.names(qual.DT.S5)==row.names(bact.DT.T5)
row.names(metafile.DT.T5)==row.names(bact.DT.T5)
row.names(qual.DT.S5)==row.names(metafile.DT.T5)  
row.names(qual.DT.S6)==row.names(bact.DT.T6)
row.names(metafile.DT.T6)==row.names(bact.DT.T6)
row.names(qual.DT.S6)==row.names(metafile.DT.T6)   
row.names(qual.DT.S7)==row.names(bact.DT.T7)
row.names(metafile.DT.T7)==row.names(bact.DT.T7)
row.names(qual.DT.S7)==row.names(metafile.DT.T7)   


row.names(qual.DT.S1)==row.names(ITS.DT.T1)
row.names(metafile.DT.T1)==row.names(ITS.DT.T1)
row.names(qual.DT.S1)==row.names(metafile.DT.T1)
row.names(qual.DT.S2)==row.names(ITS.DT.T2)
row.names(metafile.DT.T2)==row.names(ITS.DT.T2)
row.names(qual.DT.S2)==row.names(metafile.DT.T2)
row.names(qual.DT.S3)==row.names(ITS.DT.T3)
row.names(metafile.DT.T3)==row.names(ITS.DT.T3)
row.names(qual.DT.S3)==row.names(metafile.DT.T3)
row.names(qual.DT.S4)==row.names(ITS.DT.T4)
row.names(metafile.DT.T4)==row.names(ITS.DT.T4)
row.names(qual.DT.S4)==row.names(metafile.DT.T4) 
row.names(qual.DT.S5)==row.names(ITS.DT.T5)
row.names(metafile.DT.T5)==row.names(ITS.DT.T5)
row.names(qual.DT.S5)==row.names(metafile.DT.T5)  
row.names(qual.DT.S6)==row.names(ITS.DT.T6)
row.names(metafile.DT.T6)==row.names(ITS.DT.T6)
row.names(qual.DT.S6)==row.names(metafile.DT.T6)   
row.names(qual.DT.S7)==row.names(ITS.DT.T7)
row.names(metafile.DT.T7)==row.names(ITS.DT.T7)
row.names(qual.DT.S7)==row.names(metafile.DT.T7)  


#Sanity double check
row.names(qual.DS.S1)==row.names(bact.DS.T1)
row.names(metafile.DS.T1)==row.names(bact.DS.T1)
row.names(qual.DS.S1)==row.names(metafile.DS.T1)
row.names(qual.DS.S2)==row.names(bact.DS.T2)
row.names(metafile.DS.T2)==row.names(bact.DS.T2)
row.names(qual.DS.S2)==row.names(metafile.DS.T2)
row.names(qual.DS.S3)==row.names(bact.DS.T3)
row.names(metafile.DS.T3)==row.names(bact.DS.T3)
row.names(qual.DS.S3)==row.names(metafile.DS.T3)
row.names(qual.DS.S4)==row.names(bact.DS.T4)
row.names(metafile.DS.T4)==row.names(bact.DS.T4)
row.names(qual.DS.S4)==row.names(metafile.DS.T4)
row.names(qual.DS.S5)==row.names(bact.DS.T5)
row.names(metafile.DS.T5)==row.names(bact.DS.T5)
row.names(qual.DS.S5)==row.names(metafile.DS.T5)
row.names(qual.DS.S6)==row.names(bact.DS.T6)
row.names(metafile.DS.T6)==row.names(bact.DS.T6)
row.names(qual.DS.S6)==row.names(metafile.DS.T6)
row.names(qual.DS.S7)==row.names(bact.DS.T7)
row.names(metafile.DS.T7)==row.names(bact.DS.T7)
row.names(qual.DS.S7)==row.names(metafile.DS.T7)


row.names(qual.DS.S1)==row.names(ITS.DS.T1)
row.names(metafile.DS.T1)==row.names(ITS.DS.T1)
row.names(qual.DS.S1)==row.names(metafile.DS.T1)
row.names(qual.DS.S2)==row.names(ITS.DS.T2)
row.names(metafile.DS.T2)==row.names(ITS.DS.T2)
row.names(qual.DS.S2)==row.names(metafile.DS.T2)
row.names(qual.DS.S3)==row.names(ITS.DS.T3)
row.names(metafile.DS.T3)==row.names(ITS.DS.T3)
row.names(qual.DS.S3)==row.names(metafile.DS.T3)
row.names(qual.DS.S4)==row.names(ITS.DS.T4)
row.names(metafile.DS.T4)==row.names(ITS.DS.T4)
row.names(qual.DS.S4)==row.names(metafile.DS.T4)
row.names(qual.DS.S5)==row.names(ITS.DS.T5)
row.names(metafile.DS.T5)==row.names(ITS.DS.T5)
row.names(qual.DS.S5)==row.names(metafile.DS.T5)
row.names(qual.DS.S6)==row.names(ITS.DS.T6)
row.names(metafile.DS.T6)==row.names(ITS.DS.T6)
row.names(qual.DS.S6)==row.names(metafile.DS.T6)
row.names(qual.DS.S7)==row.names(ITS.DS.T7)
row.names(metafile.DS.T7)==row.names(ITS.DS.T7)
row.names(qual.DS.S7)==row.names(metafile.DS.T7)

#Data Preparation for regression model
library(vegan)
reg.data.DT.T1=cbind(qual.DT.S1,bact.DT.T1,ITS.DT.T1,pcoa.point.DT.T1,pcoa.point.ITS.DT.T1, div.16S.DT.T1,div.ITS.DT.T1,biolog.DT.T1,qpcr.T1.DT)
reg.data.DT.T1=na.omit(reg.data.DT.T1) #reg.data.DT.T1=decostand(reg.data.DT.T1,"normalize")
reg.data.DT.T2=cbind(metafile.DT.T2,qual.DT.S2,bact.DT.T2,ITS.DT.T2,pcoa.point.DT.T2,pcoa.point.ITS.DT.T2, div.16S.DT.T2,div.ITS.DT.T2,biolog.DT.T2,qpcr.T2.DT)
reg.data.DT.T2=na.omit(reg.data.DT.T2)
reg.data.DT.T3=cbind(metafile.DT.T3,qual.DT.S3,bact.DT.T3,ITS.DT.T3,pcoa.point.DT.T3,pcoa.point.ITS.DT.T3, div.16S.DT.T3,div.ITS.DT.T3,biolog.DT.T3,qpcr.T3.DT)
reg.data.DT.T3=na.omit(reg.data.DT.T3)
reg.data.DT.T4=cbind(metafile.DT.T4,qual.DT.S4,bact.DT.T4,ITS.DT.T4,pcoa.point.DT.T4,pcoa.point.ITS.DT.T4, div.16S.DT.T4,div.ITS.DT.T4,biolog.DT.T4,qpcr.T4.DT)
reg.data.DT.T4=na.omit(reg.data.DT.T4)
reg.data.DT.T5=cbind(metafile.DT.T5,qual.DT.S5,bact.DT.T5,ITS.DT.T5,pcoa.point.DT.T5,pcoa.point.ITS.DT.T5, div.16S.DT.T5,div.ITS.DT.T5,biolog.DT.T5,qpcr.T5.DT)
reg.data.DT.T5=na.omit(reg.data.DT.T5)
reg.data.DT.T6=cbind(metafile.DT.T6,qual.DT.S6,bact.DT.T6,ITS.DT.T6,pcoa.point.DT.T6,pcoa.point.ITS.DT.T6, div.16S.DT.T6,div.ITS.DT.T6,biolog.DT.T6,qpcr.T6.DT)
reg.data.DT.T6=na.omit(reg.data.DT.T6)
reg.data.DT.T7=cbind(metafile.DT.T7,qual.DT.S7,bact.DT.T7,ITS.DT.T7,pcoa.point.DT.T7,pcoa.point.ITS.DT.T7, div.16S.DT.T7,div.ITS.DT.T7,biolog.DT.T7,qpcr.T7.DT)
reg.data.DT.T7=na.omit(reg.data.DT.T7)

#Data Preparation for regression model
reg.data.DS.T1=cbind(metafile.DS.T1,qual.DS.S1,bact.DS.T1,ITS.DS.T1,pcoa.point.DS.T1,pcoa.point.ITS.DS.T1,div.16S.DS.T1,div.ITS.DS.T1,biolog.DS.T1,qpcr.T1.DS)
reg.data.DS.T1=na.omit(reg.data.DS.T1)
reg.data.DS.T2=cbind(metafile.DS.T2,qual.DS.S2,bact.DS.T2,ITS.DS.T2,pcoa.point.DS.T2,pcoa.point.ITS.DS.T2,div.16S.DS.T2,div.ITS.DS.T2,biolog.DS.T2,qpcr.T2.DS)
reg.data.DS.T2=na.omit(reg.data.DS.T2)
reg.data.DS.T3=cbind(metafile.DS.T3,qual.DS.S3,bact.DS.T3,ITS.DS.T3,pcoa.point.DS.T3,pcoa.point.ITS.DS.T3,div.16S.DS.T3,div.ITS.DS.T3,biolog.DS.T3,qpcr.T3.DS)
reg.data.DS.T3=na.omit(reg.data.DS.T3)
reg.data.DS.T4=cbind(metafile.DS.T4,qual.DS.S4,bact.DS.T4,ITS.DS.T4,pcoa.point.DS.T4,pcoa.point.ITS.DS.T4,div.16S.DS.T4,div.ITS.DS.T4,biolog.DS.T4,qpcr.T4.DS)
reg.data.DS.T4=na.omit(reg.data.DS.T4)
reg.data.DS.T5=cbind(metafile.DS.T5,qual.DS.S5,bact.DS.T5,ITS.DS.T5,pcoa.point.DS.T5,pcoa.point.ITS.DS.T5,div.16S.DS.T5,div.ITS.DS.T5,biolog.DS.T5,qpcr.T5.DS)
reg.data.DS.T5=na.omit(reg.data.DS.T5)
reg.data.DS.T6=cbind(metafile.DS.T6,qual.DS.S6,bact.DS.T6,ITS.DS.T6,pcoa.point.DS.T6,pcoa.point.ITS.DS.T6,div.16S.DS.T6,div.ITS.DS.T6,biolog.DS.T6,qpcr.T6.DS)
reg.data.DS.T6=na.omit(reg.data.DS.T6)
reg.data.DS.T7=cbind(metafile.DS.T7,qual.DS.S7,bact.DS.T7,ITS.DS.T7,pcoa.point.DS.T7,pcoa.point.ITS.DS.T7,div.16S.DS.T7,div.ITS.DS.T7,biolog.DS.T7,qpcr.T7.DS)
reg.data.DS.T7=na.omit(reg.data.DS.T7)



#Modelling for drought tolerant cultivar

## Gluten model for May 10
data.gluten.DT.T1=reg.data.DT.T1[!is.na(reg.data.DT.T1$Gluten),]

model.null = lm(Gluten ~ 1, data=data.gluten.DT.T1)
model.full = lm(Gluten ~OTU_111+OTU_1882+OTU_10146+OTU_33589+OTU_34132+OTU_34227+OTU_34381+OTU_34480+OTU_34516+OTU_41461+OTU_6761+OTU_1032+OTU_1498+OTU_1592+OTU_1936+OTU_2077+OTU_2238+OTU_2717+OTU_46+OTU_466+ 
                        bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio+D.Galacturonic.Acid  
                +y.Amino.Butyric.Acid  
                +Pyruvic.Acid.Methyl.Ester  
                +Putrescine  
                +D.Galactonic.Acid.y.Lactone  
                +Tween.80  
                +L.Serine  
                ,data=data.gluten.DT.T1)

step(model.null,scope = list(upper=model.full), direction="both",  data=data.gluten.DT.T1)  

model.gluten.DT.T1 = lm(Gluten ~ OTU_1882 + OTU_41461 + nirk + AOB + OTU_6761 + 
                          Tween.80 + ACE.fun  ,  data  = data.gluten.DT.T1)

#model.gluten.DT.T1 = lm(Gluten ~ Simpson + OTU_41461 + Simpson.fun + a.Keto.Butyric.Acid + 
                          #ACE ,  data  = data.gluten.DT.T1)

summary(model.gluten.DT.T1)

data.gluten.DT.T1$GlutenPred=predict(model.gluten.DT.T1)

p1 <- ggplot(data.gluten.DT.T1, aes(Gluten, GlutenPred)) + 
  geom_point(size = 1, alpha = 5.0) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic()+
  scale_y_log10("Predicted gluten content")+
  #scale_y_continuous("Gluten") +
  xlab("Observed gluten content") +
  ggtitle(paste("Prediction of gluten"))

p1

gridExtra::grid.arrange(p1, nrow = 1)

# fit with two strongly correlated variables
summary(model.gluten.DT.T1) %>%
  broom::tidy() %>%
  filter(term %in% c("OTU_1882" ,"OTU_41461" ,"nirk"  , "AOB" ,"OTU_6761" , 
                       "Tween.80", "ACE.fun" ))





## Protein model for May 10
data.protein.DT.T1=reg.data.DT.T1[!is.na(reg.data.DT.T1$Protein.grain),]

model.null = lm(Protein.grain ~ 1, data=data.protein.DT.T1)
model.full = lm(Protein.grain ~OTU_111+OTU_1882+OTU_10146+OTU_33589+OTU_34132+OTU_34227+OTU_34381+OTU_34480+OTU_34516+OTU_41461+OTU_6761+OTU_1032+OTU_1498+
                               OTU_1592+OTU_1936+OTU_2077+OTU_2238+OTU_2717+OTU_46+OTU_466+bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio 
                +D.Galacturonic.Acid  
                +y.Amino.Butyric.Acid  
                +Pyruvic.Acid.Methyl.Ester  
                +Putrescine  
                +D.Galactonic.Acid.y.Lactone  
                +Tween.80  
                +L.Serine  
                , data=data.protein.DT.T1)

step(model.null, scope = list(upper=model.full), direction="both",  data=data.protein.DT.T1)  

model.protein.DT.T1 = lm(Protein.grain ~OTU_41461 + OTU_111 + D.Galacturonic.Acid + 
                           y.Amino.Butyric.Acid + D.Galactonic.Acid.y.Lactone ,  data = data.protein.DT.T1)

summary(model.protein.DT.T1)

## PMT model for May 10
data.PMT.DT.T1=reg.data.DT.T1[!is.na(reg.data.DT.T1$PMT),]

model.null = lm(PMT ~ 1, data=data.PMT.DT.T1)
model.full = lm(PMT ~OTU_111+OTU_1882+OTU_10146+OTU_33589+OTU_34132+OTU_34227+OTU_34381+OTU_34480+OTU_34516+OTU_41461+OTU_6761+OTU_1032+OTU_1498+
                OTU_1592+OTU_1936+OTU_2077+OTU_2238+OTU_2717+OTU_46+OTU_466+bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+ 
                +D.Galacturonic.Acid  
                +y.Amino.Butyric.Acid  
                +Pyruvic.Acid.Methyl.Ester  
                +Putrescine  
                +D.Galactonic.Acid.y.Lactone  
                +Tween.80  
                +L.Serine  
                , data=data.PMT.DT.T1)

step(model.null, scope = list(upper=model.full),K=4,direction="both", data=data.PMT.DT.T1)  

model.PMT.DT.T1 =lm(PMT ~ y.Amino.Butyric.Acid + OTU_41461  
   ,  data = data.PMT.DT.T1)

summary(model.PMT.DT.T1)


## BEM model for May 10
data.BEM.DT.T1=reg.data.DT.T1[!is.na(reg.data.DT.T1$BEM),]

model.null = lm(BEM ~ 1, data=data.BEM.DT.T1)
model.full = lm(BEM ~OTU_111+OTU_1882+OTU_10146+OTU_33589+OTU_34132+OTU_34227+OTU_34381+OTU_34480+OTU_34516+OTU_41461+OTU_6761+OTU_1032+OTU_1498+
                  OTU_1592+OTU_1936+OTU_2077+OTU_2238+OTU_2717+OTU_46+OTU_466+bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                +D.Galacturonic.Acid  
                +y.Amino.Butyric.Acid  
                +Pyruvic.Acid.Methyl.Ester  
                +Putrescine  
                +D.Galactonic.Acid.y.Lactone  
                +Tween.80  
                +L.Serine  
                , data=data.BEM.DT.T1)

step(model.null, scope = list(upper=model.full,lower=model.null), k=4,direction="both", data=data.BEM.DT.T1)  

model.BEM.DT.T1 =lm(BEM ~ OTU_6761 + OTU_2717 + y.Amino.Butyric.Acid + 
                      Simpson.fun + AOB + OTU_46 + D.Galactonic.Acid.y.Lactone
                    ,  data = data.BEM.DT.T1)

summary(model.BEM.DT.T1)




##Protein model for May 24

data.protein.DT.T2=reg.data.DT.T2[!is.na(reg.data.DT.T2$Protein.grain),]

model.null = lm(Protein.grain ~ 1, data=data.protein.DT.T2)

model.full = lm(Protein.grain ~ OTU_111+OTU_6321+OTU_1000+OTU_10060+OTU_10454+OTU_10488+OTU_1651+OTU_17999+OTU_11759+OTU_380+OTU_17+OTU_87+OTU_6894+OTU_1831+OTU_2184+OTU_2339 +OTU_2399 +OTU_257+OTU_2641+OTU_455+bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+
                +AOA+AOB+nirk+nosZ+F.B.ratio
                +Pyruvic.Acid.Methyl.Ester  
                +D.Xylose  
                +Tween.40  
                +D.Mannitol  
                +N.Acetyl.D.Glucosamine  
                +D.Cellobiose  
                +a.Keto.Butyric.Acid  
                , data=data.protein.DT.T2)


step(model.null, scope = list(upper=model.full, lower=model.null), direction="both",k=6,  data=data.protein.DT.T2)  



model.protein.DT.T2 = lm(Protein.grain ~ OTU_111 + OTU_87 + D.Cellobiose + 
                           nosZ + OTU_380  , data = data.protein.DT.T2)


summary(model.protein.DT.T2)


library(leaps)
## PMT model for May 24
data.PMT.DT.T2=reg.data.DT.T2[!is.na(reg.data.DT.T2$PMT),]

model.null = lm(PMT ~ 1, data=data.PMT.DT.T2)
model.full = lm(PMT ~ OTU_111+OTU_6321+OTU_1000+OTU_10060+OTU_10454+OTU_10488+OTU_1651+OTU_17999+OTU_11759+OTU_380+OTU_17+OTU_87+OTU_6894+OTU_1831+OTU_2184+OTU_2339 +OTU_2399 +OTU_257+OTU_2641+OTU_455+bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+
                +AOA+AOB+nirk+nosZ+F.B.ratio
                +Pyruvic.Acid.Methyl.Ester  
                +D.Xylose  
                +Tween.40  
                +D.Mannitol  
                +N.Acetyl.D.Glucosamine  
                +D.Cellobiose  
                +a.Keto.Butyric.Acid   
                , data=data.PMT.DT.T2)


models <- leaps::regsubsets(PMT~ OTU_111+OTU_6321+OTU_1000+OTU_10060+OTU_10454+OTU_10488+OTU_1651+OTU_17999+OTU_11759+OTU_380+OTU_17+OTU_87+OTU_6894+OTU_1831+OTU_2184+OTU_2339 +OTU_2399 +OTU_257+OTU_2641+OTU_455+bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+
                              +AOA+AOB+nirk+nosZ+F.B.ratio
                            +Pyruvic.Acid.Methyl.Ester  
                            +D.Xylose  
                            +Tween.40  
                            +D.Mannitol  
                            +N.Acetyl.D.Glucosamine  
                            +D.Cellobiose  
                            +a.Keto.Butyric.Acid, data=data.PMT.DT.T2,
                            nvmax = 11, nbest=3,
                            method = "forward")



step(model.null, scope = list(upper=model.full,lower=model.null), k=3,direction="both", data=data.PMT.DT.T2)  


step(model.null, scope=list(lower=model.null, upper=model.full),direction="both", trace = 0, k=log(nrow(data.PMT.DT.T2)))



model.PMT.DT.T2 =lm(PMT ~ N.Acetyl.D.Glucosamine + fun.axis1 + D.Xylose + 
                      Pyruvic.Acid.Methyl.Ester 
                    ,  data = data.PMT.DT.T2)

summary(model.PMT.DT.T2)


# model validation
op <- par(oma=c(5,7,1,1))
par(op)
plot(models, scale="bic")



smodels = summary(models)
head(smodels$which[,1:5])


nvar <- apply(smodels$which,1,sum)-1
par(mar = c(2, 2, 2, 2))
plot(nvar, smodels$bic, xlab = "Number of Variables", ylab = "BIC")
min.bic <- which.min(smodels$bic)
points(nvar[min.bic], smodels$bic[min.bic], pch = 20, col = "red")
abline(h = smodels$bic[min.bic]+2, lty=2)


k <- apply(smodels$which,1,sum)+1
mod.aicc <- smodels$bic+k*(2+(2*k+2)/(23-k-1))-log(23)*k
mod.aic <- smodels$bic+k*2-log(23)*k
smodels 

#To find the best model, find the row of the models matrix where AIC c is the smallest. 
rmin <- which(mod.aicc==min(mod.aicc))
colnames(smodels$which)[smodels$which[rmin,]]

#In comparison, the best model with AIC is larger.
rmin <- which(mod.aic==min(mod.aic))
colnames(smodels$which)[smodels$which[rmin,]]

library(olsrr)
ols_step_best_subset(model.PMT.DT.T2)



## BEM model for May 24
data.BEM.DT.T2=reg.data.DT.T2[!is.na(reg.data.DT.T2$BEM),]

model.null = lm(BEM ~ 1, data=data.PMT.DT.T2)
model.full = lm(BEM ~OTU_111+OTU_6321+OTU_1000+OTU_10060+OTU_10454+OTU_10488+OTU_1651+OTU_17999+OTU_11759+OTU_380+OTU_17+OTU_87+OTU_6894+OTU_1831+OTU_2184+OTU_2339 +OTU_2399 +OTU_257+OTU_2641+OTU_455+bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+
                +AOA+AOB+nirk+nosZ+F.B.ratio
                +Pyruvic.Acid.Methyl.Ester  
                +D.Xylose  
                +Tween.40  
                +D.Mannitol  
                +N.Acetyl.D.Glucosamine  
                +D.Cellobiose  
                +a.Keto.Butyric.Acid  
                , data=data.BEM.DT.T2)

step(model.null, scope = list(upper=model.full,lower=model.null), k=4,direction="both", data=data.BEM.DT.T2)  

model.BEM.DT.T2 =lm(BEM ~ OTU_6321 + Pyruvic.Acid.Methyl.Ester + AOB + 
                      nirk + OTU_1000 ,data = data.BEM.DT.T2)

summary(model.BEM.DT.T2)



## Gluten for May 24
data.gluten.DT.T2=reg.data.DT.T2[!is.na(reg.data.DT.T2$Gluten),]

#Model with ASVS for gluten-MAY24
model.null = lm(Gluten ~ 1, data=data.gluten.DT.T2)
model.full = lm(Gluten ~ OTU_111+OTU_6321+OTU_1000+OTU_10060+OTU_10454+OTU_10488+OTU_1651+OTU_17999+OTU_11759+OTU_380+OTU_17+OTU_87+OTU_6894+OTU_1831+OTU_2184+OTU_2339 +OTU_2399 +OTU_257+OTU_2641+OTU_455+bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+F.B.ratio+
                        AOA+AOB+nirk+nosZ
                        +Pyruvic.Acid.Methyl.Ester  
                        +D.Xylose  
                        +Tween.40  
                        +D.Mannitol  
                        +N.Acetyl.D.Glucosamine  
                        +D.Cellobiose  
                        +a.Keto.Butyric.Acid  
                       ,data=data.gluten.DT.T2)

step(model.null, scope = list(upper=model.full), direction="both",K=7, data=data.gluten.DT.T2)  

model.gluten.DT.T2 = lm( Gluten ~ OTU_455 + D.Xylose + nosZ + Pyruvic.Acid.Methyl.Ester + 
                           AOB , data = data.gluten.DT.T2)
summary(model.gluten.DT.T2)

# Gluten for June7
data.gluten.DT.T3=reg.data.DT.T3[!is.na(reg.data.DT.T3$Gluten),]

#Model with ASVS for gluten-June7
model.null = lm(Gluten ~ 1, data=data.gluten.DT.T3)
model.full = lm(Gluten ~ OTU_10712+OTU_10782+OTU_10866+OTU_10957+OTU_2+OTU_2317+OTU_2341+OTU_2348+OTU_23745+OTU_43433+OTU_226+
                         OTU_2278+OTU_1029+OTU_1032+OTU_1040+OTU_1123+OTU_1148+OTU_1246+OTU_1278+OTU_1313+
                         bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                +B.Methyl.D.Glucoside 
                +y.Amino.Butyric.Acid 
                +Glycyl.L.Glutamic.Acid 
                +D.Cellobiose 
                +N.Acetyl.D.Glucosamine 
                +L.Threonine 
                +Tween.40 
                , data=data.gluten.DT.T3)

step(model.null, scope = list(upper=model.full), direction="both",   data=data.gluten.DT.T3)  

model.gluten.DT.T3 = lm(Gluten ~ Simpson.fun + ACE + N.Acetyl.D.Glucosamine + 
                          fun.axis2 + y.Amino.Butyric.Acid + OTU_43433 + bact.axis2, data = data.gluten.DT.T3)




summary(model.gluten.DT.T3)


#Model with ASVS for gluten-June21

data.protein.DT.T3=reg.data.DT.T3[!is.na(reg.data.DT.T3$Protein.grain),]


model.null = lm(Protein.grain ~ 1, data=data.protein.DT.T3)

model.full = lm(Protein.grain ~ OTU_10712+OTU_10782+OTU_10866+OTU_10957+OTU_2+OTU_2317+OTU_2341+OTU_2348+OTU_23745+OTU_43433+OTU_226+
                  OTU_2278+OTU_1029+OTU_1032+OTU_1040+OTU_1123+OTU_1148+OTU_1246+OTU_1278+OTU_1313+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                 +B.Methyl.D.Glucoside 
                 +y.Amino.Butyric.Acid 
                 +Glycyl.L.Glutamic.Acid 
                 +D.Cellobiose 
                 +N.Acetyl.D.Glucosamine 
                 +L.Threonine 
                 +Tween.40 
                , data=data.protein.DT.T3)

step(model.null, scope = list(upper=model.full), direction="both", K=3, data=data.protein.DT.T3)  


model.protein.DT.T3 = lm(Protein.grain ~ OTU_2 + Chao1 + ACE.fun + OTU_2278 + 
                           D.Cellobiose + Shannon.fun,  
                         data = data.protein.DT.T3)
summary(model.protein.DT.T3)


#Model for PMT-June21

data.PMT.DT.T3=reg.data.DT.T3[!is.na(reg.data.DT.T3$PMT),]


model.null = lm(PMT ~ 1, data=data.PMT.DT.T3)

model.full = lm(PMT ~ OTU_10712+OTU_10782+OTU_10866+OTU_10957+OTU_2+OTU_2317+OTU_2341+OTU_2348+OTU_23745+OTU_43433+OTU_226+
                  OTU_2278+OTU_1029+OTU_1032+OTU_1040+OTU_1123+OTU_1148+OTU_1246+OTU_1278+OTU_1313+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                +B.Methyl.D.Glucoside 
                +y.Amino.Butyric.Acid 
                +Glycyl.L.Glutamic.Acid 
                +D.Cellobiose 
                +N.Acetyl.D.Glucosamine 
                +L.Threonine 
                +Tween.40 
                , data=data.PMT.DT.T3)

step(model.null, scope = list(upper=model.full), direction="both", data=data.PMT.DT.T3)  


model.PMT.DT.T3 = lm(PMT ~ OTU_10957 + fun.axis2 + y.Amino.Butyric.Acid,  
                         data = data.PMT.DT.T3)
summary(model.PMT.DT.T3)

#Model for BEM-June21

data.BEM.DT.T3=reg.data.DT.T3[!is.na(reg.data.DT.T3$BEM),]


model.null = lm(BEM ~ 1, data=data.BEM.DT.T3)

model.full = lm(BEM ~ OTU_10712+OTU_10782+OTU_10866+OTU_10957+OTU_2+OTU_2317+OTU_2341+OTU_2348+OTU_23745+OTU_43433+OTU_226+
                  OTU_2278+OTU_1029+OTU_1032+OTU_1040+OTU_1123+OTU_1148+OTU_1246+OTU_1278+OTU_1313+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                  +B.Methyl.D.Glucoside 
                +y.Amino.Butyric.Acid 
                +Glycyl.L.Glutamic.Acid 
                +D.Cellobiose 
                +N.Acetyl.D.Glucosamine 
                +L.Threonine 
                +Tween.40 
                , data=data.BEM.DT.T3)

step(model.null, scope = list(upper=model.full), direction="both",K=4, data=data.BEM.DT.T3)  


model.BEM.DT.T3 = lm(PMT ~ OTU_10957+ OTU_2278+y.Amino.Butyric.Acid+D.Cellobiose+N.Acetyl.D.Glucosamine+
                     L.Threonine+Tween.40  ,  
                     data = data.BEM.DT.T3)
summary(model.BEM.DT.T3)

models.BEM <- leaps::regsubsets(BEM~OTU_10712+OTU_10782+OTU_10866+OTU_10957+OTU_2+OTU_2317+OTU_2341+OTU_2348+OTU_23745+OTU_43433+OTU_226+
                              OTU_2278+OTU_1029+OTU_1032+OTU_1040+OTU_1123+OTU_1148+OTU_1246+OTU_1278+OTU_1313+
                              bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                            +B.Methyl.D.Glucoside 
                            +y.Amino.Butyric.Acid 
                            +Glycyl.L.Glutamic.Acid 
                            +D.Cellobiose 
                            +N.Acetyl.D.Glucosamine 
                            +L.Threonine 
                            +Tween.40  , data=data.BEM.DT.T3,
                            nvmax = 11, nbest=3,
                            method = "forward")
smodels = summary(models.BEM)
k <- apply(smodels$which,1,sum)+1
mod.aicc <- smodels$bic+k*(2+(2*k+2)/(23-k-1))-log(23)*k
mod.aic <- smodels$bic+k*2-log(23)*k

#To find the best model, find the row of the models matrix where AIC c is the smallest. 
rmin <- which(mod.aicc==min(mod.aicc))
colnames(smodels$which)[smodels$which[rmin,]]

#In comparison, the best model with AIC is larger.
rmin <- which(mod.aic==min(mod.aic))
colnames(smodels$which)[smodels$which[rmin,]]

library(olsrr)
ols_step_best_subset(model.BEM.DT.T3)


# Gluten for July 5
data.gluten.DT.T4=reg.data.DT.T4[!is.na(reg.data.DT.T4$Gluten),]

#Model with ASVS for gluten-June21
model.null = lm(Gluten ~ 1, data=data.gluten.DT.T4)
model.full = lm(Gluten ~OTU_10000+OTU_25804+OTU_192+OTU_430+OTU_122+OTU_14774+OTU_14816+OTU_14926+OTU_1554+OTU_1605+OTU_10000+OTU_25804+
                       OTU_192+OTU_430+OTU_122+OTU_14774+OTU_14816+OTU_14926+OTU_1554+OTU_1605+bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+
                       AOA+AOB+nirk+nosZ+F.B.ratio+
                       Glucose.1..Phosphate+ 
                       Tween.80+ 
                       D.Mannitol+ 
                       X4.Hydroxy.Benzoic.Acid+ 
                       N.Acetyl.D.Glucosamine+ 
                       Itaconic.Acid+ 
                       L.Arginine 
                     ,data=data.gluten.DT.T4)

step(model.null, scope = list(upper=model.full), direction="both",  data=data.gluten.DT.T4)  

model.gluten.DT.T4 = lm(Gluten ~ OTU_10000 +nosZ , data = data.gluten.DT.T4)
summary(model.gluten.DT.T4)

# Gluten for July 5
data.protein.DT.T4=reg.data.DT.T4[!is.na(reg.data.DT.T4$Protein.grain),]

#Model with ASVS for protein
model.null = lm(Protein.grain ~ 1, data=data.protein.DT.T4)
model.full = lm(Protein.grain ~OTU_10000+OTU_25804+OTU_192+OTU_430+OTU_122+OTU_14774+OTU_14816+OTU_14926+OTU_1554+OTU_1605+OTU_10000+OTU_25804+
                  OTU_192+OTU_430+OTU_122+OTU_14774+OTU_14816+OTU_14926+OTU_1554+OTU_1605+bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun
                  +AOA+AOB+nirk+nosZ+F.B.ratio+
                  Glucose.1..Phosphate+ 
                  Tween.80+ 
                  D.Mannitol+ 
                  X4.Hydroxy.Benzoic.Acid+ 
                  N.Acetyl.D.Glucosamine+ 
                  Itaconic.Acid+ 
                  L.Arginine 
                , data=data.protein.DT.T4)

step(model.null, scope = list(upper=model.full), direction="both", K=6, data=data.protein.DT.T4)  

model.protein.DT.T4 = lm(Protein.grain ~ OTU_192 + OTU_25804 + Simpson + 
                           Shannon + Chao1 + OTU_10000 , data = data.protein.DT.T4)

summary(model.protein.DT.T4)


# Gluten for July 5
data.PMT.DT.T4=reg.data.DT.T4[!is.na(reg.data.DT.T4$PMT),]

#Model with ASVS for protein
model.null = lm(PMT ~ 1, data=data.PMT.DT.T4)
model.full = lm(PMT ~OTU_10000+OTU_25804+OTU_192+OTU_430+OTU_122+OTU_14774+OTU_14816+OTU_14926+OTU_1554+OTU_1605+OTU_10000+OTU_25804+
                  OTU_192+OTU_430+OTU_122+OTU_14774+OTU_14816+OTU_14926+OTU_1554+OTU_1605+bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun
                +AOA+AOB+nirk+nosZ+F.B.ratio+
                  Glucose.1..Phosphate+ 
                  Tween.80+ 
                  D.Mannitol+ 
                  X4.Hydroxy.Benzoic.Acid+ 
                  N.Acetyl.D.Glucosamine+ 
                  Itaconic.Acid+ 
                  L.Arginine 
                , data=data.PMT.DT.T4)

step(model.null, scope = list(upper=model.full), direction="both", K=3, data=data.PMT.DT.T4)  

model.PMT.DT.T4 = lm(PMT~ OTU_192 + OTU_25804 + Simpson + 
                           Shannon + Chao1 + OTU_10000 , data = data.PMT.DT.T4)

summary(model.PMT.DT.T4)

# Gluten for July 5
data.BEM.DT.T4=reg.data.DT.T4[!is.na(reg.data.DT.T4$BEM),]

#Model with ASVS for protein
model.null = lm(BEM ~ 1, data=data.BEM.DT.T4)
model.full = lm(BEM ~OTU_10000+OTU_25804+OTU_192+OTU_430+OTU_122+OTU_14774+OTU_14816+OTU_14926+OTU_1554+OTU_1605+OTU_10000+OTU_25804+
                  OTU_192+OTU_430+OTU_122+OTU_14774+OTU_14816+OTU_14926+OTU_1554+OTU_1605+bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun
                +AOA+AOB+nirk+nosZ+F.B.ratio+
                  Glucose.1..Phosphate+ 
                  Tween.80+ 
                  D.Mannitol+ 
                  X4.Hydroxy.Benzoic.Acid+ 
                  N.Acetyl.D.Glucosamine+ 
                  Itaconic.Acid+ 
                  L.Arginine 
                , data=data.BEM.DT.T4)

step(model.null, scope = list(upper=model.full), direction="both", data=data.BEM.DT.T4)  

model.BEM.DT.T4 = lm(BEM ~ OTU_25804 + OTU_10000 + OTU_122 + nosZ , data = data.BEM.DT.T4)

summary(model.BEM.DT.T4)


# Gluten for July5
data.gluten.DT.T5=reg.data.DT.T5[!is.na(reg.data.DT.T5$Gluten),]

#Model with ASVS for gluten-July5
model.null = lm(Gluten ~ 1, data=data.gluten.DT.T5)
model.full = lm(Gluten ~ OTU_1043+OTU_10631+OTU_12550+OTU_12560+OTU_11670+OTU_1304+OTU_1680+OTU_17011+OTU_18380+OTU_30617+OTU_16+OTU_2+OTU_27+OTU_338+OTU_4+OTU_5+OTU_52+OTU_62+OTU_6408+OTU_6576+
                bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                +D.Galacturonic.Acid 
                +L.Asparagine 
                +I.Erythritol 
                +LPhenylalanine 
                +Tween.80 
                +L.Arginine 
                +Pyruvic.Acid.Methyl.Ester 
                , data=data.gluten.DT.T5)

step(model.null, scope = list(upper=model.full), direction="both",K=6,  data=data.gluten.DT.T5)  

model.gluten.DT.T5 = lm(Gluten ~ OTU_52 + Simpson.fun + OTU_27 + OTU_338, data = data.gluten.DT.T5)
summary(model.gluten.DT.T5)


# Gluten for July5

data.protein.DT.T5=reg.data.DT.T5[!is.na(reg.data.DT.T5$Protein.grain),]

#Model with ASVS for protein
model.null = lm(Protein.grain ~ 1, data=data.protein.DT.T5)
model.full = lm(Protein.grain ~OTU_1043+OTU_10631+OTU_12550+OTU_12560+OTU_11670+OTU_1304+OTU_1680+OTU_17011+OTU_18380+OTU_30617+OTU_16+OTU_2+OTU_27+OTU_338+OTU_4+OTU_5+OTU_52+OTU_62+OTU_6408+OTU_6576+
                 bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                +D.Galacturonic.Acid 
                +L.Asparagine 
                +I.Erythritol 
                +LPhenylalanine 
                +Tween.80 
                +L.Arginine 
                +Pyruvic.Acid.Methyl.Ester 
                ,data=data.protein.DT.T5)

step(model.null, scope = list(upper=model.full), direction="both", K=6,  data=data.protein.DT.T5)  

model.protein.DT.T5 = lm(Protein.grain ~ I.Erythritol + D.Galacturonic.Acid + 
                           fun.axis2 + AOB, data = data.protein.DT.T5)


summary(model.protein.DT.T5)

# PMT for July5

data.PMT.DT.T5=reg.data.DT.T5[!is.na(reg.data.DT.T5$PMT),]

#Model with ASVS for protein
model.null = lm(PMT ~ 1, data=data.PMT.DT.T5)
model.full = lm(PMT ~OTU_1043+OTU_10631+OTU_12550+OTU_12560+OTU_11670+OTU_1304+OTU_1680+OTU_17011+OTU_18380+OTU_30617+OTU_16+OTU_2+OTU_27+OTU_338+OTU_4+OTU_5+OTU_52+OTU_62+OTU_6408+OTU_6576+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                  +D.Galacturonic.Acid 
                +L.Asparagine 
                +I.Erythritol 
                +LPhenylalanine 
                +Tween.80 
                +L.Arginine 
                +Pyruvic.Acid.Methyl.Ester 
                , data=data.PMT.DT.T5)

step(model.null, scope = list(upper=model.full), direction="both", K=6,  data=data.PMT.DT.T5)  

#model.PMT.DT.T5 = lm(PMT ~ , data = data.PMT.DT.T5)


#summary(model.PMT.DT.T5)

# PMT for July5

data.BEM.DT.T5=reg.data.DT.T5[!is.na(reg.data.DT.T5$BEM),]

#Model with ASVS for protein
model.null = lm(BEM ~ 1, data=data.BEM.DT.T5)
model.full = lm(BEM ~OTU_1043+OTU_10631+OTU_12550+OTU_12560+OTU_11670+OTU_1304+OTU_1680+OTU_17011+OTU_18380+OTU_30617+OTU_16+OTU_2+OTU_27+OTU_338+OTU_4+OTU_5+OTU_52+OTU_62+OTU_6408+OTU_6576+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                +D.Galacturonic.Acid 
                +L.Asparagine 
                +I.Erythritol 
                +LPhenylalanine 
                +Tween.80 
                +L.Arginine 
                +Pyruvic.Acid.Methyl.Ester 
                , data=data.BEM.DT.T5)

step(model.null, scope = list(upper=model.full), direction="both", K=4,  data=data.BEM.DT.T5)  

#model.BEM.DT.T5 = lm(BEM ~ , data = data.BEM.DT.T5)


#summary(model.BEM.DT.T5)



# Gluten for July19
data.gluten.DT.T6=reg.data.DT.T6[!is.na(reg.data.DT.T6$Gluten),]

#Model with ASVS for gluten-July19
model.null = lm(Gluten ~ 1, data=data.gluten.DT.T6)
model.full = lm(Gluten ~OTU_10365+OTU_1037+OTU_1044+OTU_43710+OTU_10687+OTU_11003+OTU_11321+OTU_1150+OTU_1078+OTU_1148+OTU_136+OTU_5749+OTU_1375+OTU_142+OTU_150+OTU_1568+OTU_181+OTU_276+OTU_2995+OTU_3537+
                bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                +B.Methyl.D.Glucoside 
                +L.Arginine 
                +D.Galacturonic.Acid 
                +L.Threonine 
                +a.Keto.Butyric.Acid 
                +Phenylethylamine 
                +L.Serine 
                , data=data.gluten.DT.T6)

step(model.null, scope = list(upper=model.full), direction="both", K=7, data=data.gluten.DT.T6)  

model.gluten.DT.T6 = lm(Gluten ~ OTU_276 + AOA + N.Acetyl.D.Glucosamine + 
                           y.Amino.Butyric.Acid + Chao1.fun, data = data.gluten.DT.T6)
summary(model.gluten.DT.T6)



# Gluten for July19

data.protein.DT.T6=reg.data.DT.T6[!is.na(reg.data.DT.T6$Protein.grain),]

#Model with ASVS for protein
model.null = lm(Protein.grain ~ 1, data=data.protein.DT.T6)
model.full = lm(Protein.grain ~OTU_10365+OTU_1037+OTU_1044+OTU_43710+OTU_10687+OTU_11003+OTU_11321+OTU_1150+OTU_1078+OTU_1148+OTU_136+OTU_5749+OTU_1375+OTU_142+OTU_150+OTU_1568+OTU_181+OTU_276+OTU_2995+OTU_3537+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                +B.Methyl.D.Glucoside 
                +L.Arginine 
                +D.Galacturonic.Acid 
                +L.Threonine 
                +a.Keto.Butyric.Acid 
                +Phenylethylamine 
                +L.Serine 
                , data=data.protein.DT.T6)

step(model.null, scope = list(upper=model.full), direction="both",k=6, data=data.protein.DT.T6)  

model.protein.DT.T6 = lm(Protein.grain ~ OTU_43710 + D.Galacturonic.Acid + 
                           fun.axis2 + L.Serine + Chao1, data = data.protein.DT.T6)

summary(model.protein.DT.T6)


# PMT for July19

data.PMT.DT.T6=reg.data.DT.T6[!is.na(reg.data.DT.T6$PMT),]

#Model with ASVS for protein
model.null = lm(PMT ~ 1, data=data.PMT.DT.T6)
model.full = lm(PMT ~OTU_10365+OTU_1037+OTU_1044+OTU_43710+OTU_10687+OTU_11003+OTU_11321+OTU_1150+OTU_1078+OTU_1148+OTU_136+OTU_5749+OTU_1375+OTU_142+OTU_150+OTU_1568+OTU_181+OTU_276+OTU_2995+OTU_3537+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                  +B.Methyl.D.Glucoside 
                +L.Arginine 
                +D.Galacturonic.Acid 
                +L.Threonine 
                +a.Keto.Butyric.Acid 
                +Phenylethylamine 
                +L.Serine 
                , data=data.PMT.DT.T6)

step(model.null, scope = list(upper=model.full), direction="both", K=6, data=data.PMT.DT.T6)  

#model.PMT.DT.T6 = lm(PMT ~ , data = data.PMT.DT.T6)

#summary(model.protein.DT.T6)

# BEM for July19

data.BEM.DT.T6=reg.data.DT.T6[!is.na(reg.data.DT.T6$BEM),]

#Model with ASVS for protein
model.null = lm(BEM ~ 1, data=data.BEM.DT.T6)
model.full = lm(BEM ~OTU_10365+OTU_1037+OTU_1044+OTU_43710+OTU_10687+OTU_11003+OTU_11321+OTU_1150+OTU_1078+OTU_1148+OTU_136+OTU_5749+OTU_1375+OTU_142+OTU_150+OTU_1568+OTU_181+OTU_276+OTU_2995+OTU_3537+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                  +B.Methyl.D.Glucoside 
                +L.Arginine 
                +D.Galacturonic.Acid 
                +L.Threonine 
                +a.Keto.Butyric.Acid 
                +Phenylethylamine 
                +L.Serine 
                , data=data.BEM.DT.T6)

step(model.null, scope = list(upper=model.full), direction="both", K=4, data=data.BEM.DT.T6)  

model.BEM.DT.T6 = lm(BEM ~ , data = data.BEM.DT.T6)

summary(model.BEM.DT.T6)


# Gluten for August
data.gluten.DT.T7=reg.data.DT.T7[!is.na(reg.data.DT.T7$Gluten),]

#Model with ASVS for August
model.null = lm(Gluten ~ 1, data=data.gluten.DT.T7)
model.full = lm(Gluten ~OTU_21928+OTU_1038+OTU_10633+OTU_23744+OTU_40119+OTU_620+OTU_11394+OTU_13293+OTU_14325+OTU_14336+OTU_6811+OTU_590+OTU_4+OTU_276+OTU_140+OTU_1491+OTU_1515+OTU_161+OTU_1640+OTU_1924+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                  +Pyruvic.Acid.Methyl.Ester 
                  +LPhenylalanine 
                  +Glycogen 
                  +Glycyl.L.Glutamic.Acid 
                  +a.Keto.Butyric.Acid 
                  +Putrescine 
                , data=data.gluten.DT.T7)

step(model.null, scope = list(upper=model.full), direction="both",  data=data.gluten.DT.T7)  

model.gluten.DT.T7 = lm( Gluten ~  OTU_1038 + bact.axis1 + Shannon, data = data.gluten.DT.T7)
summary(model.gluten.DT.T7)

# Protein for August
data.protein.DT.T7=reg.data.DT.T7[!is.na(reg.data.DT.T7$Protein.grain),]

#Model with ASVS for August
model.null = lm(Protein.grain~ 1, data=data.protein.DT.T7)
model.full = lm(Protein.grain~OTU_21928+OTU_1038+OTU_10633+OTU_23744+OTU_40119+OTU_620+OTU_11394+OTU_13293+OTU_14325+OTU_14336+OTU_6811+OTU_590+OTU_4+OTU_276+OTU_140+OTU_1491+OTU_1515+OTU_161+OTU_1640+OTU_1924+
                bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                +Pyruvic.Acid.Methyl.Ester 
                +LPhenylalanine 
                +Glycogen 
                +Glycyl.L.Glutamic.Acid 
                +a.Keto.Butyric.Acid 
                +Putrescine  , data=data.protein.DT.T7)

step(model.null, scope = list(upper=model.full), direction="both",K=4, data=data.protein.DT.T7)  

model.protein.DT.T7 = lm( Protein.grain ~ OTU_40119 + OTU_6811 + Glycyl.L.Glutamic.Acid + 
                            OTU_4 + Pyruvic.Acid.Methyl.Ester  , data = data.protein.DT.T7)
summary(model.protein.DT.T7)

# Protein for August
data.PMT.DT.T7=reg.data.DT.T7[!is.na(reg.data.DT.T7$PMT),]

#Model with ASVS for August
model.null = lm(PMT~ 1, data=data.PMT.DT.T7)
model.full = lm(PMT~OTU_21928+OTU_1038+OTU_10633+OTU_23744+OTU_40119+OTU_620+OTU_11394+OTU_13293+OTU_14325+OTU_14336+OTU_6811+OTU_590+OTU_4+OTU_276+OTU_140+OTU_1491+OTU_1515+OTU_161+OTU_1640+OTU_1924+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                  +Pyruvic.Acid.Methyl.Ester 
                +LPhenylalanine 
                +Glycogen 
                +Glycyl.L.Glutamic.Acid 
                +a.Keto.Butyric.Acid 
                +Putrescine  , data=data.PMT.DT.T7)

step(model.null, scope = list(upper=model.full), direction="both", data=data.PMT.DT.T7)  

model.PMT.DT.T7 = lm(PMT ~OTU_40119 + nosZ + ACE.fun + AOA + bact.axis1  , data = data.PMT.DT.T7)
summary(model.PMT.DT.T7)

# Protein for August
data.BEM.DT.T7=reg.data.DT.T7[!is.na(reg.data.DT.T7$BEM),]

#Model with ASVS for August
model.null = lm(BEM~ 1, data=data.PMT.DT.T7)
model.full = lm(BEM~OTU_21928+OTU_1038+OTU_10633+OTU_23744+OTU_40119+OTU_620+OTU_11394+OTU_13293+OTU_14325+OTU_14336+OTU_6811+OTU_590+OTU_4+OTU_276+OTU_140+OTU_1491+OTU_1515+OTU_161+OTU_1640+OTU_1924+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                +Pyruvic.Acid.Methyl.Ester 
                +LPhenylalanine 
                +Glycogen 
                +Glycyl.L.Glutamic.Acid 
                +a.Keto.Butyric.Acid 
                +Putrescine  , data=data.BEM.DT.T7)

step(model.null, scope = list(upper=model.full), direction="both",K=4, data=data.BEM.DT.T7)  

model.BEM.DT.T7 = lm(BEM ~Glycyl.L.Glutamic.Acid + OTU_140 + Putrescine + 
                       a.Keto.Butyric.Acid  , data = data.BEM.DT.T7)
summary(model.BEM.DT.T7)


## Modelling for Drought resistant cultivar
## Gluten for May 10
data.gluten.DS.T1=reg.data.DS.T1[!is.na(reg.data.DS.T1$Gluten),]

#Model with ASVS for gluten-MAY10
model.null = lm(Gluten ~ 1, data=data.gluten.DS.T1)
model.full = lm(Gluten~OTU_10214+OTU_10278+OTU_10378+OTU_11143+OTU_11397+OTU_11475+OTU_1183+OTU_12150+OTU_12183+OTU_12189+OTU_18+OTU_1077+OTU_1198+OTU_133+OTU_151+OTU_2199+OTU_330+OTU_3622+OTU_3898+OTU_638+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+LPhenylalanine +F.B.ratio 
                +Tween.80  
                +D.Mannitol  
                +Itaconic.Acid  
                +D.Cellobiose  
                +a.Keto.Butyric.Acid  
                +Phenylethylamine  
                , data=data.gluten.DS.T1)

step(model.null, scope = list(upper=model.full), direction="both", K=5,  data=data.gluten.DS.T1)  

model.gluten.DS.T1 = lm(Gluten ~ nirk + Tween.80 + D.Cellobiose + LPhenylalanine + 
                          Itaconic.Acid + Simpson.fun, data = data.gluten.DS.T1)
summary(model.gluten.DS.T1)

## Protein for May 10
data.protein.DS.T1=reg.data.DS.T1[!is.na(reg.data.DS.T1$Protein.grain),]

model.null = lm(Protein.grain ~ 1, data=data.protein.DS.T1)
model.full = lm(Protein.grain~OTU_10214+OTU_10278+OTU_10378+OTU_11143+OTU_11397+OTU_11475+OTU_1183+OTU_12150+OTU_12183+OTU_12189+OTU_18+OTU_1077+OTU_1198+OTU_133+OTU_151+OTU_2199+OTU_330+OTU_3622+OTU_3898+OTU_638+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio+LPhenylalanine+Tween.80  
                +D.Mannitol  
                +Itaconic.Acid  
                +D.Cellobiose  
                +a.Keto.Butyric.Acid  
                +Phenylethylamine  
                , data=data.protein.DS.T1)

step(model.null, scope = list(upper=model.full), direction="both",K=6,  data=data.protein.DS.T1)  

model.protein.DS.T1 = lm(Protein.grain ~ OTU_18 + LPhenylalanine + Itaconic.Acid + 
                           AOB + a.Keto.Butyric.Acid, data = data.protein.DS.T1)
summary(model.protein.DS.T1)


## Protein for May 10
data.PMT.DS.T1=reg.data.DS.T1[!is.na(reg.data.DS.T1$PMT),]

model.null = lm(PMT ~ 1, data=data.PMT.DS.T1)
model.full = lm(PMT~OTU_10214+OTU_10278+OTU_10378+OTU_11143+OTU_11397+OTU_11475+OTU_1183+OTU_12150+OTU_12183+OTU_12189+OTU_18+OTU_1077+OTU_1198+OTU_133+OTU_151+OTU_2199+OTU_330+OTU_3622+OTU_3898+OTU_638+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+LPhenylalanine+Tween.80+F.B.ratio  
                +D.Mannitol  
                +Itaconic.Acid  
                +D.Cellobiose  
                +a.Keto.Butyric.Acid  
                +Phenylethylamine  
                , data=data.PMT.DS.T1)

step(model.null, scope = list(upper=model.full), direction="both",K=5,  data=data.PMT.DS.T1)  

model.PMT.DS.T1 = lm(PMT~ OTU_18 + Itaconic.Acid + OTU_3898 + nosZ + 
                       nirk + LPhenylalanine + D.Cellobiose+F.B.ratio  , data = data.PMT.DS.T1)

summary(model.PMT.DS.T1)


## BEM for May 10
data.BEM.DS.T1=reg.data.DS.T1[!is.na(reg.data.DS.T1$BEM),]

model.null = lm(BEM ~ 1, data=data.BEM.DS.T1)
model.full = lm(BEM~OTU_10214+OTU_10278+OTU_10378+OTU_11143+OTU_11397+OTU_11475+OTU_1183+OTU_12150+OTU_12183+OTU_12189+OTU_18+OTU_1077+OTU_1198+OTU_133+OTU_151+OTU_2199+OTU_330+OTU_3622+OTU_3898+OTU_638+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+LPhenylalanine+Tween.80+F.B.ratio  
                +D.Mannitol  
                +Itaconic.Acid  
                +D.Cellobiose  
                +a.Keto.Butyric.Acid  
                +Phenylethylamine  
                , data=data.BEM.DS.T1)

step(model.null, scope = list(upper=model.full), direction="both",K=7,  data=data.BEM.DS.T1)  

#model.BEM.DS.T1 = lm(BEM~ , data = data.BEM.DS.T1)

#summary(model.BEM.DS.T1)



## Gluten for May 24
data.gluten.DS.T2=reg.data.DS.T2[!is.na(reg.data.DS.T2$Gluten),]

#Model with ASVS for gluten-MAY24
model.null = lm(Gluten ~ 1, data=data.gluten.DS.T2)
model.full = lm(Gluten~ OTU_10002+OTU_10233+OTU_10425+OTU_10948+OTU_1270+OTU_1325+OTU_1343+OTU_14681+OTU_16054+OTU_516+OTU_1184+OTU_1277+OTU_1331+OTU_1474+OTU_161+OTU_2760+OTU_2956+OTU_3382+OTU_3499+OTU_379+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                +B.Methyl.D.Glucoside  
                +L.Arginine  
                +X2.Hydroxy.Benzoic.Acid  
                +D.Mannitol  
                +Glycogen  
                +D.Malic.Acid  
                , data=data.gluten.DS.T2)

step(model.null, scope = list(upper=model.full), direction="both",  data=data.gluten.DS.T2)  

model.gluten.DS.T2 = lm(Gluten ~ OTU_516 + Glycogen + B.Methyl.D.Glucoside + 
                          Chao1.fun + OTU_161+ L.Arginine  , data = data.gluten.DS.T2)
summary(model.gluten.DS.T2)

## Protein for May 24
data.protein.DS.T2=reg.data.DS.T2[!is.na(reg.data.DS.T2$Protein.grain),]

#Model with ASVS for Protein-MAY24
model.null = lm(Protein.grain ~ 1, data=data.protein.DS.T2)
model.full = lm(Protein.grain~OTU_10002+OTU_10233+OTU_10425+OTU_10948+OTU_1270+OTU_1325+OTU_1343+OTU_14681+OTU_16054+OTU_516+OTU_1184+OTU_1277+OTU_1331+OTU_1474+OTU_161+OTU_2760+OTU_2956+OTU_3382+OTU_3499+OTU_379+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+D.Mannitol+F.B.ratio 	
                +B.Methyl.D.Glucoside  
                +L.Arginine  
                +X2.Hydroxy.Benzoic.Acid  
                +D.Mannitol  
                +Glycogen  
                +D.Malic.Acid, data=data.protein.DS.T2)

step(model.null, scope = list(upper=model.full), direction="both",K=7,  data=data.protein.DS.T2)  

model.protein.DS.T2 = lm(Protein.grain ~ OTU_516 + Chao1, data = data.protein.DS.T2)
summary(model.protein.DS.T2)

## PMT for May 24
data.PMT.DS.T2=reg.data.DS.T2[!is.na(reg.data.DS.T2$PMT),]

#Model with ASVS for Protein-MAY24
model.null = lm(PMT ~ 1, data=data.PMT.DS.T2)
model.full = lm(PMT~OTU_10002+OTU_10233+OTU_10425+OTU_10948+OTU_1270+OTU_1325+OTU_1343+OTU_14681+OTU_16054+OTU_516+OTU_1184+OTU_1277+OTU_1331+OTU_1474+OTU_161+OTU_2760+OTU_2956+OTU_3382+OTU_3499+OTU_379+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+D.Mannitol+F.B.ratio 	
                +B.Methyl.D.Glucoside  
                +L.Arginine  
                +X2.Hydroxy.Benzoic.Acid  
                +D.Mannitol  
                +Glycogen  
                +D.Malic.Acid, data=data.PMT.DS.T2)

step(model.null, scope = list(upper=model.full), direction="both",  data=data.PMT.DS.T2)  

model.PMT.DS.T2 = lm(PMT ~ OTU_1474 + X2.Hydroxy.Benzoic.Acid + fun.axis2 + 
                       F.B.ratio, data = data.PMT.DS.T2)
summary(model.PMT.DS.T2)

## PMT for May 24
data.BEM.DS.T2=reg.data.DS.T2[!is.na(reg.data.DS.T2$BEM),]

#Model with ASVS for Protein-MAY24
model.null = lm(BEM ~ 1, data=data.BEM.DS.T2)
model.full = lm(BEM~OTU_10002+OTU_10233+OTU_10425+OTU_10948+OTU_1270+OTU_1325+OTU_1343+OTU_14681+OTU_16054+OTU_516+OTU_1184+OTU_1277+OTU_1331+OTU_1474+OTU_161+OTU_2760+OTU_2956+OTU_3382+OTU_3499+OTU_379+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+D.Mannitol +F.B.ratio	
                +B.Methyl.D.Glucoside  
                +L.Arginine  
                +X2.Hydroxy.Benzoic.Acid  
                +D.Mannitol  
                +Glycogen  
                +D.Malic.Acid, data=data.BEM.DS.T2)

step(model.null, scope = list(upper=model.full), direction="both",K=4,  data=data.BEM.DS.T2)  

model.BEM.DS.T2 = lm(BEM ~ OTU_1474 + X2.Hydroxy.Benzoic.Acid + fun.axis2, data = data.BEM.DS.T2)
summary(model.BEM.DS.T2)


## Gluten for June 7
data.gluten.DS.T3=reg.data.DS.T3[!is.na(reg.data.DS.T3$Gluten),]

#Model with ASVS for gluten-June7
model.null = lm(Gluten ~ 1, data=data.gluten.DS.T3)
model.full = lm(Gluten~ OTU_10331+ OTU_10632+OTU_10742+OTU_13085+OTU_9449+OTU_11673+OTU_1174+OTU_1187+OTU_11878+OTU_12307+OTU_1088+OTU_1142+OTU_1330+OTU_1541+OTU_156+OTU_1703+OTU_2361+OTU_2557+OTU_3008+OTU_3172+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio 	
                +D.Galactonic.Acid.y.Lactone 
                +Pyruvic.Acid.Methyl.Ester 
                +L.Asparagine 
                +I.Erythritol 
                , data=data.gluten.DS.T3)

step(model.null, scope = list(upper=model.full), direction="both", k=4,  data=data.gluten.DS.T3)  

model.gluten.DS.T3 = lm(Gluten ~ OTU_156 + I.Erythritol + Pyruvic.Acid.Methyl.Ester + 
                          Simpson + fun.axis2 + nosZ
                          , data = data.gluten.DS.T3)
summary(model.gluten.DS.T3)

## Protein for June 7
data.protein.DS.T3=reg.data.DS.T3[!is.na(reg.data.DS.T3$Protein.grain),]

#Model with ASVS for Protein-MAY24
model.null = lm(Protein.grain ~ 1, data=data.protein.DS.T3)
model.full = lm(Protein.grain~OTU_10331+ OTU_10632+OTU_10742+OTU_13085+OTU_9449+OTU_11673+OTU_1174+OTU_1187+OTU_11878+OTU_12307+OTU_1088+OTU_1142+OTU_1330+OTU_1541+OTU_156+OTU_1703+OTU_2361+OTU_2557+OTU_3008+OTU_3172+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio	
                +D.Galactonic.Acid.y.Lactone 
                +Pyruvic.Acid.Methyl.Ester 
                +L.Asparagine 
                +I.Erythritol, data=data.protein.DS.T3)

step(model.null, scope = list(upper=model.full), direction="both",K=8,  data=data.protein.DS.T3)  

model.protein.DS.T3= lm(Protein.grain ~ Chao1 + OTU_156 + D.Galactonic.Acid.y.Lactone + 
                          ACE + Simpson + Shannon  + Simpson.fun + 
                          bact.axis1 + nirk + Pyruvic.Acid.Methyl.Ester + fun.axis2+ 
                          OTU_10331 + L.Asparagine + ACE.fun  
                          , data = data.protein.DS.T3)
summary(model.protein.DS.T3)

## Protein for June 7
data.PMT.DS.T3=reg.data.DS.T3[!is.na(reg.data.DS.T3$PMT),]

#Model with ASVS for Protein-MAY24
model.null = lm(PMT ~ 1, data=data.PMT.DS.T3)
model.full = lm(PMT~OTU_10331+ OTU_10632+OTU_10742+OTU_13085+OTU_9449+OTU_11673+OTU_1174+OTU_1187+OTU_11878+OTU_12307+OTU_1088+OTU_1142+OTU_1330+OTU_1541+OTU_156+OTU_1703+OTU_2361+OTU_2557+OTU_3008+OTU_3172+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio 	
                +D.Galactonic.Acid.y.Lactone 
                +Pyruvic.Acid.Methyl.Ester 
                +L.Asparagine 
                +I.Erythritol, data=data.PMT.DS.T3)

step(model.null, scope = list(upper=model.full), direction="both",  data=data.PMT.DS.T3)  

model.PMT.DS.T3= lm(PMT ~ OTU_10331 + bact.axis2 + ACE.fun + nirk, data = data.PMT.DS.T3)
summary(model.PMT.DS.T3)

## Protein for June 7
data.BEM.DS.T3=reg.data.DS.T3[!is.na(reg.data.DS.T3$BEM),]

#Model with ASVS for Protein-MAY24
model.null = lm(BEM ~ 1, data=data.BEM.DS.T3)
model.full = lm(BEM~OTU_10331+ OTU_10632+OTU_10742+OTU_13085+OTU_9449+OTU_11673+OTU_1174+OTU_1187+OTU_11878+OTU_12307+OTU_1088+OTU_1142+OTU_1330+OTU_1541+OTU_156+OTU_1703+OTU_2361+OTU_2557+OTU_3008+OTU_3172+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ +F.B.ratio	
                +D.Galactonic.Acid.y.Lactone 
                +Pyruvic.Acid.Methyl.Ester 
                +L.Asparagine 
                +I.Erythritol, data=data.BEM.DS.T3)

step(model.null, scope = list(upper=model.full), direction="both",K=4,  data=data.BEM.DS.T3)  

model.BEM.DS.T3= lm(BEM ~ OTU_1088 + fun.axis2 + Pyruvic.Acid.Methyl.Ester + 
                      ACE.fun, data = data.BEM.DS.T3)
summary(model.BEM.DS.T3)


#Model with ASVS for gluten-June21

data.gluten.DS.T4=reg.data.DS.T4[!is.na(reg.data.DS.T4$Gluten),]

#Model with ASVS for gluten-June7
model.null = lm(Gluten ~ 1, data=data.gluten.DS.T4)
model.full = lm(Gluten~OTU_10024+OTU_10214+OTU_10258+OTU_10298+OTU_10470+OTU_10497+OTU_10644+OTU_10676+OTU_10740+OTU_10981+OTU_105+OTU_1090+OTU_132+OTU_2597+OTU_4728+OTU_5266+OTU_589+OTU_6120+OTU_7030+OTU_813+
                bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio 	
                +D.Mannitol 
                +X4.Hydroxy.Benzoic.Acid 
                +D.Galactonic.Acid.y.Lactone 
                +L.Threonine 
                +a.Keto.Butyric.Acid 
                , data=data.gluten.DS.T4)

step(model.null, scope = list(upper=model.full), direction="both",  data=data.gluten.DS.T4)  

model.gluten.DS.T4 = lm(Gluten ~ OTU_105 + OTU_589 + Chao1.fun + bact.axis2, data = data.gluten.DS.T4)

summary(model.gluten.DS.T4)

#Model with ASVS for protein-June21

data.protein.DS.T4=reg.data.DS.T4[!is.na(reg.data.DS.T4$Protein.grain),]

#Model with ASVS for Protein
model.null = lm(Protein.grain ~ 1, data=data.protein.DS.T4)
model.full = lm(Protein.grain~OTU_10024+OTU_10214+OTU_10258+OTU_10298+OTU_10470+OTU_10497+OTU_10644+OTU_10676+OTU_10740+OTU_10981+OTU_105+OTU_1090+OTU_132+OTU_2597+OTU_4728+OTU_5266+OTU_589+OTU_6120+OTU_7030+OTU_813+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                +D.Mannitol 
                +X4.Hydroxy.Benzoic.Acid 
                +D.Galactonic.Acid.y.Lactone 
                +L.Threonine 
                +a.Keto.Butyric.Acid, data=data.protein.DS.T4)

step(model.null, scope = list(upper=model.full), direction="both",  data=data.protein.DS.T4)  

model.protein.DS.T4= lm( Protein.grain ~ D.Mannitol + D.Galactonic.Acid.y.Lactone + 
                           bact.axis1 + Simpson + fun.axis1 + OTU_589, data = data.protein.DS.T4)

summary(model.protein.DS.T4)


#Model for PMT-June21

data.PMT.DS.T4=reg.data.DS.T4[!is.na(reg.data.DS.T4$PMT),]

#Model with ASVS for Protein
model.null = lm(PMT ~ 1, data=data.PMT.DS.T4)
model.full = lm(PMT~OTU_10024+OTU_10214+OTU_10258+OTU_10298+OTU_10470+OTU_10497+OTU_10644+OTU_10676+OTU_10740+OTU_10981+OTU_105+OTU_1090+OTU_132+OTU_2597+OTU_4728+OTU_5266+OTU_589+OTU_6120+OTU_7030+OTU_813+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio 	
                +D.Mannitol 
                +X4.Hydroxy.Benzoic.Acid 
                +D.Galactonic.Acid.y.Lactone 
                +L.Threonine 
                +a.Keto.Butyric.Acid, data=data.PMT.DS.T4)

step(model.null, scope = list(upper=model.full), direction="both",K=6,  data=data.PMT.DS.T4)  



model.PMT.DS.T4= lm( PMT ~ D.Mannitol + D.Galactonic.Acid.y.Lactone + 
                           Chao1.fun + ACE.fun + bact.axis2 , data = data.PMT.DS.T4)

summary(model.PMT.DS.T4)


#Model for PMT-June21

data.BEM.DS.T4=reg.data.DS.T4[!is.na(reg.data.DS.T4$BEM),]

#Model with ASVS for Protein
model.null = lm(BEM ~ 1, data=data.BEM.DS.T4)
model.full = lm(BEM~OTU_10024+OTU_10214+OTU_10258+OTU_10298+OTU_10470+OTU_10497+OTU_10644+OTU_10676+OTU_10740+OTU_10981+OTU_105+OTU_1090+OTU_132+OTU_2597+OTU_4728+OTU_5266+OTU_589+OTU_6120+OTU_7030+OTU_813+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+ 	
                +D.Mannitol 
                +X4.Hydroxy.Benzoic.Acid 
                +D.Galactonic.Acid.y.Lactone 
                +L.Threonine 
                +a.Keto.Butyric.Acid, data=data.BEM.DS.T4)

step(model.null, scope = list(upper=model.full), direction="both",K=7,  data=data.BEM.DS.T4)  


model.BEM.DS.T4= lm(BEM ~ OTU_589 + Chao1.fun + bact.axis2 + ACE.fun + 
                      AOB + D.Galactonic.Acid.y.Lactone + Shannon + D.Mannitol  
                      , data = data.BEM.DS.T4)

summary(model.BEM.DS.T4)



#Model with ASVS for gluten-July5

data.gluten.DS.T5=reg.data.DS.T5[!is.na(reg.data.DS.T5$Gluten),]

#Model with ASVS for gluten-July5
model.null = lm(Gluten ~ 1, data=data.gluten.DS.T5)
model.full = lm(Gluten~OTU_10038+OTU_10282+OTU_10553+OTU_11058+OTU_11242+OTU_1143+OTU_11844+OTU_1217+OTU_1242+OTU_13629+OTU_1032+OTU_948+OTU_1416+OTU_1585+OTU_1617+OTU_1659+OTU_1965+OTU_2062+OTU_2133+OTU_2788+
                       bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+
                       AOA+AOB+nirk+nosZ + F.B.ratio	
                       +Pyruvic.Acid.Methyl.Ester 
                       +LPhenylalanine 
                       +a.Cyclodextrin 
                       +Glycogen 
                       +D.L.a.Glycerol.Phosphate 
                       +L.Threonine 
                       +L.Arginine 
                , data=data.gluten.DS.T5)

step(model.null, scope = list(upper=model.full), direction="both", k=3, data=data.gluten.DS.T5)  

model.gluten.DS.T5 = lm(Gluten ~ D.L.a.Glycerol.Phosphate + L.Arginine + 
                          Glycogen + F.B.ratio + Chao1.fun + fun.axis2 + fun.axis1 , data = data.gluten.DS.T5)

summary(model.gluten.DS.T5)

#Model with ASVS for protein-June21

data.protein.DS.T5=reg.data.DS.T5[!is.na(reg.data.DS.T5$Protein.grain),]

#Model with ASVS for Protein
model.null = lm(Protein.grain ~ 1, data=data.protein.DS.T5)
model.full = lm(Protein.grain~OTU_10038+OTU_10282+OTU_10553+OTU_11058+OTU_11242+OTU_1143+OTU_11844+OTU_1217+OTU_1242+OTU_13629+OTU_1032+OTU_948+OTU_1416+OTU_1585+OTU_1617+OTU_1659+OTU_1965+OTU_2062+OTU_2133+OTU_2788+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+
                  AOA+AOB+nirk+nosZ +	F.B.ratio
                  +Pyruvic.Acid.Methyl.Ester 
                 +LPhenylalanine 
                 +a.Cyclodextrin 
                 +Glycogen 
                 +D.L.a.Glycerol.Phosphate 
                 +L.Threonine 
                 +L.Arginine , data=data.protein.DS.T5)

step(model.null, scope = list(upper=model.full), direction="both",  data=data.protein.DS.T5)  

model.protein.DS.T5= lm( Protein.grain ~ D.L.a.Glycerol.Phosphate + LPhenylalanine + 
                           Glycogen + AOA + ACE.fun + Simpson.fun + fun.axis1 + bact.axis1 , data = data.protein.DS.T5)

summary(model.protein.DS.T5)


#Model with ASVS for PMT-June21

data.PMT.DS.T5=reg.data.DS.T5[!is.na(reg.data.DS.T5$PMT),]

#Model with ASVS for Protein
model.null = lm(PMT ~ 1, data=data.PMT.DS.T5)
model.full = lm(PMT~OTU_10038+OTU_10282+OTU_10553+OTU_11058+OTU_11242+OTU_1143+OTU_11844+OTU_1217+OTU_1242+OTU_13629+OTU_1032+OTU_948+OTU_1416+OTU_1585+OTU_1617+OTU_1659+OTU_1965+OTU_2062+OTU_2133+OTU_2788+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+
                  AOA+AOB+nirk+nosZ +	F.B.ratio
                  +Pyruvic.Acid.Methyl.Ester 
                +LPhenylalanine 
                +a.Cyclodextrin 
                +Glycogen 
                +D.L.a.Glycerol.Phosphate 
                +L.Threonine 
                +L.Arginine , data=data.PMT.DS.T5)

step(model.null, scope = list(upper=model.full), direction="both", k=3,  data=data.PMT.DS.T5)  

model.PMT.DS.T5= lm( PMT ~OTU_1965 + LPhenylalanine + L.Threonine  , data = data.PMT.DS.T5)

summary(model.PMT.DS.T5)


#Model with ASVS for BEM-June21

data.BEM.DS.T5=reg.data.DS.T5[!is.na(reg.data.DS.T5$BEM),]

#Model with ASVS for Protein
model.null = lm(BEM ~ 1, data=data.BEM.DS.T5)
model.full = lm(BEM~OTU_10038+OTU_10282+OTU_10553+OTU_11058+OTU_11242+OTU_1143+OTU_11844+OTU_1217+OTU_1242+OTU_13629+OTU_1032+OTU_948+OTU_1416+OTU_1585+OTU_1617+OTU_1659+OTU_1965+OTU_2062+OTU_2133+OTU_2788+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+
                  AOA+AOB+nirk+nosZ + F.B.ratio
                +Pyruvic.Acid.Methyl.Ester 
                +LPhenylalanine 
                +a.Cyclodextrin 
                +Glycogen 
                +D.L.a.Glycerol.Phosphate 
                +L.Threonine 
                +L.Arginine , data=data.BEM.DS.T5)

step(model.null, scope = list(upper=model.full), direction="both",  data=data.BEM.DS.T5)  

model.BEM.DS.T5= lm( BEM ~ OTU_1032 + a.Cyclodextrin + fun.axis1 + bact.axis1 + 
                       D.L.a.Glycerol.Phosphate , data = data.BEM.DS.T5)

summary(model.BEM.DS.T5)





#Model with ASVS for gluten-July19

data.gluten.DS.T6=reg.data.DS.T6[!is.na(reg.data.DS.T6$Gluten),]

#Model with ASVS for gluten-July19
model.null = lm(Gluten ~ 1, data=data.gluten.DS.T6)
model.full = lm(Gluten~ OTU_10050+OTU_1007+OTU_10416+OTU_1053+OTU_1067+OTU_1079+OTU_1082+OTU_1091+OTU_11025+OTU_1106+OTU_1236+OTU_1247+OTU_1250+OTU_303+OTU_94+OTU_1667+OTU_281+OTU_5084+OTU_5873+
                        bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                +D.Galacturonic.Acid 
                +LPhenylalanine 
                +X4.Hydroxy.Benzoic.Acid 
                +y.Amino.Butyric.Acid 
                +Glycogen 
                +Glucose.1..Phosphate 
                +Putrescine 
                , data=data.gluten.DS.T6)

step(model.null, scope = list(upper=model.full), direction="both", data=data.gluten.DS.T6)  

model.gluten.DS.T6 = lm(Gluten ~ Glucose.1..Phosphate + Glycogen + ACE + 
                          Shannon + fun.axis2 + F.B.ratio + nosZ + Chao1 + OTU_10050
                             , data = data.gluten.DS.T6)

summary(model.gluten.DS.T6)


#Model with ASVS for protein-Jul19

data.protein.DS.T6=reg.data.DS.T6[!is.na(reg.data.DS.T6$Protein.grain),]

#Model with ASVS for Protein
model.null = lm(Protein.grain ~ 1, data=data.protein.DS.T6)
model.full = lm(Protein.grain~OTU_10050+OTU_1007+OTU_10416+OTU_1053+OTU_1067+OTU_1079+OTU_1082+OTU_1091+OTU_11025+OTU_1106+OTU_1236+OTU_1247+OTU_1250+OTU_303+OTU_94+OTU_1667+OTU_281+OTU_5084+OTU_5873+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio	
                +D.Galacturonic.Acid 
                +LPhenylalanine 
                +X4.Hydroxy.Benzoic.Acid 
                +y.Amino.Butyric.Acid 
                +Glycogen 
                +Glucose.1..Phosphate 
                +Putrescine , data=data.protein.DS.T6)

step(model.null, scope = list(upper=model.full), direction="both",k=3, data=data.protein.DS.T6)  

model.protein.DS.T6= lm(  Protein.grain ~ bact.axis2   , data = data.protein.DS.T6)

summary(model.protein.DS.T6)


#Model with ASVS for PMT-Jul19

data.PMT.DS.T6=reg.data.DS.T6[!is.na(reg.data.DS.T6$PMT),]

#Model with ASVS for PMT
model.null = lm(PMT ~ 1, data=data.PMT.DS.T6)
model.full = lm(PMT~OTU_10050+OTU_1007+OTU_10416+OTU_1053+OTU_1067+OTU_1079+OTU_1082+OTU_1091+OTU_11025+OTU_1106+OTU_1236+OTU_1247+OTU_1250+OTU_303+OTU_94+OTU_1667+OTU_281+OTU_5084+OTU_5873+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio	
                  +D.Galacturonic.Acid 
                +LPhenylalanine 
                +X4.Hydroxy.Benzoic.Acid 
                +y.Amino.Butyric.Acid 
                +Glycogen 
                +Glucose.1..Phosphate 
                +Putrescine , data=data.PMT.DS.T6)

step(model.null, scope = list(upper=model.full), direction="both",  data=data.PMT.DS.T6)  

model.PMT.DS.T6= lm( PMT ~ OTU_10050 + Simpson.fun + Chao1.fun + ACE.fun + 
                       AOA + fun.axis2 + F.B.ratio + bact.axis2 + AOB   , data = data.PMT.DS.T6)

summary(model.PMT.DS.T6)


#Model with ASVS for BEM-Jul19

data.BEM.DS.T6=reg.data.DS.T6[!is.na(reg.data.DS.T6$BEM),]

#Model with ASVS for BEM
model.null = lm(BEM ~ 1, data=data.BEM.DS.T6)
model.full = lm(BEM~OTU_10050+OTU_1007+OTU_10416+OTU_1053+OTU_1067+OTU_1079+OTU_1082+OTU_1091+OTU_11025+OTU_1106+OTU_1236+OTU_1247+OTU_1250+OTU_303+OTU_94+OTU_1667+OTU_281+OTU_5084+OTU_5873+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio	
                  +D.Galacturonic.Acid 
                +LPhenylalanine 
                +X4.Hydroxy.Benzoic.Acid 
                +y.Amino.Butyric.Acid 
                +Glycogen 
                +Glucose.1..Phosphate 
                +Putrescine , data=data.BEM.DS.T6)

step(model.null, scope = list(upper=model.full), direction="both",  data=data.BEM.DS.T6)  

model.BEM.DS.T6= lm( BEM ~ OTU_10050 + ACE + fun.axis2 + F.B.ratio + 
                       bact.axis1 + Chao1 + AOA   , data = data.BEM.DS.T6)

summary(model.BEM.DS.T6)


#Model with ASVS for gluten-August1

data.gluten.DS.T7=reg.data.DS.T7[!is.na(reg.data.DS.T7$Gluten),]

#Model with ASVS for gluten-August1
model.null = lm(Gluten ~ 1, data=data.gluten.DS.T7)
model.full = lm(Gluten~ OTU_10371+OTU_10391+OTU_10409+OTU_12375+OTU_12465+OTU_13582+OTU_14604+OTU_15282+OTU_18314+OTU_18504+OTU_1694+OTU_2165+OTU_330+OTU_393+OTU_4025+OTU_4059+OTU_4372+OTU_5879+OTU_77+OTU_791+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOA+AOB+nirk+nosZ+F.B.ratio
                  +Phenylethylamine 
                  +a.Keto.Butyric.Acid 
                  +X4.Hydroxy.Benzoic.Acid 
                  +B.Methyl.D.Glucoside 
                  +Tween.40 
                 , data=data.gluten.DS.T7)

step(model.null, scope = list(upper=model.full), direction="both", K=6,  data=data.gluten.DS.T7)  

model.gluten.DS.T7 = lm( Gluten ~ AOA + Simpson + B.Methyl.D.Glucoside + 
                          nirk + fun.axis1 + OTU_5879, data = data.gluten.DS.T7)

summary(model.gluten.DS.T7)


#Model with ASVS for Protein-August1

data.protein.DS.T7=reg.data.DS.T7[!is.na(reg.data.DS.T7$Protein.grain),]

#Model with ASVS for Protein-August1
model.null = lm(Protein.grain ~ 1, data=data.protein.DS.T7)
model.full = lm(Protein.grain~ OTU_10371+OTU_10391+OTU_10409+OTU_12375+OTU_12465+OTU_13582+OTU_14604+OTU_15282+OTU_18314+OTU_18504+OTU_1694+OTU_2165+OTU_330+OTU_393+OTU_4025+OTU_4059+OTU_4372+OTU_5879+OTU_77+OTU_791+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOB+nirk+nosZ+F.B.ratio
                +Phenylethylamine 
                +a.Keto.Butyric.Acid 
                +X4.Hydroxy.Benzoic.Acid 
                +B.Methyl.D.Glucoside 
                +Tween.40 
                 , data=data.protein.DS.T7)

step(model.null, scope = list(upper=model.full), direction="both",  data=data.protein.DS.T7)  

model.protein.DS.T7 = lm(Protein.grain ~ OTU_12465 + OTU_393 + OTU_5879 + 
                           B.Methyl.D.Glucoside + a.Keto.Butyric.Acid  , data = data.protein.DS.T7)

summary(model.protein.DS.T7)

#Model with ASVS for PMT-August1

data.PMT.DS.T7=reg.data.DS.T7[!is.na(reg.data.DS.T7$PMT),]

#Model with ASVS for Protein-August1
model.null = lm(PMT ~ 1, data=data.PMT.DS.T7)
model.full = lm(PMT~ OTU_10371+OTU_10391+OTU_10409+OTU_12375+OTU_12465+OTU_13582+OTU_14604+OTU_15282+OTU_18314+OTU_18504+OTU_1694+OTU_2165+OTU_330+OTU_393+OTU_4025+OTU_4059+OTU_4372+OTU_5879+OTU_77+OTU_791+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOB+nirk+nosZ+F.B.ratio
                  +Phenylethylamine 
                +a.Keto.Butyric.Acid 
                +X4.Hydroxy.Benzoic.Acid 
                +B.Methyl.D.Glucoside 
                +Tween.40 
                , data=data.PMT.DS.T7)

step(model.null, scope = list(upper=model.full), direction="both",  data=data.PMT.DS.T7)  

model.PMT.DS.T7 = lm( PMT ~  OTU_10371 + OTU_5879 + Phenylethylamine + 
                        Chao1.fun + fun.axis2, data = data.PMT.DS.T7)

summary(model.PMT.DS.T7)


#Model with ASVS for BEM-August1

data.BEM.DS.T7=reg.data.DS.T7[!is.na(reg.data.DS.T7$BEM),]

#Model with ASVS for Protein-August1
model.null = lm(BEM~ 1, data=data.BEM.DS.T7)
model.full = lm(BEM~ OTU_10371+OTU_10391+OTU_10409+OTU_12375+OTU_12465+OTU_13582+OTU_14604+OTU_15282+OTU_18314+OTU_18504+OTU_1694+OTU_2165+OTU_330+OTU_393+OTU_4025+OTU_4059+OTU_4372+OTU_5879+OTU_77+OTU_791+
                  bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOB+nirk+nosZ+F.B.ratio
                  +Phenylethylamine 
                +a.Keto.Butyric.Acid 
                +X4.Hydroxy.Benzoic.Acid 
                +B.Methyl.D.Glucoside 
                +Tween.40 
                , data=data.BEM.DS.T7)

step(model.null, scope = list(upper=model.full), direction="both",k=3,  data=data.BEM.DS.T7)  

model.BEM.DS.T7 = lm( BEM ~OTU_12465 + nirk + OTU_393 + Chao1 + OTU_1694 , data = data.BEM.DS.T7)

summary(model.BEM.DS.T7)

car::vif(model.BEM.DS.T7)
olsrr::ols_vif_tol(model.BEM.DS.T7)

#Cross validation
library(leaps)
library(Hmisc)



predict.regsubsets <- function(object, newdata, id, ...) {
  form <- as.formula(object$call[[2]])
  mat <- model.matrix(form, newdata)
  coefi <- leaps:::coef.regsubsets(object, id = id)
  mat[, names(coefi)] %*% coefi
}



nfolds <- 5
nreps <- 10
folds <- matrix(NA, nreps, nrow(data.BEM.DS.T7))
for(i in 1:nreps) 
  folds[i,] <- sample(rep(1:nfolds, length = nrow(data.BEM.DS.T7)))


nvmax <- 7
cv.errors <- matrix(0, nreps*nfolds, nvmax)



for(r in 1:nreps){
  for (k in 1:nfolds) {
    traindat <- data.BEM.DS.T7[folds[r,]!=k,]
    testdat <- data.BEM.DS.T7[folds[r,]==k,]
    best.fit <- leaps::regsubsets(BEM ~OTU_10371+OTU_10391+OTU_10409+OTU_12375+OTU_12465+OTU_13582+OTU_14604+OTU_15282+OTU_18314+OTU_18504+OTU_1694+OTU_2165+OTU_330+OTU_393+OTU_4025+OTU_4059+OTU_4372+OTU_5879+OTU_77+OTU_791+
                                    bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOB+nirk+nosZ+
                                    +Phenylethylamine 
                                  +a.Keto.Butyric.Acid 
                                  +X4.Hydroxy.Benzoic.Acid 
                                  +B.Methyl.D.Glucoside 
                                  +Tween.40, data=traindat, nvmax = nvmax, method = "seqrep")
    for (i in 1:nvmax) {
      pred <- predict.regsubsets(best.fit, testdat, id = i)
      cv.errors[r+(k-1)*nreps, i] <- 
        mean((testdat$BEM - pred)^2)
    }
  }
}
rmse.cv <- sqrt(apply(cv.errors, 2, mean, na.rm=TRUE))
graphics.off()
par("mar")
par(mar=c(5,5,5,5))
plot(1:nvmax, rmse.cv, pch = 19, type = "b",xlab="Number of Variables", ylab="RMSE")


# 23-fold cross-validation (Leave One Out CV)

best.fit <- leaps::regsubsets(BEM ~ OTU_10371+OTU_10391+OTU_10409+OTU_12375+OTU_12465+OTU_13582+OTU_14604+OTU_15282+OTU_18314+OTU_18504+OTU_1694+OTU_2165+OTU_330+OTU_393+OTU_4025+OTU_4059+OTU_4372+OTU_5879+OTU_77+OTU_791+
                                bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOB+nirk+nosZ+
                              +Phenylethylamine 
                              +a.Keto.Butyric.Acid 
                              +X4.Hydroxy.Benzoic.Acid 
                              +B.Methyl.D.Glucoside 
                              +Tween.40, data=traindat, nvmax = 5, method = "seqrep")
tmp <- summary(best.fit)$which
colnames(tmp)[tmp[5,]]

#Cross-validation with caret package
library(caret)
# Set up repeated k-fold cross-validation
train.control <- trainControl(method = "repeatedcv", number=5, repeats=20)
# Train the model
step.model <- train(BEM ~ OTU_10371+OTU_10391+OTU_10409+OTU_12375+OTU_12465+OTU_13582+OTU_14604+OTU_15282+OTU_18314+OTU_18504+OTU_1694+OTU_2165+OTU_330+OTU_393+OTU_4025+OTU_4059+OTU_4372+OTU_5879+OTU_77+OTU_791+
                      bact.axis1+bact.axis2+fun.axis1+fun.axis2+Shannon+Simpson+Chao1+ACE+Shannon.fun+Simpson.fun+Chao1.fun+ACE.fun+AOB+nirk+nosZ+
                      +Phenylethylamine 
                    +a.Keto.Butyric.Acid 
                    +X4.Hydroxy.Benzoic.Acid 
                    +B.Methyl.D.Glucoside 
                    +Tween.40, data=reg.data.DS.T7,
                    method = "leapSeq",
                    tuneGrid = data.frame(nvmax = 1:nvmax),
                    trControl = train.control)
plot(step.model$results$RMSE, pch = 19, type = "b", ylab="RMSE")

summary(warnings())

step.model$results

coef(step.model$finalModel, id=2)


dev.off()

