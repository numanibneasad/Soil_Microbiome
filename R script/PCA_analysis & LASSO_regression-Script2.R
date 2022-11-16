#Author: Numman Ibne Asad
#2022/10/01
#Project (Microbial based predictive modeling)


#PCA analysis of microbial OTUs among all the sample dates
#Principle component analysis
library("FactoMineR")
library("factoextra")
library(tidyverse)
library(caret)
library(glmnet)
library(vip)
library(tidyverse)
library(modelr)
library(broom)
library(vegan)
library(corrplot)
library(psych)
require(ggpubr)
require(tidyverse)
require(Hmisc)
library(ggplot2)
library(scales)
library(corrplot)
library(psych)
library("corrplot")

#Bacterial OTU transformation by Hellinger

#Cultivar DT (Bacteria)
otu.16S.DT.T1.hel=decostand(otu.16S.DT.T1, method="hellinger",na.rm = TRUE)
otu.16S.DT.T2.hel=decostand(otu.16S.DT.T2, method="hellinger",na.rm = TRUE)
otu.16S.DT.T3.hel=decostand(otu.16S.DT.T3, method="hellinger",na.rm = TRUE)
otu.16S.DT.T4.hel=decostand(otu.16S.DT.T4, method="hellinger",na.rm = TRUE)
otu.16S.DT.T5.hel=decostand(otu.16S.DT.T5, method="hellinger",na.rm = TRUE)
otu.16S.DT.T6.hel=decostand(otu.16S.DT.T6, method="hellinger",na.rm = TRUE)
otu.16S.DT.T7.hel=decostand(otu.16S.DT.T7, method="hellinger",na.rm = TRUE)


#Cultivar DS (Bacteria)
otu.16S.DS.T1.hel=decostand(otu.16S.DS.T1, method="hellinger",na.rm = TRUE)
otu.16S.DS.T2.hel=decostand(otu.16S.DS.T2, method="hellinger",na.rm = TRUE)
otu.16S.DS.T3.hel=decostand(otu.16S.DS.T3, method="hellinger",na.rm = TRUE)
otu.16S.DS.T4.hel=decostand(otu.16S.DS.T4, method="hellinger",na.rm = TRUE)
otu.16S.DS.T5.hel=decostand(otu.16S.DS.T5, method="hellinger",na.rm = TRUE)
otu.16S.DS.T6.hel=decostand(otu.16S.DS.T6, method="hellinger",na.rm = TRUE)
otu.16S.DS.T7.hel=decostand(otu.16S.DS.T7, method="hellinger",na.rm = TRUE)

# Fungal OTUs transformation by Hellinger

# Cultivar DT (Fungi)
otu.ITS.DT.T1.hel=decostand(otu.ITS.DT.T1, method="hellinger",na.rm = TRUE)
otu.ITS.DT.T2.hel=decostand(otu.ITS.DT.T2, method="hellinger",na.rm = TRUE)
otu.ITS.DT.T3.hel=decostand(otu.ITS.DT.T3, method="hellinger",na.rm = TRUE)
otu.ITS.DT.T4.hel=decostand(otu.ITS.DT.T4, method="hellinger",na.rm = TRUE)
otu.ITS.DT.T5.hel=decostand(otu.ITS.DT.T5, method="hellinger",na.rm = TRUE)
otu.ITS.DT.T6.hel=decostand(otu.ITS.DT.T6, method="hellinger",na.rm = TRUE)
otu.ITS.DT.T7.hel=decostand(otu.ITS.DT.T7, method="hellinger",na.rm = TRUE)

#Cultivar DS (Fungi)
otu.ITS.DS.T1.hel=decostand(otu.ITS.DS.T1, method="hellinger",na.rm = TRUE)
otu.ITS.DS.T2.hel=decostand(otu.ITS.DS.T2, method="hellinger",na.rm = TRUE)
otu.ITS.DS.T3.hel=decostand(otu.ITS.DS.T3, method="hellinger",na.rm = TRUE)
otu.ITS.DS.T4.hel=decostand(otu.ITS.DS.T4, method="hellinger",na.rm = TRUE)
otu.ITS.DS.T5.hel=decostand(otu.ITS.DS.T5, method="hellinger",na.rm = TRUE)
otu.ITS.DS.T6.hel=decostand(otu.ITS.DS.T6, method="hellinger",na.rm = TRUE)
otu.ITS.DS.T7.hel=decostand(otu.ITS.DS.T7, method="hellinger",na.rm = TRUE)

#Cultivar DT (Biolog)
biolog.DT.T1.hel<-decostand(biolog.DT.T1,method="hellinger", na.rm=TRUE)
biolog.DT.T2.hel<-decostand(biolog.DT.T2,method="hellinger", na.rm=TRUE)

#Biolog.DT.T2.hel<-na.omit(biolog.DT.T2.hel)
biolog.DT.T3.hel<-decostand(biolog.DT.T3,method="hellinger", na.rm=TRUE)
biolog.DT.T4.hel<-decostand(biolog.DT.T4,method="hellinger", na.rm=TRUE)
biolog.DT.T5.hel<-decostand(biolog.DT.T5,method="hellinger", na.rm=TRUE)
biolog.DT.T6.hel<-decostand(biolog.DT.T6,method="hellinger", na.rm=TRUE)
biolog.DT.T7.hel<-decostand(biolog.DT.T7,method="hellinger", na.rm=TRUE)

#Cultivar DS(Biolog)
biolog.DS.T1.hel<-decostand(biolog.DS.T1,method="hellinger", na.rm=TRUE)
biolog.DS.T2.hel<-decostand(biolog.DS.T2,method="hellinger", na.rm=TRUE)
biolog.DS.T3.hel<-decostand(biolog.DS.T3,method="hellinger", na.rm=TRUE)
biolog.DS.T4.hel<-decostand(biolog.DS.T4,method="hellinger", na.rm=TRUE)
biolog.DS.T5.hel<-decostand(biolog.DS.T5,method="hellinger", na.rm=TRUE)
biolog.DS.T6.hel<-decostand(biolog.DS.T6,method="hellinger", na.rm=TRUE)
biolog.DS.T7.hel<-decostand(biolog.DT.T7,method="hellinger", na.rm=TRUE)




#PCA analysis for 10th May
pca.16S.DT.T1 <- PCA(otu.16S.DT.T1.hel, graph = FALSE)
print(pca.16S.DT.T1)

##Cumulative percentage of eigenvalues
fviz_eig(pca.16S.DT.T1, addlabels = TRUE, ylim = c(0,7.5))
var.16S.DT.T1 <- get_pca_var(pca.16S.DT.T1)
var.16S.DT.T1$contrib
# Extraction of PCA points
pca.point.16S.DT.T1=data.frame(pca.16S.DT.T1$ind)
pca.data.point.16S.DT.T1=pca.point.16S.DT.T1[,which(colnames(pca.point.16S.DT.T1) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.16S.DT.T1)[which(colnames(pca.data.point.16S.DT.T1) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("bact.DT.T1.axis1", "bact.DT.T1.axis2","bact.DT.T1.axis3","bact.DT.T1.axis4","bact.DT.T1.axis5") 

#Analysis of significant dimensions

res.desc.DT.T1 <- dimdesc(pca.16S.DT.T1, axes = c(1,2,3,4,5), proba = 0.05)
# Description of dimension 1
res.desc.DT.T1$Dim.1
# Description of dimension 2
res.desc.DT.T1$Dim.2
# Contributions of variables to PC1
fviz_contrib(pca.16S.DT.T1, choice = "var", axes = 1, top = 5)
# Contributions of variables to PC2
fviz_contrib(pca.16S.DT.T1, choice = "var", axes = 2, top = 5)
#Total contribution to PC1 and PC2
fviz_contrib(pca.16S.DT.T1, choice = "var", axes = 1:2, top = 10)

# Description of dimension 2
pca.axis1.DT.T1<-data.frame(res.desc.DT.T1$Dim.1)
pca.axis2.DT.T1<-data.frame(res.desc.DT.T1$Dim.2)
pca.axis3.DT.T1<-data.frame(res.desc.DT.T1$Dim.3)
pca.axis4.DT.T1<-data.frame(res.desc.DT.T1$Dim.4)
pca.axis5.DT.T1<-data.frame(res.desc.DT.T1$Dim.5)

res.desc.DT.T1 <- dimdesc(pca.16S.DT.T1, axes = c(1,2,3,4,5), proba = 0.05)
write.table(pca.axis1.DT.T1,here("output/tables","PCA.axis1.DT.T1.csv"),  sep="\t")
write.table(pca.axis2.DT.T1,here("output/tables","PCA.axis2.DT.T1.csv"),  sep="\t")
write.table(pca.axis3.DT.T1,here("output/tables","PCA.axis3.DT.T1.csv"),  sep="\t")
write.table(pca.axis4.DT.T1,here("output/tables","PCA.axis4.DT.T1.csv"),  sep="\t")
write.table(pca.axis5.DT.T1,here("output/tables","PCA.axis5.DT.T1.csv"),  sep="\t")



bact.DT.T1.cor <- corr.test(qual.DT.S1,pca.data.point.16S.DT.T1, method = 'spearman')

myfun <- function(bact.DT.T1.cor){
  corrplot(bact.DT.T1.cor$r, 
           p.mat = round(as.matrix(bact.DT.T1.cor$p),3),
           method = 'circle',
           type = 'lower',
           sig.level = c(.001, .01, .05), 
           tl.pos="lt", 
           tl.col="black", tl.cex=1.0, 
           tl.offset=0.2,
           cl.pos="r",
           insig = "label_sig",
           pch.cex = 1.0,
           pch.col="red",
           cl.cex = 1.0)
  corrplot(bact.DT.T1.cor$r,  type="upper", method="number",
           col="coral4",  tl.pos="n", cl.pos="n", number.cex = 0.9, add=T,diag=F)
  recordPlot() # record the latest plot
}


myfun(bact.DT.T1.cor) 




#PCA analysis for 24th May
pca.16S.DT.T2 <- PCA(otu.16S.DT.T2.hel, graph = FALSE)
print(pca.16S.DT.T2)

##Cumulative percentage of eigenvalues
fviz_eig(pca.16S.DT.T2, addlabels = TRUE, ylim = c(0, 10))
var.16S.DT.T2 <- get_pca_var(pca.16S.DT.T2)

# Extraction of PCA points
pca.point.16S.DT.T2=data.frame(pca.16S.DT.T2$ind)
pca.data.point.16S.DT.T2=pca.point.16S.DT.T2[,which(colnames(pca.point.16S.DT.T2) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.16S.DT.T2)[which(colnames(pca.data.point.16S.DT.T2) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("bact.DT.T2.axis1", "bact.DT.T2.axis2","bact.DT.T2.axis3","bact.DT.T2.axis4","bact.DT.T2.axis5") 

#Analysis of significant dimensions

res.desc.DT.T2 <- dimdesc(pca.16S.DT.T2, axes = c(2), proba = 0.005)

# Description of dimension 3
pca.axis2.DT.T2<-data.frame(res.desc.DT.T2$Dim.2)

write.table(pca.axis2.DT.T2,here("output/tables","PCA.axis2.DT.T2.txt"),sep="\t")


#PCA analysis for 7th June
pca.16S.DT.T3 <- PCA(otu.16S.DT.T3.hel, graph = FALSE)
print(pca.16S.DT.T3)

##Cumulative percentage of eigenvalues
fviz_eig(pca.16S.DT.T3, addlabels = TRUE, ylim = c(0, 7.5))
var.16S.DT.T3 <- get_pca_var(pca.16S.DT.T3)

# Extraction of PCA points
pca.point.16S.DT.T3=data.frame(pca.16S.DT.T3$ind)
pca.data.point.16S.DT.T3=pca.point.16S.DT.T3[,which(colnames(pca.point.16S.DT.T3) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.16S.DT.T3)[which(colnames(pca.data.point.16S.DT.T3) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("bact.DT.T3.axis1", "bact.DT.T3.axis2","bact.DT.T3.axis3","bact.DT.T3.axis4","bact.DT.T3.axis5") 


# Correlation  between quality and bacterial PCA

bact.DT.T3.cor <- corr.test(qual.DT.S3,pca.data.point.16S.DT.T3, method = 'spearman')

myfun <- function(bact.DT.T3.cor){
  corrplot(bact.DT.T3.cor$r, 
           p.mat = round(as.matrix(bact.DT.T3.cor$p),3),
           method = 'circle',
           type = 'lower',
           sig.level = c(.001, .01, .05), 
           tl.pos="lt", 
           tl.col="black", tl.cex=1.0, 
           tl.offset=0.2,
           cl.pos="r",
           insig = "label_sig",
           pch.cex = 1.0,
           pch.col="red",
           cl.cex = 1.0)
  corrplot(bact.DT.T3.cor$r,  type="upper", method="number",
           col="coral4",  tl.pos="n", cl.pos="n", number.cex = 0.9, add=T,diag=F)
  recordPlot() # record the latest plot
}


myfun (bact.DT.T3.cor )

#Analysis of significant dimensions

res.desc.DT.T3 <- dimdesc(pca.16S.DT.T3, axes = c(1,2,3,4,5), proba = 0.05)

# Description of dimension 3
pca.axis1.DT.T3<-data.frame(res.desc.DT.T3$Dim.1)
pca.axis2.DT.T3<-data.frame(res.desc.DT.T3$Dim.2)
pca.axis3.DT.T3<-data.frame(res.desc.DT.T3$Dim.3)
pca.axis4.DT.T3<-data.frame(res.desc.DT.T3$Dim.4)
pca.axis5.DT.T3<-data.frame(res.desc.DT.T3$Dim.5)

write.table(pca.axis1.DT.T3,here("output/tables","PCA.axis1.16S.DT.T3.txt"),sep="\t")
write.table(pca.axis2.DT.T3,here("output/tables","PCA.axis2.16S.DT.T3.txt"),sep="\t")
write.table(pca.axis3.DT.T3,here("output/tables","PCA.axis3.16S.DT.T3.txt"),sep="\t")
write.table(pca.axis4.DT.T3,here("output/tables","PCA.axis4.16S.DT.T3.txt"),sep="\t")
write.table(pca.axis5.DT.T3,here("output/tables","PCA.axis5.16S.DT.T3.txt"),sep="\t")


# Description of dimension 2
res.desc.DT.T3$Dim.2
# Contributions of variables to PC1
fviz_contrib(pca.16S.DT.T3, choice = "var", axes = 1, top = 5)
# Contributions of variables to PC2
fviz_contrib(pca.16S.DT.T3, choice = "var", axes = 3, top = 25)
#Total contribution to PC1 and PC2
fviz_contrib(pca.16S.DT.T3, choice = "var", axes = 1:2, top = 10)



#PCA analysis for 24th June
pca.16S.DT.T4 <- PCA(otu.16S.DT.T4.hel, graph = FALSE)
print(pca.16S.DT.T4)

##Cumulative percentage of eigenvalues
fviz_eig(pca.16S.DT.T4, addlabels = TRUE, ylim = c(0, 8))
var.16S.DT.T4 <- get_pca_var(pca.16S.DT.T4)

# Extraction of PCA points
pca.point.16S.DT.T4=data.frame(pca.16S.DT.T4$ind)
pca.data.point.16S.DT.T4=pca.point.16S.DT.T4[,which(colnames(pca.point.16S.DT.T4) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.16S.DT.T4)[which(colnames(pca.data.point.16S.DT.T4) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("bact.DT.T4.axis1", "bact.DT.T4.axis2","bact.DT.T4.axis3","bact.DT.T4.axis4","bact.DT.T4.axis5") 


bact.DT.T4.cor <- corr.test(qual.DT.S4,pca.data.point.16S.DT.T4, method = 'spearman')
myfun <- function(bact.DT.T4.cor){
  corrplot(bact.DT.T4.cor$r, 
           p.mat = round(as.matrix(bact.DT.T4.cor$p),3),
           method = 'circle',
           type = 'lower',
           sig.level = c(.001, .01, .05), 
           tl.pos="lt", 
           tl.col="black", tl.cex=1.0, 
           tl.offset=0.2,
           cl.pos="r",
           insig = "label_sig",
           pch.cex = 1.0,
           pch.col="red",
           cl.cex = 1.0)
  corrplot(bact.DT.T4.cor$r,  type="upper", method="number",
           col="coral4",  tl.pos="n", cl.pos="n", number.cex = 0.9, add=T,diag=F)
  recordPlot() # record the latest plot
}


myfun(bact.DT.T4.cor) 




#PCA analysis for 7th July
pca.16S.DT.T5 <- PCA(otu.16S.DT.T5.hel, graph = FALSE)
print(pca.16S.DT.T5)

##Cumulative percentage of eigenvalues
fviz_eig(pca.16S.DT.T5, addlabels = TRUE, ylim = c(0, 10))
var.16S.DT.T5 <- get_pca_var(pca.16S.DT.T5)

# Extraction of PCA points
pca.point.16S.DT.T5=data.frame(pca.16S.DT.T5$ind)
pca.data.point.16S.DT.T5=pca.point.16S.DT.T5[,which(colnames(pca.point.16S.DT.T5) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.16S.DT.T5)[which(colnames(pca.data.point.16S.DT.T5) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("bact.DT.T5.axis1", "bact.DT.T5.axis2","bact.DT.T5.axis3","bact.DT.T5.axis4","bact.DT.T5.axis5") 

#PCA analysis for 19th July
pca.16S.DT.T6 <- PCA(otu.16S.DT.T6.hel, graph = FALSE)
print(pca.16S.DT.T5)

##Cumulative percentage of eigenvalues
fviz_eig(pca.16S.DT.T6, addlabels = TRUE, ylim = c(0, 10))
var.16S.DT.T6 <- get_pca_var(pca.16S.DT.T6)

# Extraction of PCA points
pca.point.16S.DT.T6=data.frame(pca.16S.DT.T6$ind)
pca.data.point.16S.DT.T6=pca.point.16S.DT.T6[,which(colnames(pca.point.16S.DT.T6) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.16S.DT.T6)[which(colnames(pca.data.point.16S.DT.T6) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("bact.DT.T6.axis1", "bact.DT.T6.axis2","bact.DT.T6.axis3","bact.DT.T6.axis4","bact.DT.T6.axis5") 

#PCA analysis for 1Aug
pca.16S.DT.T7 <- PCA(otu.16S.DT.T7.hel, graph = FALSE)
print(pca.16S.DT.T7)

##Cumulative percentage of eigenvalues
fviz_eig(pca.16S.DT.T7, addlabels = TRUE, ylim = c(0, 10))
var.16S.DT.T7 <- get_pca_var(pca.16S.DT.T7)

# Extraction of PCA points
pca.point.16S.DT.T7=data.frame(pca.16S.DT.T7$ind)
pca.data.point.16S.DT.T7=pca.point.16S.DT.T7[,which(colnames(pca.point.16S.DT.T7) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.16S.DT.T7)[which(colnames(pca.data.point.16S.DT.T7) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("bact.DT.T7.axis1", "bact.DT.T7.axis2","bact.DT.T7.axis3","bact.DT.T7.axis4","bact.DT.T7.axis5") 

res.desc.DT.T7 <- dimdesc(pca.16S.DT.T7, axes = c(1,2,3,4,5), proba = 0.05)

# Description of dimension 2
pca.axis1.DT.T7<-data.frame(res.desc.DT.T7$Dim.1)
pca.axis4.DT.T7<-data.frame(res.desc.DT.T7$Dim.4)
pca.axis5.DT.T7<-data.frame(res.desc.DT.T7$Dim.5)


#write.table(pca.axis1.DT.T7,"PCA.axis1.DT.T7.csv",  sep="\t")
#write.table(pca.axis4.DT.T7,"PCA.axis4.DT.T7.csv",  sep="\t")
#write.table(pca.axis5.DT.T7,"PCA.axis5.DT.T7.csv",  sep="\t")



# Fungal OTUs for PCA analysis

#PCA run for 10 MAY
pca.ITS.DT.T1 <- PCA(otu.ITS.DT.T1.hel, graph = FALSE)
print(pca.ITS.DT.T1)

#Visualization and interpretation of eigenvalues
eig.val.ITS.DT.T1 <- get_eigenvalue(pca.ITS.DT.T1)
eig.val.ITS.DT.T1

#Cumulative percentage of eigenvalues
fviz_eig(pca.ITS.DT.T1, addlabels = TRUE, ylim = c(0, 8))
var.ITS.DT.T1<- get_pca_var(pca.ITS.DT.T1)


pca.point.DT.ITS.T1<-data.frame(pca.ITS.DT.T1$ind)

pca.data.point.ITS.DT.T1<-pca.point.DT.ITS.T1[,which(colnames(pca.point.DT.ITS.T1) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.ITS.DT.T1)[which(colnames(pca.data.point.ITS.DT.T1) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("fun.DT.T1.axis1", "fun.DT.T1.axis2","fun.DT.T1.axis3","fun.DT.T1.axis4","fun.DT.T1.axis5") 


#Analysis of significant dimensions

res.desc.ITS.DT.T1 <- dimdesc(pca.ITS.DT.T1, axes = c(1,2,3,4,5), proba = 0.05)

# Description of dimension 
pca.axis1.ITS.DT.T1<-data.frame(res.desc.ITS.DT.T1$Dim.1)
pca.axis2.ITS.DT.T1<-data.frame(res.desc.ITS.DT.T1$Dim.2)
pca.axis3.ITS.DT.T1<-data.frame(res.desc.ITS.DT.T1$Dim.3)
pca.axis4.ITS.DT.T1<-data.frame(res.desc.ITS.DT.T1$Dim.4)
pca.axis5.ITS.DT.T1<-data.frame(res.desc.ITS.DT.T1$Dim.5)

write.table(pca.axis1.ITS.DT.T1,here("output/tables","PCA.axis1.ITS.DT.T1.txt"),sep="\t")
write.table(pca.axis2.ITS.DT.T1,here("output/tables","PCA.axis2.ITS.DT.T1.txt"),sep="\t")
write.table(pca.axis3.ITS.DT.T1,here("output/tables","PCA.axis3.ITS.DT.T1.txt"),sep="\t")
write.table(pca.axis4.ITS.DT.T1,here("output/tables","PCA.axis4.ITS.DT.T1.txt"),sep="\t")
write.table(pca.axis5.ITS.DT.T1,here("output/tables","PCA.axis5.ITS.DT.T1.txt"),sep="\t")

fungi.DT.T1.cor <- corr.test(qual.DT.S1,pca.data.point.ITS.DT.T1, method = 'spearman')

myfun <- function(fungi.DT.T1.cor){
  corrplot(fungi.DT.T1.cor$r, 
           p.mat = round(as.matrix(fungi.DT.T1.cor$p),3),
           method = 'circle',
           type = 'lower',
           sig.level = c(.001, .01, .05), 
           tl.pos="lt", 
           tl.col="black", tl.cex=1.0, 
           tl.offset=0.2,
           cl.pos="r",
           insig = "label_sig",
           pch.cex = 1.0,
           pch.col="red",
           cl.cex = 1.0)
  corrplot(fungi.DT.T1.cor$r,  type="upper", method="number",
           col="coral4",  tl.pos="n", cl.pos="n", number.cex = 0.9, add=T,diag=F)
  recordPlot() # record the latest plot
}


myfun(fungi.DT.T1.cor) 



#PCA run for 24 MAY
pca.ITS.DT.T2 <- PCA(otu.ITS.DT.T2.hel, graph = FALSE)
print(pca.ITS.DT.T2)

#Visualization and interpretation of eigenvalues
eig.val.ITS.DT.T2 <- get_eigenvalue(pca.ITS.DT.T2)
eig.val.ITS.DT.T2

#Cumulative percentage of eigenvalues
fviz_eig(pca.ITS.DT.T2, addlabels = TRUE, ylim = c(0, 10))
var.ITS.DT.T2<- get_pca_var(pca.ITS.DT.T2)


pca.point.DT.ITS.T2<-data.frame(pca.ITS.DT.T2$ind)

pca.data.point.ITS.DT.T2<-pca.point.DT.ITS.T2[,which(colnames(pca.point.DT.ITS.T2) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.ITS.DT.T2)[which(colnames(pca.data.point.ITS.DT.T2) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("fun.DT.T2.axis1", "fun.DT.T2.axis2","fun.DT.T1.axis3","fun.DT.T2.axis4","fun.DT.T2.axis5") 





#PCA run for 07 th June
pca.ITS.DT.T3 <- PCA(otu.ITS.DT.T3.hel, graph = FALSE)
print(pca.ITS.DT.T3)

#Visualization and interpretation of eigenvalues
eig.val.ITS.DT.T3 <- get_eigenvalue(pca.ITS.DT.T3)
eig.val.ITS.DT.T3

#Cumulative percentage of eigenvalues
fviz_eig(pca.ITS.DT.T3, addlabels = TRUE, ylim = c(0, 10))
var.ITS.DT.T3<- get_pca_var(pca.ITS.DT.T3)


pca.point.DT.ITS.T3<-data.frame(pca.ITS.DT.T3$ind)

pca.data.point.ITS.DT.T3<-pca.point.DT.ITS.T3[,which(colnames(pca.point.DT.ITS.T3) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.ITS.DT.T3)[which(colnames(pca.data.point.ITS.DT.T3) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("fun.DT.T3.axis1", "fun.DT.T3.axis2","fun.DT.T3.axis3","fun.DT.T3.axis4","fun.DT.T3.axis5") 



res.desc.ITS.DT.T3 <- dimdesc(pca.ITS.DT.T3, axes = c(1,2,3,4,5), proba = 0.05)

PCA.axis1.ITS.DT.T3<-data.frame(res.desc.ITS.DT.T3$Dim.1)
PCA.axis2.ITS.DT.T3<-data.frame(res.desc.ITS.DT.T3$Dim.2)
PCA.axis3.ITS.DT.T3<-data.frame(res.desc.ITS.DT.T3$Dim.3)
PCA.axis4.ITS.DT.T3<-data.frame(res.desc.ITS.DT.T3$Dim.4)
PCA.axis5.ITS.DT.T3<-data.frame(res.desc.ITS.DT.T3$Dim.5)


write.table(PCA.axis1.ITS.DT.T3,here("output/tables","PCA.axis1.ITS.DT.T3.txt"),sep="\t")
write.table(PCA.axis2.ITS.DT.T3,here("output/tables","PCA.axis2.ITS.DT.T3.txt"),sep="\t")
write.table(PCA.axis3.ITS.DT.T3,here("output/tables","PCA.axis3.ITS.DT.T3.txt"),sep="\t")
write.table(PCA.axis4.ITS.DT.T3,here("output/tables","PCA.axis4.ITS.DT.T3.txt"),sep="\t")
write.table(PCA.axis5.ITS.DT.T3,here("output/tables","PCA.axis5.ITS.DT.T3.txt"),sep="\t")

################################################################################
#Corr-plot visualizations

fungi.DT.T3.cor <- corr.test(qual.DT.S3,pca.data.point.ITS.DT.T3, method = 'spearman')

myfun <- function(fungi.DT.T3.cor){
  corrplot(fungi.DT.T3.cor$r, 
           p.mat = round(as.matrix(fungi.DT.T3.cor$p),3),
           method = 'circle',
           type = 'lower',
           sig.level = c(.001, .01, .05), 
           tl.pos="lt", 
           tl.col="black", tl.cex=1.0, 
           tl.offset=0.2,
           cl.pos="r",
           insig = "label_sig",
           pch.cex = 1.0,
           pch.col="red",
           cl.cex = 1.0)
  corrplot(fungi.DT.T3.cor$r,  type="upper", method="number",
           col="coral4",  tl.pos="n", cl.pos="n", number.cex = 0.9, add=T,diag=F)
  recordPlot() # record the latest plot
}

myfun(fungi.DT.T3.cor) 

#PCA run for 27th June
pca.ITS.DT.T4 <- PCA(otu.ITS.DT.T4.hel, graph = FALSE)
print(pca.ITS.DT.T4)

#Visualization and interpretation of eigenvalues
eig.val.ITS.DT.T4 <- get_eigenvalue(pca.ITS.DT.T4)
eig.val.ITS.DT.T4

#Cumulative percentage of eigenvalues
fviz_eig(pca.ITS.DT.T4, addlabels = TRUE, ylim = c(0, 10))
var.ITS.DT.T4<- get_pca_var(pca.ITS.DT.T4)


pca.point.DT.ITS.T4<-data.frame(pca.ITS.DT.T4$ind)

pca.data.point.ITS.DT.T4<-pca.point.DT.ITS.T4[,which(colnames(pca.point.DT.ITS.T4) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.ITS.DT.T4)[which(colnames(pca.data.point.ITS.DT.T4) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("fun.DT.T4.axis1", "fun.DT.T4.axis2","fun.DT.T4.axis3","fun.DT.T4.axis4","fun.DT.T4.axis5") 


res.desc.ITS.DT.T4 <- dimdesc(pca.ITS.DT.T4, axes = c(3), proba = 0.05)
PCA.axis3.ITS.DT.T4<-data.frame(res.desc.ITS.DT.T4$Dim.3)


write.table(PCA.axis3.ITS.DT.T4,here("output/tables","PCA.axis3.ITS.DT.T4.txt"),sep="\t")

fungi.DT.T4.cor <- corr.test(qual.DT.S4,pca.data.point.ITS.DT.T4, method = 'spearman')

myfun <- function(fungi.DT.T4.cor){
  corrplot(fungi.DT.T4.cor$r, 
           p.mat = round(as.matrix(fungi.DT.T4.cor$p),3),
           method = 'circle',
           type = 'lower',
           sig.level = c(.001, .01, .05), 
           tl.pos="lt", 
           tl.col="black", tl.cex=1.0, 
           tl.offset=0.2,
           cl.pos="r",
           insig = "label_sig",
           pch.cex = 1.0,
           pch.col="red",
           cl.cex = 1.0)
  corrplot(fungi.DT.T4.cor$r,  type="upper", method="number",
           col="coral4",  tl.pos="n", cl.pos="n", number.cex = 0.9, add=T,diag=F)
  recordPlot() # record the latest plot
}


myfun(fungi.DT.T4.cor) 




#PCA run for 27th June
pca.ITS.DT.T5 <- PCA(otu.ITS.DT.T5.hel, graph = FALSE)
print(pca.ITS.DT.T5)

#Visualization and interpretation of eigenvalues
eig.val.ITS.DT.T5 <- get_eigenvalue(pca.ITS.DT.T5)
eig.val.ITS.DT.T5

#Cumulative percentage of eigenvalues
fviz_eig(pca.ITS.DT.T5, addlabels = TRUE, ylim = c(0, 10))
var.ITS.DT.T5<- get_pca_var(pca.ITS.DT.T5)


pca.point.DT.ITS.T5<-data.frame(pca.ITS.DT.T5$ind)

pca.data.point.ITS.DT.T5<-pca.point.DT.ITS.T5[,which(colnames(pca.point.DT.ITS.T5) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.ITS.DT.T5)[which(colnames(pca.data.point.ITS.DT.T5) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("fun.DT.T5.axis1", "fun.DT.T5.axis2","fun.DT.T5.axis3","fun.DT.T4.axis4","fun.DT.T5.axis5") 


#PCA run for 19th JUly
pca.ITS.DT.T6 <- PCA(otu.ITS.DT.T6.hel, graph = FALSE)
print(pca.ITS.DT.T6)

#Visualization and interpretation of eigenvalues
eig.val.ITS.DT.T6 <- get_eigenvalue(pca.ITS.DT.T6)
eig.val.ITS.DT.T6

#Cumulative percentage of eigenvalues
fviz_eig(pca.ITS.DT.T6, addlabels = TRUE, ylim = c(0, 10))
var.ITS.DT.T6<- get_pca_var(pca.ITS.DT.T6)


pca.point.DT.ITS.T6<-data.frame(pca.ITS.DT.T6$ind)

pca.data.point.ITS.DT.T6<-pca.point.DT.ITS.T6[,which(colnames(pca.point.DT.ITS.T6) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.ITS.DT.T6)[which(colnames(pca.data.point.ITS.DT.T6) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("fun.DT.T6.axis1", "fun.DT.T6.axis2","fun.DT.T6.axis3","fun.DT.T6.axis4","fun.DT.T6.axis5") 

#PCA run for 1st August
pca.ITS.DT.T7 <- PCA(otu.ITS.DT.T7.hel, graph = FALSE)
print(pca.ITS.DT.T7)

#Visualization and interpretation of eigenvalues
eig.val.ITS.DT.T7 <- get_eigenvalue(pca.ITS.DT.T7)
eig.val.ITS.DT.T7

#Cumulative percentage of eigenvalues
fviz_eig(pca.ITS.DT.T7, addlabels = TRUE, ylim = c(0, 10))
var.ITS.DT.T7<- get_pca_var(pca.ITS.DT.T7)


pca.point.DT.ITS.T7<-data.frame(pca.ITS.DT.T7$ind)

pca.data.point.ITS.DT.T7<-pca.point.DT.ITS.T7[,which(colnames(pca.point.DT.ITS.T7) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.ITS.DT.T7)[which(colnames(pca.data.point.ITS.DT.T7) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("fun.DT.T7.axis1", "fun.DT.T7.axis2","fun.DT.T7.axis3","fun.DT.T7.axis4","fun.DT.T7.axis5") 


res.desc.ITS.DT.T7 <- dimdesc(pca.ITS.DT.T7, axes = c(1,2,3,4,5), proba = 0.05)

# Description of dimension 2
pca.axis3.ITS.DT.T7<-data.frame(res.desc.ITS.DT.T7$Dim.3)
pca.axis4.ITS.DT.T7<-data.frame(res.desc.ITS.DT.T7$Dim.4)
pca.axis5.ITS.DT.T7<-data.frame(res.desc.ITS.DT.T7$Dim.5)


write.table(pca.axis3.ITS.DT.T7,here("output/tables","PCA.axis3.ITS.DT.T7.csv"),  sep="\t")
write.table(pca.axis4.ITS.DT.T7,here("output/tables","PCA.axis4.ITS.DT.T7.csv"),  sep="\t")
write.table(pca.axis5.ITS.DT.T7,here("output/tables","PCA.axis5.ITS.DT.T7.csv"),  sep="\t")




#Biolog PCA

#10 May data (BIOLOG)


#PCA run for T1
pca.biolog.DT.T1 <- PCA(biolog.DT.T1.hel, graph = FALSE)
print(pca.biolog.DT.T1)

#Visualization and interpretation of eigenvalues
eig.val.biolog.DT.T1 <- get_eigenvalue(pca.biolog.DT.T1)
eig.val.biolog.DT.T1

#Cumulative percentage of eigenvalues
fviz_eig(pca.biolog.DT.T1, addlabels = TRUE, ylim = c(0, 40))
var.biolog.DT.T1<- get_pca_var(pca.biolog.DT.T1)

pca.point.biolog.DT.T1<-data.frame(pca.biolog.DT.T1$ind)

pca.data.point.biolog.DT.T1<-pca.point.biolog.DT.T1[,which(colnames(pca.point.biolog.DT.T1) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.biolog.DT.T1)[which(colnames(pca.data.point.biolog.DT.T1) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("biolog.DT.T1.axis1", "biolog.DT.T1.axis2","biolog.DT.T1.axis3","biolog.DT.T1.axis4","biolog.DT.T1.axis5") 


res.desc.biolog.DT.T1 <- dimdesc(pca.biolog.DT.T1, axes = c(1,2,3,4,5), proba = 0.05)
res.desc.biolog.DT.T1$Dim.1
PCA.axis1.biolog.DT.T1<-data.frame(res.desc.biolog.DT.T1$Dim.1)
PCA.axis2.biolog.DT.T1<-data.frame(res.desc.biolog.DT.T1$Dim.2)
PCA.axis3.biolog.DT.T1<-data.frame(res.desc.biolog.DT.T1$Dim.3)
PCA.axis4.biolog.DT.T1<-data.frame(res.desc.biolog.DT.T1$Dim.4)
PCA.axis5.biolog.DT.T1<-data.frame(res.desc.biolog.DT.T1$Dim.5)

write.table(PCA.axis1.biolog.DT.T1,here("output/tables","PCA.axis1.biolog.DT.T1.txt"),sep="\t")
write.table(PCA.axis2.biolog.DT.T1,here("output/tables","PCA.axis2.biolog.DT.T1.txt"),sep="\t")
write.table(PCA.axis3.biolog.DT.T1,here("output/tables","PCA.axis3.biolog.DT.T1.txt"),sep="\t")
write.table(PCA.axis4.biolog.DT.T1,here("output/tables","PCA.axis4.biolog.DT.T1.txt"),sep="\t")
write.table(PCA.axis5.biolog.DT.T1,here("output/tables","PCA.axis5.biolog.DT.T1.txt"),sep="\t")

library(corrplot)
library(psych)

biolog.DT.T1.cor <- corr.test(qual.DT.S1,pca.data.point.biolog.DT.T1, method = 'spearman')

myfun <- function(biolog.DT.T1.cor){
  corrplot(biolog.DT.T1.cor$r, 
           p.mat = round(as.matrix(biolog.DT.T1.cor$p),3),
           method = 'circle',
           type = 'lower',
           sig.level = c(.001, .01, .05), 
           tl.pos="lt", 
           tl.col="black", tl.cex=1.0, 
           tl.offset=0.2,
           cl.pos="r",
           insig = "label_sig",
           pch.cex = 1.0,
           pch.col="red",
           cl.cex = 1.0)
  corrplot(biolog.DT.T1.cor$r,  type="upper", method="number",
           col="coral4",  tl.pos="n", cl.pos="n", number.cex = 0.9, add=T,diag=F)
  recordPlot() # record the latest plot
}


myfun(biolog.DT.T1.cor) # shows the plot

p <- myfun(biolog.DT.T1.cor) 



library(missMDA)


#PCA run FOR T2
pca.biolog.DT.T2 <- PCA(biolog.DT.T2.hel, graph = FALSE)
print(pca.biolog.DT.T2)

#Visualization and interpretation of eigenvalues
eig.val.biolog.DT.T2 <- get_eigenvalue(pca.biolog.DT.T2)
eig.val.biolog.DT.T2

#Cumulative percentage of eigenvalues
fviz_eig(pca.biolog.DT.T2, addlabels = TRUE, ylim = c(0, 25))
var.biolog.DT.T2<- get_pca_var(pca.biolog.DT.T2)

pca.point.biolog.DT.T2<-data.frame(pca.biolog.DT.T2$ind)

pca.data.point.biolog.DT.T2<-pca.point.biolog.DT.T2[,which(colnames(pca.point.biolog.DT.T2) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.biolog.DT.T2)[which(colnames(pca.data.point.biolog.DT.T2) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("biolog.DT.T2.axis1", "biolog.DT.T2.axis2","biolog.DT.T2.axis3","biolog.DT.T2.axis4","biolog.DT.T2.axis5") 



res.desc.biolog.DT.T2 <- dimdesc(pca.biolog.DT.T2, axes = c(2,3,4), proba = 0.05)
PCA.axis2.biolog.DT.T2<-data.frame(res.desc.biolog.DT.T2$Dim.2)
PCA.axis2.biolog.DT.T2<-na.omit(PCA.axis2.biolog.DT.T2)

write.table(PCA.axis2.biolog.DT.T2,here("output/tables","PCA.axis2.biolog.DT.T2.txt"),sep="\t")




#PCA run FOR T3

pca.biolog.DT.T3 <- PCA(biolog.DT.T3.hel, graph = FALSE)
print(pca.biolog.DT.T3)

#Visualization and interpretation of eigenvalues
eig.val.biolog.DT.T3 <- get_eigenvalue(pca.biolog.DT.T3)
eig.val.biolog.DT.T3

#Cumulative percentage of eigenvalues
fviz_eig(pca.biolog.DT.T3, addlabels = TRUE, ylim = c(0, 60))
var.biolog.DT.T3<- get_pca_var(pca.biolog.DT.T3)

pca.point.biolog.DT.T3<-data.frame(pca.biolog.DT.T3$ind)

pca.data.point.biolog.DT.T3<-pca.point.biolog.DT.T3[,which(colnames(pca.point.biolog.DT.T3) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.biolog.DT.T3)[which(colnames(pca.data.point.biolog.DT.T3) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("biolog.DT.T3.axis1", "biolog.DT.T3.axis2","biolog.DT.T3.axis3","biolog.DT.T3.axis4","biolog.DT.T3.axis5") 



res.desc.biolog.DT.T3 <- dimdesc(pca.biolog.DT.T3, axes = c(2,3,4), proba = 0.05)
PCA.axis2.biolog.DT.T3<-data.frame(res.desc.biolog.DT.T3$Dim.2)
PCA.axis3.biolog.DT.T3<-data.frame(res.desc.biolog.DT.T3$Dim.3)
PCA.axis4.biolog.DT.T3<-data.frame(res.desc.biolog.DT.T3$Dim.4)

write.table(PCA.axis2.biolog.DT.T3,here("output/tables","PCA.axis2.biolog.DT.T3.txt"),sep="\t")
write.table(PCA.axis3.biolog.DT.T3,here("output/tables","PCA.axis3.biolog.DT.T3.txt"),sep="\t")
write.table(PCA.axis4.biolog.DT.T3,here("output/tables","PCA.axis4.biolog.DT.T3.txt"),sep="\t")



biolog.DT.T3.cor <- corr.test(qual.DT.S3,pca.data.point.biolog.DT.T3, method = 'spearman')

myfun <- function(biolog.DT.T3.cor){
  corrplot(biolog.DT.T3.cor$r, 
           p.mat = round(as.matrix(biolog.DT.T3.cor$p),3),
           method = 'circle',
           type = 'lower',
           sig.level = c(.001, .01, .05), 
           tl.pos="lt", 
           tl.col="black", tl.cex=1.0, 
           tl.offset=0.2,
           cl.pos="r",
           insig = "label_sig",
           pch.cex = 1.0,
           pch.col="red",
           cl.cex = 1.0)
  corrplot(biolog.DT.T3.cor$r,  type="upper", method="number",
           col="coral4",  tl.pos="n", cl.pos="n", number.cex = 0.9, add=T,diag=F)
  recordPlot() # record the latest plot
}


myfun(biolog.DT.T3.cor)







#PCA run FOR T4

pca.biolog.DT.T4 <- PCA(biolog.DT.T4.hel, graph = FALSE)
print(pca.biolog.DT.T4)

#Visualization and interpretation of eigenvalues
eig.val.biolog.DT.T4 <- get_eigenvalue(pca.biolog.DT.T4)
eig.val.biolog.DT.T4

#Cumulative percentage of eigenvalues
fviz_eig(pca.biolog.DT.T4, addlabels = TRUE, ylim = c(0, 35))
var.biolog.DT.T4<- get_pca_var(pca.biolog.DT.T4)

pca.point.biolog.DT.T4<-data.frame(pca.biolog.DT.T4$ind)

pca.data.point.biolog.DT.T4<-pca.point.biolog.DT.T4[,which(colnames(pca.point.biolog.DT.T4) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.biolog.DT.T4)[which(colnames(pca.data.point.biolog.DT.T4) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("biolog.DT.T4.axis1", "biolog.DT.T4.axis2","biolog.DT.T4.axis3","biolog.DT.T4.axis4","biolog.DT.T4.axis5") 

res.desc.biolog.DT.T4 <- dimdesc(pca.biolog.DT.T4, axes = c(1,4,5), proba = 0.05)
PCA.axis1.biolog.DT.T4<-data.frame(res.desc.biolog.DT.T4$Dim.1)
PCA.axis4.biolog.DT.T4<-data.frame(res.desc.biolog.DT.T4$Dim.4)
PCA.axis5.biolog.DT.T4<-data.frame(res.desc.biolog.DT.T4$Dim.5)

write.table(PCA.axis1.biolog.DT.T4,here("output/tables","PCA.axis1.biolog.DT.T4.txt"),sep="\t")
write.table(PCA.axis4.biolog.DT.T4,here("output/tables","PCA.axis4.biolog.DT.T4.txt"),sep="\t")
write.table(PCA.axis5.biolog.DT.T4,here("output/tables","PCA.axis5.biolog.DT.T4.txt"),sep="\t")

biolog.DT.T4.cor <- corr.test(qual.DT.S4,pca.data.point.biolog.DT.T4, method = 'spearman')

myfun <- function(biolog.DT.T4.cor){
  corrplot(biolog.DT.T4.cor$r, 
           p.mat = round(as.matrix(biolog.DT.T4.cor$p),3),
           method = 'circle',
           type = 'lower',
           sig.level = c(.001, .01, .05), 
           tl.pos="lt", 
           tl.col="black", tl.cex=1.0, 
           tl.offset=0.2,
           cl.pos="r",
           insig = "label_sig",
           pch.cex = 1.0,
           pch.col="red",
           cl.cex = 1.0)
  corrplot(biolog.DT.T4.cor$r,  type="upper", method="number",
           col="coral4",  tl.pos="n", cl.pos="n", number.cex = 0.9, add=T,diag=F)
  recordPlot() # record the latest plot
}


myfun(biolog.DT.T4.cor)



#PCA run FOR T5

pca.biolog.DT.T5 <- PCA(biolog.DT.T5.hel, graph = FALSE)
print(pca.biolog.DT.T5)

#Visualization and interpretation of eigenvalues
eig.val.biolog.DT.T5 <- get_eigenvalue(pca.biolog.DT.T5)
eig.val.biolog.DT.T5

#Cumulative percentage of eigenvalues
fviz_eig(pca.biolog.DT.T5, addlabels = TRUE, ylim = c(0, 25))
var.biolog.DT.T5<- get_pca_var(pca.biolog.DT.T5)

pca.point.biolog.DT.T5<-data.frame(pca.biolog.DT.T5$ind)

pca.data.point.biolog.DT.T5<-pca.point.biolog.DT.T5[,which(colnames(pca.point.biolog.DT.T5) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.biolog.DT.T5)[which(colnames(pca.data.point.biolog.DT.T5) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("biolog.DT.T5.axis1", "biolog.DT.T5.axis2","biolog.DT.T5.axis3","biolog.DT.T5.axis4","biolog.DT.T5.axis5") 


#PCA run FOR T6
pca.biolog.DT.T6 <- PCA(biolog.DT.T6.hel, graph = FALSE)
print(pca.biolog.DT.T6)

#Visualization and interpretation of eigenvalues
eig.val.biolog.DT.T6 <- get_eigenvalue(pca.biolog.DT.T6)
eig.val.biolog.DT.T6

#Cumulative percentage of eigenvalues
fviz_eig(pca.biolog.DT.T6, addlabels = TRUE, ylim = c(0, 25))
var.biolog.DT.T6<- get_pca_var(pca.biolog.DT.T6)

pca.point.biolog.DT.T6<-data.frame(pca.biolog.DT.T6$ind)

pca.data.point.biolog.DT.T6<-pca.point.biolog.DT.T6[,which(colnames(pca.point.biolog.DT.T6) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.biolog.DT.T6)[which(colnames(pca.data.point.biolog.DT.T6) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("biolog.DT.T6.axis1", "biolog.DT.T6.axis2","biolog.DT.T6.axis3","biolog.DT.T6.axis4","biolog.DT.T6.axis5") 

#PCA run for T7
pca.biolog.DT.T7 <- PCA(biolog.DT.T7.hel, graph = FALSE)
print(pca.biolog.DT.T7)

#Visualization and interpretation of eigenvalues
eig.val.biolog.DT.T7 <- get_eigenvalue(pca.biolog.DT.T7)
eig.val.biolog.DT.T7

#Cumulative percentage of eigenvalues
fviz_eig(pca.biolog.DT.T7, addlabels = TRUE, ylim = c(0, 25))
var.biolog.DT.T7<- get_pca_var(pca.biolog.DT.T7)

pca.point.biolog.DT.T7<-data.frame(pca.biolog.DT.T7$ind)

pca.data.point.biolog.DT.T7<-pca.point.biolog.DT.T7[,which(colnames(pca.point.biolog.DT.T7) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.biolog.DT.T7)[which(colnames(pca.data.point.biolog.DT.T7) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("biolog.DT.T7.axis1", "biolog.DT.T7.axis2","biolog.DT.T7.axis3","biolog.DT.T7.axis4","biolog.DT.T7.axis5") 


res.desc.biolog.DT.T7 <- dimdesc(pca.biolog.DT.T7, axes = c(1,2,3,4,5), proba = 0.05)

PCA.axis4.biolog.DT.T7<-data.frame(res.desc.biolog.DT.T7$Dim.4)
PCA.axis5.biolog.DT.T7<-data.frame(res.desc.biolog.DT.T7$Dim.5)


write.table(PCA.axis4.biolog.DT.T7,here("output/tables","PCA.axis4.biolog.DT.T7.txt"),sep="\t")
write.table(PCA.axis5.biolog.DT.T7,here("output/tables","PCA.axis5.biolog.DT.T7.txt"),sep="\t")






# LASSO regression analysis

# Preparing data for LASSO regression for T1

reg.lasso.DT.T1<-cbind(pca.data.point.16S.DT.T1,pca.data.point.ITS.DT.T1,pca.data.point.biolog.DT.T1,div.16S.DT.T1,div.ITS.DT.T1,qpcr.T1.DT,qual.DT.S1)
row.names(reg.lasso.DT.T1)==row.names(qual.DT.S1)
#reg.lasso.DT.T1<-decostand(reg.lasso.DT.T1,method="standardize")

#qpcr.T1.DT<-log(qpcr.T1.DT)

## Protein model for May 10
data.glut.DT.T1=reg.lasso.DT.T1[!is.na(reg.lasso.DT.T1$Gluten),]
#data.glut.DT.T1=data.glut.DT.T1[-c(5,18),]

data.glut.DT.T1=data.glut.DT.T1[-c(5,18),]
data.glut.DT.T1=data.glut.DT.T1 [,-which(names(data.glut.DT.T1) %in% c("F.B.ratio"))]

data.prot.DT.T1=reg.lasso.DT.T1[!is.na(reg.lasso.DT.T1$Protein.grain),]
data.pmt.DT.T1=reg.lasso.DT.T1[!is.na(reg.lasso.DT.T1$PMT),]
data.bem.DT.T1=reg.lasso.DT.T1[!is.na(reg.lasso.DT.T1$BEM),]



data.explain.glut.DT.T1=data.glut.DT.T1[,-which(names(data.glut.DT.T1) %in% c("Protein.grain", "Proteine.flour","Humidity","Ash","Gluten", "IDC","PMT","BEM"))]
data.explain.prot.DT.T1=data.prot.DT.T1[,-which(names(data.prot.DT.T1) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.pmt.DT.T1=data.pmt.DT.T1[,-which(names(data.pmt.DT.T1) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.bem.DT.T1=data.bem.DT.T1[,-which(names(data.bem.DT.T1) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]





# X-Y data prep
X.glut.DT.T1=as.matrix(scale(data.explain.glut.DT.T1,center = T, scale = T))
Y.glut.DT.T1=as.matrix(scale(data.glut.DT.T1$Gluten,center = T, scale = T))
#Y.glut.DT.T1=as.matrix(data.glut.DT.T1$Gluten)

X.prot.DT.T1=as.matrix(scale(data.explain.prot.DT.T1,center = T, scale = T))
#Y.reg.prot.DT.T1<-as.matrix(data.prot.DT.T1$Protein.grain)
Y.prot.DT.T1=as.matrix(scale(data.prot.DT.T1$Protein.grain,center = T, scale = T))

X.pmt.DT.T1=as.matrix(scale(data.explain.pmt.DT.T1,center = T, scale = T))
Y.pmt.DT.T1=as.matrix(scale(data.pmt.DT.T1$PMT,center = T, scale = T))

X.bem.DT.T1=as.matrix(scale(data.explain.bem.DT.T1,center = T, scale = T))
Y.bem.DT.T1=as.matrix(scale(data.bem.DT.T1$BEM,center = T, scale = T))

#LASSO regression
fit.lasso.glut <-glmnet(X.glut.DT.T1,Y.glut.DT.T1, family="gaussian", alpha=1)
fit.lasso.prot<- glmnet(X.prot.DT.T1,Y.prot.DT.T1, family="gaussian", alpha=1)
fit.lasso.pmt <- glmnet(X.pmt.DT.T1,Y.pmt.DT.T1, family="gaussian", alpha=1)
fit.lasso.bem <- glmnet(X.bem.DT.T1,Y.bem.DT.T1, family="gaussian", alpha=1)


#Gluten
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.glut.DT.T1, Y.glut.DT.T1, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.glut, s=s.best.lasso))[-1]
s.best.lasso

lasso.glut.DT.T1<-coef(fit.lasso.glut, s=s.best.lasso)
lasso.glut.DT.T1

lasso_test.DT.T1.glut<-glmnet(X.glut.DT.T1, Y.glut.DT.T1, alpha = 1, standardize=TRUE)
plot(lasso_test.DT.T1.glut, xvar = "lambda")
#perform k-fold cross-validation to find optimal lambda value
##grid=10^seq(10,-2, length=100)# We can apply this extra grid function for tuning lambdas. This "grid" search begins with null model up-to 10 fold repetitive search for optimal lambda values
# 10 fold random division of data without specifying the row number may yield different lambda values)
#X= number of row of the input data e.g. X.glut.DT.T1 (to obtain the constant lambda value during cross validation, as it contains only 23 sample.
#If the objective X is not working in the nrow function, then use directly X data matrix, like nfolds=nrow(X.glut.DT.T1).

#One can change the cross validation threshold to 3 or 5 folds, but in my analysis 10 folds produced the most optimal lambda.

cv_model.glut.DT.T1 <- cv.glmnet(X.glut.DT.T1, Y.glut.DT.T1, k=10, nfolds=nrow(X.glut.DT.T1), alpha = 1, standardize=TRUE, grouped=FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.glut.DT.T1 <- cv_model.glut.DT.T1$lambda.min
best_lambda.glut.DT.T1

#produce plot of test MSE by lambda value
plot(cv_model.glut.DT.T1) 

#find coefficients of best model
best_model.glut.DT.T1 <- glmnet(X.glut.DT.T1, Y.glut.DT.T1, alpha = 1, lambda = best_lambda.glut.DT.T1, standardize = TRUE)
coef(best_model.glut.DT.T1)


#use fitted best model to make predictions
y_predicted.glut.DT.T1 <- predict(best_model.glut.DT.T1, s = best_lambda.glut.DT.T1, newx = X.glut.DT.T1)

predicted.value.glut.DT.T1<-data.frame(cbind(y_predicted.glut.DT.T1,Y.glut.DT.T1))
#find SST and SSE
sst <- sum((Y.glut.DT.T1 - mean(Y.glut.DT.T1))^2)
sse <- sum((y_predicted.glut.DT.T1 - Y.glut.DT.T1)^2)

#find R-Squared
rsq <- 1 - sse/sst
rsq

#Evaluation of the accuracy
Rsquare <- cor(Y.glut.DT.T1,y_predicted.glut.DT.T1)^2
Rsquare

MSE <- mean((Y.glut.DT.T1 - y_predicted.glut.DT.T1)^2)
MSE

plot(Y.glut.DT.T1, y_predicted.glut.DT.T1,col="#00000050",pch = 19,main = "Test result for lasso")
abline(lm(y_predicted.glut.DT.T1~Y.glut.DT.T1),lty = 2,lwd = 2,col = "gray")


plot(best_model.glut.DT.T1, xvar = "lambda")

# plot lasso model
lam <- best_model.glut.DT.T1$lambda %>% 
  as.data.frame() %>%
  mutate(penalty = best_model.glut.DT.T1$a0 %>% names()) %>%
  rename(lambda = ".")

results.DT.glut.T1 <- best_model.glut.DT.T1$beta %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(penalty, coefficients, -rowname) %>%
  left_join(lam)

results.DT.glut.T1.filt<-filter_if(results.DT.glut.T1, is.numeric, all_vars((.) != 0))
results.DT.glut.T1.filt



# Visualize the predictive plot
write.table(results.DT.glut.T1.filt,here("output/tables", "results_labels_glut_DT.T1.txt"),sep="\t" )

plot.glut.DT.T1<-ggplot(predicted.value.glut.DT.T1, aes(x=V2, y=s1)) + 
  geom_point(color="#00000050")+
  geom_smooth(method=lm, color="lightblue")+
  #geom_text(x = 25, y = 31, label = eq(predicted.value.glut.DT.T3$Y.glut.DT.T3,predicted.value.glut.DT.T3$s1), parse = TRUE) +
  
  stat_regline_equation(label.y =2, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 1.7, aes(label = ..rr.label..)) +
  labs(title="Gluten (May10)  ",
       y=" Predicted gluten content", x = " Observed gluten content")+
  theme_classic() 

plot.glut.DT.T1

ggsave(file=here("output/photo","model.DT.T1.glut.tiff"), plot.glut.DT.T1, height=3.5, width=4.0, units="in", dpi=600)




#Protein
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.prot.DT.T1,Y.prot.DT.T1, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.prot, s=s.best.lasso))[-1]
s.best.lasso

lasso.prot.DT.T1<-coef(fit.lasso.prot, s=s.best.lasso)
lasso.prot.DT.T1


lasso_test.DT.T1.prot<-glmnet(X.prot.DT.T1, Y.prot.DT.T1, alpha = 1, standardize=TRUE)
plot(lasso_test.DT.T1.glut, xvar = "lambda")

#perform k-fold cross-validation to find optimal lambda value
cv_model.prot.DT.T1 <- cv.glmnet(X.prot.DT.T1, Y.prot.DT.T1, k=10,nfold=nrow(X.prot.DT.T1), alpha = 1, standardize=TRUE, grouped=FALSE)


#find optimal lambda value that minimizes test MSE
best_lambda.prot.DT.T1 <- cv_model.prot.DT.T1$lambda.min
best_lambda.prot.DT.T1

#produce plot of test MSE by lambda value
plot(cv_model.prot.DT.T1) 

#find coefficients of best model
best_model.prot.DT.T1 <- glmnet(X.prot.DT.T1, Y.prot.DT.T1, alpha = 1, lambda = best_lambda.prot.DT.T1, standardize = TRUE)
coef(best_model.prot.DT.T1)

#use fitted best model to make predictions

y_predicted.prot.DT.T1 <- predict(best_model.prot.DT.T1, s = best_lambda.prot.DT.T1, newx = X.prot.DT.T1)
#y_predicted.prot.DT.T1<-scale(y_predicted.prot.DT.T1, center=TRUE)


#r1 <- y_predicted.prot.DT.T1* attr(y_predicted.prot.DT.T1, 'scaled:scale')[col(y_predicted.prot.DT.T1)] + attr(y_predicted.prot.DT.T1, 'scaled:center')[col(y_predicted.prot.DT.T1)]
#all.equal(as.data.frame(r1), y_predicted.prot.DT.T1)

predicted.value.prot.DT.T1<-data.frame(cbind(y_predicted.prot.DT.T1,Y.prot.DT.T1))

#Evaluation of the accuracy
Rsquare <- cor(Y.prot.DT.T1,y_predicted.prot.DT.T1)^2
Rsquare

MSE <- mean((Y.prot.DT.T1 - y_predicted.prot.DT.T1)^2)
MSE


# plot lasso model
lam <- best_model.prot.DT.T1$lambda %>% 
  as.data.frame() %>%
  mutate(penalty = best_model.prot.DT.T1$a0 %>% names()) %>%
  rename(lambda = ".")

results.DT.prot.T1 <- best_model.prot.DT.T1$beta %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(penalty, coefficients, -rowname) %>%
  left_join(lam)       

results.DT.prot.T1.filt<-filter_if(results.DT.prot.T1, is.numeric, all_vars((.) != 0))
results.DT.prot.T1.filt



# Visualize the predictive plot
write.table(results.DT.prot.T1.filt,file = here("output/tables", "results_labels_protein_DT.T1.txt"),sep="\t" )

eq <- function(x,y) {
  m <- lm(y ~ x)
  as.character(
    as.expression(
      substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                 list(a = format(coef(m)[1], digits = 4),
                      b = format(coef(m)[2], digits = 4),
                      r2 = format(summary(m)$r.squared, digits = 3)))
    )
  )
}





plot.prot.DT.T1<-ggplot(predicted.value.prot.DT.T1, aes(x=V2, y=s1)) + 
  geom_point(color="#00000050")+
  geom_smooth(method=lm, color="lightblue")+
  #geom_text(x = 25, y = 31, label = eq(predicted.value.glut.DT.T3$Y.glut.DT.T3,predicted.value.glut.DT.T3$s1), parse = TRUE) +
  
  stat_regline_equation(label.y =2, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 1.7, aes(label = ..rr.label..)) +
  labs(title="Protein (May10)  ",
       y=" Predicted protein content ", x = " Observed protein content ")+
  theme_classic() 

plot.prot.DT.T1

ggsave(file=here("output/photo","model.DT.T1.prot.tiff"), plot.prot.DT.T1, height=3.5, width=4.0, units="in", dpi=600)





#PMT
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.pmt.DT.T1, Y.pmt.DT.T1, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.pmt, s=s.best.lasso))[-1]
s.best.lasso

lasso.pmt.DT.T1<-coef(fit.lasso.pmt, s=s.best.lasso)
lasso.pmt.DT.T1

#perform k-fold cross-validation to find optimal lambda value
cv_model.pmt.DT.T1 <- cv.glmnet(X.pmt.DT.T1, Y.pmt.DT.T1, K=10, nfolds=nrow(X.pmt.DT.T1), alpha = 1, standardize=TRUE, grouped =FALSE )

#find optimal lambda value that minimizes test MSE
best_lambda.pmt.DT.T1 <- cv_model.pmt.DT.T1$lambda.min
best_lambda.pmt.DT.T1

#produce plot of test MSE by lambda value
plot(cv_model.pmt.DT.T1) 

#find coefficients of best model
best_model.pmt.DT.T1 <- glmnet(X.pmt.DT.T1, Y.pmt.DT.T1, alpha = 1, lambda = best_lambda.pmt.DT.T1, standardize = TRUE)
coef(best_model.pmt.DT.T1)




#use fitted best model to make predictions
y_predicted.pmt.DT.T1 <- predict(best_model.pmt.DT.T1, s = best_lambda.pmt.DT.T1, newx = X.pmt.DT.T1)

#Evaluation of the accuracy
Rsquare <- cor(Y.pmt.DT.T1,y_predicted.pmt.DT.T1)^2
Rsquare

MSE <- mean((Y.pmt.DT.T1 - y_predicted.pmt.DT.T1)^2)
MSE




#BEM
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.bem.DT.T1, Y.bem.DT.T1,nfolds = 3, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.bem, s=s.best.lasso))[-1]
s.best.lasso

lasso.bem.DT.T1<-coef(fit.lasso.bem, s=s.best.lasso)
lasso.bem.DT.T1

#perform k-fold cross-validation to find optimal lambda value
cv_model.bem.DT.T1 <- cv.glmnet(X.bem.DT.T1, Y.bem.DT.T1, k=10,nfolds=nrow(X.bem.DT.T1), alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.bem.DT.T1 <- cv_model.bem.DT.T1$lambda.min
best_lambda.bem.DT.T1

#produce plot of test MSE by lambda value
plot(cv_model.bem.DT.T1) 

#find coefficients of best model
best_model.bem.DT.T1 <- glmnet(X.bem.DT.T1, Y.bem.DT.T1, alpha = 1, lambda = best_lambda.bem.DT.T1, standardize = TRUE)
coef(best_model.bem.DT.T1)




#use fitted best model to make predictions
y_predicted.bem.DT.T1 <- predict(best_model.bem.DT.T1, s = best_lambda.bem.DT.T1, newx = X.bem.DT.T1)



#Evaluation of the accuracy
Rsquare <- cor(Y.bem.DT.T1,y_predicted.bem.DT.T1)^2
Rsquare

MSE <- mean((Y.bem.DT.T1 - y_predicted.bem.DT.T1)^2)
MSE

# plot lasso model
lam <- best_model.bem.DT.T1$lambda %>% 
  as.data.frame() %>%
  mutate(penalty = best_model.bem.DT.T1$a0 %>% names()) %>%
  rename(lambda = ".")

results.DT.bem.T1 <- best_model.bem.DT.T1$beta %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(penalty, coefficients, -rowname) %>%
  left_join(lam)

results.DT.bem.T1.filt<-filter_if(results.DT.bem.T1, is.numeric, all_vars((.) != 0))
results.DT.bem.T1.filt


#lasso.DT.T1=as.matrix(lasso.DT.T1)
#lasso.df.DT.T1=as.data.frame(lasso.DT.T1)
#write.table(lasso.df.DT.T1,"lasso.df.protien.txt",sep="\t")


model.glut.DT.T1=lm(Gluten~bact.DT.T1.axis2+fun.DT.T1.axis2+fun.DT.T1.axis3,
                    data=data.glut.DT.T1)


model.protein.DT.T1=lm(Protein.grain~bact.DT.T1.axis2+fun.DT.T1.axis1+fun.DT.T1.axis3+ biolog.DT.T1.axis4+biolog.DT.T1.axis5+nirk
                       ,data=data.prot.DT.T1)

model.pmt.DT.T1=lm(PMT~biolog.DT.T1.axis1, data=data.pmt.DT.T1)

model.bem.DT.T1=lm(BEM~biolog.DT.T1.axis1, data=data.bem.DT.T1)





vif(model.protein.DT.T1)

summary(model.glut.DT.T1)
summary(model.protein.DT.T1)
summary(model.pmt.DT.T1)
summary(model.bem.DT.T1)






# LASSO regression analysis

# Preparing data for LASSO regression for T2

reg.lasso.DT.T2<-cbind(pca.data.point.16S.DT.T2,pca.data.point.ITS.DT.T2,pca.data.point.biolog.DT.T2,div.16S.DT.T2,div.ITS.DT.T2,qpcr.T2.DT,qual.DT.S2)
row.names(reg.lasso.DT.T2)==row.names(qual.DT.S2)
#reg.lasso.DT.T2<-decostand(reg.lasso.DT.T2,method="standardize")

#qpcr.T1.DT<-log(qpcr.T1.DT)

## Protein model for May 10
data.glut.DT.T2=reg.lasso.DT.T2[!is.na(reg.lasso.DT.T2$Gluten),]
data.glut.DT.T2=data.glut.DT.T2[-c(5,18),]

data.prot.DT.T2=reg.lasso.DT.T2[!is.na(reg.lasso.DT.T2$Protein.grain),]
data.pmt.DT.T2=reg.lasso.DT.T2[!is.na(reg.lasso.DT.T2$PMT),]
data.bem.DT.T2=reg.lasso.DT.T2[!is.na(reg.lasso.DT.T2$BEM),]

data.explain.glut.DT.T2=data.glut.DT.T2[,-which(names(data.glut.DT.T2) %in% c("Protein.grain", "Proteine.flour","Humidity","Ash","Gluten", "IDC","PMT","BEM"))]
data.explain.prot.DT.T2=data.prot.DT.T2[,-which(names(data.prot.DT.T2) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.pmt.DT.T2=data.pmt.DT.T2[,-which(names(data.pmt.DT.T2) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.bem.DT.T2=data.prot.DT.T2[,-which(names(data.glut.DT.T2) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]

# X-Y data prep
X.glut.DT.T2=as.matrix(scale(data.explain.glut.DT.T2,center = T,scale = T))
Y.glut.DT.T2=as.matrix(scale(data.glut.DT.T2$Gluten,center = T,scale = T))

X.prot.DT.T2=as.matrix(scale(data.explain.prot.DT.T2,center = T,scale = T))
Y.prot.DT.T2=as.matrix(scale(data.prot.DT.T2$Protein.grain,center = T,scale = T))

X.pmt.DT.T2=as.matrix(scale(data.explain.pmt.DT.T2,center = T,scale = T))
Y.pmt.DT.T2=as.matrix(scale(data.pmt.DT.T2$PMT,center = T,scale = T))

X.bem.DT.T2=as.matrix(scale(data.explain.bem.DT.T2,center = T,scale = T))
Y.bem.DT.T2=as.matrix(scale(data.bem.DT.T2$BEM,center = T,scale = T))

#LASSO regression
fit.lasso.glut.DT.T2 <-glmnet(X.glut.DT.T2,Y.glut.DT.T2, family="gaussian", alpha=1)
fit.lasso.prot.DT.T2<- glmnet(X.prot.DT.T2,Y.prot.DT.T2, family="gaussian", alpha=1)
fit.lasso.pmt.DT.T2 <- glmnet(X.pmt.DT.T2,Y.pmt.DT.T2, family="gaussian", alpha=1)
fit.lasso.bem.DT.T2 <- glmnet(X.bem.DT.T2,Y.bem.DT.T2, family="gaussian", alpha=1)


#Gluten
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.glut.DT.T2, Y.glut.DT.T2, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.glut.DT.T2, s=s.best.lasso))[-1]
s.best.lasso

lasso.glut.DT.T2<-coef(fit.lasso.glut.DT.T2, s=s.best.lasso)
lasso.glut.DT.T2


#perform k-fold cross-validation to find optimal lambda value
cv_model.glut.DT.T2<- cv.glmnet(X.glut.DT.T2, Y.glut.DT.T2, k=10, nfolds = nrow(X.glut.DT.T2), alpha = 1, standardize=TRUE, grouped = FALSE)



#find optimal lambda value that minimizes test MSE
best_lambda.glut.DT.T2 <- cv_model.glut.DT.T2$lambda.min
best_lambda.glut.DT.T2

#produce plot of test MSE by lambda value
plot(cv_model.glut.DT.T2) 

#find coefficients of best model
best_model.glut.DT.T2 <- glmnet(X.glut.DT.T2, Y.glut.DT.T2, alpha = 1, lambda = best_lambda.glut.DT.T2, standardize = TRUE)
coef(best_model.glut.DT.T2)

#use fitted best model to make predictions
y_predicted.glut.DT.T2 <- predict(best_model.glut.DT.T2, s = best_lambda.glut.DT.T2, newx = X.glut.DT.T2)



#Evaluation of the accuracy
Rsquare <- cor(Y.glut.DT.T2,y_predicted.glut.DT.T2)^2
Rsquare

MSE <- mean((Y.glut.DT.T2 - y_predicted.glut.DT.T2)^2)
MSE

plot(Y.glut.DT.T2, y_predicted.glut.DT.T2,col="#00000050",pch = 19,main = "Train result for lasso")
abline(lm(y_predicted.glut.DT.T2~Y.glut.DT.T2),lty = 2,lwd = 2,col = "gray")





#Protein
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.prot.DT.T2, Y.prot.DT.T2, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.prot.DT.T2, s=s.best.lasso))[-1]
s.best.lasso

lasso.prot.DT.T2<-coef(fit.lasso.prot.DT.T2, s=s.best.lasso)
lasso.prot.DT.T2

#perform k-fold cross-validation to find optimal lambda value
cv_model.prot.DT.T2<- cv.glmnet(X.prot.DT.T2, Y.prot.DT.T2, k=10, nfolds = nrow(X.prot.DT.T2), alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.prot.DT.T2 <- cv_model.prot.DT.T2$lambda.min
best_lambda.prot.DT.T2

#produce plot of test MSE by lambda value
plot(cv_model.prot.DT.T2) 

#find coefficients of best model
best_model.prot.DT.T2 <- glmnet(X.prot.DT.T2, Y.prot.DT.T2, alpha = 1, lambda = best_lambda.prot.DT.T2, standardize = TRUE)
coef(best_model.prot.DT.T2)

#use fitted best model to make predictions
y_predicted.prot.DT.T2 <- predict(best_model.prot.DT.T2, s = best_lambda.prot.DT.T2, newx = X.prot.DT.T2)


#Evaluation of the accuracy
Rsquare <- cor(Y.prot.DT.T2,y_predicted.prot.DT.T2)^2
Rsquare

MSE <- mean((Y.prot.DT.T2 - y_predicted.prot.DT.T2)^2)
MSE

plot(Y.prot.DT.T2, y_predicted.prot.DT.T2,col="#00000050",pch = 19,main = "Test result for lasso")
abline(lm(y_predicted.prot.DT.T2~Y.prot.DT.T2),lty = 2,lwd = 2,col = "gray")



# plot lasso model
lam <- best_model.prot.DT.T2$lambda %>% 
  as.data.frame() %>%
  mutate(penalty = best_model.prot.DT.T2$a0 %>% names()) %>%
  rename(lambda = ".")

results.DT.prot.T2 <- best_model.prot.DT.T2$beta %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(penalty, coefficients, -rowname) %>%
  left_join(lam)

results.DT.prot.T2.filt<-filter_if(results.DT.prot.T2, is.numeric, all_vars((.) != 0))
results.DT.prot.T2.filt



# Visualize the predictive plot
write.table(results.DT.prot.T2.filt,here("output/tables", "results_labels_prot_DT.T2.txt"),sep="\t" )





#PMT
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.pmt.DT.T2, Y.pmt.DT.T2, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.pmt.DT.T2, s=s.best.lasso))[-1]
s.best.lasso

lasso.pmt.DT.T2<-coef(fit.lasso.pmt.DT.T2, s=s.best.lasso)
lasso.pmt.DT.T2

#perform k-fold cross-validation to find optimal lambda value
cv_model.pmt.DT.T2<- cv.glmnet(X.pmt.DT.T2, Y.pmt.DT.T2, k=10,nfolds=nrow(X.pmt.DT.T2), alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.pmt.DT.T2 <- cv_model.pmt.DT.T2$lambda.min
best_lambda.pmt.DT.T2

#produce plot of test MSE by lambda value
plot(cv_model.pmt.DT.T2) 

#find coefficients of best model
best_model.pmt.DT.T2 <- glmnet(X.pmt.DT.T2, Y.pmt.DT.T2, alpha = 1, lambda = best_lambda.pmt.DT.T2, standardize = TRUE)
coef(best_model.pmt.DT.T2)

#use fitted best model to make predictions
y_predicted.pmt.DT.T2 <- predict(best_model.pmt.DT.T2, s = best_lambda.pmt.DT.T2, newx = X.pmt.DT.T2)

predicted.value.pmt.DT.T2<-data.frame(cbind(y_predicted.pmt.DT.T2,Y.pmt.DT.T2))

#Evaluation of the accuracy
Rsquare<- cor(Y.pmt.DT.T2,y_predicted.pmt.DT.T2)^2
Rsquare

MSE<-mean((Y.pmt.DT.T2 - y_predicted.pmt.DT.T2)^2)
MSE




plot.pmt.DT.T2<-ggplot(predicted.value.pmt.DT.T2, aes(x=V2, y=s1)) + 
  geom_point(color="#00000050")+
  geom_smooth(method=lm, color="lightblue")+
  #geom_text(x = 25, y = 31, label = eq(predicted.value.glut.DT.T3$Y.glut.DT.T3,predicted.value.glut.DT.T3$s1), parse = TRUE) +
  
  stat_regline_equation(label.y =2, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 1.7, aes(label = ..rr.label..)) +
  labs(title="PMT (May 21)  ",
       y=" Predicted peak maximum time ", x = " Observed peak maximum time ")+
  theme_classic() 

plot.pmt.DT.T2

ggsave(file=here("output/photo","model.DT.T2.PMT.tiff"), plot.pmt.DT.T2, height=3.5, width=4.0, units="in", dpi=600)





#BEM
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.bem.DT.T2, Y.bem.DT.T2, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.bem.DT.T2, s=s.best.lasso))[-1]
s.best.lasso

lasso.bem.DT.T2<-coef(fit.lasso.bem.DT.T2, s=s.best.lasso)
lasso.bem.DT.T2


#perform k-fold cross-validation to find optimal lambda value
cv_model.bem.DT.T2<- cv.glmnet(X.bem.DT.T2, Y.bem.DT.T2, k=10, nfolds=nrow(X.bem.DT.T2), alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.bem.DT.T2 <- cv_model.bem.DT.T2$lambda.min
best_lambda.bem.DT.T2

#produce plot of test MSE by lambda value
plot(cv_model.bem.DT.T2) 

#find coefficients of best model
best_model.bem.DT.T2 <- glmnet(X.bem.DT.T2, Y.bem.DT.T2, alpha = 1, lambda = best_lambda.bem.DT.T2, standardize = TRUE)
coef(best_model.pmt.DT.T2)

#use fitted best model to make predictions
y_predicted.bem.DT.T2 <- predict(best_model.bem.DT.T2, s = best_lambda.bem.DT.T2, newx = X.bem.DT.T2)


#Evaluation of the accuracy
Rsquare <- cor(Y.bem.DT.T2, y_predicted.bem.DT.T2)^2
Rsquare

MSE <- mean((Y.bem.DT.T2 - y_predicted.bem.DT.T2)^2)
MSE





#lasso.DT.T1=as.matrix(lasso.DT.T1)
#lasso.df.DT.T1=as.data.frame(lasso.DT.T1)
#write.table(lasso.df.DT.T1,"lasso.df.protien.txt",sep="\t")


model.glut.DT.T2=lm(Gluten~biolog.DT.T2.axis4+PD.fun ,
                    data=data.glut.DT.T2)


model.protein.DT.T2=lm(Protein.grain~bact.DT.T2.axis5+biolog.DT.T2.axis2+Simpson.fun+Chao1.fun+AOA+F.B.ratio 
                       ,data=data.prot.DT.T2)

model.pmt.DT.T2=lm(PMT~biolog.DT.T2.axis2+ACE.fun  , data=data.pmt.DT.T2)

model.bem.DT.T2=lm(BEM~biolog.DT.T2.axis5+Simpson.fun , data=data.bem.DT.T2)


vif(model.protein.DT.T2)

summary(model.glut.DT.T2)
summary(model.protein.DT.T2)
summary(model.pmt.DT.T2)
summary(model.bem.DT.T2)


# LASSO regression analysis

# Preparing data for LASSO regression for T3
reg.lasso.DT.T3<-cbind(pca.data.point.16S.DT.T3,pca.data.point.ITS.DT.T3,pca.data.point.biolog.DT.T3,div.16S.DT.T3,div.ITS.DT.T3,qpcr.T3.DT,qual.DT.S3)
row.names(reg.lasso.DT.T3)==row.names(qual.DT.S3)


#reg.lasso.DT.T3<-decostand(reg.lasso.DT.T3,method="standardize")

#qpcr.T1.DT<-log(qpcr.T1.DT)

## Protein model for May 10
data.glut.DT.T3=reg.lasso.DT.T3[!is.na(reg.lasso.DT.T3$Gluten),]
data.glut.DT.T3<-data.glut.DT.T3[-c(5,18),]

data.prot.DT.T3=reg.lasso.DT.T3[!is.na(reg.lasso.DT.T3$Protein.grain),]
data.pmt.DT.T3=reg.lasso.DT.T3[!is.na(reg.lasso.DT.T3$PMT),]
data.bem.DT.T3=reg.lasso.DT.T3[!is.na(reg.lasso.DT.T3$BEM),]

data.explain.glut.DT.T3=data.glut.DT.T3[,-which(names(data.glut.DT.T3) %in% c("Protein.grain", "Proteine.flour","Humidity","Ash","Gluten", "IDC","PMT","BEM"))]
data.explain.prot.DT.T3=data.prot.DT.T3[,-which(names(data.prot.DT.T3) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.pmt.DT.T3=data.pmt.DT.T3[,-which(names(data.pmt.DT.T3) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.bem.DT.T3=data.bem.DT.T3[,-which(names(data.bem.DT.T3) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]

# X-Y data prep
X.glut.DT.T3=as.matrix(scale(data.explain.glut.DT.T3,center = T,scale = T))
#Y.glut.DT.T3=as.matrix(data.glut.DT.T3$Gluten)
Y.glut.DT.T3=as.matrix(scale(data.glut.DT.T3$Gluten,center = T,scale = T))

X.prot.DT.T3=as.matrix(scale(data.explain.prot.DT.T3,center = T,scale = T))
Y.prot.DT.T3=as.matrix(scale(data.prot.DT.T3$Protein.grain,center = T,scale = T))

X.pmt.DT.T3=as.matrix(scale(data.explain.pmt.DT.T3,center = T,scale = T))
Y.pmt.DT.T3=as.matrix(scale(data.pmt.DT.T3$PMT,center = T,scale = T))

X.bem.DT.T3=as.matrix(scale(data.explain.bem.DT.T3,center = T,scale = T))
Y.bem.DT.T3=as.matrix(scale(data.bem.DT.T3$BEM,center = T,scale = T))

#LASSO regression
fit.lasso.glut.DT.T3 <-glmnet(X.glut.DT.T3,Y.glut.DT.T3, family="gaussian", alpha=1)
fit.lasso.prot.DT.T3<- glmnet(X.prot.DT.T3,Y.prot.DT.T3, family="gaussian", alpha=1)
fit.lasso.pmt.DT.T3 <- glmnet(X.pmt.DT.T3,Y.pmt.DT.T3, family="gaussian", alpha=1)
fit.lasso.bem.DT.T3 <- glmnet(X.bem.DT.T3,Y.bem.DT.T3, family="gaussian", alpha=1)


#Gluten
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.glut.DT.T3, Y.glut.DT.T3, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.glut.DT.T3, s=s.best.lasso))[-1]
s.best.lasso


lasso.glut.DT.T3<-coef(fit.lasso.glut.DT.T3, s=s.best.lasso)
lasso.glut.DT.T3



plot(fit.lasso.glut.DT.T3,xvar="lambda",label=TRUE)

#perform k-fold cross-validation to find optimal lambda value
cv_model.glut.DT.T3<- cv.glmnet(X.glut.DT.T3, Y.glut.DT.T3, k=10,nfolds = nrow(X.glut.DT.T3), alpha = 1,standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.glut.DT.T3 <- cv_model.glut.DT.T3$lambda.min
best_lambda.glut.DT.T3

#produce plot of test MSE by lambda value
plot(cv_model.glut.DT.T3) 

#find coefficients of best model
best_model.glut.DT.T3 <- glmnet(X.glut.DT.T3, Y.glut.DT.T3, alpha = 1, lambda = best_lambda.glut.DT.T3, standardize = TRUE)
coef(best_model.glut.DT.T3)

#use fitted best model to make predictions
y_predicted.glut.DT.T3 <- predict(best_model.glut.DT.T3, s = best_lambda.glut.DT.T3, newx = X.glut.DT.T3)

predicted.value.glut.DT.T3<-data.frame(y_predicted.glut.DT.T3)
predicted.value.glut.DT.T3<-cbind(Y.glut.DT.T3,predicted.value.glut.DT.T3)
#Evaluation of the accuracy
Rsquare <- cor(Y.glut.DT.T3,y_predicted.glut.DT.T3)^2
Rsquare

MSE <- mean((Y.glut.DT.T3 - y_predicted.glut.DT.T3)^2)
MSE

plot(Y.glut.DT.T3, y_predicted.glut.DT.T3,col="#00000050",pch = 19,main = "Test result for lasso")
abline(lm(y_predicted.glut.DT.T3~Y.glut.DT.T3),lty = 2,lwd = 2,col = "gray")






# Visualize the predictive plot


eq <- function(x,y) {
  m <- lm(y ~ x)
  as.character(
    as.expression(
      substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                 list(a = format(coef(m)[1], digits = 4),
                      b = format(coef(m)[2], digits = 4),
                      r2 = format(summary(m)$r.squared, digits = 3)))
    )
  )
}




#Protein
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.prot.DT.T3, Y.prot.DT.T3, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.prot.DT.T3, s=s.best.lasso))[-1]
s.best.lasso

lasso.prot.DT.T3<-coef(fit.lasso.prot.DT.T3, s=s.best.lasso)
lasso.prot.DT.T3

set.seed(1)

#perform k-fold cross-validation to find optimal lambda value
cv_model.prot.DT.T3<- cv.glmnet(X.prot.DT.T3, Y.prot.DT.T3, nfolds=nrow(X.prot.DT.T3), alpha =1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.prot.DT.T3 <- cv_model.prot.DT.T3$lambda.min
best_lambda.prot.DT.T3

#produce plot of test MSE by lambda value
plot(cv_model.prot.DT.T3) 

#find coefficients of best model
best_model.prot.DT.T3 <- glmnet(X.prot.DT.T3, Y.prot.DT.T3, alpha = 1, lambda = best_lambda.prot.DT.T3, standardize = TRUE)
coef(best_model.prot.DT.T3)


#use fitted best model to make predictions
y_predicted.prot.DT.T3 <- predict(best_model.prot.DT.T3, s = best_lambda.prot.DT.T3, newx = X.prot.DT.T3)


#Evaluation of the accuracy
Rsquare <- cor(Y.prot.DT.T3,y_predicted.prot.DT.T3)^2
Rsquare

MSE <- mean((Y.prot.DT.T3 - y_predicted.prot.DT.T3)^2)
MSE

plot(Y.prot.DT.T3, y_predicted.prot.DT.T3,col="#00000050",pch = 19,main = "Test result for lasso")
abline(lm(y_predicted.prot.DT.T3~Y.prot.DT.T3),lty = 2,lwd = 2,col = "gray")

# plot lasso model
lam <- best_model.prot.DT.T3$lambda %>% 
  as.data.frame() %>%
  mutate(penalty = best_model.prot.DT.T3$a0 %>% names()) %>%
  rename(lambda = ".")

results.DT.prot.T3 <- best_model.prot.DT.T3$beta %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(penalty, coefficients, -rowname) %>%
  left_join(lam)

results.DT.prot.T3.filt<-filter_if(results.DT.prot.T3, is.numeric, all_vars((.) != 0))
results.DT.prot.T3.filt



# Visualize the predictive plot
write.table(results.DT.prot.T3.filt,file= here("output/tables","results_labels_prot_DT.T3.txt"),sep="\t" )






#PMT
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.pmt.DT.T3, Y.pmt.DT.T3, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.pmt.DT.T3, s=s.best.lasso))[-1]
s.best.lasso

lasso.pmt.DT.T3<-coef(fit.lasso.pmt.DT.T3, s=s.best.lasso)
lasso.pmt.DT.T3

#perform k-fold cross-validation to find optimal lambda value
cv_model.pmt.DT.T3<- cv.glmnet(X.pmt.DT.T3, Y.pmt.DT.T3, k=10,nfolds = nrow(X.pmt.DT.T3), alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.pmt.DT.T3 <- cv_model.pmt.DT.T3$lambda.min
best_lambda.pmt.DT.T3

#produce plot of test MSE by lambda value
plot(cv_model.pmt.DT.T3) 

#find coefficients of best model
best_model.pmt.DT.T3 <- glmnet(X.pmt.DT.T3, Y.pmt.DT.T3, alpha = 1, lambda = best_lambda.pmt.DT.T3, standardize = TRUE)
coef(best_model.pmt.DT.T3)


#use fitted best model to make predictions
y_predicted.pmt.DT.T3 <- predict(best_model.pmt.DT.T3, s = best_lambda.pmt.DT.T3, newx = X.pmt.DT.T3)


#Evaluation of the accuracy
Rsquare <- cor(Y.pmt.DT.T3,y_predicted.pmt.DT.T3)^2
Rsquare

MSE <- mean((Y.prot.DT.T3 - y_predicted.prot.DT.T3)^2)
MSE

#BEM
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.bem.DT.T3, Y.bem.DT.T3, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.bem.DT.T3, s=s.best.lasso))[-1]
s.best.lasso

lasso.bem.DT.T3<-coef(fit.lasso.bem.DT.T3, s=s.best.lasso)
lasso.bem.DT.T3

#perform k-fold cross-validation to find optimal lambda value
cv_model.bem.DT.T3<- cv.glmnet(X.bem.DT.T3, Y.bem.DT.T3,k=10, nfolds = nrow(X.bem.DT.T3), alpha = 1, standardize=TRUE, grouped=FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.bem.DT.T3 <- cv_model.bem.DT.T3$lambda.min
best_lambda.bem.DT.T3

#produce plot of test MSE by lambda value
plot(cv_model.bem.DT.T3) 

#find coefficients of best model
best_model.bem.DT.T3 <- glmnet(X.bem.DT.T3, Y.bem.DT.T3, alpha = 1, lambda = best_lambda.bem.DT.T3, standardize = TRUE)
coef(best_model.bem.DT.T3)


#use fitted best model to make predictions
y_predicted.bem.DT.T3 <- predict(best_model.bem.DT.T3, s = best_lambda.bem.DT.T3, newx = X.bem.DT.T3)
predicted.value.bem.DT.T3<-data.frame(cbind(y_predicted.bem.DT.T3, Y.bem.DT.T3))

#Evaluation of the accuracy
Rsquare <- cor(Y.bem.DT.T3,y_predicted.bem.DT.T3)^2
Rsquare

MSE <- mean((Y.bem.DT.T3 - y_predicted.bem.DT.T3)^2)
MSE


# plot lasso model
lam <- best_model.bem.DT.T3$lambda %>% 
  as.data.frame() %>%
  mutate(penalty = best_model.bem.DT.T3$a0 %>% names()) %>%
  rename(lambda = ".")

results.DT.bem.T3 <- best_model.bem.DT.T3$beta %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(penalty, coefficients, -rowname) %>%
  left_join(lam)       

results.DT.bem.T3.filt<-filter_if(results.DT.bem.T3, is.numeric, all_vars((.) != 0))
results.DT.bem.T3.filt

# Visualize the predictive plot
write.table(results.DT.bem.T3.filt,file=here("output/tables", "result_labels_bem_DT.T3.txt"),sep="\t" )

plot.bem.DT.T3<-ggplot(predicted.value.bem.DT.T3, aes(x=V2, y=s1)) + 
  geom_point(color="#00000050")+
  geom_smooth(method=lm, color="lightblue")+
  #geom_text(x = 25, y = 31, label = eq(predicted.value.glut.DT.T3$Y.glut.DT.T3,predicted.value.glut.DT.T3$s1), parse = TRUE) +
  
  stat_regline_equation(label.y =2, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y =1.7, aes(label = ..rr.label..)) +
  labs(title="BEM (June 07)",
       y="Predicted maximum torque", x = "Observed maximum torque")+
  theme_classic() 

plot.bem.DT.T3


ggsave(file = here("output/photo","model.DT.T3.BEM.tiff"), plot.bem.DT.T3, height=3.5, width=4.0, units="in", dpi=600)








#lasso.DT.T1=as.matrix(lasso.DT.T1)
#lasso.df.DT.T1=as.data.frame(lasso.DT.T1)
#write.table(lasso.df.DT.T1,"lasso.df.protien.txt",sep="\t")


model.glut.DT.T3=lm(Gluten~bact.DT.T3.axis3 +fun.DT.T3.axis2+fun.DT.T3.axis3+fun.DT.T3.axis4+ACE+Shannon.fun+PD.fun ,
                    data=data.glut.DT.T3)


model.protein.DT.T3=lm(Protein.grain~bact.DT.T3.axis3+biolog.DT.T3.axis4+biolog.DT.T3.axis5+ACE+PD.fun+AOA+F.B.ratio
                       ,data=data.prot.DT.T3)

#model.pmt.DT.T3=lm(PMT~......., data=data.pmt.DT.T3)

model.bem.DT.T3=lm(BEM~bact.DT.T3.axis2+bact.DT.T3.axis3+fun.DT.T3.axis2+fun.DT.T3.axis4+fun.DT.T3.axis5+biolog.DT.T3.axis2+biolog.DT.T3.axis4+Simpson+F.B.ratio   , data=data.bem.DT.T3)


vif(model.protein.DT.T3)

summary(model.glut.DT.T3)
summary(model.protein.DT.T3)
summary(model.bem.DT.T3)

# Preparing data for LASSO regression for T4

reg.lasso.DT.T4<-cbind(pca.data.point.16S.DT.T4,pca.data.point.ITS.DT.T4,pca.data.point.biolog.DT.T4,div.16S.DT.T4,div.ITS.DT.T4,qpcr.T4.DT,qual.DT.S4)
row.names(reg.lasso.DT.T4)==row.names(qual.DT.S4)

#reg.lasso.DT.T4<-decostand(reg.lasso.DT.T4,method="standardize")

#qpcr.T1.DT<-log(qpcr.T1.DT)

## Protein model for May 10
data.glut.DT.T4=reg.lasso.DT.T4[!is.na(reg.lasso.DT.T4$Gluten),]
data.glut.DT.T4=data.glut.DT.T4[-c(5.18),]

data.prot.DT.T4=reg.lasso.DT.T4[!is.na(reg.lasso.DT.T4$Protein.grain),]
data.pmt.DT.T4=reg.lasso.DT.T4[!is.na(reg.lasso.DT.T4$PMT),]
data.bem.DT.T4=reg.lasso.DT.T4[!is.na(reg.lasso.DT.T4$BEM),]

data.explain.glut.DT.T4=data.glut.DT.T4[,-which(names(data.glut.DT.T4) %in% c("Protein.grain", "Proteine.flour","Humidity","Ash","Gluten", "IDC","PMT","BEM"))]
data.explain.prot.DT.T4=data.prot.DT.T4[,-which(names(data.prot.DT.T4) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.pmt.DT.T4=data.pmt.DT.T4[,-which(names(data.pmt.DT.T4) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.bem.DT.T4=data.bem.DT.T4[,-which(names(data.bem.DT.T4) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]

# X-Y data prep
X.glut.DT.T4=as.matrix(scale(data.explain.glut.DT.T4,center = T,scale = T))
Y.glut.DT.T4=as.matrix(scale(data.glut.DT.T4$Gluten,center = T,scale = T))

X.prot.DT.T4=as.matrix(scale(data.explain.prot.DT.T4,center = T,scale = T))
Y.prot.DT.T4=as.matrix(scale(data.prot.DT.T4$Protein.grain,center = T,scale = T))

#X.pmt.DT.T4=as.matrix(data.explain.pmt.DT.T4)
X.pmt.DT.T4=as.matrix(scale(data.explain.pmt.DT.T4,center = T,scale = T))
Y.pmt.DT.T4=as.matrix(scale(data.pmt.DT.T4$PMT,center = T,scale = T))

X.bem.DT.T4=as.matrix(scale(data.explain.bem.DT.T4,center = T,scale = T))
Y.bem.DT.T4=as.matrix(scale(data.bem.DT.T4$BEM,center = T,scale = T))

#LASSO regression
fit.lasso.glut.DT.T4 <-glmnet(X.glut.DT.T4,Y.glut.DT.T4, family="gaussian", alpha=1)
fit.lasso.prot.DT.T4<- glmnet(X.prot.DT.T4,Y.prot.DT.T4, family="gaussian", alpha=1)
fit.lasso.pmt.DT.T4 <- glmnet(X.pmt.DT.T4,Y.pmt.DT.T4, family="gaussian", alpha=1)
fit.lasso.bem.DT.T4 <- glmnet(X.bem.DT.T4,Y.bem.DT.T4, family="gaussian", alpha=1)


#Gluten
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.glut.DT.T4, Y.glut.DT.T4, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.glut.DT.T4, s=s.best.lasso))[-1]
s.best.lasso

lasso.glut.DT.T4<-coef(fit.lasso.glut.DT.T4, s=s.best.lasso)
lasso.glut.DT.T4

#perform k-fold cross-validation to find optimal lambda value



cv_model.glut.DT.T4<- cv.glmnet(X.glut.DT.T4, Y.glut.DT.T4, k=10,  alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.glut.DT.T4 <- cv_model.glut.DT.T4$lambda.min
best_lambda.glut.DT.T4

#produce plot of test MSE by lambda value
plot(cv_model.glut.DT.T4) 

#find coefficients of best model
best_model.glut.DT.T4 <- glmnet(X.glut.DT.T4, Y.glut.DT.T4, alpha = 1, lambda = best_lambda.glut.DT.T4, standardize = TRUE)
coef(best_model.glut.DT.T4)


#use fitted best model to make predictions
y_predicted.glut.DT.T4<- predict(best_model.glut.DT.T4, s = best_lambda.glut.DT.T4, newx = X.glut.DT.T4)


#Evaluation of the accuracy
Rsquare <- cor(Y.glut.DT.T4,y_predicted.glut.DT.T4)^2
Rsquare

MSE <- mean((Y.glut.DT.T4 - y_predicted.glut.DT.T4)^2)
MSE






#Protein
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.prot.DT.T4, Y.prot.DT.T4, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.prot.DT.T4, s=s.best.lasso))[-1]
s.best.lasso

lasso.prot.DT.T4<-coef(fit.lasso.prot.DT.T4, s=s.best.lasso)
lasso.prot.DT.T4

set.seed (1)
#perform k-fold cross-validation to find optimal lambda value
cv_model.prot.DT.T4<- cv.glmnet(X.prot.DT.T4, Y.prot.DT.T4, nfolds=nrow(X.prot.DT.T4),alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.prot.DT.T4 <- cv_model.prot.DT.T4$lambda.min
best_lambda.prot.DT.T4

#produce plot of test MSE by lambda value
plot(cv_model.prot.DT.T4) 

#find coefficients of best model
best_model.prot.DT.T4 <- glmnet(X.prot.DT.T4, Y.prot.DT.T4, alpha = 1, lambda = best_lambda.prot.DT.T4, standardize = TRUE)
coef(best_model.prot.DT.T4)


#use fitted best model to make predictions
y_predicted.prot.DT.T4<- predict(best_model.prot.DT.T4, s = best_lambda.prot.DT.T4, newx = X.prot.DT.T4)


#Evaluation of the accuracy
Rsquare <- cor(Y.prot.DT.T4,y_predicted.prot.DT.T4)^2
Rsquare

MSE <- mean((Y.prot.DT.T4 - y_predicted.prot.DT.T4)^2)
MSE





#PMT
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.pmt.DT.T4, Y.pmt.DT.T4, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.pmt.DT.T4, s=s.best.lasso))[-1]
s.best.lasso

lasso.pmt.DT.T4<-coef(fit.lasso.pmt.DT.T4, s=s.best.lasso)
lasso.pmt.DT.T4

#perform k-fold cross-validation to find optimal lambda value
cv_model.pmt.DT.T4<- cv.glmnet(X.pmt.DT.T4, Y.pmt.DT.T4, k=10,nfolds=nrow(X.pmt.DT.T4), alpha = 1, standardize=TRUE, grouped =FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.pmt.DT.T4 <- cv_model.pmt.DT.T4$lambda.min
best_lambda.pmt.DT.T4

#produce plot of test MSE by lambda value
plot(cv_model.pmt.DT.T4) 

#find coefficients of best model
best_model.pmt.DT.T4 <- glmnet(X.pmt.DT.T4, Y.pmt.DT.T4, alpha = 1, lambda = best_lambda.pmt.DT.T4, standardize = TRUE)
coef(best_model.pmt.DT.T4)




#use fitted best model to make predictions
y_predicted.pmt.DT.T4<- predict(best_model.pmt.DT.T4, s = best_lambda.pmt.DT.T4, newx = X.pmt.DT.T4)
predicted.value.pmt.DT.T4<-data.frame(cbind(y_predicted.pmt.DT.T4,Y.pmt.DT.T4))

#Evaluation of the accuracy
Rsquare <- cor(Y.pmt.DT.T4,y_predicted.pmt.DT.T4)^2
Rsquare

MSE <- mean((Y.pmt.DT.T4 - y_predicted.pmt.DT.T4)^2)
MSE


# plot lasso model
lam <- best_model.pmt.DT.T4$lambda %>% 
  as.data.frame() %>%
  mutate(penalty = best_model.pmt.DT.T4$a0 %>% names()) %>%
  rename(lambda = ".")

results.DT.pmt.T4 <- best_model.pmt.DT.T4$beta %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(penalty, coefficients, -rowname) %>%
  left_join(lam)       

results.DT.pmt.T4.filt<-filter_if(results.DT.pmt.T4, is.numeric, all_vars((.) != 0))
results.DT.pmt.T4.filt

# Visualize the predictive plot
write.table(results.DT.pmt.T4.filt, file = here("output/tables","result_labels_pmt_DT.T4.txt"),sep="\t" )

plot.pmt.DT.T4<-ggplot(predicted.value.pmt.DT.T4, aes(x=V2, y=s1)) + 
  geom_point(color="#00000050")+
  geom_smooth(method=lm, color="lightblue")+
  #geom_text(x = 25, y = 31, label = eq(predicted.value.glut.DT.T3$Y.glut.DT.T3,predicted.value.glut.DT.T3$s1), parse = TRUE) +
  
  stat_regline_equation(label.y =1.9, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y =1.7, aes(label = ..rr.label..)) +
  labs(title="PMT (June 21)",
       x="Observed peak maximum time",y = "Predicted peak maximum time" )+
  theme_classic() 

plot.pmt.DT.T4


ggsave(file=here("output/photo","model.DT.T4.pmt.tiff"), plot.pmt.DT.T4, height=3.5, width=4.0, units="in", dpi=600)







#BEM
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.bem.DT.T4, Y.bem.DT.T4, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.bem.DT.T4, s=s.best.lasso))[-1]
s.best.lasso

lasso.bem.DT.T4<-coef(fit.lasso.bem.DT.T4, s=s.best.lasso)
lasso.bem.DT.T4

#perform k-fold cross-validation to find optimal lambda value
train.control <- trainControl(method = "repeatedcv", 
                              number = 10, repeats = 3)


cv_model.bem.DT.T4<- cv.glmnet(X.bem.DT.T4, Y.bem.DT.T4,k=10, nfolds  =nrow(X.bem.DT.T4), alpha = 1,standardize=TRUE, grouped = FALSE)



#find optimal lambda value that minimizes test MSE
best_lambda.bem.DT.T4 <- cv_model.bem.DT.T4$lambda.min
best_lambda.bem.DT.T4

#produce plot of test MSE by lambda value
plot(cv_model.bem.DT.T4) 

#find coefficients of best model
best_model.bem.DT.T4 <- glmnet(X.bem.DT.T4, Y.bem.DT.T4, alpha = 1, lambda = best_lambda.bem.DT.T4, standardize = TRUE)
coef(best_model.bem.DT.T4)


#use fitted best model to make predictions
y_predicted.bem.DT.T4<- predict(best_model.bem.DT.T4, s = best_lambda.bem.DT.T4, newx = X.bem.DT.T4)


#Evaluation of the accuracy
Rsquare <- cor(Y.bem.DT.T4,y_predicted.bem.DT.T4)^2
Rsquare

MSE <- mean((Y.bem.DT.T4 - y_predicted.bem.DT.T4)^2)
MSE



#lasso.DT.T1=as.matrix(lasso.DT.T1)
#lasso.df.DT.T1=as.data.frame(lasso.DT.T1)
#write.table(lasso.df.DT.T1,"lasso.df.protien.txt",sep="\t")


model.glut.DT.T4=lm(Gluten~bact.DT.T4.axis3+Simpson.fun  ,
                    data=data.glut.DT.T4)


model.protein.DT.T4=lm(Protein.grain~AOA+nirk
                       ,data=data.prot.DT.T4)

model.pmt.DT.T4=lm(PMT~fun.DT.T4.axis3+biolog.DT.T4.axis1+biolog.DT.T4.axis4, data=data.pmt.DT.T4)

model.bem.DT.T4=lm(BEM~Simpson.fun , data=data.bem.DT.T4)


vif(model.protein.DT.T3)

summary(model.glut.DT.T4)
summary(model.protein.DT.T4)
summary(model.pmt.DT.T4)
summary(model.bem.DT.T4)



# Preparing data for LASSO regression for T5

reg.lasso.DT.T5<-cbind(pca.data.point.16S.DT.T5,pca.data.point.ITS.DT.T5,pca.data.point.biolog.DT.T5,div.16S.DT.T5,div.ITS.DT.T5,qpcr.T5.DT,qual.DT.S5)
row.names(reg.lasso.DT.T5)==row.names(qual.DT.S5)

#reg.lasso.DT.T5<-decostand(reg.lasso.DT.T5,method="standardize")

#qpcr.T1.DT<-log(qpcr.T1.DT)

## lasso model for May 10
data.glut.DT.T5=reg.lasso.DT.T5[!is.na(reg.lasso.DT.T5$Gluten),]
data.glut.DT.T5=data.glut.DT.T5[-c(5.18),]

data.prot.DT.T5=reg.lasso.DT.T5[!is.na(reg.lasso.DT.T5$Protein.grain),]
data.pmt.DT.T5=reg.lasso.DT.T5[!is.na(reg.lasso.DT.T5$PMT),]
data.bem.DT.T5=reg.lasso.DT.T5[!is.na(reg.lasso.DT.T5$BEM),]

data.explain.glut.DT.T5=data.glut.DT.T5[,-which(names(data.glut.DT.T5) %in% c("Protein.grain", "Proteine.flour","Humidity","Ash","Gluten", "IDC","PMT","BEM"))]
data.explain.prot.DT.T5=data.prot.DT.T5[,-which(names(data.prot.DT.T5) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.pmt.DT.T5=data.pmt.DT.T5[,-which(names(data.pmt.DT.T5) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.bem.DT.T5=data.bem.DT.T5[,-which(names(data.bem.DT.T5) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]

# X-Y data prep
X.glut.DT.T5=as.matrix(scale(data.explain.glut.DT.T5,center = T,scale = T))
Y.glut.DT.T5=as.matrix(scale(data.glut.DT.T5$Gluten,center = T,scale = T))

X.prot.DT.T5=as.matrix(scale(data.explain.prot.DT.T5,center = T,scale = T))
Y.prot.DT.T5=as.matrix(scale(data.prot.DT.T5$Protein.grain,center = T,scale = T))

X.pmt.DT.T5=as.matrix(scale(data.explain.pmt.DT.T5,center = T,scale = T))
Y.pmt.DT.T5=as.matrix(scale(data.pmt.DT.T5$PMT,center = T,scale = T))

X.bem.DT.T5=as.matrix(scale(data.explain.bem.DT.T5,center = T,scale = T))
Y.bem.DT.T5=as.matrix(scale(data.bem.DT.T5$BEM,center = T,scale = T))

#LASSO regression
fit.lasso.glut.DT.T5 <-glmnet(X.glut.DT.T5,Y.glut.DT.T5, family="gaussian", alpha=1)
fit.lasso.prot.DT.T5<- glmnet(X.prot.DT.T5,Y.prot.DT.T5, family="gaussian", alpha=1)
fit.lasso.pmt.DT.T5 <- glmnet(X.pmt.DT.T5,Y.pmt.DT.T5, family="gaussian", alpha=1)
fit.lasso.bem.DT.T5 <- glmnet(X.bem.DT.T5,Y.bem.DT.T5, family="gaussian", alpha=1)


#Gluten
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.glut.DT.T5, Y.glut.DT.T5, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.glut.DT.T5, s=s.best.lasso))[-1]
s.best.lasso

lasso.glut.DT.T5<-coef(fit.lasso.glut.DT.T5, s=s.best.lasso)
lasso.glut.DT.T5


#perform k-fold cross-validation to find optimal lambda value

cv_model.glut.DT.T5<- cv.glmnet(X.glut.DT.T5, Y.glut.DT.T5, k=10, trControl=train.control, alpha = 1, standardize=TRUE,grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.glut.DT.T5 <- cv_model.glut.DT.T5$lambda.min
best_lambda.glut.DT.T5

#produce plot of test MSE by lambda value
plot(cv_model.glut.DT.T5) 

#find coefficients of best model
best_model.glut.DT.T5 <- glmnet(X.glut.DT.T5, Y.glut.DT.T5, alpha = 1, lambda = best_lambda.glut.DT.T5, standardize = TRUE)
coef(best_model.glut.DT.T5)


#use fitted best model to make predictions
y_predicted.glut.DT.T5<- predict(best_model.glut.DT.T5, s = best_lambda.glut.DT.T5, newx = X.glut.DT.T5)


#Evaluation of the accuracy
Rsquare <- cor(Y.glut.DT.T5,y_predicted.glut.DT.T5)^2
Rsquare

MSE <- mean((Y.glut.DT.T5 - y_predicted.glut.DT.T5)^2)
MSE



#Protein
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.prot.DT.T5, Y.prot.DT.T5, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.prot.DT.T5, s=s.best.lasso))[-1]
s.best.lasso

lasso.prot.DT.T5<-coef(fit.lasso.prot.DT.T5, s=s.best.lasso)
lasso.prot.DT.T5


#perform k-fold cross-validation to find optimal lambda value

cv_model.prot.DT.T5<- cv.glmnet(X.prot.DT.T5, Y.prot.DT.T5, k=10,nfolds=nrow(X.prot.DT.T5), trControl=train.control, alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.prot.DT.T5 <- cv_model.prot.DT.T5$lambda.min
best_lambda.prot.DT.T5

#produce plot of test MSE by lambda value
plot(cv_model.prot.DT.T5) 

#find coefficients of best model
best_model.prot.DT.T5 <- glmnet(X.prot.DT.T5, Y.prot.DT.T5, alpha = 1, lambda = best_lambda.prot.DT.T5, standardize = TRUE)
coef(best_model.prot.DT.T5)


#use fitted best model to make predictions
y_predicted.prot.DT.T5<- predict(best_model.prot.DT.T5, best_lambda.prot.DT.T5, newx = X.prot.DT.T5)


#Evaluation of the accuracy
Rsquare <- cor(Y.prot.DT.T5,y_predicted.prot.DT.T5)^2
Rsquare

MSE <- mean((Y.prot.DT.T5 - y_predicted.prot.DT.T5)^2)
MSE



#PMT
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.pmt.DT.T5, Y.pmt.DT.T5, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.pmt.DT.T5, s=s.best.lasso))[-1]
s.best.lasso

lasso.pmt.DT.T5<-coef(fit.lasso.pmt.DT.T5, s=s.best.lasso)
lasso.pmt.DT.T5


#perform k-fold cross-validation to find optimal lambda value

cv_model.pmt.DT.T5<- cv.glmnet(X.pmt.DT.T5, Y.pmt.DT.T5, k=10,nfolds = nrow(X.pmt.DT.T5), alpha = 1, standardize=TRUE,grouped=FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.pmt.DT.T5 <- cv_model.pmt.DT.T5$lambda.min
best_lambda.pmt.DT.T5

#produce plot of test MSE by lambda value
plot(cv_model.pmt.DT.T5) 

#find coefficients of best model
best_model.pmt.DT.T5 <- glmnet(X.pmt.DT.T5, Y.pmt.DT.T5, alpha = 1, lambda = best_lambda.pmt.DT.T5, standardize = TRUE)
coef(best_model.pmt.DT.T5)


#use fitted best model to make predictions
y_predicted.pmt.DT.T5<- predict(best_model.pmt.DT.T5, best_lambda.pmt.DT.T5, newx = X.pmt.DT.T5)


#Evaluation of the accuracy
Rsquare <- cor(Y.pmt.DT.T5,y_predicted.pmt.DT.T5)^2
Rsquare

MSE <- mean((Y.pmt.DT.T5 - y_predicted.pmt.DT.T5)^2)
MSE



#BEM
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.bem.DT.T5, Y.bem.DT.T5, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.bem.DT.T5, s=s.best.lasso))[-1]
s.best.lasso

lasso.bem.DT.T5<-coef(fit.lasso.bem.DT.T5, s=s.best.lasso)
lasso.bem.DT.T5


#perform k-fold cross-validation to find optimal lambda value

cv_model.bem.DT.T5<- cv.glmnet(X.bem.DT.T5, Y.bem.DT.T5, k=10, nfolds = nrow(X.bem.DT.T5), alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.bem.DT.T5 <- cv_model.bem.DT.T5$lambda.min
best_lambda.bem.DT.T5

#produce plot of test MSE by lambda value
plot(cv_model.bem.DT.T5) 

#find coefficients of best model
best_model.bem.DT.T5 <- glmnet(X.bem.DT.T5, Y.bem.DT.T5, alpha = 1, lambda = best_lambda.bem.DT.T5, standardize = TRUE)
coef(best_model.bem.DT.T5)


#use fitted best model to make predictions
y_predicted.bem.DT.T5<- predict(best_model.bem.DT.T5, best_lambda.bem.DT.T5, newx = X.bem.DT.T5)


#Evaluation of the accuracy
Rsquare <- cor(Y.bem.DT.T5,y_predicted.bem.DT.T5)^2
Rsquare

MSE <- mean((Y.bem.DT.T5 - y_predicted.bem.DT.T5)^2)
MSE



#lasso.DT.T1=as.matrix(lasso.DT.T1)
#lasso.df.DT.T1=as.data.frame(lasso.DT.T1)
#write.table(lasso.df.DT.T1,"lasso.df.protien.txt",sep="\t")


model.glut.DT.T5=lm(Gluten~OTU.richness.fun  ,
                    data=data.glut.DT.T5)


model.protein.DT.T5=lm(Protein.grain~fun.DT.T5.axis1+biolog.DT.T5.axis3
                       ,data=data.prot.DT.T5)

model.pmt.DT.T5=lm(PMT~bact.DT.T5.axis4+nosZ, data=data.pmt.DT.T5)

model.bem.DT.T5=lm(BEM~ biolog.DT.T5.axis5, data=data.bem.DT.T5)


vif(model.protein.DT.T3)

summary(model.glut.DT.T5)
summary(model.protein.DT.T5)
summary(model.pmt.DT.T5)
summary(model.bem.DT.T5)



# Preparing data for LASSO regression for T6

reg.lasso.DT.T6<-cbind(pca.data.point.16S.DT.T6,pca.data.point.ITS.DT.T6,pca.data.point.biolog.DT.T6,div.16S.DT.T6,div.ITS.DT.T6,qpcr.T6.DT,qual.DT.S6)
row.names(reg.lasso.DT.T6)==row.names(qual.DT.S6)

#reg.lasso.DT.T6<-decostand(reg.lasso.DT.T6,method="standardize")

#qpcr.T1.DT<-log(qpcr.T1.DT)

## lasso model for May 10
data.glut.DT.T6=reg.lasso.DT.T6[!is.na(reg.lasso.DT.T6$Gluten),]
data.glut.DT.T6=data.glut.DT.T6[-c(5,18),]


data.prot.DT.T6=reg.lasso.DT.T6[!is.na(reg.lasso.DT.T6$Protein.grain),]
data.pmt.DT.T6=reg.lasso.DT.T6[!is.na(reg.lasso.DT.T6$PMT),]
data.bem.DT.T6=reg.lasso.DT.T6[!is.na(reg.lasso.DT.T6$BEM),]

data.explain.glut.DT.T6=data.glut.DT.T6[,-which(names(data.glut.DT.T6) %in% c("Protein.grain", "Proteine.flour","Humidity","Ash","Gluten", "IDC","PMT","BEM"))]
data.explain.prot.DT.T6=data.prot.DT.T6[,-which(names(data.prot.DT.T6) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.pmt.DT.T6=data.pmt.DT.T6[,-which(names(data.pmt.DT.T6) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.bem.DT.T6=data.bem.DT.T6[,-which(names(data.bem.DT.T6) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]

# X-Y data prep
X.glut.DT.T6=as.matrix(scale(data.explain.glut.DT.T6,center=T,scale=T))
Y.glut.DT.T6=as.matrix(scale(data.glut.DT.T6$Gluten,center=T,scale=T))

X.prot.DT.T6=as.matrix(scale(data.explain.prot.DT.T6,center=T,scale=T))
Y.prot.DT.T6=as.matrix(scale(data.prot.DT.T6$Protein.grain,center=T,scale=T))

X.pmt.DT.T6=as.matrix(scale(data.explain.pmt.DT.T6,center=T,scale=T))
Y.pmt.DT.T6=as.matrix(scale(data.pmt.DT.T6$PMT,center=T,scale=T))

X.bem.DT.T6=as.matrix(scale(data.explain.bem.DT.T6,center=T,scale=T))
Y.bem.DT.T6=as.matrix(scale(data.bem.DT.T6$BEM,center=T,scale=T))

#LASSO regression
fit.lasso.glut.DT.T6 <-glmnet(X.glut.DT.T6,Y.glut.DT.T6, family="gaussian", alpha=1)
fit.lasso.prot.DT.T6<- glmnet(X.prot.DT.T6,Y.prot.DT.T6, family="gaussian", alpha=1)
fit.lasso.pmt.DT.T6 <- glmnet(X.pmt.DT.T6,Y.pmt.DT.T6, family="gaussian", alpha=1)
fit.lasso.bem.DT.T6 <- glmnet(X.bem.DT.T6,Y.bem.DT.T6, family="gaussian", alpha=1)


#Gluten
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.glut.DT.T6, Y.glut.DT.T6, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.glut.DT.T6, s=s.best.lasso))[-1]
s.best.lasso

lasso.glut.DT.T6<-coef(fit.lasso.glut.DT.T6, s=s.best.lasso)
lasso.glut.DT.T6

#perform k-fold cross-validation to find optimal lambda value

cv_model.glut.DT.T6<- cv.glmnet(X.glut.DT.T6, Y.glut.DT.T6, k=10, trControl=train.control, alpha = 1, standardize=TRUE,grouped =FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.glut.DT.T6 <- cv_model.glut.DT.T6$lambda.min
best_lambda.glut.DT.T6

#produce plot of test MSE by lambda value
plot(cv_model.glut.DT.T6) 

#find coefficients of best model
best_model.glut.DT.T6 <- glmnet(X.glut.DT.T6, Y.glut.DT.T6, alpha = 1, lambda = best_lambda.glut.DT.T6, standardize = TRUE)
coef(best_model.glut.DT.T6)


#use fitted best model to make predictions
y_predicted.glut.DT.T6<- predict(best_model.glut.DT.T6, s = best_lambda.glut.DT.T6, newx = X.glut.DT.T6)


#Evaluation of the accuracy
Rsquare <- cor(Y.glut.DT.T6,y_predicted.glut.DT.T6)^2
Rsquare

MSE <- mean((Y.glut.DT.T6 - y_predicted.glut.DT.T6)^2)
MSE

set.seed(1)

#Protein
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.prot.DT.T6, Y.prot.DT.T6, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.prot.DT.T6, s=s.best.lasso))[-1]
s.best.lasso

lasso.prot.DT.T6<-coef(fit.lasso.prot.DT.T6, s=s.best.lasso)
lasso.prot.DT.T6

#perform k-fold cross-validation to find optimal lambda value

cv_model.prot.DT.T6<- cv.glmnet(X.prot.DT.T6, Y.prot.DT.T6, nfolds=10, alpha = 1, standardize=TRUE, grouped=FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.prot.DT.T6 <- cv_model.prot.DT.T6$lambda.min
best_lambda.prot.DT.T6

#produce plot of test MSE by lambda value
plot(cv_model.prot.DT.T6) 

#find coefficients of best model
best_model.prot.DT.T6 <- glmnet(X.prot.DT.T6, Y.prot.DT.T6, alpha = 1, lambda = best_lambda.prot.DT.T6, standardize = TRUE)
coef(best_model.prot.DT.T6)


#use fitted best model to make predictions
y_predicted.prot.DT.T6<- predict(best_model.prot.DT.T6, s = best_lambda.prot.DT.T6, newx = X.prot.DT.T6)


#Evaluation of the accuracy
Rsquare <- cor(Y.prot.DT.T6,y_predicted.prot.DT.T6)^2
Rsquare

MSE <- mean((Y.prot.DT.T6 - y_predicted.prot.DT.T6)^2)
MSE


#PMT
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.pmt.DT.T6, Y.pmt.DT.T6, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.pmt.DT.T6, s=s.best.lasso))[-1]
s.best.lasso

lasso.pmt.DT.T6<-coef(fit.lasso.pmt.DT.T6, s=s.best.lasso)
lasso.pmt.DT.T6


#perform k-fold cross-validation to find optimal lambda value

cv_model.pmt.DT.T6<- cv.glmnet(X.pmt.DT.T6, Y.pmt.DT.T6, nfolds = 10,  alpha = 1, standardize=TRUE, grouped=FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.pmt.DT.T6 <- cv_model.pmt.DT.T6$lambda.min
best_lambda.pmt.DT.T6

#produce plot of test MSE by lambda value
plot(cv_model.pmt.DT.T6) 

#find coefficients of best model
best_model.pmt.DT.T6 <- glmnet(X.pmt.DT.T6, Y.pmt.DT.T6, alpha = 1, lambda = best_lambda.pmt.DT.T6, standardize = TRUE)
coef(best_model.pmt.DT.T6)


#use fitted best model to make predictions
y_predicted.pmt.DT.T6<- predict(best_model.pmt.DT.T6, s = best_lambda.pmt.DT.T6, newx = X.pmt.DT.T6)


#Evaluation of the accuracy
Rsquare <- cor(Y.pmt.DT.T6,y_predicted.pmt.DT.T6)^2
Rsquare

MSE <- mean((Y.pmt.DT.T6 - y_predicted.pmt.DT.T6)^2)
MSE



#BEM
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.bem.DT.T6, Y.bem.DT.T6, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.bem.DT.T6, s=s.best.lasso))[-1]
s.best.lasso

lasso.bem.DT.T6<-coef(fit.lasso.bem.DT.T6, s=s.best.lasso)
lasso.bem.DT.T6

#perform k-fold cross-validation to find optimal lambda value

cv_model.bem.DT.T6<- cv.glmnet(X.bem.DT.T6, Y.bem.DT.T6, nfolds =nrow(X.bem.DT.T6), trControl=train.control, alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.bem.DT.T6 <- cv_model.bem.DT.T6$lambda.min
best_lambda.bem.DT.T6

#produce plot of test MSE by lambda value
plot(cv_model.bem.DT.T6) 

#find coefficients of best model
best_model.bem.DT.T6 <- glmnet(X.bem.DT.T6, Y.bem.DT.T6, alpha = 1, lambda = best_lambda.bem.DT.T6, standardize = TRUE)
coef(best_model.bem.DT.T6)


#use fitted best model to make predictions
y_predicted.bem.DT.T6<- predict(best_model.bem.DT.T6, s = best_lambda.bem.DT.T6, newx = X.bem.DT.T6)


#Evaluation of the accuracy
Rsquare <- cor(Y.bem.DT.T6,y_predicted.bem.DT.T6)^2
Rsquare

MSE <- mean((Y.bem.DT.T6 - y_predicted.bem.DT.T6)^2)
MSE


# plot lasso model
lam <- best_model.bem.DT.T6$lambda %>% 
  as.data.frame() %>%
  mutate(penalty = best_model.bem.DT.T6$a0 %>% names()) %>%
  rename(lambda = ".")

results.DT.bem.T6 <- best_model.bem.DT.T6$beta %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(penalty, coefficients, -rowname) %>%
  left_join(lam)

results.DT.bem.T6.filt<-filter_if(results.DT.bem.T6, is.numeric, all_vars((.) != 0))
results.DT.bem.T6.filt



# Visualize the predictive plot
write.table(results.DT.bem.T6.filt, file=here("output/tables","results_labels_bem_DT.T6.txt"),sep="\t" )


#lasso.DT.T1=as.matrix(lasso.DT.T1)
#lasso.df.DT.T1=as.data.frame(lasso.DT.T1)
#write.table(lasso.df.DT.T1,"lasso.df.protien.txt",sep="\t")


model.glut.DT.T6=lm(Gluten~bact.DT.T6.axis5+biolog.DT.T6.axis3  ,
                    data=data.glut.DT.T6)


model.protein.DT.T6=lm(Protein.grain~nirk 
                       ,data=data.prot.DT.T6)

model.pmt.DT.T6=lm(PMT~fun.DT.T6.axis4+biolog.DT.T6.axis2+Simpson.fun+nirk  , data=data.pmt.DT.T6)

model.bem.DT.T6=lm(BEM~bact.DT.T6.axis4+bact.DT.T6.axis5+Simpson+Simpson.fun    , data=data.bem.DT.T6)


vif(model.protein.DT.T3)

summary(model.glut.DT.T6)
summary(model.protein.DT.T6)
summary(model.pmt.DT.T6)
summary(model.bem.DT.T6)



# Preparing data for LASSO regression for T7

reg.lasso.DT.T7<-cbind(pca.data.point.16S.DT.T7,pca.data.point.ITS.DT.T7,pca.data.point.biolog.DT.T7,div.16S.DT.T7,div.ITS.DT.T7,qpcr.T7.DT,qual.DT.S7)
row.names(reg.lasso.DT.T7)==row.names(qual.DT.S7)

#reg.lasso.DT.T7<-decostand(reg.lasso.DT.T7,method="standardize")

#qpcr.T1.DT<-log(qpcr.T1.DT)

## lasso model for May 10
data.glut.DT.T7=reg.lasso.DT.T7[!is.na(reg.lasso.DT.T7$Gluten),]
data.glut.DT.T7=data.glut.DT.T7[-c(5,18),]

data.prot.DT.T7=reg.lasso.DT.T7[!is.na(reg.lasso.DT.T7$Protein.grain),]
data.pmt.DT.T7=reg.lasso.DT.T7[!is.na(reg.lasso.DT.T7$PMT),]
data.bem.DT.T7=reg.lasso.DT.T7[!is.na(reg.lasso.DT.T7$BEM),]

data.explain.glut.DT.T7=data.glut.DT.T7[,-which(names(data.glut.DT.T7) %in% c("Protein.grain", "Proteine.flour","Humidity","Ash","Gluten", "IDC","PMT","BEM"))]
data.explain.prot.DT.T7=data.prot.DT.T7[,-which(names(data.glut.DT.T7) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.pmt.DT.T7=data.pmt.DT.T7[,-which(names(data.pmt.DT.T7) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.bem.DT.T7=data.bem.DT.T7[,-which(names(data.bem.DT.T7) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]

# X-Y data prep
X.glut.DT.T7=as.matrix(scale(data.explain.glut.DT.T7,center =T ,scale = T))
Y.glut.DT.T7=as.matrix(scale(data.glut.DT.T7$Gluten,center =T ,scale = T))

X.prot.DT.T7=as.matrix(scale(data.explain.prot.DT.T7,center =T ,scale = T))
Y.prot.DT.T7=as.matrix(scale(data.prot.DT.T7$Protein.grain,center =T ,scale = T))

X.pmt.DT.T7=as.matrix(scale(data.explain.pmt.DT.T7,center =T ,scale = T))
Y.pmt.DT.T7=as.matrix(scale(data.pmt.DT.T7$PMT,center =T ,scale = T))

X.bem.DT.T7=as.matrix(scale(data.explain.bem.DT.T7,center =T ,scale = T))
Y.bem.DT.T7=as.matrix(scale(data.bem.DT.T7$BEM,center =T ,scale = T))

#LASSO regression
fit.lasso.glut.DT.T7 <-glmnet(X.glut.DT.T7,Y.glut.DT.T7, family="gaussian", alpha=1)
fit.lasso.prot.DT.T7<- glmnet(X.prot.DT.T7,Y.prot.DT.T7, family="gaussian", alpha=1)
fit.lasso.pmt.DT.T7 <- glmnet(X.pmt.DT.T7,Y.pmt.DT.T7, family="gaussian", alpha=1)
fit.lasso.bem.DT.T7 <- glmnet(X.bem.DT.T7,Y.bem.DT.T7, family="gaussian", alpha=1)


#Gluten
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.glut.DT.T7, Y.glut.DT.T7, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.glut.DT.T7, s=s.best.lasso))[-1]
s.best.lasso

lasso.glut.DT.T7<-coef(fit.lasso.glut.DT.T7, s=s.best.lasso)
lasso.glut.DT.T7


#perform k-fold cross-validation to find optimal lambda value

cv_model.glut.DT.T7<- cv.glmnet(X.glut.DT.T7, Y.glut.DT.T7, nfolds  =10, trControl=train.control, alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.glut.DT.T7 <- cv_model.glut.DT.T7$lambda.min
best_lambda.glut.DT.T7

#produce plot of test MSE by lambda value
plot(cv_model.glut.DT.T7) 

#find coefficients of best model
best_model.glut.DT.T7 <- glmnet(X.glut.DT.T7, Y.glut.DT.T7, alpha = 1, lambda = best_lambda.glut.DT.T7, standardize = TRUE)
coef(best_model.glut.DT.T7)


#use fitted best model to make predictions
y_predicted.glut.DT.T7<- predict(best_model.glut.DT.T7, best_lambda.glut.DT.T7, newx = X.glut.DT.T7)


#Evaluation of the accuracy
Rsquare <- cor(Y.glut.DT.T7,y_predicted.glut.DT.T7)^2
Rsquare

MSE <- mean((Y.glut.DT.T7 - y_predicted.glut.DT.T7)^2)
MSE






#Protein
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.prot.DT.T7, Y.prot.DT.T7, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.prot.DT.T7, s=s.best.lasso))[-1]
s.best.lasso

lasso.prot.DT.T7<-coef(fit.lasso.prot.DT.T7, s=s.best.lasso)
lasso.prot.DT.T7

#perform k-fold cross-validation to find optimal lambda value

cv_model.prot.DT.T7<- cv.glmnet(X.prot.DT.T7, Y.prot.DT.T7, nfolds =nrow(X.prot.DT.T7),  alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.prot.DT.T7 <- cv_model.prot.DT.T7$lambda.min
best_lambda.prot.DT.T7

#produce plot of test MSE by lambda value
plot(cv_model.prot.DT.T7) 

#find coefficients of best model
best_model.prot.DT.T7 <- glmnet(X.prot.DT.T7, Y.prot.DT.T7, alpha = 1, lambda = best_lambda.prot.DT.T7, standardize = TRUE)
coef(best_model.prot.DT.T7)


#use fitted best model to make predictions
y_predicted.prot.DT.T7<- predict(best_model.prot.DT.T7, best_lambda.prot.DT.T7, newx = X.prot.DT.T7)


#Evaluation of the accuracy
Rsquare <- cor(Y.prot.DT.T7,y_predicted.prot.DT.T7)^2
Rsquare

MSE <- mean((Y.prot.DT.T7 - y_predicted.prot.DT.T7)^2)
MSE




#PMT
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.pmt.DT.T7, Y.pmt.DT.T7, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.pmt.DT.T7, s=s.best.lasso))[-1]
s.best.lasso

lasso.pmt.DT.T7<-coef(fit.lasso.pmt.DT.T7, s=s.best.lasso)
lasso.pmt.DT.T7

#perform k-fold cross-validation to find optimal lambda value

cv_model.pmt.DT.T7<- cv.glmnet(X.pmt.DT.T7, Y.pmt.DT.T7,k=10, nfolds=nrow(X.pmt.DT.T7), alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.pmt.DT.T7 <- cv_model.pmt.DT.T7$lambda.min
best_lambda.pmt.DT.T7

#produce plot of test MSE by lambda value
plot(cv_model.pmt.DT.T7) 

#find coefficients of best model
best_model.pmt.DT.T7 <- glmnet(X.pmt.DT.T7, Y.pmt.DT.T7, alpha = 1, lambda = best_lambda.pmt.DT.T7, standardize = TRUE)
coef(best_model.pmt.DT.T7)


#use fitted best model to make predictions
y_predicted.pmt.DT.T7<- predict(best_model.pmt.DT.T7, best_lambda.pmt.DT.T7, newx = X.pmt.DT.T7)


#Evaluation of the accuracy
Rsquare <- cor(Y.pmt.DT.T7,y_predicted.pmt.DT.T7)^2
Rsquare

MSE <- mean((Y.pmt.DT.T7 - y_predicted.pmt.DT.T7)^2)
MSE

#BEM
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.bem.DT.T7, Y.bem.DT.T7, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.bem.DT.T7, s=s.best.lasso))[-1]
s.best.lasso

lasso.bem.DT.T7<-coef(fit.lasso.bem.DT.T7, s=s.best.lasso)
lasso.bem.DT.T7


#perform k-fold cross-validation to find optimal lambda value

cv_model.bem.DT.T7<- cv.glmnet(X.bem.DT.T7, Y.bem.DT.T7, k=10, nfolds = nrow(X.bem.DT.T7), alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.bem.DT.T7 <- cv_model.bem.DT.T7$lambda.min
best_lambda.bem.DT.T7

#produce plot of test MSE by lambda value
plot(cv_model.bem.DT.T7) 

#find coefficients of best model
best_model.bem.DT.T7 <- glmnet(X.bem.DT.T7, Y.bem.DT.T7, alpha = 1, lambda = best_lambda.bem.DT.T7, standardize = TRUE)
coef(best_model.bem.DT.T7)


#use fitted best model to make predictions
y_predicted.bem.DT.T7<- predict(best_model.bem.DT.T7, best_lambda.bem.DT.T7, newx = X.bem.DT.T7)


#Evaluation of the accuracy
Rsquare <- cor(Y.bem.DT.T7,y_predicted.bem.DT.T7)^2
Rsquare

MSE <- mean((Y.bem.DT.T7 - y_predicted.bem.DT.T7)^2)
MSE




#lasso.DT.T1=as.matrix(lasso.DT.T1)
#lasso.df.DT.T1=as.data.frame(lasso.DT.T1)
#write.table(lasso.df.DT.T1,"lasso.df.protien.txt",sep="\t")


model.glut.DT.T7=lm(Gluten~fun.DT.T7.axis2  ,
                    data=data.glut.DT.T7)


model.protein.DT.T7=lm(Protein.grain~bact.DT.T7.axis4+fun.DT.T7.axis1 +Simpson+Simpson.fun  
                       ,data=data.prot.DT.T7)

model.pmt.DT.T7=lm(PMT~nosZ, data=data.pmt.DT.T7)

model.bem.DT.T7=lm(BEM~bact.DT.T7.axis5+fun.DT.T7.axis2+fun.DT.T7.axis5+biolog.DT.T7.axis3 , data=data.bem.DT.T7)


vif(model.protein.DT.T7)

summary(model.glut.DT.T7)
summary(model.protein.DT.T7)
summary(model.pmt.DT.T7)
summary(model.bem.DT.T7)


################################################################################
################################################################################


#Lasso modelling for Drought sensitive  Cultivar

#PCA analysis for 10th May
pca.16S.DS.T1 <- PCA(otu.16S.DS.T1.hel, graph = FALSE)
print(pca.16S.DS.T1)

##Cumulative percentage of eigenvalues
fviz_eig(pca.16S.DS.T1, addlabels = TRUE, ylim = c(0, 7.5))
var.16S.DS.T1 <- get_pca_var(pca.16S.DS.T1)

# Extraction of PCA points
pca.point.16S.DS.T1=data.frame(pca.16S.DS.T1$ind)
pca.data.point.16S.DS.T1=pca.point.16S.DS.T1[,which(colnames(pca.point.16S.DS.T1) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.16S.DS.T1)[which(colnames(pca.data.point.16S.DS.T1) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("bact.DS.T1.axis1", "bact.DS.T1.axis2","bact.DS.T1.axis3","bact.DS.T1.axis4","bact.DS.T1.axis5") 

#Analysis of significant dimensions

res.desc.16S.DS.T1 <- dimdesc(pca.16S.DS.T1, axes = c(1,2,3,4,5), proba = 0.05)

# Contributions of variables to PC1
fviz_contrib(pca.16S.DS.T1, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(pca.16S.DS.T1, choice = "var", axes = 2, top = 10)
#Total contribution to PC1 and PC2
fviz_contrib(pca.16S.DS.T1, choice = "var", axes = 1:2, top = 10)

#Axis4

PCA.axis1.16S.DS.T1<-data.frame(res.desc.16S.DS.T1$Dim.1)
PCA.axis2.16S.DS.T1<-data.frame(res.desc.16S.DS.T1$Dim.2)
PCA.axis3.16S.DS.T1<-data.frame(res.desc.16S.DS.T1$Dim.3)
PCA.axis4.16S.DS.T1<-data.frame(res.desc.16S.DS.T1$Dim.4)
PCA.axis5.16S.DS.T1<-data.frame(res.desc.16S.DS.T1$Dim.5)


write.table(PCA.axis1.16S.DS.T1,file=here("output/tables","PCA.axis1.16S.DS.T1.txt"),sep="\t")
write.table(PCA.axis2.16S.DS.T1,file=here("output/tables","PCA.axis2.16S.DS.T1.txt"),sep="\t")
write.table(PCA.axis3.16S.DS.T1,file=here("output/tables","PCA.axis3.16S.DS.T1.txt"),sep="\t")
write.table(PCA.axis4.16S.DS.T1,file=here("output/tables","PCA.axis4.16S.DS.T1.txt"),sep="\t")
write.table(PCA.axis5.16S.DS.T1,file=here("output/tables","PCA.axis5.16S.DS.T1.txt"),sep="\t")

bact.DS.T1.cor <- corr.test(qual.DS.S1,pca.data.point.16S.DS.T1, method = 'spearman')
myfun <- function(bact.DS.T1.cor){
  corrplot(bact.DS.T1.cor$r, 
           p.mat = round(as.matrix(bact.DS.T1.cor$p),3),
           method = 'circle',
           type = 'lower',
           sig.level = c(.001, .01, .05), 
           tl.pos="lt", 
           tl.col="black", tl.cex=1.0, 
           tl.offset=0.2,
           cl.pos="r",
           insig = "label_sig",
           pch.cex = 1.0,
           pch.col="red",
           cl.cex = 1.0)
  corrplot(bact.DS.T1.cor$r,  type="upper", method="number",
           col="coral4",  tl.pos="n", cl.pos="n", number.cex = 0.9, add=T,diag=F)
  recordPlot() # record the latest plot
}


myfun(bact.DS.T1.cor) 



#PCA analysis for 24th May
pca.16S.DS.T2 <- PCA(otu.16S.DS.T2.hel, graph = FALSE)
print(pca.16S.DS.T2)

##Cumulative percentage of eigenvalues
fviz_eig(pca.16S.DS.T2, addlabels = TRUE, ylim = c(0, 10))
var.16S.DS.T2 <- get_pca_var(pca.16S.DS.T2)

# Extraction of PCA points
pca.point.16S.DS.T2=data.frame(pca.16S.DS.T2$ind)
pca.data.point.16S.DS.T2=pca.point.16S.DS.T2[,which(colnames(pca.point.16S.DS.T2) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.16S.DS.T2)[which(colnames(pca.data.point.16S.DS.T2) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("bact.DS.T2.axis1", "bact.DS.T2.axis2","bact.DS.T2.axis3","bact.DS.T2.axis4","bact.DS.T2.axis5") 



#PCA analysis for 7th June
pca.16S.DS.T3 <- PCA(otu.16S.DS.T3.hel, graph = FALSE)
print(pca.16S.DS.T3)

##Cumulative percentage of eigenvalues
fviz_eig(pca.16S.DS.T3, addlabels = TRUE, ylim = c(0, 10))
var.16S.DS.T3 <- get_pca_var(pca.16S.DS.T3)

# Extraction of PCA points
pca.point.16S.DS.T3=data.frame(pca.16S.DS.T3$ind)
pca.data.point.16S.DS.T3=pca.point.16S.DS.T3[,which(colnames(pca.point.16S.DS.T3) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.16S.DS.T3)[which(colnames(pca.data.point.16S.DS.T3) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("bact.DS.T3.axis1", "bact.DS.T3.axis2","bact.DS.T3.axis3","bact.DS.T3.axis4","bact.DS.T3.axis5") 


res.desc.16S.DS.T3 <- dimdesc(pca.16S.DS.T3, axes = c(1,2,3,4,5), proba = 0.05)

PCA.axis1.16S.DS.T3<-data.frame(res.desc.16S.DS.T3$Dim.1)
PCA.axis2.16S.DS.T3<-data.frame(res.desc.16S.DS.T3$Dim.2)
PCA.axis3.16S.DS.T3<-data.frame(res.desc.16S.DS.T3$Dim.3)
PCA.axis4.16S.DS.T3<-data.frame(res.desc.16S.DS.T3$Dim.4)
PCA.axis5.16S.DS.T3<-data.frame(res.desc.16S.DS.T3$Dim.5)

write.table(PCA.axis1.16S.DS.T3,file=here("output/tables","PCA.axis1.16S.DS.T3.txt"),sep="\t")
write.table(PCA.axis2.16S.DS.T3,file=here("output/tables","PCA.axis2.16S.DS.T3.txt"),sep="\t")
write.table(PCA.axis3.16S.DS.T3,file=here("output/tables","PCA.axis3.16S.DS.T3.txt"),sep="\t")
write.table(PCA.axis4.16S.DS.T3,file=here("output/tables","PCA.axis4.16S.DS.T3.txt"),sep="\t")
write.table(PCA.axis5.16S.DS.T3,file=here("output/tables","PCA.axis5.16S.DS.T3.txt"),sep="\t")

bact.DS.T3.cor <- corr.test(qual.DS.S3,pca.data.point.16S.DS.T3, method = 'spearman')
myfun <- function(bact.DS.T3.cor){
  corrplot(bact.DS.T3.cor$r, 
           p.mat = round(as.matrix(bact.DS.T3.cor$p),3),
           method = 'circle',
           type = 'lower',
           sig.level = c(.001, .01, .05), 
           tl.pos="lt", 
           tl.col="black", tl.cex=1.0, 
           tl.offset=0.2,
           cl.pos="r",
           insig = "label_sig",
           pch.cex = 1.0,
           pch.col="red",
           cl.cex = 1.0)
  corrplot(bact.DS.T3.cor$r,  type="upper", method="number",
           col="coral4",  tl.pos="n", cl.pos="n", number.cex = 0.9, add=T,diag=F)
  recordPlot() # record the latest plot
}


myfun(bact.DS.T3.cor) 



#PCA analysis for 24th June
pca.16S.DS.T4 <- PCA(otu.16S.DS.T4.hel, graph = FALSE)
print(pca.16S.DS.T4)

##Cumulative percentage of eigenvalues
fviz_eig(pca.16S.DS.T4, addlabels = TRUE, ylim = c(0,8))
var.16S.DS.T4 <- get_pca_var(pca.16S.DS.T4)

# Extraction of PCA points
pca.point.16S.DS.T4=data.frame(pca.16S.DS.T4$ind)
pca.data.point.16S.DS.T4=pca.point.16S.DS.T4[,which(colnames(pca.point.16S.DS.T4) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.16S.DS.T4)[which(colnames(pca.data.point.16S.DS.T4) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("bact.DS.T4.axis1", "bact.DS.T4.axis2","bact.DS.T4.axis3","bact.DS.T4.axis4","bact.DS.T4.axis5") 








bact.DS.T4.cor <- corr.test(qual.DS.S4,pca.data.point.16S.DS.T4, method = 'spearman')
myfun <- function(bact.DS.T4.cor){
  corrplot(bact.DS.T4.cor$r, 
           p.mat = round(as.matrix(bact.DS.T4.cor$p),3),
           method = 'circle',
           type = 'lower',
           sig.level = c(.001, .01, .05), 
           tl.pos="lt", 
           tl.col="black", tl.cex=1.0, 
           tl.offset=0.2,
           cl.pos="r",
           insig = "label_sig",
           pch.cex = 1.0,
           pch.col="red",
           cl.cex = 1.0)
  corrplot(bact.DS.T4.cor$r,  type="upper", method="number",
           col="coral4",  tl.pos="n", cl.pos="n", number.cex = 0.9, add=T,diag=F)
  recordPlot() # record the latest plot
}


myfun(bact.DS.T4.cor) 


res.desc.16S.DS.T4 <- dimdesc(pca.16S.DS.T4, axes = c(1,2,3,4,5), proba = 0.05)

PCA.axis2.16S.DS.T4<-data.frame(res.desc.16S.DS.T4$Dim.2)
PCA.axis4.16S.DS.T4<-data.frame(res.desc.16S.DS.T4$Dim.4)

write.table(PCA.axis2.16S.DS.T4,file=here("output/tables","PCA.axis2.16S.DS.T4.txt"),sep="\t")
write.table(PCA.axis4.16S.DS.T4,file=here("output/tables","PCA.axis4.16S.DS.T4.txt"),sep="\t")


#PCA analysis for 7th July
pca.16S.DS.T5 <- PCA(otu.16S.DS.T5.hel, graph = FALSE)
print(pca.16S.DS.T5)

##Cumulative percentage of eigenvalues
fviz_eig(pca.16S.DS.T5, addlabels = TRUE, ylim = c(0, 10))
var.16S.DS.T5 <- get_pca_var(pca.16S.DS.T5)

# Extraction of PCA points
pca.point.16S.DS.T5=data.frame(pca.16S.DS.T5$ind)
pca.data.point.16S.DS.T5=pca.point.16S.DS.T5[,which(colnames(pca.point.16S.DS.T5) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.16S.DS.T5)[which(colnames(pca.data.point.16S.DS.T5) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("bact.DS.T5.axis1", "bact.DS.T5.axis2","bact.DS.T5.axis3","bact.DS.T5.axis4","bact.DS.T5.axis5") 

#PCA analysis for 19th July
pca.16S.DS.T6 <- PCA(otu.16S.DS.T6.hel, graph = FALSE)
print(pca.16S.DS.T5)

##Cumulative percentage of eigenvalues
fviz_eig(pca.16S.DS.T6, addlabels = TRUE, ylim = c(0, 10))
var.16S.DS.T6 <- get_pca_var(pca.16S.DS.T6)

# Extraction of PCA points
pca.point.16S.DS.T6=data.frame(pca.16S.DS.T6$ind)
pca.data.point.16S.DS.T6=pca.point.16S.DS.T6[,which(colnames(pca.point.16S.DS.T6) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.16S.DS.T6)[which(colnames(pca.data.point.16S.DS.T6) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("bact.DS.T6.axis1", "bact.DS.T6.axis2","bact.DS.T6.axis3","bact.DS.T6.axis4","bact.DS.T6.axis5") 

#PCA analysis for 1Aug
pca.16S.DS.T7 <- PCA(otu.16S.DS.T7.hel, graph = FALSE)
print(pca.16S.DT.T7)

##Cumulative percentage of eigenvalues
fviz_eig(pca.16S.DS.T7, addlabels = TRUE, ylim = c(0, 10))
var.16S.DS.T7 <- get_pca_var(pca.16S.DS.T7)

# Extraction of PCA points
pca.point.16S.DS.T7=data.frame(pca.16S.DS.T7$ind)
pca.data.point.16S.DS.T7=pca.point.16S.DS.T7[,which(colnames(pca.point.16S.DS.T7) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.16S.DS.T7)[which(colnames(pca.data.point.16S.DS.T7) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("bact.DS.T7.axis1", "bact.DS.T7.axis2","bact.DS.T7.axis3","bact.DS.T7.axis4","bact.DS.T7.axis5") 

# Fungal OTUs for PCA analysis

#PCA run for 10 MAY
pca.ITS.DS.T1 <- PCA(otu.ITS.DS.T1.hel, graph = FALSE)
print(pca.ITS.DS.T1)

#Visualization and interpretation of eigenvalues
eig.val.ITS.DS.T1 <- get_eigenvalue(pca.ITS.DS.T1)
eig.val.ITS.DS.T1

#Cumulative percentage of eigenvalues
fviz_eig(pca.ITS.DS.T1, addlabels = TRUE, ylim = c(0, 8))
var.ITS.DS.T1<- get_pca_var(pca.ITS.DS.T1)


pca.point.DS.ITS.T1<-data.frame(pca.ITS.DS.T1$ind)

pca.data.point.ITS.DS.T1<-pca.point.DS.ITS.T1[,which(colnames(pca.point.DS.ITS.T1) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.ITS.DS.T1)[which(colnames(pca.data.point.ITS.DS.T1) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("fun.DS.T1.axis1", "fun.DS.T1.axis2","fun.DS.T1.axis3","fun.DS.T1.axis4","fun.DS.T1.axis5") 


res.desc.ITS.DS.T1 <- dimdesc(pca.ITS.DS.T1, axes = c(1,2,3,4,5), proba = 0.05)
PCA.axis1.ITS.DS.T1<-data.frame(res.desc.ITS.DS.T1$Dim.1)
PCA.axis2.ITS.DS.T1<-data.frame(res.desc.ITS.DS.T1$Dim.2)
PCA.axis3.ITS.DS.T1<-data.frame(res.desc.ITS.DS.T1$Dim.3)
PCA.axis4.ITS.DS.T1<-data.frame(res.desc.ITS.DS.T1$Dim.4)
PCA.axis5.ITS.DS.T1<-data.frame(res.desc.ITS.DS.T1$Dim.5)

write.table(PCA.axis1.ITS.DS.T1,file = here("output/tables","PCA.axis1.ITS.DS.T1.txt"),sep="\t")
write.table(PCA.axis2.ITS.DS.T1,file = here("output/tables","PCA.axis2.ITS.DS.T1.txt"),sep="\t")
write.table(PCA.axis3.ITS.DS.T1,file = here("output/tables","PCA.axis3.ITS.DS.T1.txt"),sep="\t")
write.table(PCA.axis4.ITS.DS.T1,file = here("output/tables","PCA.axis4.ITS.DS.T1.txt"),sep="\t")
write.table(PCA.axis5.ITS.DS.T1,file = here("output/tables","PCA.axis5.ITS.DS.T1.txt"),sep="\t")

fun.DS.T1.cor <- corr.test(qual.DS.S1,pca.data.point.ITS.DS.T1, method = 'spearman')
myfun <- function(fun.DS.T1.cor){
  corrplot(fun.DS.T1.cor$r, 
           p.mat = round(as.matrix(fun.DS.T1.cor$p),3),
           method = 'circle',
           type = 'lower',
           sig.level = c(.001, .01, .05), 
           tl.pos="lt", 
           tl.col="black", tl.cex=1.0, 
           tl.offset=0.2,
           cl.pos="r",
           insig = "label_sig",
           pch.cex = 1.0,
           pch.col="red",
           cl.cex = 1.0)
  corrplot(fun.DS.T1.cor$r,  type="upper", method="number",
           col="coral4",  tl.pos="n", cl.pos="n", number.cex = 0.9, add=T,diag=F)
  recordPlot() # record the latest plot
}


myfun(fun.DS.T1.cor) 


#PCA run for 24 MAY
pca.ITS.DS.T2 <- PCA(otu.ITS.DS.T2.hel, graph = FALSE)
print(pca.ITS.DS.T2)

#Visualization and interpretation of eigenvalues
eig.val.ITS.DS.T2 <- get_eigenvalue(pca.ITS.DS.T2)
eig.val.ITS.DS.T2

#Cumulative percentage of eigenvalues
fviz_eig(pca.ITS.DS.T2, addlabels = TRUE, ylim = c(0, 10))
var.ITS.DS.T2<- get_pca_var(pca.ITS.DS.T2)


pca.point.DS.ITS.T2<-data.frame(pca.ITS.DS.T2$ind)

pca.data.point.ITS.DS.T2<-pca.point.DS.ITS.T2[,which(colnames(pca.point.DS.ITS.T2) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.ITS.DS.T2)[which(colnames(pca.data.point.ITS.DS.T2) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("fun.DS.T2.axis1", "fun.DS.T2.axis2","fun.DS.T1.axis3","fun.DS.T2.axis4","fun.DS.T2.axis5") 



#PCA run for 07 th June
pca.ITS.DS.T3 <- PCA(otu.ITS.DS.T3.hel, graph = FALSE)
print(pca.ITS.DS.T3)

#Visualization and interpretation of eigenvalues
eig.val.ITS.DS.T3 <- get_eigenvalue(pca.ITS.DS.T3)
eig.val.ITS.DS.T3

#Cumulative percentage of eigenvalues
fviz_eig(pca.ITS.DS.T3, addlabels = TRUE, ylim = c(0, 8))
var.ITS.DS.T3<- get_pca_var(pca.ITS.DS.T3)


pca.point.DS.ITS.T3<-data.frame(pca.ITS.DS.T3$ind)

pca.data.point.ITS.DS.T3<-pca.point.DS.ITS.T3[,which(colnames(pca.point.DS.ITS.T3) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.ITS.DS.T3)[which(colnames(pca.data.point.ITS.DS.T3) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("fun.DS.T3.axis1", "fun.DS.T3.axis2","fun.DS.T3.axis3","fun.DS.T3.axis4","fun.DS.T3.axis5") 



res.desc.ITS.DS.T3 <- dimdesc(pca.ITS.DS.T3, axes = c(1,2,3,4,5), proba = 0.05)

PCA.axis1.ITS.DS.T3<-data.frame(res.desc.ITS.DS.T3$Dim.1)
PCA.axis2.ITS.DS.T3<-data.frame(res.desc.ITS.DS.T3$Dim.2)
PCA.axis3.ITS.DS.T3<-data.frame(res.desc.ITS.DS.T3$Dim.3)
PCA.axis4.ITS.DS.T3<-data.frame(res.desc.ITS.DS.T3$Dim.4)
PCA.axis5.ITS.DS.T3<-data.frame(res.desc.ITS.DS.T3$Dim.5)

write.table(PCA.axis1.ITS.DS.T3,file = here("output/tables","PCA.axis1.ITS.DS.T3.txt"),sep="\t")
write.table(PCA.axis2.ITS.DS.T3,file = here("output/tables","PCA.axis2.ITS.DS.T3.txt"),sep="\t")
write.table(PCA.axis3.ITS.DS.T3,file = here("output/tables","PCA.axis3.ITS.DS.T3.txt"),sep="\t")
write.table(PCA.axis4.ITS.DS.T3,file = here("output/tables","PCA.axis4.ITS.DS.T3.txt"),sep="\t")
write.table(PCA.axis5.ITS.DS.T3,file = here("output/tables","PCA.axis5.ITS.DS.T3.txt"),sep="\t")

fun.DS.T3.cor <- corr.test(qual.DS.S3,pca.data.point.ITS.DS.T3, method = 'spearman')
myfun <- function(fun.DS.T3.cor){
  corrplot(fun.DS.T3.cor$r, 
           p.mat = round(as.matrix(fun.DS.T3.cor$p),3),
           method = 'circle',
           type = 'lower',
           sig.level = c(.001, .01, .05), 
           tl.pos="lt", 
           tl.col="black", tl.cex=1.0, 
           tl.offset=0.2,
           cl.pos="r",
           insig = "label_sig",
           pch.cex = 1.0,
           pch.col="red",
           cl.cex = 1.0)
  corrplot(fun.DS.T3.cor$r,  type="upper", method="number",
           col="coral4",  tl.pos="n", cl.pos="n", number.cex = 0.9, add=T,diag=F)
  recordPlot() # record the latest plot
}


myfun(fun.DS.T3.cor) 








# Extraction of PCA points from DT-T3

res.desc.ITS.DT.T3 <- dimdesc(pca.ITS.DS.T3, axes = c(1,2,4,5), proba = 0.05)

PCA.axis1.ITS.DT.T3<-data.frame(res.desc.ITS.DT.T3$Dim.1)
PCA.axis2.ITS.DT.T3<-data.frame(res.desc.ITS.DT.T3$Dim.2)
PCA.axis4.ITS.DT.T3<-data.frame(res.desc.ITS.DT.T3$Dim.4)
PCA.axis5.ITS.DT.T3<-data.frame(res.desc.ITS.DT.T3$Dim.5)

write.table(PCA.axis1.ITS.DT.T3,file = here("output/tables","PCA.axis1.ITS.DT.T3.txt"),sep="\t")
write.table(PCA.axis2.ITS.DT.T3,file = here("output/tables","PCA.axis2.ITS.DT.T3.txt"),sep="\t")
write.table(PCA.axis4.ITS.DT.T3,file = here("output/tables","PCA.axis4.ITS.DT.T3.txt"),sep="\t")
write.table(PCA.axis5.ITS.DT.T3,file = here("output/tables","PCA.axis5.ITS.DT.T3.txt"),sep="\t")



#PCA run for 27th June
pca.ITS.DS.T4 <- PCA(otu.ITS.DS.T4.hel, graph = FALSE)
print(pca.ITS.DS.T4)

#Visualization and interpretation of eigenvalues
eig.val.ITS.DS.T4 <- get_eigenvalue(pca.ITS.DS.T4)
eig.val.ITS.DS.T4

#Cumulative percentage of eigenvalues
fviz_eig(pca.ITS.DS.T4, addlabels = TRUE, ylim = c(0, 10))
var.ITS.DS.T4<- get_pca_var(pca.ITS.DS.T4)


pca.point.DS.ITS.T4<-data.frame(pca.ITS.DS.T4$ind)

pca.data.point.ITS.DS.T4<-pca.point.DS.ITS.T4[,which(colnames(pca.point.DS.ITS.T4) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.ITS.DS.T4)[which(colnames(pca.data.point.ITS.DS.T4) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("fun.DS.T4.axis1", "fun.DS.T4.axis2","fun.DS.T4.axis3","fun.DS.T4.axis4","fun.DS.T4.axis5") 



res.desc.ITS.DS.T4 <- dimdesc(pca.ITS.DS.T4, axes = c(1,2,3,4,5), proba = 0.05)

PCA.axis2.ITS.DS.T4<-data.frame(res.desc.ITS.DS.T4$Dim.2)

write.table(PCA.axis2.ITS.DS.T4,file = here("output/tables","PCA.axis2.ITS.DS.T4.txt"),sep="\t")




#PCA run for 27th June
pca.ITS.DS.T5 <- PCA(otu.ITS.DS.T5.hel, graph = FALSE)
print(pca.ITS.DS.T5)

#Visualization and interpretation of eigenvalues
eig.val.ITS.DS.T5 <- get_eigenvalue(pca.ITS.DS.T5)
eig.val.ITS.DS.T5

#Cumulative percentage of eigenvalues
fviz_eig(pca.ITS.DS.T5, addlabels = TRUE, ylim = c(0, 10))
var.ITS.DS.T5<- get_pca_var(pca.ITS.DS.T5)


pca.point.DS.ITS.T5<-data.frame(pca.ITS.DS.T5$ind)

pca.data.point.ITS.DS.T5<-pca.point.DS.ITS.T5[,which(colnames(pca.point.DS.ITS.T5) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.ITS.DS.T5)[which(colnames(pca.data.point.ITS.DS.T5) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("fun.DS.T5.axis1", "fun.DS.T5.axis2","fun.DS.T5.axis3","fun.DS.T4.axis4","fun.DS.T5.axis5") 


#PCA run for 19th JUly
pca.ITS.DS.T6 <- PCA(otu.ITS.DS.T6.hel, graph = FALSE)
print(pca.ITS.DS.T6)

#Visualization and interpretation of eigenvalues
eig.val.ITS.DS.T6 <- get_eigenvalue(pca.ITS.DS.T6)
eig.val.ITS.DS.T6

#Cumulative percentage of eigenvalues
fviz_eig(pca.ITS.DS.T6, addlabels = TRUE, ylim = c(0, 10))
var.ITS.DS.T6<- get_pca_var(pca.ITS.DS.T6)


pca.point.DS.ITS.T6<-data.frame(pca.ITS.DS.T6$ind)

pca.data.point.ITS.DS.T6<-pca.point.DS.ITS.T6[,which(colnames(pca.point.DS.ITS.T6) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.ITS.DS.T6)[which(colnames(pca.data.point.ITS.DS.T6) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("fun.DS.T6.axis1", "fun.DS.T6.axis2","fun.DS.T6.axis3","fun.DT.T6.axis4","fun.DS.T6.axis5") 

#PCA run for 1st August
pca.ITS.DS.T7 <- PCA(otu.ITS.DS.T7.hel, graph = FALSE)
print(pca.ITS.DS.T7)

#Visualization and interpretation of eigenvalues
eig.val.ITS.DS.T7 <- get_eigenvalue(pca.ITS.DS.T7)
eig.val.ITS.DS.T7

#Cumulative percentage of eigenvalues
fviz_eig(pca.ITS.DS.T7, addlabels = TRUE, ylim = c(0, 10))
var.ITS.DS.T7<- get_pca_var(pca.ITS.DS.T7)


pca.point.DS.ITS.T7<-data.frame(pca.ITS.DS.T7$ind)

pca.data.point.ITS.DS.T7<-pca.point.DS.ITS.T7[,which(colnames(pca.point.DS.ITS.T7) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.ITS.DS.T7)[which(colnames(pca.data.point.ITS.DS.T7) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("fun.DS.T7.axis1", "fun.DS.T7.axis2","fun.DS.T7.axis3","fun.DS.T7.axis4","fun.DS.T7.axis5") 


#Biolog PCA

# 5 May data (BIOLOG)


#PCA run for T1
pca.biolog.DS.T1 <- PCA(biolog.DS.T1.hel, graph = FALSE)
print(pca.biolog.DS.T1)

#Visualization and interpretation of eigenvalues
eig.val.biolog.DS.T1 <- get_eigenvalue(pca.biolog.DS.T1)
eig.val.biolog.DS.T1

#Cumulative percentage of eigenvalues
fviz_eig(pca.biolog.DS.T1, addlabels = TRUE, ylim = c(0, 35))
var.biolog.DS.T1<- get_pca_var(pca.biolog.DS.T1)

pca.point.biolog.DS.T1<-data.frame(pca.biolog.DS.T1$ind)

pca.data.point.biolog.DS.T1<-pca.point.biolog.DS.T1[,which(colnames(pca.point.biolog.DS.T1) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.biolog.DS.T1)[which(colnames(pca.data.point.biolog.DS.T1) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("biolog.DS.T1.axis1", "biolog.DS.T1.axis2","biolog.DS.T1.axis3","biolog.DS.T1.axis4","biolog.DS.T1.axis5") 

#Analysis of significant dimensions

res.desc.biolog.DS.T1 <- dimdesc(pca.biolog.DS.T1, axes = c(3,4,5), proba = 0.05)
# Description of dimension 1
res.desc.biolog.DS.T1$Dim.1

# Description of dimension 2
res.desc.biolog.DS.T1$Dim.2

#Clustering
set.seed(123)
res.km.bio.DS.T1 <- kmeans(var.biolog.DS.T1$coord, centers = 3, nstart = 25)
grp.bio.DS.T1 <- as.factor(res.km.bio.DS.T1$cluster)


# Color variables by groups
fviz_pca_var(pca.biolog.DS.T1, col.var = grp.bio.DS.T1, 
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),
             legend.title = "Cluster")



# Contributions of variables to PC1
fviz_contrib(pca.biolog.DS.T1, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(pca.biolog.DS.T1, choice = "var", axes = 3, top = 10)

## Contributions of variables to PC2
fviz_contrib(pca.biolog.DS.T1, choice = "var", axes = 4, top = 10)


#Total contribution to PC1 and PC2
fviz_contrib(pca.biolog.DS.T1, choice = "var", axes = 1:2, top = 10)

# Create a random continuous variable of length 10
set.seed(123)
my.cont.var <- rnorm(10)


# Color variables by the continuous variable
# Color individuals by the continuous variable
# Color variables by groups


fviz_pca_ind(pca.biolog.DS.T1, col.ind = "cos2", pointsize = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)

# Color individuals by the continuous variable

fviz_pca_ind(pca.biolog.DS.T1,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = factor(metafile.DS.T1$Treatment), # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07","#D55E00"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)



#PCA biplot

fviz_pca_biplot(pca.biolog.DS.T1, 
                # Fill individuals by groups
                geom.ind = "point",
                pointshape = 21,
                pointsize = 2.5,
                fill.ind = factor(metafile.DS.T1$Treatment),
                col.ind = "black",
                # Color variable by groups
                #col.var = factor(c("100", "75", "25", "50")),
                
                legend.title = list(fill = "Treatment", color = "Cultivar"),
                repel = TRUE        # Avoid label over plotting
)+
  ggpubr::fill_palette("jco")+      # Individual fill color
  ggpubr::color_palette("npg")  


fviz_pca_biplot(pca.biolog.DS.T1, 
                # Individuals
                geom.ind = "point",
                fill.ind = factor(metafile.DS.T1$Treatment), col.ind = "black",
                pointshape = 21, pointsize = 2,
                palette = "jco",
                addEllipses = TRUE,
                # Variables
                alpha.var ="contrib", col.var = "contrib",
                gradient.cols = "RdYlBu",
                
                legend.title = list(fill = "Treatment", color = "Contrib",
                                    alpha = "Contrib")
)



res.desc.biolog.DS.T1 <- dimdesc(pca.biolog.DS.T1, axes = c(3,4,5), proba = 0.05)


PCA.axis3.biolog.DS.T1<-data.frame(res.desc.biolog.DS.T1$Dim.3)
PCA.axis4.biolog.DS.T1<-data.frame(res.desc.biolog.DS.T1$Dim.4)
PCA.axis5.biolog.DS.T1<-data.frame(res.desc.biolog.DS.T1$Dim.5)

write.table(PCA.axis3.biolog.DS.T1,file = here("output/tables","PCA.axis3.biolog.DS.T1.txt"),sep="\t")
write.table(PCA.axis4.biolog.DS.T1,file = here("output/tables","PCA.axis4.biolog.DS.T1.txt"),sep="\t")
write.table(PCA.axis5.biolog.DS.T1,file = here("output/tables","PCA.axis5.biolog.DS.T1.txt"),sep="\t")


biolog.DS.T1.cor <- corr.test(qual.DS.S1,pca.data.point.biolog.DS.T1, method = 'spearman')
myfun <- function(biolog.DS.T1.cor){
  corrplot(biolog.DS.T1.cor$r, 
           p.mat = round(as.matrix(biolog.DS.T1.cor$p),3),
           method = 'circle',
           type = 'lower',
           sig.level = c(.001, .01, .05), 
           tl.pos="lt", 
           tl.col="black", tl.cex=1.0, 
           tl.offset=0.2,
           cl.pos="r",
           insig = "label_sig",
           pch.cex = 1.0,
           pch.col="red",
           cl.cex = 1.0)
  corrplot(biolog.DS.T1.cor$r,  type="upper", method="number",
           col="coral4",  tl.pos="n", cl.pos="n", number.cex = 0.9, add=T,diag=F)
  recordPlot() # record the latest plot
}


myfun(biolog.DS.T1.cor) 





#PCA run FOR T2
pca.biolog.DS.T2 <- PCA(biolog.DS.T2.hel, graph = FALSE)
print(pca.biolog.DS.T2)

#Visualization and interpretation of eigenvalues
eig.val.biolog.DS.T2 <- get_eigenvalue(pca.biolog.DS.T2)
eig.val.biolog.DS.T2

#Cumulative percentage of eigenvalues
fviz_eig(pca.biolog.DS.T2, addlabels = TRUE, ylim = c(0, 25))
var.biolog.DS.T2<- get_pca_var(pca.biolog.DS.T2)

pca.point.biolog.DS.T2<-data.frame(pca.biolog.DS.T2$ind)

pca.data.point.biolog.DS.T2<-pca.point.biolog.DS.T2[,which(colnames(pca.point.biolog.DS.T2) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.biolog.DS.T2)[which(colnames(pca.data.point.biolog.DS.T2) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("biolog.DS.T2.axis1", "biolog.DS.T2.axis2","biolog.DS.T2.axis3","biolog.DS.T2.axis4","biolog.DS.T2.axis5") 


#PCA run FOR T3
pca.biolog.DS.T3 <- PCA(biolog.DS.T3.hel, graph = FALSE)
print(pca.biolog.DS.T3)

#Visualization and interpretation of eigenvalues
eig.val.biolog.DS.T3 <- get_eigenvalue(pca.biolog.DS.T3)
eig.val.biolog.DS.T3

#Cumulative percentage of eigenvalues
fviz_eig(pca.biolog.DS.T3, addlabels = TRUE, ylim = c(0, 25))
var.biolog.DS.T3<- get_pca_var(pca.biolog.DS.T3)

pca.point.biolog.DS.T3<-data.frame(pca.biolog.DS.T3$ind)

pca.data.point.biolog.DS.T3<-pca.point.biolog.DS.T3[,which(colnames(pca.point.biolog.DS.T3) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]
colnames(pca.data.point.biolog.DS.T3)[which(colnames(pca.data.point.biolog.DS.T3) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("biolog.DS.T3.axis1", "biolog.DS.T3.axis2","biolog.DS.T3.axis3","biolog.DS.T3.axis4","biolog.DS.T3.axis5") 


pca.data.point.biolog.DS.T3.cor<-pca.point.biolog.DS.T3[,which(colnames(pca.point.biolog.DS.T3) %in% c("coord.Dim.1","coord.Dim.3", "coord.Dim.5"))]
colnames(pca.data.point.biolog.DS.T3.cor)[which(colnames(pca.data.point.biolog.DS.T3.cor) %in% c("coord.Dim.1","coord.Dim.3", "coord.Dim.5"))]<- c("biolog.DS.T3.axis1", "biolog.DS.T3.axis3","biolog.DS.T3.axis5") 





corr.bio.DS.T3<- corr.test(qual.DS.S3,pca.data.point.biolog.DS.T3, method='spearman')



myfun <- function(corr.bio.DS.T3){
  corrplot(corr.bio.DS.T3$r, 
           p.mat = round(as.matrix(corr.bio.DS.T3$p),3),
           method = 'circle',
           type = 'lower',
           sig.level = c(.001, .01, .05), 
           tl.pos="lt", 
           tl.col="black", tl.cex=1, 
           tl.offset=0.09,
           cl.pos="r",
           insig = "label_sig",
           pch.cex = 1,
           pch.col="red",
           cl.cex = 1)
  corrplot(corr.bio.DS.T3$r,  type="upper", method="number",
           col="coral4",  tl.pos="n", cl.pos="n", number.cex = 1, add=T,diag=F)
  recordPlot() # record the latest plot
}

myfun (corr.bio.DS.T3)

# Extraction of PCA points

res.desc.biolog.DS.T3 <- dimdesc(pca.biolog.DS.T3, axes = c(1,3,5), proba = 0.05)


PCA.axis1.biolog.DS.T3<-data.frame(res.desc.biolog.DS.T3$Dim.1)
PCA.axis3.biolog.DS.T3<-data.frame(res.desc.biolog.DS.T3$Dim.3)
PCA.axis5.biolog.DS.T3<-data.frame(res.desc.biolog.DS.T3$Dim.5)

write.table(PCA.axis1.biolog.DS.T3,file = here("output/tables","PCA.axis1.biolog.DS.T3.txt"),sep="\t")
write.table(PCA.axis3.biolog.DS.T3,file = here("output/tables","PCA.axis3.biolog.DS.T3.txt"),sep="\t")
write.table(PCA.axis5.biolog.DS.T3,file = here("output/tables","PCA.axis5.biolog.DS.T3.txt"),sep="\t")





#Clustering
set.seed(123)
res.km.bio.DS.T3 <- kmeans(var.biolog.DS.T3$coord, centers = 5, nstart = 25)
grp.bio.DS.T3 <- as.factor(res.km.bio.DS.T3$cluster)


# Color variables by groups
fviz_pca_var(pca.biolog.DS.T3, col.var = grp.bio.DS.T3, 
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF","darkgreen", "deeppink"),
             legend.title = "Dimensions")

# Contributions of variables to PC1
fviz_contrib(pca.biolog.DS.T3, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(pca.biolog.DS.T3, choice = "var", axes = 3, top = 10)

## Contributions of variables to PC2
fviz_contrib(pca.biolog.DS.T3, choice = "var", axes = 4, top = 10)

#Total contribution to PC1 and PC2
fviz_contrib(pca.biolog.DS.T3, choice = "var", axes = 1:2, top = 10)

# Create a random continuous variable of length 10
set.seed(123)
my.cont.var <- rnorm(10)


# Color variables by the continuous variable
# Color individuals by the continuous variable
# Color variables by groups


fviz_pca_ind(pca.biolog.DS.T3, col.ind = "cos2", pointsize = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)

# Color individuals by the continuous variable

fviz_pca_ind(pca.biolog.DS.T3,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = factor(metafile.DS.T3$Treatment), # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07","#D55E00"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)



#PCA biplot

fviz_pca_biplot(pca.biolog.DS.T3, 
                # Fill individuals by groups
                geom.ind = "point",
                pointshape = 21,
                pointsize = 2.5,
                fill.ind = factor(metafile.DS.T3$Treatment),
                col.ind = "black",
                # Color variable by groups
                #col.var = factor(c("100", "75", "25", "50")),
                
                legend.title = list(fill = "Treatment", color = "Cultivar"),
                repel = TRUE        # Avoid label over plotting
)+
  ggpubr::fill_palette("jco")+      # Individual fill color
  ggpubr::color_palette("npg")  


fviz_pca_biplot(pca.biolog.DS.T3, 
                # Individuals
                geom.ind = "point",
                fill.ind = factor(metafile.DS.T3$Treatment), col.ind = "black",
                pointshape = 21, pointsize = 2,
                palette = "jco",
                addEllipses = TRUE,
                # Variables
                alpha.var ="contrib", col.var = "contrib",
                gradient.cols = "RdYlBu",
                
                legend.title = list(fill = "Treatment", color = "Contrib",
                                    alpha = "Contrib")
)




#PCA run FOR T4

pca.biolog.DS.T4 <- PCA(biolog.DS.T4.hel, graph = FALSE)
print(pca.biolog.DS.T4)

#Visualization and interpretation of eigenvalues
eig.val.biolog.DS.T4 <- get_eigenvalue(pca.biolog.DS.T4)
eig.val.biolog.DS.T4

#Cumulative percentage of eigenvalues
fviz_eig(pca.biolog.DS.T4, addlabels = TRUE, ylim = c(0, 35))
var.biolog.DS.T4<- get_pca_var(pca.biolog.DS.T4)

pca.point.biolog.DS.T4<-data.frame(pca.biolog.DS.T4$ind)

pca.data.point.biolog.DS.T4<-pca.point.biolog.DS.T4[,which(colnames(pca.point.biolog.DS.T4) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.biolog.DS.T4)[which(colnames(pca.data.point.biolog.DS.T4) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("biolog.DS.T4.axis1", "biolog.DS.T4.axis2","biolog.DS.T4.axis3","biolog.DS.T4.axis4","biolog.DS.T4.axis5") 


corr.bio.DS.T4<- corr.test(qual.DS.S4,pca.data.point.biolog.DS.T4, method='spearman')




myfun <- function(corr.bio.DS.T4){
  corrplot(corr.bio.DS.T4$r, 
           p.mat = round(as.matrix(corr.bio.DS.T4$p),3),
           method = 'circle',
           type = 'lower',
           sig.level = c(.001, .01, .05), 
           tl.pos="lt", 
           tl.col="black", tl.cex=1, 
           tl.offset=0.09,
           cl.pos="r",
           insig = "label_sig",
           pch.cex = 1,
           pch.col="red",
           cl.cex = 1)
  corrplot(corr.bio.DS.T4$r,  type="upper", method="number",
           col="coral4",  tl.pos="n", cl.pos="n", number.cex = 1, add=T,diag=F)
  recordPlot() # record the latest plot
}

myfun (corr.bio.DS.T4)


# Extraction of PCA points

res.desc.biolog.DS.T4 <- dimdesc(pca.biolog.DS.T3, axes = c(1,2,3,4), proba = 0.05)

PCA.axis2.biolog.DS.T4<-data.frame(res.desc.biolog.DS.T4$Dim.2)
PCA.axis4.biolog.DS.T4<-data.frame(res.desc.biolog.DS.T4$Dim.4)

write.table(PCA.axis4.biolog.DS.T4,file = here("output/tables","PCA.axis4.biolog.DS.T4.txt"),sep="\t")
write.table(PCA.axis2.biolog.DS.T4,file = here("output/tables","PCA.axis2.biolog.DS.T4.txt"),sep="\t")



#PCA run FOR T5

pca.biolog.DS.T5 <- PCA(biolog.DS.T5.hel, graph = FALSE)
print(pca.biolog.DS.T5)

#Visualization and interpretation of eigenvalues
eig.val.biolog.DS.T5 <- get_eigenvalue(pca.biolog.DS.T5)
eig.val.biolog.DS.T5

#Cumulative percentage of eigenvalues
fviz_eig(pca.biolog.DS.T5, addlabels = TRUE, ylim = c(0, 25))
var.biolog.DS.T5<- get_pca_var(pca.biolog.DS.T5)

pca.point.biolog.DS.T5<-data.frame(pca.biolog.DS.T5$ind)

pca.data.point.biolog.DS.T5<-pca.point.biolog.DS.T5[,which(colnames(pca.point.biolog.DS.T5) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.biolog.DS.T5)[which(colnames(pca.data.point.biolog.DS.T5) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("biolog.DS.T5.axis1", "biolog.DS.T5.axis2","biolog.DS.T5.axis3","biolog.DS.T5.axis4","biolog.DS.T5.axis5") 


#PCA run FOR T6
pca.biolog.DS.T6 <- PCA(biolog.DS.T6.hel, graph = FALSE)
print(pca.biolog.DS.T6)

#Visualization and interpretation of eigenvalues
eig.val.biolog.DS.T6 <- get_eigenvalue(pca.biolog.DS.T6)
eig.val.biolog.DS.T6

#Cumulative percentage of eigenvalues
fviz_eig(pca.biolog.DS.T6, addlabels = TRUE, ylim = c(0, 25))
var.biolog.DS.T6<- get_pca_var(pca.biolog.DS.T6)

pca.point.biolog.DS.T6<-data.frame(pca.biolog.DS.T6$ind)

pca.data.point.biolog.DS.T6<-pca.point.biolog.DS.T6[,which(colnames(pca.point.biolog.DS.T6) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.biolog.DS.T6)[which(colnames(pca.data.point.biolog.DS.T6) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("biolog.DS.T6.axis1", "biolog.DS.T6.axis2","biolog.DS.T6.axis3","biolog.DS.T6.axis4","biolog.DS.T6.axis5") 

#PCA run for T7
pca.biolog.DS.T7 <- PCA(biolog.DS.T7.hel, graph = FALSE)
print(pca.biolog.DS.T7)

#Visualization and interpretation of eigenvalues
eig.val.biolog.DS.T7 <- get_eigenvalue(pca.biolog.DS.T7)
eig.val.biolog.DS.T7

#Cumulative percentage of eigenvalues
fviz_eig(pca.biolog.DS.T7, addlabels = TRUE, ylim = c(0, 25))
var.biolog.DS.T7<- get_pca_var(pca.biolog.DS.T7)

pca.point.biolog.DS.T7<-data.frame(pca.biolog.DS.T7$ind)

pca.data.point.biolog.DS.T7<-pca.point.biolog.DS.T7[,which(colnames(pca.point.biolog.DS.T7) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]

colnames(pca.data.point.biolog.DS.T7)[which(colnames(pca.data.point.biolog.DS.T7) %in% c("coord.Dim.1","coord.Dim.2","coord.Dim.3","coord.Dim.4", "coord.Dim.5"))]<- c("biolog.DS.T7.axis1", "biolog.DS.T7.axis2","biolog.DS.T7.axis3","biolog.DS.T7.axis4","biolog.DS.T7.axis5") 



# LASSO regression analysis

# Preparing data for LASSO regression for T1

reg.lasso.DS.T1<-cbind(pca.data.point.16S.DS.T1,pca.data.point.ITS.DS.T1,pca.data.point.biolog.DS.T1,div.16S.DS.T1,div.ITS.DS.T1,qpcr.T1.DS,qual.DS.S1)
row.names(reg.lasso.DS.T1)==row.names(qual.DS.S1)
#reg.lasso.DS.T1<-decostand(reg.lasso.DS.T1,method="standardize")

#qpcr.T1.DT<-log(qpcr.T1.DT)

## Protein model for May 10
data.glut.DS.T1=reg.lasso.DS.T1[!is.na(reg.lasso.DS.T1$Gluten),]
data.prot.DS.T1=reg.lasso.DS.T1[!is.na(reg.lasso.DS.T1$Protein.grain),]
data.pmt.DS.T1=reg.lasso.DS.T1[!is.na(reg.lasso.DS.T1$PMT),]
data.bem.DS.T1=reg.lasso.DS.T1[!is.na(reg.lasso.DS.T1$BEM),]
data.bem.DS.T1=data.bem.DS.T1[-c(22),]



data.explain.glut.DS.T1=data.glut.DS.T1[,-which(names(data.glut.DS.T1) %in% c("Protein.grain", "Proteine.flour","Humidity","Ash","Gluten", "IDC","PMT","BEM"))]
data.explain.prot.DS.T1=data.prot.DS.T1[,-which(names(data.prot.DS.T1) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.pmt.DS.T1=data.pmt.DS.T1[,-which(names(data.pmt.DS.T1) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.bem.DS.T1=data.bem.DS.T1[,-which(names(data.bem.DS.T1) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]

# X-Y data prep
X.glut.DS.T1=as.matrix(scale(data.explain.glut.DS.T1,center = T,scale = T))
#X.glut.DS.T1=as.matrix(data.explain.glut.DS.T1,center = T,scale = T)
#Y.reg.DS.T1=as.matrix(data.glut.DS.T1$Gluten)

Y.glut.DS.T1=as.matrix(scale(data.glut.DS.T1$Gluten,center = T,scale = T))

X.prot.DS.T1=as.matrix(scale(data.explain.prot.DS.T1,center = T,scale = T))
Y.prot.DS.T1=as.matrix(scale(data.prot.DS.T1$Protein.grain,center = T,scale = T))

X.pmt.DS.T1=as.matrix(scale(data.explain.pmt.DS.T1,center = T,scale = T))
Y.pmt.DS.T1=as.matrix(scale(data.pmt.DS.T1$PMT, center = T,scale = T))

X.bem.DS.T1=as.matrix(scale(data.explain.bem.DS.T1,center = T,scale = T))
Y.bem.DS.T1=as.matrix(scale(data.bem.DS.T1$BEM,center = T,scale = T))

#LASSO regression
fit.lasso.glut.DS.T1 <-glmnet(X.glut.DS.T1,Y.glut.DS.T1, family="gaussian", alpha=1)
fit.lasso.prot.DS.T1<- glmnet(X.prot.DS.T1,Y.prot.DS.T1, family="gaussian", alpha=1)
fit.lasso.pmt.DS.T1 <- glmnet(X.pmt.DS.T1,Y.pmt.DS.T1, family="gaussian", alpha=1)
fit.lasso.bem.DS.T1 <- glmnet(X.bem.DS.T1,Y.bem.DS.T1, family="gaussian", alpha=1)


#Gluten

#perform k-fold cross-validation to find optimal lambda value
n <- 100; s <- 0
for(i in 1:n) 
  
  s <- s + cv.glmnet(X.glut.DS.T1, Y.glut.DS.T1, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.glut.DS.T1, s=s.best.lasso))[-1]
s.best.lasso

lasso.glut.DS.T1<-coef(fit.lasso.glut, s=s.best.lasso)
lasso.glut.DS.T1


#perform k-fold cross-validation to find optimal lambda value
cv_model.glut.DS.T1 <- cv.glmnet(X.glut.DS.T1, Y.glut.DS.T1, k=10, nfolds=nrow(X.glut.DS.T1),alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.glut.DS.T1 <- cv_model.glut.DS.T1$lambda.min
best_lambda.glut.DS.T1

#produce plot of test MSE by lambda value
plot(cv_model.glut.DS.T1) 

#find coefficients of best model
best_model.glut.DS.T1 <- glmnet(X.glut.DS.T1, Y.glut.DS.T1, alpha = 1, lambda = best_lambda.glut.DS.T1, standardize = TRUE)
coef(best_model.glut.DS.T1)


#use fitted best model to make predictions
y_predicted.glut.DS.T1 <- predict(best_model.glut.DS.T1, s = best_lambda.glut.DS.T1, newx = X.glut.DS.T1)

#find SST and SSE
sst <- sum((Y.glut.DS.T1 - mean(Y.glut.DS.T1))^2)
sse <- sum((y_predicted.glut.DS.T1 - Y.glut.DS.T1)^2)

#find R-Squared
rsq <- 1 - sse/sst
rsq

#Evaluation of the accuracy
Rsquare <- cor(Y.glut.DS.T1,y_predicted.glut.DS.T1)^2
Rsquare

MSE <- mean((Y.glut.DS.T1 - y_predicted.glut.DS.T1)^2)
MSE

plot(Y.glut.DS.T1, y_predicted.glut.DS.T1,col="#00000050",pch = 19,main = "Test result for lasso")
abline(lm(y_predicted.glut.DS.T1~Y.glut.DS.T1),lty = 2,lwd = 2,col = "gray")


# plot lasso model
lam <- best_model.glut.DS.T1$lambda %>% 
  as.data.frame() %>%
  mutate(penalty = best_model.glut.DS.T1$a0 %>% names()) %>%
  rename(lambda = ".")

results.DS.glut.T1 <- best_model.glut.DS.T1$beta %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(penalty, coefficients, -rowname) %>%
  left_join(lam)       

results.DS.glut.T1.filt<-filter_if(results.DS.glut.T1, is.numeric, all_vars((.) != 0))
results.DS.glut.T1.filt

# Visualize the predictive plot
write.table(results.DS.glut.T1.filt,file = here("output/tables", "result_labels_glut_DS.T1.txt"),sep="\t" )



#Protein
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.prot.DS.T1, Y.prot.DS.T1, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.prot.DS.T1, s=s.best.lasso))[-1]
s.best.lasso

lasso.prot.DS.T1<-coef(fit.lasso.prot, s=s.best.lasso)
lasso.prot.DS.T1


#perform k-fold cross-validation to find optimal lambda value
cv_model.prot.DS.T1 <- cv.glmnet(X.prot.DS.T1, Y.prot.DS.T1,K=10, nfolds=nrow(X.prot.DS.T1), alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.prot.DS.T1 <- cv_model.prot.DS.T1$lambda.min
best_lambda.prot.DS.T1

#produce plot of test MSE by lambda value
plot(cv_model.prot.DS.T1) 

#find coefficients of best model
best_model.prot.DS.T1 <- glmnet(X.prot.DS.T1, Y.prot.DS.T1, alpha = 1, lambda = best_lambda.prot.DS.T1, standardize = TRUE )
coef(best_model.prot.DS.T1)


#use fitted best model to make predictions
y_predicted.prot.DS.T1 <- predict(best_model.prot.DS.T1, s = best_lambda.prot.DS.T1, newx = X.prot.DS.T1)
predicted.value.prot.DS.T1<-data.frame(cbind(y_predicted.prot.DS.T1,Y.prot.DS.T1))
#find SST and SSE
sst <- sum((Y.prot.DS.T1 - mean(Y.prot.DS.T1))^2)
sse <- sum((y_predicted.prot.DS.T1 - Y.prot.DS.T1)^2)

#find R-Squared
rsq <- 1 - sse/sst
rsq

Rsquare <- cor(Y.prot.DS.T1,y_predicted.prot.DS.T1)^2
Rsquare

MSE <- mean((Y.prot.DS.T1 - y_predicted.prot.DS.T1)^2)
MSE

plot(Y.prot.DS.T1, y_predicted.prot.DS.T1,col="#00000050",pch = 19,main = "Test result for lasso")
abline(lm(y_predicted.prot.DS.T1~Y.prot.DS.T1),lty = 2,lwd = 2,col = "gray")

#text(x = 16,y = 17,labels = paste0("MSE = ",round(MSE,4),"\n","R-sequared = ",round(Rsquare,4)))

# plot lasso model
lam <- best_model.prot.DS.T1$lambda %>% 
  as.data.frame() %>%
  mutate(penalty = best_model.prot.DS.T1$a0 %>% names()) %>%
  rename(lambda = ".")

results.DS.prot.T1 <- best_model.prot.DS.T1$beta %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(penalty, coefficients, -rowname) %>%
  left_join(lam)       

results.DS.prot.T1.filt<-filter_if(results.DS.prot.T1, is.numeric, all_vars((.) != 0))
results.DS.prot.T1.filt

# Visualize the predictive plot
write.table(results.DS.prot.T1.filt,file = here("output/tables","result_labels_prot_DS.T1.txt"),sep="\t" )

plot.prot.DS.T1<-ggplot(predicted.value.prot.DS.T1, aes(x=V2, y=s1)) + 
  geom_point(color="#00000050")+
  geom_smooth(method=lm, color="lightblue")+
  #geom_text(x = 25, y = 31, label = eq(predicted.value.glut.DT.T3$Y.glut.DT.T3,predicted.value.glut.DT.T3$s1), parse = TRUE) +
  
  stat_regline_equation(label.y =2, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y =1.7, aes(label = ..rr.label..)) +
  labs(title="Protein (May10)",
       y=" Predicted protein content ", x = " Observed protein content ")+
  theme_classic() 

plot.prot.DS.T1


ggsave(file=here("output/photo","model.DS.T1.prot.tiff"), plot.prot.DS.T1, height=3.5, width=4.0, units="in", dpi=600)





#PMT
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.pmt.DS.T1, Y.pmt.DS.T1, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.pmt.DS.T1, s=s.best.lasso))[-1]
s.best.lasso

lasso.pmt.DS.T1<-coef(fit.lasso.pmt, s=s.best.lasso)
lasso.pmt.DS.T1


#perform k-fold cross-validation to find optimal lambda value
cv_model.pmt.DS.T1 <- cv.glmnet(X.pmt.DS.T1, Y.pmt.DS.T1, k=10, nfolds=nrow(X.pmt.DS.T1), alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.pmt.DS.T1 <- cv_model.pmt.DS.T1$lambda.min
best_lambda.pmt.DS.T1

#produce plot of test MSE by lambda value
plot(cv_model.pmt.DS.T1) 

#find coefficients of best model
best_model.pmt.DS.T1 <- glmnet(X.pmt.DS.T1, Y.pmt.DS.T1, alpha = 1, lambda = best_lambda.pmt.DS.T1, standardize=TRUE)
coef(best_model.pmt.DS.T1)


#use fitted best model to make predictions
y_predicted.pmt.DS.T1 <- predict(best_model.pmt.DS.T1, s = best_lambda.pmt.DS.T1, newx = X.pmt.DS.T1)

#find SST and SSE
sst <- sum((Y.pmt.DS.T1 - mean(Y.pmt.DS.T1))^2)
sse <- sum((y_predicted.pmt.DS.T1 - Y.pmt.DS.T1)^2)

#find R-Squared
rsq <- 1 - sse/sst
rsq

Rsquare <- cor(Y.pmt.DS.T1,y_predicted.pmt.DS.T1)^2
Rsquare

MSE <- mean((Y.pmt.DS.T1 - y_predicted.pmt.DS.T1)^2)
MSE

plot(Y.pmt.DS.T1, y_predicted.pmt.DS.T1,col="#00000050",pch = 19,main = "Test result for lasso")
abline(lm(y_predicted.pmt.DS.T1~Y.pmt.DS.T1),lty = 2,lwd = 2,col = "gray")



#BEM
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.bem.DS.T1, Y.bem.DS.T1, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.bem.DS.T1, s=s.best.lasso))[-1]
s.best.lasso

lasso.bem.DS.T1<-coef(fit.lasso.bem.DS.T1, s=s.best.lasso)
lasso.bem.DS.T1

#perform k-fold cross-validation to find optimal lambda value
cv_model.bem.DS.T1 <- cv.glmnet(X.bem.DS.T1, Y.bem.DS.T1, k=10, nfolds=nrow(X.bem.DS.T1),alpha = 1,standardize=1,grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.bem.DS.T1 <- cv_model.bem.DS.T1$lambda.min
best_lambda.bem.DS.T1

#produce plot of test MSE by lambda value
plot(cv_model.bem.DS.T1) 

#find coefficients of best model
best_model.bem.DS.T1 <- glmnet(X.bem.DS.T1, Y.bem.DS.T1, alpha = 1, lambda = best_lambda.bem.DS.T1,standardize = TRUE)
coef(best_model.bem.DS.T1)


#use fitted best model to make predictions
y_predicted.bem.DS.T1 <- predict(best_model.bem.DS.T1, s = best_lambda.bem.DS.T1, newx = X.bem.DS.T1)
predicted.value.bem.DS.T1<-data.frame(cbind(y_predicted.bem.DS.T1,Y.bem.DS.T1))


#find SST and SSE
sst <- sum((Y.bem.DS.T1 - mean(Y.bem.DS.T1))^2)
sse <- sum((y_predicted.bem.DS.T1 - Y.bem.DS.T1)^2)

#find R-Squared
rsq <- 1 - sse/sst
rsq

Rsquare <- cor(Y.bem.DS.T1,y_predicted.bem.DS.T1)^2
Rsquare

MSE <- mean((Y.bem.DS.T1 - y_predicted.bem.DS.T1)^2)
MSE

plot(Y.bem.DS.T1, y_predicted.bem.DS.T1,col="#00000050",pch = 19,main = "Test result for lasso")
abline(lm(y_predicted.bem.DS.T1~Y.bem.DS.T1),lty = 2,lwd = 2,col = "gray")


# plot lasso model
lam <- best_model.bem.DS.T1$lambda %>% 
  as.data.frame() %>%
  mutate(penalty = best_model.bem.DS.T1$a0 %>% names()) %>%
  rename(lambda = ".")

results.DS.bem.T1 <- best_model.bem.DS.T1$beta %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(penalty, coefficients, -rowname) %>%
  left_join(lam)       

results.DS.bem.T1.filt<-filter_if(results.DS.bem.T1, is.numeric, all_vars((.) != 0))
results.DS.bem.T1.filt

write.table(results.DS.bem.T1.filt, here("output/tables","result_labels_bem_DS.T1.txt"),sep="\t" )

plot.bem.DS.T1<-ggplot(predicted.value.bem.DS.T1, aes(x=V2, y=s1)) + 
  geom_point(color="#00000050")+
  geom_smooth(method=lm, color="lightblue")+
  #geom_text(x = 25, y = 31, label = eq(predicted.value.glut.DT.T3$Y.glut.DT.T3,predicted.value.glut.DT.T3$s1), parse = TRUE) +
  
  stat_regline_equation(label.y =2, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y =1.7, aes(label = ..rr.label..)) +
  labs(title="D.Maximum Torque BEM (DS-10 May)",
       y=" Predicted maximum torque", x = " Observed maximum torque")+
  theme_classic() 

plot.bem.DS.T1


ggsave(file=here("output/photo","model.DS.T1.bem.tiff"), plot.bem.DS.T1, height=3.5, width=4.0, units="in", dpi=600)









#lasso.DT.T1=as.matrix(lasso.DT.T1)
#lasso.df.DT.T1=as.data.frame(lasso.DT.T1)
#write.table(lasso.df.DT.T1,"lasso.df.protien.txt",sep="\t")


model.glut.DS.T1=lm(Gluten~ bact.DS.T1.axis4+bact.DS.T1.axis4+biolog.DS.T1.axis1+biolog.DS.T1.axis4+PD+PD.fun, data=data.glut.DS.T1)


model.protein.DS.T1=lm(Protein.grain~bact.DS.T1.axis4+fun.DS.T1.axis1+fun.DS.T1.axis4+ biolog.DS.T1.axis3+Chao1+Chao1.fun+PD.fun
                       ,data=data.prot.DS.T1)

model.pmt.DS.T1=lm(PMT~biolog.DS.T1.axis3+OTU.richness.fun+PD.fun+F.B.ratio , data=data.pmt.DS.T1)






vif(model.protein.DS.T1)

summary(model.glut.DS.T1)
summary(model.protein.DS.T1)
summary(model.pmt.DS.T1)





# LASSO regression analysis

# Preparing data for LASSO regression for T2

reg.lasso.DS.T2<-cbind(pca.data.point.16S.DS.T2,pca.data.point.ITS.DS.T2,pca.data.point.biolog.DS.T2,div.16S.DS.T2,div.ITS.DS.T2,qpcr.T2.DS,qual.DS.S2)
row.names(reg.lasso.DS.T2)==row.names(qual.DS.S2)
#reg.lasso.DS.T2<-deconstand(reg.lasso.DS.T2,method="standardize")

#qpcr.T1.DT<-log(qpcr.T1.DT)

## Protein model for May 10
data.glut.DS.T2=reg.lasso.DS.T2[!is.na(reg.lasso.DS.T2$Gluten),]
data.prot.DS.T2=reg.lasso.DS.T2[!is.na(reg.lasso.DS.T2$Protein.grain),]
data.pmt.DS.T2=reg.lasso.DS.T2[!is.na(reg.lasso.DS.T2$PMT),]
data.bem.DS.T2=reg.lasso.DS.T2[!is.na(reg.lasso.DS.T2$BEM),]
data.bem.DS.T2=data.bem.DS.T2[-c(22),]


data.explain.glut.DS.T2=data.glut.DS.T2[,-which(names(data.glut.DS.T2) %in% c("Protein.grain", "Proteine.flour","Humidity","Ash","Gluten", "IDC","PMT","BEM"))]
data.explain.prot.DS.T2=data.prot.DS.T2[,-which(names(data.prot.DS.T2) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.pmt.DS.T2=data.pmt.DS.T2[,-which(names(data.pmt.DS.T2) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.bem.DS.T2=data.bem.DS.T2[,-which(names(data.bem.DS.T2) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]

# X-Y data prep
X.glut.DS.T2=as.matrix(scale(data.explain.glut.DS.T2,center = T,scale = T))
#Y.glut.DS.T2=as.matrix(data.glut.DS.T2$Gluten)
Y.glut.DS.T2=as.matrix(scale(data.glut.DS.T2$Gluten,center = T,scale = T))

X.prot.DS.T2=as.matrix(scale(data.explain.prot.DS.T2,center = T,scale = T))
Y.prot.DS.T2=as.matrix(scale(data.prot.DS.T2$Protein.grain,center = T,scale = T))

X.pmt.DS.T2=as.matrix(scale(data.explain.pmt.DS.T2,center = T,scale = T))
Y.pmt.DS.T2=as.matrix(scale(data.pmt.DS.T2$PMT,center = T,scale = T))

X.bem.DS.T2=as.matrix(scale(data.explain.bem.DS.T2,center = T,scale = T))
#Y.bem.DS.T2=as.matrix(data.bem.DS.T2$BEM)
Y.bem.DS.T2=as.matrix(scale(data.bem.DS.T2$BEM,center = T,scale = T))

#LASSO regression
fit.lasso.glut.DS.T2 <-glmnet(X.glut.DS.T2,Y.glut.DS.T2, family="gaussian", alpha=1)
fit.lasso.prot.DS.T2<- glmnet(X.prot.DS.T2,Y.prot.DS.T2, family="gaussian", alpha=1)
fit.lasso.pmt.DS.T2 <- glmnet(X.pmt.DS.T2,Y.pmt.DS.T2, family="gaussian", alpha=1)
fit.lasso.bem.DS.T2 <- glmnet(X.bem.DS.T2,Y.bem.DS.T2, family="gaussian", alpha=1)


#Gluten
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.glut.DS.T2, Y.glut.DS.T2, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.glut.DS.T2, s=s.best.lasso))[-1]
s.best.lasso

lasso.glut.DS.T2<-coef(fit.lasso.glut.DS.T2, s=s.best.lasso)
lasso.glut.DS.T2


#perform k-fold cross-validation to find optimal lambda value
cv_model.glut.DS.T2 <- cv.glmnet(X.glut.DS.T2, Y.glut.DS.T2, k=10,nfolds = nrow(X.glut.DS.T2), alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.glut.DS.T2 <- cv_model.glut.DS.T2$lambda.min
best_lambda.glut.DS.T2

#produce plot of test MSE by lambda value
plot(cv_model.glut.DS.T2) 

#find coefficients of best model
best_model.glut.DS.T2 <- glmnet(X.glut.DS.T2, Y.glut.DS.T2, alpha = 1, lambda = best_lambda.glut.DS.T2, standardize = TRUE)
coef(best_model.glut.DS.T2)


#use fitted best model to make predictions
y_predicted.glut.DS.T2 <- predict(best_model.glut.DS.T2, s = best_lambda.glut.DS.T2, newx = X.glut.DS.T2)

#find SST and SSE
sst <- sum((Y.glut.DS.T2 - mean(Y.glut.DS.T2))^2)
sse <- sum((y_predicted.glut.DS.T2 - Y.glut.DS.T2)^2)

#find R-Squared
rsq <- 1 - sse/sst
rsq

#Evaluation of the accuracy
Rsquare <- cor(Y.glut.DS.T2,y_predicted.glut.DS.T2)^2
Rsquare

MSE <- mean((Y.glut.DS.T2 - y_predicted.glut.DS.T2)^2)
MSE

plot(Y.glut.DS.T2, y_predicted.glut.DS.T2,col="#00000050",pch = 19,main = "Test result for lasso")
abline(lm(y_predicted.glut.DS.T2~Y.glut.DS.T2),lty = 2,lwd = 2,col = "gray")




#Protein
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.prot.DS.T2, Y.prot.DS.T2, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.prot.DS.T2, s=s.best.lasso))[-1]
s.best.lasso

lasso.prot.DS.T2<-coef(fit.lasso.prot.DS.T2, s=s.best.lasso)
lasso.prot.DS.T2

#perform k-fold cross-validation to find optimal lambda value
cv_model.prot.DS.T2 <- cv.glmnet(X.prot.DS.T2, Y.prot.DS.T2, k=10,nfolds = nrow(X.prot.DS.T2), alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.prot.DS.T2 <- cv_model.prot.DS.T2$lambda.min
best_lambda.prot.DS.T2

#produce plot of test MSE by lambda value
plot(cv_model.prot.DS.T2) 

#find coefficients of best model
best_model.prot.DS.T2 <- glmnet(X.prot.DS.T2, Y.prot.DS.T2, alpha = 1, lambda = best_lambda.prot.DS.T2, standardize = TRUE )
coef(best_model.prot.DS.T2)


#use fitted best model to make predictions
y_predicted.prot.DS.T2<- predict(best_model.prot.DS.T2, s = best_lambda.prot.DS.T2, newx = X.prot.DS.T2)


Rsquare <- cor(Y.prot.DS.T2,y_predicted.prot.DS.T2)^2
Rsquare

MSE <- mean((Y.prot.DS.T2 - y_predicted.prot.DS.T2)^2)
MSE

plot(Y.prot.DS.T1, y_predicted.prot.DS.T1,col="#00000050",pch = 19,main = "Test result for lasso")
abline(lm(y_predicted.prot.DS.T1~Y.prot.DS.T1),lty = 2,lwd = 2,col = "gray")

#PMT
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.pmt.DS.T2, Y.pmt.DS.T2, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.pmt.DS.T2, s=s.best.lasso))[-1]
s.best.lasso

lasso.pmt.DS.T2<-coef(fit.lasso.pmt.DS.T2, s=s.best.lasso)
lasso.pmt.DS.T2


#perform k-fold cross-validation to find optimal lambda value
cv_model.pmt.DS.T2 <- cv.glmnet(X.pmt.DS.T2, Y.pmt.DS.T2, k=10,nfolds = nrow(X.pmt.DS.T2), alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.pmt.DS.T2 <- cv_model.pmt.DS.T2$lambda.min
best_lambda.pmt.DS.T2

#produce plot of test MSE by lambda value
plot(cv_model.pmt.DS.T2) 

#find coefficients of best model
best_model.pmt.DS.T2 <- glmnet(X.pmt.DS.T2, Y.pmt.DS.T2, alpha = 1, lambda = best_lambda.pmt.DS.T2, standardize = TRUE )
coef(best_model.pmt.DS.T2)


#use fitted best model to make predictions
y_predicted.pmt.DS.T2<- predict(best_model.pmt.DS.T2, s = best_lambda.pmt.DS.T2, newx = X.pmt.DS.T2)


Rsquare <- cor(Y.pmt.DS.T2,y_predicted.pmt.DS.T2)^2
Rsquare

MSE <- mean((Y.pmt.DS.T2 - y_predicted.pmt.DS.T2)^2)
MSE

plot(Y.pmt.DS.T2, y_predicted.pmt.DS.T2,col="#00000050",pch = 19,main = "Test result for lasso")
abline(lm(y_predicted.pmt.DS.T2~Y.pmt.DS.T2),lty = 2,lwd = 2,col = "gray")




#BEM
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.bem.DS.T2, Y.bem.DS.T2, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.bem.DS.T2, s=s.best.lasso))[-1]
s.best.lasso

lasso.bem.DS.T2<-coef(fit.lasso.bem.DS.T2, s=s.best.lasso)
lasso.bem.DS.T2

#perform k-fold cross-validation to find optimal lambda value
cv_model.bem.DS.T2 <- cv.glmnet(X.bem.DS.T2, Y.bem.DS.T2, k=10, nfolds = nrow(X.bem.DS.T2),alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.bem.DS.T2 <- cv_model.bem.DS.T2$lambda.min
best_lambda.bem.DS.T2

#produce plot of test MSE by lambda value
plot(cv_model.bem.DS.T2) 

#find coefficients of best model
best_model.bem.DS.T2 <- glmnet(X.bem.DS.T2, Y.bem.DS.T2, alpha = 1, lambda = best_lambda.bem.DS.T2, standardize = TRUE )
coef(best_model.bem.DS.T2)


#use fitted best model to make predictions
y_predicted.bem.DS.T2<- predict(best_model.bem.DS.T2, s = best_lambda.bem.DS.T2, newx = X.bem.DS.T2)
predicted.value.bem.DS.T2<-data.frame(cbind(y_predicted.bem.DS.T2, Y.bem.DS.T2))


ssr_cv <- t(Y.bem.DS.T2 - y_predicted.bem.DS.T2) %*% (Y.bem.DS.T2- y_predicted.bem.DS.T2)

Rsquare <- cor(Y.bem.DS.T2,y_predicted.bem.DS.T2)^2
Rsquare

MSE <- mean((Y.bem.DS.T2 - y_predicted.bem.DS.T2)^2)
MSE

#Visualization

plot.bem.DS.T2<-ggplot(predicted.value.bem.DS.T2, aes(x=V2, y=s1)) + 
  geom_point(color="#00000050")+
  geom_smooth(method=lm, color="lightblue")+
  #geom_text(x = 25, y = 31, label = eq(predicted.value.glut.DT.T3$Y.glut.DT.T3,predicted.value.glut.DT.T3$s1), parse = TRUE) +
  
  stat_regline_equation(label.y =2, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y =1.6, aes(label = ..rr.label..)) +
  labs(title="C.BEM maximum torque (DS-24 May)",
       y="Predicted maximum torque(BU)", x = "Observed maximum torque (BU)")+
  theme_classic() 

plot.bem.DS.T2


ggsave(file=here("output/photo","model.DS.T2.BEM.tiff"), plot.bem.DS.T2, height=3.5, width=4.0, units="in", dpi=300)





plot(Y.bem.DS.T2, y_predicted.bem.DS.T2,col="#00000050",pch = 19,main = "Test result for lasso")
abline(lm(y_predicted.bem.DS.T2~Y.bem.DS.T2),lty = 2,lwd = 2,col = "gray")


assess.glmnet(best_model.bem.DS.T2,           #in this case, we are evaluating the model
              newx = X.bem.DS.T2,              #in the same data used to fit the model
              newy = Y.bem.DS.T2 )




# LASSO regression analysis

# Preparing data for LASSO regression for T3

reg.lasso.DS.T3<-cbind(pca.data.point.16S.DS.T3,pca.data.point.ITS.DS.T3,pca.data.point.biolog.DS.T3,div.16S.DS.T3,div.ITS.DS.T3,qpcr.T3.DS,qual.DS.S3)
row.names(reg.lasso.DS.T3)==row.names(qual.DS.S3)
#reg.lasso.DS.T3<-decostand(reg.lasso.DS.T3,method="standardize")

#qpcr.T1.DS<-log(qpcr.T1.DT)



## Protein model for May 10
data.glut.DS.T3=reg.lasso.DS.T3[!is.na(reg.lasso.DS.T3$Gluten),]
data.prot.DS.T3=reg.lasso.DS.T3[!is.na(reg.lasso.DS.T3$Protein.grain),]
data.pmt.DS.T3=reg.lasso.DS.T3[!is.na(reg.lasso.DS.T3$PMT),]
data.bem.DS.T3=reg.lasso.DS.T3[!is.na(reg.lasso.DS.T3$BEM),]

data.bem.DS.T3=data.bem.DS.T3[-c(23),]



data.explain.glut.DS.T3=data.glut.DS.T3[,-which(names(data.glut.DS.T3) %in% c("Protein.grain", "Proteine.flour","Humidity","Ash","Gluten", "IDC","PMT","BEM"))]
data.explain.prot.DS.T3=data.prot.DS.T3[,-which(names(data.prot.DS.T3) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.pmt.DS.T3=data.pmt.DS.T3[,-which(names(data.pmt.DS.T3) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.bem.DS.T3=data.bem.DS.T3[,-which(names(data.bem.DS.T3) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]

# X-Y data prep
X.glut.DS.T3=as.matrix(scale(data.explain.glut.DS.T3,center = T,scale = T))
#Y.glut.DS.T3=as.matrix(data.glut.DS.T3$Gluten)
Y.glut.DS.T3=as.matrix(scale(data.glut.DS.T3$Gluten,center = T,scale = T))

X.prot.DS.T3=as.matrix(scale(data.explain.prot.DS.T3,center = T,scale = T))
Y.prot.DS.T3=as.matrix(scale(data.prot.DS.T3$Protein.grain,center = T,scale = T))

X.pmt.DS.T3=as.matrix(scale(data.explain.pmt.DS.T3,center = T,scale = T))
Y.pmt.DS.T3=as.matrix(scale(data.pmt.DS.T3$PMT,center = T,scale = T))

X.bem.DS.T3=as.matrix(scale(data.explain.bem.DS.T3,center = T,scale = T))
Y.bem.DS.T3=as.matrix(scale(data.bem.DS.T3$BEM,center = T,scale = T))




#LASSO regression
fit.lasso.glut.DS.T3 <-glmnet(X.glut.DS.T3,Y.glut.DS.T3, family="gaussian", alpha=1)
fit.lasso.prot.DS.T3<- glmnet(X.prot.DS.T3,Y.prot.DS.T3, family="gaussian", alpha=1)
fit.lasso.pmt.DS.T3 <- glmnet(X.pmt.DS.T3,Y.pmt.DS.T3, family="gaussian", alpha=1)
fit.lasso.bem.DS.T3 <- glmnet(X.bem.DS.T3,Y.bem.DS.T3, family="gaussian", alpha=1)


#Gluten
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.glut.DS.T3, Y.glut.DS.T3, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.glut.DS.T3, s=s.best.lasso))[-1]
s.best.lasso

lasso.glut.DS.T3<-coef(fit.lasso.glut.DS.T3, s=s.best.lasso)
lasso.glut.DS.T3

#perform k-fold cross-validation to find optimal lambda value
cv_model.glut.DS.T3 <- cv.glmnet(X.glut.DS.T3, Y.glut.DS.T3,  nfolds=nrow(X.glut.DS.T3), alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.glut.DS.T3 <- cv_model.glut.DS.T3$lambda.min
best_lambda.glut.DS.T3

#produce plot of test MSE by lambda value
plot(cv_model.glut.DS.T3) 

#find coefficients of best model
best_model.glut.DS.T3 <- glmnet(X.glut.DS.T3, Y.glut.DS.T3, alpha = 1, lambda = best_lambda.glut.DS.T3, standardize = TRUE)
coef(best_model.glut.DS.T3)


#use fitted best model to make predictions
y_predicted.glut.DS.T3 <- predict(best_model.glut.DS.T3, s = best_lambda.glut.DS.T3, newx = X.glut.DS.T3)
predicted.value.glut.DS.T3<-data.frame(cbind(Y.glut.DS.T3,y_predicted.glut.DS.T3))

#Evaluation of the accuracy
Rsquare <- cor(Y.glut.DS.T3,y_predicted.glut.DS.T3)^2
Rsquare

MSE <- mean((Y.glut.DS.T3 - y_predicted.glut.DS.T3)^2)
MSE

plot(Y.glut.DS.T3, y_predicted.glut.DS.T3,col="#00000050",pch = 19,main = "Test result for lasso")
abline(lm(y_predicted.glut.DS.T3~Y.glut.DS.T3),lty = 2,lwd = 2,col = "gray")




# plot lasso model
lam <- best_model.glut.DS.T3$lambda %>% 
  as.data.frame() %>%
  mutate(penalty = best_model.glut.DS.T3$a0 %>% names()) %>%
  rename(lambda = ".")

results.DS.glut.T3 <- best_model.glut.DS.T3$beta %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(penalty, coefficients, -rowname) %>%
  left_join(lam)       

results.DS.glut.T3.filt<-filter_if(results.DS.glut.T3, is.numeric, all_vars((.) != 0))
results.DS.glut.T3.filt

result_labels_gluten_DS.T3 <- results.DS.glut.T3 %>%
  #colSums(result_labels_gluten_DS.T3 != 0)%>%
  group_by(rowname) %>%
  filter(lambda == min(lambda)) %>%
  ungroup() %>%
  top_n(8, wt = abs(coefficients))%>% 
  mutate(var = paste0("x", 1:8))


# Visualize the predictive plot
write.table(results.DS.glut.T3.filt,file = here("output/tables", "result_labels_gluten_DS.T3.txt"),sep="\t" )

eq <- function(x,y) {
  m <- lm(y ~ x)
  as.character(
    as.expression(
      substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                 list(a = format(coef(m)[1], digits = 4),
                      b = format(coef(m)[2], digits = 4),
                      r2 = format(summary(m)$r.squared, digits = 3)))
    )
  )
}


plot.glut.DS.T3<-ggplot(predicted.value.glut.DS.T3, aes(x=V1, y=s1)) + 
  geom_point(color="#00000050")+
  geom_smooth(method=lm, color="lightblue")+
  #geom_text(x = 25, y = 31, label = eq(predicted.value.glut.DT.T3$Y.glut.DT.T3,predicted.value.glut.DT.T3$s1), parse = TRUE) +
  
  stat_regline_equation(label.y =2, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y =1.7, aes(label = ..rr.label..)) +
  labs(title="Gluten (June 07)",
       y=" Predicted gluten content ", x = " Observed gluten content ")+
  theme_classic() 

plot.glut.DS.T3

ggsave(file = here("output/photo","model.DS.T3.glut.tiff"), plot.glut.DS.T3, height=3.5, width=4.0, units="in", dpi=600)



#Protein
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.prot.DS.T3, Y.prot.DS.T3, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.prot.DS.T3, s=s.best.lasso))[-1]
s.best.lasso

lasso.prot.DS.T3<-coef(fit.lasso.prot.DS.T3, s=s.best.lasso)
lasso.prot.DS.T3

#perform k-fold cross-validation to find optimal lambda value
cv_model.prot.DS.T3 <- cv.glmnet(X.prot.DS.T3, Y.prot.DS.T3, k=10, nfolds = nrow(X.prot.DS.T3), alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.prot.DS.T3 <- cv_model.prot.DS.T3$lambda.min
best_lambda.prot.DS.T3

#produce plot of test MSE by lambda value
plot(cv_model.prot.DS.T3) 

#find coefficients of best model
best_model.prot.DS.T3 <- glmnet(X.prot.DS.T3, Y.prot.DS.T3, alpha = 1, lambda = best_lambda.prot.DS.T3, standardize = TRUE)
coef(best_model.prot.DS.T3)


#use fitted best model to make predictions
y_predicted.prot.DS.T3 <- predict(best_model.prot.DS.T3, s = best_lambda.prot.DS.T3, newx = X.prot.DS.T3)


#Evaluation of the accuracy
Rsquare <- cor(Y.prot.DS.T3,y_predicted.prot.DS.T3)^2
Rsquare

MSE <- mean((Y.prot.DS.T3 - y_predicted.prot.DS.T3)^2)
MSE


#PMT
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.pmt.DS.T3, Y.pmt.DS.T3, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.pmt.DS.T3, s=s.best.lasso))[-1]
s.best.lasso

lasso.pmt.DS.T3<-coef(fit.lasso.pmt.DS.T3, s=s.best.lasso)
lasso.pmt.DS.T3

#perform k-fold cross-validation to find optimal lambda value


cv_model.pmt.DS.T3 <- cv.glmnet(X.pmt.DS.T3, Y.pmt.DS.T3, K=10,  alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.pmt.DS.T3 <- cv_model.pmt.DS.T3$lambda.min
best_lambda.pmt.DS.T3
# best lambda 0.3267792# To get this value, repeatedly C.V tuning would be required. 


#produce plot of test MSE by lambda value
plot(cv_model.pmt.DS.T3) 

#find coefficients of best model
best_model.pmt.DS.T3 <- glmnet(X.pmt.DS.T3, Y.pmt.DS.T3, alpha = 1, lambda = best_lambda.pmt.DS.T3, standardize = TRUE)
coef(best_model.pmt.DS.T3)


#use fitted best model to make predictions
y_predicted.pmt.DS.T3 <- predict(best_model.pmt.DS.T3, s = best_lambda.pmt.DS.T3, newx = X.pmt.DS.T3)
predicted.value.pmt.DS.T3<-data.frame(cbind(y_predicted.pmt.DS.T3, Y.pmt.DS.T3))

#Evaluation of the accuracy
Rsquare <- cor(Y.pmt.DS.T3,y_predicted.pmt.DS.T3)^2
Rsquare

MSE <- mean((Y.pmt.DS.T3 - y_predicted.pmt.DS.T3)^2)
MSE


# plot lasso model
lam <- best_model.pmt.DS.T3$lambda %>% 
  as.data.frame() %>%
  mutate(penalty = best_model.pmt.DS.T3$a0 %>% names()) %>%
  rename(lambda = ".")

results.DS.pmt.T3 <- best_model.pmt.DS.T3$beta %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(penalty, coefficients, -rowname) %>%
  left_join(lam)       

results.DS.pmt.T3.filt<-filter_if(results.DS.pmt.T3, is.numeric, all_vars((.) != 0))
results.DS.pmt.T3.filt

# Visualize the predictive plot
write.table(results.DS.pmt.T3.filt,file = here("output/tables","result_labels_pmt_DS.T3.txt"),sep="\t" )

plot.pmt.DS.T3<-ggplot(predicted.value.pmt.DS.T3, aes(x=V2, y=s1)) + 
  geom_point(color="#00000050")+
  geom_smooth(method=lm, color="lightblue")+
  #geom_text(x = 25, y = 31, label = eq(predicted.value.glut.DT.T3$Y.glut.DT.T3,predicted.value.glut.DT.T3$s1), parse = TRUE) +
  
  stat_regline_equation(label.y =1.9, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y =1.7, aes(label = ..rr.label..)) +
  labs(title="PMT (June 07)",
       y="Predicted peak maximum time",x = "Observed peak maximum time" )+
  theme_classic() 

plot.pmt.DS.T3


ggsave(file = here("output/photo","model.DS.T3.pmt.tiff"), plot.pmt.DS.T3, height=3.5, width=4.0, units="in", dpi=600)




#BEM
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.bem.DS.T3, Y.bem.DS.T3, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.bem.DS.T3, s=s.best.lasso))[-1]
s.best.lasso

lasso.bem.DS.T3<-coef(fit.lasso.bem.DS.T3, s=s.best.lasso)
lasso.bem.DS.T3

#perform k-fold cross-validation to find optimal lambda value
cv_model.bem.DS.T3 <- cv.glmnet(X.bem.DS.T3, Y.bem.DS.T3, k=10, nfolds=nrow(X.bem.DS.T3),alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.bem.DS.T3 <- cv_model.bem.DS.T3$lambda.min
best_lambda.bem.DS.T3

#produce plot of test MSE by lambda value
plot(cv_model.bem.DS.T3) 

#find coefficients of best model
best_model.bem.DS.T3 <- glmnet(X.bem.DS.T3, Y.bem.DS.T3, alpha = 1, lambda = best_lambda.bem.DS.T3, standardize = TRUE)
coef(best_model.bem.DS.T3)


#use fitted best model to make predictions
y_predicted.bem.DS.T3 <- predict(best_model.bem.DS.T3, s = best_lambda.bem.DS.T3, newx = X.bem.DS.T3)
predicted.value.pmt.DS.T3<-data.frame(cbind(y_predicted.bem.DS.T3, Y.bem.DS.T3))

#Evaluation of the accuracy
Rsquare <- cor(Y.bem.DS.T3,y_predicted.bem.DS.T3)^2
Rsquare

MSE <- mean((Y.bem.DS.T3 - y_predicted.bem.DS.T3)^2)
MSE

#lasso.DT.T1=as.matrix(lasso.DT.T1)
#lasso.df.DT.T1=as.data.frame(lasso.DT.T1)
#write.table(lasso.df.DT.T1,"lasso.df.protien.txt",sep="\t")


model.glut.DS.T3=lm(Gluten~bact.DS.T3.axis5 +fun.DS.T3.axis2 +fun.DS.T3.axis4+biolog.DS.T3.axis1+biolog.DS.T3.axis3+
                      biolog.DS.T3.axis5+Chao1.fun +AOB   ,
                    data=data.glut.DS.T3)


model.protein.DS.T3=lm(Protein.grain~fun.DS.T3.axis5+Chao1, data=data.glut.DS.T3)

model.pmt.DS.T3=lm(PMT~bact.DS.T3.axis4+fun.DS.T3.axis4+fun.DS.T3.axis5+biolog.DS.T3.axis5, data=data.pmt.DS.T3)

model.bem.DS.T3=lm(BEM~fun.DS.T3.axis3 +fun.DS.T3.axis4  , data=data.bem.DS.T3)



vif(model.protein.DS.T3)

summary(model.glut.DS.T3)
summary(model.protein.DS.T3)
summary(model.pmt.DS.T3)
summary(model.bem.DS.T3)

plot(Y.bem.DS.T3, y_predicted.bem.DS.T3,col="#00000050",pch = 19,main = "Test result for lasso")
abline(lm(y_predicted.bem.DS.T3~Y.bem.DS.T3),lty = 2,lwd = 2,col = "gray")



# Preparing data for LASSO regression for T4

reg.lasso.DS.T4<-cbind(pca.data.point.16S.DS.T4,pca.data.point.ITS.DS.T4,pca.data.point.biolog.DS.T4,div.16S.DS.T4,div.ITS.DS.T4,qpcr.T4.DS,qual.DS.S4)
row.names(reg.lasso.DS.T4)==row.names(qual.DS.S4)

#reg.lasso.DS.T4<-decostand(reg.lasso.DS.T4,method="standardize")

#qpcr.T1.DT<-log(qpcr.T1.DT)

## Protein model for May 10
data.glut.DS.T4=reg.lasso.DS.T4[!is.na(reg.lasso.DS.T4$Gluten),]
data.prot.DS.T4=reg.lasso.DS.T4[!is.na(reg.lasso.DS.T4$Protein.grain),]
data.pmt.DS.T4=reg.lasso.DS.T4[!is.na(reg.lasso.DS.T4$PMT),]
data.bem.DS.T4=reg.lasso.DS.T4[!is.na(reg.lasso.DS.T4$BEM),]
data.bem.DS.T4=data.bem.DS.T4[-c(22),]

data.explain.glut.DS.T4=data.glut.DS.T4[,-which(names(data.glut.DS.T4) %in% c("Protein.grain", "Proteine.flour","Humidity","Ash","Gluten", "IDC","PMT","BEM"))]
data.explain.prot.DS.T4=data.prot.DS.T4[,-which(names(data.glut.DS.T4) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.pmt.DS.T4=data.pmt.DS.T4[,-which(names(data.pmt.DS.T4) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.bem.DS.T4=data.bem.DS.T4[,-which(names(data.bem.DS.T4) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]

# X-Y data prep
X.glut.DS.T4=as.matrix(scale(data.explain.glut.DS.T4,center = T,scale = T))
Y.glut.DS.T4=as.matrix(scale(data.glut.DS.T4$Gluten,center = T,scale = T))

X.prot.DS.T4=as.matrix(scale(data.explain.prot.DS.T4,center = T,scale = T))
Y.prot.DS.T4=as.matrix(data.prot.DS.T4$Protein.grain)

Y.prot.DS.T4=as.matrix(scale(data.prot.DS.T4$Protein.grain,center = T,scale = T))

X.pmt.DS.T4=as.matrix(scale(data.explain.pmt.DS.T4,center = T,scale = T))
Y.pmt.DS.T4=as.matrix(scale(data.pmt.DS.T4$PMT,center = T,scale = T))

X.bem.DS.T4=as.matrix(scale(data.explain.bem.DS.T4,center = T,scale = T))
Y.bem.DS.T4=as.matrix(scale(data.bem.DS.T4$BEM,center = T,scale = T))

#LASSO regression
fit.lasso.glut.DS.T4 <-glmnet(X.glut.DS.T4,Y.glut.DS.T4, family="gaussian", alpha=1)
fit.lasso.prot.DS.T4<- glmnet(X.prot.DS.T4,Y.prot.DS.T4, family="gaussian", alpha=1)
fit.lasso.pmt.DS.T4 <- glmnet(X.pmt.DS.T4,Y.pmt.DS.T4, family="gaussian", alpha=1)
fit.lasso.bem.DS.T4 <- glmnet(X.bem.DS.T4,Y.bem.DS.T4, family="gaussian", alpha=1)


#Gluten
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.glut.DS.T4, Y.glut.DS.T4, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.glut.DS.T4, s=s.best.lasso))[-1]
s.best.lasso

lasso.glut.DS.T4<-coef(fit.lasso.glut.DS.T4, s=s.best.lasso)
lasso.glut.DS.T4


#perform k-fold cross-validation to find optimal lambda value
cv_model.glut.DS.T4 <- cv.glmnet(X.glut.DS.T4,Y.glut.DS.T4, k=10, nfolds = nrow(X.glut.DS.T4), alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.glut.DS.T4 <- cv_model.glut.DS.T4$lambda.min
best_lambda.glut.DS.T4

#produce plot of test MSE by lambda value
plot(cv_model.glut.DS.T4) 

#find coefficients of best model
best_model.glut.DS.T4 <- glmnet(X.glut.DS.T4, Y.glut.DS.T4, alpha = 1, lambda = best_lambda.glut.DS.T4, standardize = TRUE)
coef(best_model.glut.DS.T4)


#use fitted best model to make predictions
y_predicted.glut.DS.T4 <- predict(best_model.glut.DS.T4, s = best_lambda.glut.DS.T4, newx = X.glut.DS.T4)


#Evaluation of the accuracy
Rsquare <- cor(Y.glut.DS.T4,y_predicted.glut.DS.T4)^2
Rsquare

MSE <- mean((Y.glut.DS.T4 - y_predicted.glut.DS.T4)^2)
MSE



# plot lasso model
lam <- best_model.glut.DS.T4$lambda %>% 
  as.data.frame() %>%
  mutate(penalty = best_model.glut.DS.T4$a0 %>% names()) %>%
  rename(lambda = ".")

results.DS.glut.T4 <- best_model.glut.DS.T4$beta %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(penalty, coefficients, -rowname) %>%
  left_join(lam)       

results.DS.glut.T4.filt<-filter_if(results.DS.glut.T4, is.numeric, all_vars((.) != 0))
results.DS.glut.T4.filt



# Visualize the predictive plot
write.table(results.DS.glut.T4.filt,file = here("output/tables", "result_labels_gluten_DS.T4.txt"),sep="\t" )




#Protein
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.prot.DS.T4, Y.prot.DS.T4, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.prot.DS.T4, s=s.best.lasso))[-1]
s.best.lasso

lasso.prot.DS.T4<-coef(fit.lasso.prot.DS.T4, s=s.best.lasso)
lasso.prot.DS.T4

#perform k-fold cross-validation to find optimal lambda value
cv_model.prot.DS.T4 <- cv.glmnet(X.prot.DS.T4,Y.prot.DS.T4, k=10, nfolds = nrow(X.prot.DS.T4), alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.prot.DS.T4 <- cv_model.prot.DS.T4$lambda.min
best_lambda.prot.DS.T4

#produce plot of test MSE by lambda value
plot(cv_model.prot.DS.T4) 

#find coefficients of best model
best_model.prot.DS.T4 <- glmnet(X.prot.DS.T4, Y.prot.DS.T4, alpha = 1, lambda = best_lambda.prot.DS.T4, standardize = TRUE)
coef(best_model.prot.DS.T4)


#use fitted best model to make predictions
y_predicted.prot.DS.T4 <- predict(best_model.prot.DS.T4, s = best_lambda.prot.DS.T4, newx = X.prot.DS.T4)


#Evaluation of the accuracy
Rsquare <- cor(Y.prot.DS.T4,y_predicted.prot.DS.T4)^2
Rsquare

MSE <- mean((Y.prot.DS.T4 - y_predicted.prot.DS.T4)^2)
MSE




#PMT
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.pmt.DS.T4, Y.pmt.DS.T4, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.pmt.DS.T4, s=s.best.lasso))[-1]
s.best.lasso

lasso.pmt.DS.T4<-coef(fit.lasso.pmt.DS.T4, s=s.best.lasso)
lasso.pmt.DS.T4

#perform k-fold cross-validation to find optimal lambda value
cv_model.pmt.DS.T4<- cv.glmnet(X.pmt.DS.T4,Y.pmt.DS.T4, k=10, nfolds = nrow(X.pmt.DS.T4), alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE

best_lambda.pmt.DS.T4 <-cv_model.pmt.DS.T4$lambda.min
best_lambda.pmt.DS.T4

#produce plot of test MSE by lambda value
plot(cv_model.pmt.DS.T4) 

#find coefficients of best model
best_model.pmt.DS.T4 <- glmnet(X.pmt.DS.T4, Y.pmt.DS.T4, alpha = 1, lambda = best_lambda.pmt.DS.T4, standardize = TRUE)
coef(best_model.pmt.DS.T4)


#use fitted best model to make predictions
y_predicted.pmt.DS.T4 <- predict(best_model.pmt.DS.T4, s = best_lambda.pmt.DS.T4, newx = X.pmt.DS.T4)


#Evaluation of the accuracy
Rsquare <- cor(Y.pmt.DS.T4,y_predicted.pmt.DS.T4)^2
Rsquare

MSE <- mean((Y.pmt.DS.T4 - y_predicted.pmt.DS.T4)^2)
MSE




#BEM
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.bem.DS.T4, Y.bem.DS.T4, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.bem.DS.T4, s=s.best.lasso))[-1]
s.best.lasso

lasso.bem.DS.T4<-coef(fit.lasso.bem.DS.T4, s=s.best.lasso)
lasso.bem.DS.T4

#perform k-fold cross-validation to find optimal lambda value
cv_model.bem.DS.T4<- cv.glmnet(X.bem.DS.T4,Y.bem.DS.T4, k=10,nfolds = nrow(X.bem.DS.T4), alpha = 1, standardize=TRUE, grouped=FALSE)

#find optimal lambda value that minimizes test MSE

best_lambda.bem.DS.T4 <-cv_model.bem.DS.T4$lambda.min
best_lambda.bem.DS.T4

#produce plot of test MSE by lambda value
plot(cv_model.bem.DS.T4) 

#find coefficients of best model
best_model.bem.DS.T4 <- glmnet(X.bem.DS.T4, Y.bem.DS.T4, alpha = 1,  lambda = best_lambda.bem.DS.T4, standardize = TRUE)
coef(best_model.bem.DS.T4)


#use fitted best model to make predictions
y_predicted.bem.DS.T4 <- predict(best_model.bem.DS.T4, s = best_lambda.bem.DS.T4, newx = X.bem.DS.T4)
predicted.value.bem.DS.T4<-data.frame(cbind(y_predicted.bem.DS.T4,Y.bem.DS.T4))

#Evaluation of the accuracy
Rsquare <- cor(Y.bem.DS.T4,y_predicted.bem.DS.T4)^2
Rsquare

MSE <- mean((Y.bem.DS.T4 - y_predicted.bem.DS.T4)^2)
MSE


# plot lasso model
lam <- best_model.bem.DS.T4$lambda %>% 
  as.data.frame() %>%
  mutate(penalty = best_model.bem.DS.T4$a0 %>% names()) %>%
  rename(lambda = ".")

results.DS.bem.T4 <- best_model.bem.DS.T4$beta %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(penalty, coefficients, -rowname) %>%
  left_join(lam)       

results.DS.bem.T4.filt<-filter_if(results.DS.bem.T4, is.numeric, all_vars((.) != 0))
results.DS.bem.T4.filt



# Visualize the predictive plot
write.table(results.DS.bem.T4.filt,file = here("output/tables","result_labels_bem_DS.T4.txt"),sep="\t" )



#lasso.DT.T1=as.matrix(lasso.DT.T1)
#lasso.df.DT.T1=as.data.frame(lasso.DT.T1)
#write.table(lasso.df.DT.T1,"lasso.df.protien.txt",sep="\t")


model.glut.DS.T4=lm(Gluten~bact.DS.T4.axis2 +fun.DS.T4.axis3 +biolog.DS.T4.axis2+biolog.DS.T4.axis2+
                      biolog.DS.T4.axis5+Simpson+PD+Chao1+Chao1.fun,
                    data=data.glut.DS.T4)


model.protein.DS.T4=lm(Protein.grain~PD+ACE
                       ,data=data.prot.DS.T4)



model.bem.DS.T4=lm(BEM~Chao1.fun , data=data.bem.DS.T4)


vif(model.protein.DT.T3)

summary(model.glut.DS.T4)
summary(model.protein.DS.T4)

summary(model.bem.DS.T4)


plot.bem.DS.T4<-ggplot(predicted.value.bem.DS.T4, aes(x=V2, y=s1)) + 
  geom_point(color="#00000050")+
  geom_smooth(method=lm, color="lightblue")+
  #geom_text(x = 25, y = 31, label = eq(predicted.value.glut.DT.T3$Y.glut.DT.T3,predicted.value.glut.DT.T3$s1), parse = TRUE) +
  
  stat_regline_equation(label.y =2, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y =1.7, aes(label = ..rr.label..)) +
  labs(title="BEM (June 21)",
       y=" Predicted maximum torque", x = " Observed maximum torque")+
  theme_classic() 

plot.bem.DS.T4


ggsave(here("output/photo","model.DS.T4.bem.tiff"), plot.bem.DS.T4, height=3.5, width=4.0, units="in", dpi=600)



# Preparing data for LASSO regression for T5

reg.lasso.DS.T5<-cbind(pca.data.point.16S.DS.T5,pca.data.point.ITS.DS.T5,pca.data.point.biolog.DS.T5,div.16S.DS.T5,div.ITS.DS.T5,qpcr.T5.DS,qual.DS.S5)
row.names(reg.lasso.DS.T5)==row.names(qual.DS.S5)

#reg.lasso.DS.T5<-decostand(reg.lasso.DS.T5,method="standardize")

#qpcr.T1.DT<-log(qpcr.T1.DT)

## lasso model for May 10
data.glut.DS.T5=reg.lasso.DS.T5[!is.na(reg.lasso.DS.T5$Gluten),]
data.prot.DS.T5=reg.lasso.DS.T5[!is.na(reg.lasso.DS.T5$Protein.grain),]
data.pmt.DS.T5=reg.lasso.DS.T5[!is.na(reg.lasso.DS.T5$PMT),]
data.bem.DS.T5=reg.lasso.DS.T5[!is.na(reg.lasso.DS.T5$BEM),]
data.bem.DS.T5=data.bem.DS.T5[-c(22),]

data.explain.glut.DS.T5=data.glut.DS.T5[,-which(names(data.glut.DS.T5) %in% c("Protein.grain", "Proteine.flour","Humidity","Ash","Gluten", "IDC","PMT","BEM"))]
data.explain.prot.DS.T5=data.prot.DS.T5[,-which(names(data.glut.DS.T5) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.pmt.DS.T5=data.glut.DS.T5[,-which(names(data.glut.DS.T5) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.bem.DS.T5=data.bem.DS.T5[,-which(names(data.bem.DS.T5) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]

# X-Y data prep
X.glut.DS.T5=as.matrix(scale(data.explain.glut.DS.T5,center = T,scale = T))
Y.glut.DS.T5=as.matrix(scale(data.glut.DS.T5$Gluten,center = T,scale = T))

X.prot.DS.T5=as.matrix(scale(data.explain.prot.DS.T5,center = T,scale = T))
Y.prot.DS.T5=as.matrix(scale(data.prot.DS.T5$Protein.grain,center = T,scale = T))

X.pmt.DS.T5=as.matrix(scale(data.explain.pmt.DS.T5,center = T,scale = T))
Y.pmt.DS.T5=as.matrix(scale(data.pmt.DS.T5$PMT,center = T,scale = T))

X.bem.DS.T5=as.matrix(scale(data.explain.bem.DS.T5,center = T,scale = T))
Y.bem.DS.T5=as.matrix(scale(data.bem.DS.T5$BEM,center = T,scale = T))

#LASSO regression
fit.lasso.glut.DS.T5 <-glmnet(X.glut.DS.T5,Y.glut.DS.T5, family="gaussian", alpha=1)
fit.lasso.prot.DS.T5<- glmnet(X.prot.DS.T5,Y.prot.DS.T5, family="gaussian", alpha=1)
fit.lasso.pmt.DS.T5 <- glmnet(X.pmt.DS.T5,Y.pmt.DS.T5, family="gaussian", alpha=1)
fit.lasso.bem.DS.T5 <- glmnet(X.bem.DS.T5,Y.bem.DS.T5, family="gaussian", alpha=1)


#Gluten
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.glut.DS.T5, Y.glut.DS.T5, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.glut.DS.T5, s=s.best.lasso))[-1]
s.best.lasso

lasso.glut.DS.T5<-coef(fit.lasso.glut.DS.T5, s=s.best.lasso)
lasso.glut.DS.T5

#perform k-fold cross-validation to find optimal lambda value
cv_model.glut.DS.T5 <- cv.glmnet(X.glut.DS.T5,Y.glut.DS.T5, k=10, nfolds=nrow(X.glut.DS.T5), alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.glut.DS.T5 <- cv_model.glut.DS.T5$lambda.min
best_lambda.glut.DS.T5

#produce plot of test MSE by lambda value
plot(cv_model.glut.DS.T5) 

#find coefficients of best model
best_model.glut.DS.T5<- glmnet(X.glut.DS.T5, Y.glut.DS.T5, alpha = 1, lambda = best_lambda.glut.DS.T5, standardize = TRUE)
coef(best_model.glut.DS.T5)


#use fitted best model to make predictions
y_predicted.glut.DS.T5 <- predict(best_model.glut.DS.T5, s = best_lambda.glut.DS.T5, newx = X.glut.DS.T5)


#Evaluation of the accuracy
Rsquare <- cor(Y.glut.DS.T5,y_predicted.glut.DS.T5)^2
Rsquare

MSE <- mean((Y.glut.DS.T5- y_predicted.glut.DS.T5)^2)
MSE




#Protein
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.prot.DS.T5, Y.prot.DS.T5, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.prot.DS.T5, s=s.best.lasso))[-1]
s.best.lasso

lasso.prot.DS.T5<-coef(fit.lasso.prot.DS.T5, s=s.best.lasso)
lasso.prot.DS.T5

#perform k-fold cross-validation to find optimal lambda value
cv_model.prot.DS.T5 <- cv.glmnet(X.prot.DS.T5,Y.prot.DS.T5, K=10,  nfolds=nrow(X.prot.DS.T5), alpha = 1, standardize=TRUE)

#find optimal lambda value that minimizes test MSE
best_lambda.prot.DS.T5 <- cv_model.prot.DS.T5$lambda.min
best_lambda.prot.DS.T5

#produce plot of test MSE by lambda value
plot(cv_model.prot.DS.T5) 

#find coefficients of best model
best_model.prot.DS.T5<- glmnet(X.prot.DS.T5, Y.prot.DS.T5, alpha = 1, lambda = best_lambda.prot.DS.T5, standardize = TRUE)
coef(best_model.prot.DS.T5)


#use fitted best model to make predictions
y_predicted.prot.DS.T5 <- predict(best_model.prot.DS.T5, s = best_lambda.prot.DS.T5, newx = X.prot.DS.T5)


#Evaluation of the accuracy
Rsquare <- cor(Y.prot.DS.T5,y_predicted.prot.DS.T5)^2
Rsquare

MSE <- mean((Y.prot.DS.T5- y_predicted.prot.DS.T5)^2)
MSE


# plot lasso model
lam <- best_model.prot.DS.T5$lambda %>% 
  as.data.frame() %>%
  mutate(penalty = best_model.prot.DS.T5$a0 %>% names()) %>%
  rename(lambda = ".")

results.DS.prot.T5 <- best_model.prot.DS.T5$beta %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(penalty, coefficients, -rowname) %>%
  left_join(lam)       

results.DS.prot.T5.filt<-filter_if(results.DS.prot.T5, is.numeric, all_vars((.) != 0))
results.DS.prot.T5.filt



# Visualize the predictive plot
write.table(results.DS.prot.T5.filt, file = here("output/tables","result_labels_prot_DS.T5.txt"),sep="\t" )





#PMT
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.pmt.DS.T5, Y.pmt.DS.T5, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.pmt.DS.T5, s=s.best.lasso))[-1]
s.best.lasso

lasso.pmt.DS.T5<-coef(fit.lasso.pmt.DS.T5, s=s.best.lasso)
lasso.pmt.DS.T5


#perform k-fold cross-validation to find optimal lambda value
cv_model.pmt.DS.T5 <- cv.glmnet(X.pmt.DS.T5,Y.pmt.DS.T5, k=10,nfolds=nrow(X.pmt.DS.T5), alpha = 1, standardize=TRUE, grouped=FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.pmt.DS.T5 <- cv_model.pmt.DS.T5$lambda.min
best_lambda.pmt.DS.T5

#produce plot of test MSE by lambda value
plot(cv_model.pmt.DS.T5) 

#find coefficients of best model
best_model.pmt.DS.T5<- glmnet(X.pmt.DS.T5, Y.pmt.DS.T5, alpha = 1,  lambda = best_lambda.pmt.DS.T5, standardize = TRUE)
coef(best_model.pmt.DS.T5)


#use fitted best model to make predictions
y_predicted.pmt.DS.T5 <- predict(best_model.pmt.DS.T5, s = best_lambda.pmt.DS.T5, newx = X.pmt.DS.T5)


#Evaluation of the accuracy
Rsquare <- cor(Y.pmt.DS.T5,y_predicted.pmt.DS.T5)^2
Rsquare

MSE <- mean((Y.pmt.DS.T5- y_predicted.pmt.DS.T5)^2)
MSE






#BEM
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.bem.DS.T5, Y.bem.DS.T5, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.bem.DS.T5, s=s.best.lasso))[-1]
s.best.lasso

lasso.bem.DS.T5<-coef(fit.lasso.bem.DS.T5, s=s.best.lasso)
lasso.bem.DS.T5


#perform k-fold cross-validation to find optimal lambda value
cv_model.bem.DS.T5 <- cv.glmnet(X.bem.DS.T5,Y.bem.DS.T5, k=10,nfolds=nrow(X.bem.DS.T5), alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.bem.DS.T5 <- cv_model.bem.DS.T5$lambda.min
best_lambda.bem.DS.T5

#produce plot of test MSE by lambda value
plot(cv_model.bem.DS.T5) 

#find coefficients of best model
best_model.bem.DS.T5<- glmnet(X.bem.DS.T5, Y.bem.DS.T5, alpha = 1, lambda = best_lambda.bem.DS.T5, standardize = TRUE)
coef(best_model.bem.DS.T5)


#use fitted best model to make predictions
y_predicted.bem.DS.T5 <- predict(best_model.bem.DS.T5, s = best_lambda.bem.DS.T5, newx = X.bem.DS.T5)


#Evaluation of the accuracy
Rsquare <- cor(Y.bem.DS.T5,y_predicted.bem.DS.T5)^2
Rsquare

MSE <- mean((Y.bem.DS.T5 - y_predicted.bem.DS.T5)^2)
MSE




# plot lasso model
lam <- best_model.bem.DS.T5$lambda %>% 
  as.data.frame() %>%
  mutate(penalty = best_model.bem.DS.T5$a0 %>% names()) %>%
  rename(lambda = ".")

results.DS.bem.T5 <- best_model.bem.DS.T5$beta %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(penalty, coefficients, -rowname) %>%
  left_join(lam)       

results.DS.bem.T5.filt<-filter_if(results.DS.bem.T5, is.numeric, all_vars((.) != 0))
results.DS.bem.T5.filt



# Visualize the predictive plot
write.table(results.DS.bem.T5.filt,here("output/tables","result.bem.DS.T5.txt"), sep="\t" )




#lasso.DT.T1=as.matrix(lasso.DT.T1)
#lasso.df.DT.T1=as.data.frame(lasso.DT.T1)
#write.table(lasso.df.DT.T1,"lasso.df.protien.txt",sep="\t")


model.glut.DS.T5=lm(Gluten~biolog.DS.T5.axis4,
                    data=data.glut.DS.T5)


model.protein.DS.T5=lm(Protein.grain~biolog.DS.T5.axis4
                       ,data=data.prot.DS.T5)

model.pmt.DS.T5=lm(PMT~biolog.DS.T5.axis1+biolog.DS.T5.axis3, data=data.pmt.DS.T5)

model.bem.DS.T5=lm(BEM~biolog.DS.T5.axis5 , data=data.bem.DS.T5)


vif(model.protein.DS.T3)

summary(model.glut.DS.T5)
summary(model.protein.DS.T5)
summary(model.pmt.DS.T5)
summary(model.bem.DS.T5)



# Preparing data for LASSO regression for T6

reg.lasso.DS.T6<-cbind(pca.data.point.16S.DS.T6,pca.data.point.ITS.DS.T6,pca.data.point.biolog.DS.T6,div.16S.DS.T6,div.ITS.DS.T6,qpcr.T6.DS,qual.DS.S6)
row.names(reg.lasso.DS.T6)==row.names(qual.DS.S6)

#reg.lasso.DS.T6<-decostand(reg.lasso.DS.T6,method="standardize")

#qpcr.T1.DT<-log(qpcr.T1.DT)

## lasso model for May 10
data.glut.DS.T6=reg.lasso.DS.T6[!is.na(reg.lasso.DS.T6$Gluten),]
data.prot.DS.T6=reg.lasso.DS.T6[!is.na(reg.lasso.DS.T6$Protein.grain),]
data.pmt.DS.T6=reg.lasso.DS.T6[!is.na(reg.lasso.DS.T6$PMT),]
data.bem.DS.T6=reg.lasso.DS.T6[!is.na(reg.lasso.DS.T6$BEM),]
data.bem.DS.T6=data.bem.DS.T6[-c(22),]


data.explain.glut.DS.T6=data.glut.DS.T6[,-which(names(data.glut.DS.T6) %in% c("Protein.grain", "Proteine.flour","Humidity","Ash","Gluten", "IDC","PMT","BEM"))]
data.explain.prot.DS.T6=data.prot.DS.T6[,-which(names(data.glut.DS.T6) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.pmt.DS.T6=data.glut.DS.T6[,-which(names(data.glut.DS.T6) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.bem.DS.T6=data.bem.DS.T6[,-which(names(data.bem.DS.T6) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]

# X-Y data prep
X.glut.DS.T6=as.matrix(scale(data.explain.glut.DS.T6,center = T,scale = T))
Y.glut.DS.T6=as.matrix(scale(data.glut.DS.T6$Gluten,center = T,scale = T))

X.prot.DS.T6=as.matrix(scale(data.explain.prot.DS.T6,center = T,scale = T))
Y.prot.DS.T6=as.matrix(scale(data.prot.DS.T6$Protein.grain,center = T,scale = T))

X.pmt.DS.T6=as.matrix(scale(data.explain.pmt.DS.T6,center = T,scale = T))
Y.pmt.DS.T6=as.matrix(scale(data.pmt.DS.T6$PMT,center = T,scale = T))

X.bem.DS.T6=as.matrix(scale(data.explain.bem.DS.T6,center = T,scale = T))
Y.bem.DS.T6=as.matrix(scale(data.bem.DS.T6$BEM,center = T,scale = T))



#LASSO regression
fit.lasso.glut.DS.T6 <-glmnet(X.glut.DS.T6,Y.glut.DS.T6, family="gaussian", alpha=1)
fit.lasso.prot.DS.T6<- glmnet(X.prot.DS.T6,Y.prot.DS.T6, family="gaussian", alpha=1)
fit.lasso.pmt.DS.T6 <- glmnet(X.pmt.DS.T6,Y.pmt.DS.T6, family="gaussian", alpha=1)
fit.lasso.bem.DS.T6 <- glmnet(X.bem.DS.T6,Y.bem.DS.T6, family="gaussian", alpha=1)


#Gluten
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.glut.DS.T6, Y.glut.DS.T6, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.glut.DS.T6, s=s.best.lasso))[-1]
s.best.lasso

lasso.glut.DS.T6<-coef(fit.lasso.glut.DS.T6, s=s.best.lasso)
lasso.glut.DS.T6


#perform k-fold cross-validation to find optimal lambda value
cv_model.glut.DS.T6 <- cv.glmnet(X.glut.DS.T6,Y.glut.DS.T6, k=10, nfolds=nrow(X.glut.DS.T6),alpha = 1, standardize=TRUE, grouped=FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.glut.DS.T6 <- cv_model.glut.DS.T6$lambda.min
best_lambda.glut.DS.T6

#produce plot of test MSE by lambda value
plot(cv_model.glut.DS.T6) 

#find coefficients of best model
best_model.glut.DS.T6<- glmnet(X.glut.DS.T6, Y.glut.DS.T6, alpha = 1, family="gaussian", lambda = best_lambda.glut.DS.T6, standardize = TRUE)
coef(best_model.glut.DS.T6)


#use fitted best model to make predictions
y_predicted.glut.DS.T6 <- predict(best_model.glut.DS.T6, s = best_lambda.glut.DS.T6,  newx  = X.glut.DS.T6)
y_predicted.glut.DS.T6


#Evaluation of the accuracy by linear model
Rsquare <- cor(Y.glut.DS.T6,y_predicted.glut.DS.T6)^2
Rsquare



MSE <- mean((Y.glut.DS.T6- y_predicted.glut.DS.T6)^2)
MSE


lam <- best_model.glut.DS.T6$lambda %>% 
  as.data.frame() %>%
  mutate(penalty = best_model.glut.DS.T6$a0 %>% names()) %>%
  rename(lambda = ".")

results.DS.glut.T6 <- best_model.glut.DS.T6$beta %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(penalty, coefficients, -rowname) %>%
  left_join(lam)       

results.DS.glut.T6.filt<-filter_if(results.DS.glut.T6, is.numeric, all_vars((.) != 0))
results.DS.glut.T6.filt



# Visualize the predictive plot
write.table(results.DS.glut.T6.filt,file = here("output/tables", "result_labels_glut_DS.T6.txt"),sep="\t" )







#Protein
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.prot.DS.T6, Y.prot.DS.T6, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.prot.DS.T6, s=s.best.lasso))[-1]
s.best.lasso

lasso.prot.DS.T6<-coef(fit.lasso.prot.DS.T6, s=s.best.lasso)
lasso.prot.DS.T6

#perform k-fold cross-validation to find optimal lambda value
cv_model.prot.DS.T6 <- cv.glmnet(X.prot.DS.T6,Y.prot.DS.T6, k=10,nfolds=nrow(X.prot.DS.T6), alpha = 1, standardize=TRUE, grouped = FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.prot.DS.T6 <- cv_model.prot.DS.T6$lambda.min
best_lambda.prot.DS.T6

#produce plot of test MSE by lambda value
plot(cv_model.prot.DS.T6) 

#find coefficients of best model
best_model.prot.DS.T6<- glmnet(X.prot.DS.T6, Y.prot.DS.T6, alpha = 1, lambda = best_lambda.prot.DS.T6, standardize = TRUE)
coef(best_model.prot.DS.T6)


#use fitted best model to make predictions
y_predicted.prot.DS.T6 <- predict(best_model.prot.DS.T6, s = best_lambda.prot.DS.T6, newx = X.prot.DS.T6)


#Evaluation of the accuracy
Rsquare <- cor(Y.prot.DS.T6,y_predicted.prot.DS.T6)^2
Rsquare

MSE <- mean((Y.prot.DS.T6- y_predicted.prot.DS.T6)^2)
MSE


# Extraction of regression coefficients
lam <- best_model.prot.DS.T6$lambda %>% 
  as.data.frame() %>%
  mutate(penalty = best_model.prot.DS.T6$a0 %>% names()) %>%
  rename(lambda = ".")

results.DS.prot.T6 <- best_model.prot.DS.T6$beta %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(penalty, coefficients, -rowname) %>%
  left_join(lam)       

results.DS.prot.T6.filt<-filter_if(results.DS.prot.T6, is.numeric, all_vars((.) != 0))
results.DS.prot.T6.filt



write.table(results.DS.prot.T6.filt, here("output/tables", "result_labels_prot_DS.T6.txt"),sep="\t" )


#PMT
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.pmt.DS.T6, Y.pmt.DS.T6, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.pmt.DS.T6, s=s.best.lasso))[-1]
s.best.lasso

lasso.pmt.DS.T6<-coef(fit.lasso.pmt.DS.T6, s=s.best.lasso)
lasso.pmt.DS.T6

#perform k-fold cross-validation to find optimal lambda value
cv_model.pmt.DS.T6 <- cv.glmnet(X.pmt.DS.T6,Y.pmt.DS.T6, k=10,nfolds=nrow(X.pmt.DS.T6), alpha = 1, standardize=TRUE)

#find optimal lambda value that minimizes test MSE
best_lambda.pmt.DS.T6 <- cv_model.pmt.DS.T6$lambda.min
best_lambda.pmt.DS.T6

#produce plot of test MSE by lambda value
plot(cv_model.pmt.DS.T6) 

#find coefficients of best model
best_model.pmt.DS.T6<- glmnet(X.pmt.DS.T6, Y.pmt.DS.T6, alpha = 1, lambda = best_lambda.pmt.DS.T6, standardize = TRUE)
coef(best_model.pmt.DS.T6)


#use fitted best model to make predictions
y_predicted.pmt.DS.T6 <- predict(best_model.pmt.DS.T6, s = best_lambda.pmt.DS.T6, newx = X.pmt.DS.T6)


#Evaluation of the accuracy
Rsquare <- cor(Y.pmt.DS.T6,y_predicted.pmt.DS.T6)^2
Rsquare

MSE <- mean((Y.pmt.DS.T6- y_predicted.pmt.DS.T6)^2)
MSE


#BEM
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.bem.DS.T6, Y.bem.DS.T6, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.bem.DS.T6, s=s.best.lasso))[-1]
s.best.lasso

lasso.bem.DS.T6<-coef(fit.lasso.bem.DS.T6, s=s.best.lasso)
lasso.bem.DS.T6


#perform k-fold cross-validation to find optimal lambda value
cv_model.bem.DS.T6 <- cv.glmnet(X.bem.DS.T6,Y.bem.DS.T6, k=10,nfolds=nrow(X.bem.DS.T6), alpha = 1, standardize=TRUE)

#find optimal lambda value that minimizes test MSE
best_lambda.bem.DS.T6 <- cv_model.bem.DS.T6$lambda.min
best_lambda.bem.DS.T6

#produce plot of test MSE by lambda value
plot(cv_model.bem.DS.T6) 

#find coefficients of best model
best_model.bem.DS.T6<- glmnet(X.bem.DS.T6, Y.bem.DS.T6, alpha = 1, lambda = best_lambda.bem.DS.T6, standardize = TRUE)
coef(best_model.bem.DS.T6)


#use fitted best model to make predictions
y_predicted.bem.DS.T6 <- predict(best_model.bem.DS.T6, s = best_lambda.bem.DS.T6, newx = X.bem.DS.T6)


#Evaluation of the accuracy
Rsquare <- cor(Y.bem.DS.T6,y_predicted.bem.DS.T6)^2
Rsquare

MSE <- mean((Y.bem.DS.T6- y_predicted.bem.DS.T6)^2)
MSE


#lasso.DT.T1=as.matrix(lasso.DT.T1)
#lasso.df.DT.T1=as.data.frame(lasso.DT.T1)
#write.table(lasso.df.DT.T1,"lasso.df.protien.txt",sep="\t")


model.glut.DS.T6=lm(Gluten~bact.DS.T6.axis2 +fun.DT.T6.axis4+biolog.DS.T6.axis3+biolog.DS.T6.axis4+Simpson+Chao1+PD.fun ,
                    data=data.glut.DS.T6)


model.protein.DS.T6=lm(Protein.grain~bact.DS.T6.axis2+fun.DT.T6.axis4+fun.DS.T6.axis5+biolog.DS.T6.axis1+OTU.richness+ACE  
                       ,data=data.prot.DS.T6)

model.pmt.DS.T6=lm(PMT~bact.DS.T6.axis4   , data=data.pmt.DS.T6)

model.bem.DS.T6=lm(BEM~bact.DS.T6.axis4   , data=data.bem.DS.T6)


vif(model.protein.DS.T3)

summary(model.glut.DS.T6)
summary(model.protein.DS.T6)
summary(model.pmt.DS.T6)
summary(model.bem.DS.T6)



plot(Y.bem.DS.T6, y_predicted.bem.DS.T6,col="#00000050",pch = 19,main = "Test result for lasso")
abline(lm(y_predicted.bem.DS.T6~Y.bem.DS.T6),lty = 2,lwd = 2,col = "gray")





# Preparing data for LASSO regression for T7

reg.lasso.DS.T7<-cbind(pca.data.point.16S.DS.T7,pca.data.point.ITS.DS.T7,div.16S.DS.T7,div.ITS.DS.T7,qpcr.T7.DS,qual.DS.S7)
row.names(reg.lasso.DS.T7)==row.names(qual.DS.S7)

#reg.lasso.DS.T7<-decostand(reg.lasso.DS.T7,method="standardize")

#qpcr.T1.DT<-log(qpcr.T1.DT)

## lasso model for May 10
data.glut.DS.T7=reg.lasso.DS.T7[!is.na(reg.lasso.DS.T7$Gluten),]
data.prot.DS.T7=reg.lasso.DS.T7[!is.na(reg.lasso.DS.T7$Protein.grain),]
data.pmt.DS.T7=reg.lasso.DS.T7[!is.na(reg.lasso.DS.T7$PMT),]
data.bem.DS.T7=reg.lasso.DS.T7[!is.na(reg.lasso.DS.T7$BEM),]
data.bem.DS.T7=data.bem.DS.T7[-c(22),]

data.explain.glut.DS.T7=data.glut.DS.T7[,-which(names(data.glut.DS.T7) %in% c("Protein.grain", "Proteine.flour","Humidity","Ash","Gluten", "IDC","PMT","BEM"))]
data.explain.prot.DS.T7=data.prot.DS.T7[,-which(names(data.glut.DS.T7) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.pmt.DS.T7=data.glut.DS.T7[,-which(names(data.glut.DS.T7) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]
data.explain.bem.DS.T7=data.bem.DS.T7[,-which(names(data.bem.DS.T7) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","Gluten","IDC","PMT","BEM"))]

# X-Y data prep
X.glut.DS.T7=as.matrix(scale(data.explain.glut.DS.T7,center = T,scale = T))
Y.glut.DS.T7=as.matrix(scale(data.glut.DS.T7$Gluten,center = T,scale = T))

X.prot.DS.T7=as.matrix(scale(data.explain.prot.DS.T7,center = T,scale = T))
Y.prot.DS.T7=as.matrix(scale(data.prot.DS.T7$Protein.grain,center = T,scale = T))

X.pmt.DS.T7=as.matrix(scale(data.explain.pmt.DS.T7,center = T,scale = T))
Y.pmt.DS.T7=as.matrix(scale(data.pmt.DS.T7$PMT,center = T,scale = T))

X.bem.DS.T7=as.matrix(scale(data.explain.bem.DS.T7,center = T,scale = T))
Y.bem.DS.T7=as.matrix(scale(data.bem.DS.T7$BEM,center = T,scale = T))

#LASSO regression
fit.lasso.glut.DS.T7 <-glmnet(X.glut.DS.T7,Y.glut.DS.T7, family="gaussian", alpha=1)
fit.lasso.prot.DS.T7<- glmnet(X.prot.DS.T7,Y.prot.DS.T7, family="gaussian", alpha=1)
fit.lasso.pmt.DS.T7 <- glmnet(X.pmt.DS.T7,Y.pmt.DS.T7, family="gaussian", alpha=1)
fit.lasso.bem.DS.T7 <- glmnet(X.bem.DS.T7,Y.bem.DS.T7, family="gaussian", alpha=1)


#Gluten
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.glut.DS.T7, Y.glut.DS.T7, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.glut.DS.T7, s=s.best.lasso))[-1]
s.best.lasso

lasso.glut.DS.T7<-coef(fit.lasso.glut.DS.T7, s=s.best.lasso)
lasso.glut.DS.T7


#perform k-fold cross-validation to find optimal lambda value
cv_model.glut.DS.T7 <- cv.glmnet(X.glut.DS.T7,Y.glut.DS.T7, k=10,nfolds=nrow(X.glut.DS.T7), alpha = 1, standardize=TRUE, grouped=FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.glut.DS.T7 <- cv_model.glut.DS.T7$lambda.min
best_lambda.glut.DS.T7

#produce plot of test MSE by lambda value
plot(cv_model.glut.DS.T7) 

#find coefficients of best model
best_model.glut.DS.T7<- glmnet(X.glut.DS.T7, Y.glut.DS.T7, alpha = 1, lambda = best_lambda.glut.DS.T7, standardize = TRUE)
coef(best_model.glut.DS.T7)


#use fitted best model to make predictions
y_predicted.glut.DS.T7 <- predict(best_model.glut.DS.T7, s = best_lambda.glut.DS.T7, newx = X.glut.DS.T7)


#Evaluation of the accuracy
Rsquare <- cor(Y.glut.DS.T7,y_predicted.glut.DS.T7)^2
Rsquare

MSE <- mean((Y.glut.DS.T7- y_predicted.glut.DS.T7)^2)
MSE


#Protein
n <- 100; s <- 0
for(i in 1:n) 
  s <- s + cv.glmnet(X.prot.DS.T7, Y.prot.DS.T7, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.prot.DS.T7, s=s.best.lasso))[-1]
s.best.lasso

lasso.prot.DS.T7<-coef(fit.lasso.prot.DS.T7, s=s.best.lasso)
lasso.prot.DS.T7


#perform k-fold cross-validation to find optimal lambda value
cv_model.prot.DS.T7 <- cv.glmnet(X.prot.DS.T7,Y.prot.DS.T7, k=10, nfolds = nrow(X.prot.DS.T7), alpha = 1, standardize=TRUE, grouped=FALSE)

#find optimal lambda value that minimizes test MSE
best_lambda.prot.DS.T7 <- cv_model.prot.DS.T7$lambda.min
best_lambda.prot.DS.T7

#produce plot of test MSE by lambda value
plot(cv_model.prot.DS.T7) 

#find coefficients of best model
best_model.prot.DS.T7<- glmnet(X.prot.DS.T7, Y.prot.DS.T7, alpha = 1, lambda = best_lambda.prot.DS.T7, standardize = TRUE)
coef(best_model.prot.DS.T7)


#use fitted best model to make predictions
y_predicted.prot.DS.T7 <- predict(best_model.prot.DS.T7, s = best_lambda.prot.DS.T7,newx = X.prot.DS.T7 )#type= "coefficients"
y_predicted.prot.DT.T7

#Evaluation of the accuracy
Rsquare <- cor(Y.prot.DS.T7,y_predicted.prot.DS.T7)^2
Rsquare

MSE <- mean((Y.prot.DS.T7- y_predicted.prot.DS.T7)^2)
MSE

################################################################################

lam <- best_model.prot.DS.T7$lambda %>% 
  as.data.frame() %>%
  mutate(penalty = best_model.prot.DS.T7$a0 %>% names()) %>%
  rename(lambda = ".")

results.DS.prot.T7 <- best_model.prot.DS.T7$beta %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(penalty, coefficients, -rowname) %>%
  left_join(lam)       

results.DS.prot.T7.filt<-filter_if(results.DS.prot.T7, is.numeric, all_vars((.) != 0))
results.DS.prot.T7.filt



# Visualize the predictive plot
write.table(results.DS.prot.T7.filt, here("output/tables", "result_labels_prot_DS.T7.txt"),sep="\t" )




#PMT
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.pmt.DS.T7, Y.pmt.DS.T7, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n
coef.lasso <-as.vector(coef(fit.lasso.pmt.DS.T7, s=s.best.lasso))[-1]
s.best.lasso

lasso.pmt.DS.T7<-coef(fit.lasso.pmt.DS.T7, s=s.best.lasso)
lasso.pmt.DS.T7

#perform k-fold cross-validation to find optimal lambda value
grid=10^seq(10,-2, length=100)
cv_model.pmt.DS.T7 <- cv.glmnet(X.pmt.DS.T7,Y.pmt.DS.T7, nfolds=nrow(X.pmt.DS.T7), lambda=grid, alpha = 1,  standardize=TRUE)

#find optimal lambda value that minimizes test MSE
best_lambda.pmt.DS.T7 <- cv_model.pmt.DS.T7$lambda.min
best_lambda.pmt.DS.T7



#produce plot of test MSE by lambda value
plot(cv_model.pmt.DS.T7) 

#find coefficients of best model
best_model.pmt.DS.T7<- glmnet(X.pmt.DS.T7, Y.pmt.DS.T7, alpha = 1, lambda = best_lambda.pmt.DS.T7, standardize = TRUE)
coef(best_model.pmt.DS.T7)


#use fitted best model to make predictions with the train data set
y_predicted.pmt.DS.T7 <- predict(best_model.pmt.DS.T7, s = best_lambda.pmt.DS.T7,  newx = X.pmt.DS.T7)

y_predicted.pmt.DS.T7

#Evaluation of the accuracy
Rsquare <- cor(Y.pmt.DS.T7,y_predicted.pmt.DS.T7)^2
Rsquare




MSE <- mean((Y.pmt.DS.T7- y_predicted.pmt.DS.T7)^2)
MSE



#BEM
n <- 100; s <- 0
for(i in 1:n) s <- s + cv.glmnet(X.bem.DS.T7, Y.bem.DS.T7, nfolds=10, alpha=1)$lambda.min
s.best.lasso <- s/n

coef.lasso <-as.vector(coef(fit.lasso.bem.DS.T7, s=s.best.lasso))[-1]
s.best.lasso

lasso.bem.DS.T7<-coef(fit.lasso.bem.DS.T7, s=s.best.lasso)
lasso.bem.DS.T7



#perform k-fold cross-validation to find optimal lambda value
cv_model.bem.DS.T7 <- cv.glmnet(X.bem.DS.T7,Y.bem.DS.T7, k=10,nfolds=nrow(X.bem.DS.T7), alpha = 1, standardize=TRUE)

#find optimal lambda value that minimizes test MSE
best_lambda.bem.DS.T7 <- cv_model.bem.DS.T7$lambda.min
best_lambda.bem.DS.T7

#produce plot of test MSE by lambda value
plot(cv_model.bem.DS.T7) 

#find coefficients of best model
best_model.bem.DS.T7<- glmnet(X.bem.DS.T7, Y.bem.DS.T7, alpha = 1, lambda = best_lambda.bem.DS.T7, standardize = TRUE)
coef(best_model.bem.DS.T7)


#use fitted best model to make predictions
y_predicted.bem.DS.T7 <- predict(best_model.bem.DS.T7, s = best_lambda.bem.DS.T7, newx = X.bem.DS.T7)


#Evaluation of the accuracy
Rsquare <- cor(Y.bem.DS.T7,y_predicted.bem.DS.T7)^2
Rsquare

MSE <- mean((Y.bem.DS.T7- y_predicted.bem.DS.T7)^2)
MSE






#lasso.DT.T1=as.matrix(lasso.DT.T1)
#lasso.df.DT.T1=as.data.frame(lasso.DT.T1)
#write.table(lasso.df.DT.T1,"lasso.df.protien.txt",sep="\t")


model.glut.DS.T7=lm(Gluten~Simpson ,
                    data=data.glut.DS.T7)


model.protein.DS.T7=lm(Protein.grain~bact.DS.T7.axis1   
                       ,data=data.prot.DS.T7)

model.pmt.DS.T7=lm(PMT~ACE+PD.fun , data=data.pmt.DS.T7)




summary(model.glut.DS.T7)
summary(model.protein.DS.T7)
summary(model.pmt.DS.T7)



#Model evaluation

#library(vip)

#vip(fit.lasso.glut.DS.T1, num_features = 15, geom = "point")




#################################################################################

