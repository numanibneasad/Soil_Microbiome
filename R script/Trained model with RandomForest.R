#Random Forest Regression
library(rsample)      # data splitting 
library(randomForest) # basic implementation
library(ranger)       # a faster implementation of random-forest
library(caret)        # an aggregation package for performing many machine learning models
library(h2o)
library(vegan)

#Data preparation for model-DT
data.rf.glut.DT.T1=data.glut.DT.T1[,-which(names(data.glut.DT.T1) %in% c("Protein.grain", "Proteine.flour","Humidity","Ash", "IDC","PMT","BEM"))]
data.rf.glut.DT.T2=data.glut.DT.T2[,-which(names(data.glut.DT.T2) %in% c("Protein.grain", "Proteine.flour","Humidity","Ash", "IDC","PMT","BEM"))]
data.rf.glut.DT.T3=data.glut.DT.T3[,-which(names(data.glut.DT.T3) %in% c("Protein.grain", "Proteine.flour","Humidity","Ash", "IDC","PMT","BEM"))]
data.rf.glut.DT.T3=na.omit(data.rf.glut.DT.T3)
data.rf.glut.DT.T4=data.glut.DT.T4[,-which(names(data.glut.DT.T4) %in% c("Protein.grain", "Proteine.flour","Humidity","Ash", "IDC","PMT","BEM"))]
data.rf.glut.DT.T4=na.omit(data.rf.glut.DT.T4)
data.rf.glut.DT.T5=data.glut.DT.T5[,-which(names(data.glut.DT.T5) %in% c("Protein.grain", "Proteine.flour","Humidity","Ash", "IDC","PMT","BEM"))]
data.rf.glut.DT.T5=na.omit(data.rf.glut.DT.T5)
data.rf.glut.DT.T6=data.glut.DT.T6[,-which(names(data.glut.DT.T6) %in% c("Protein.grain", "Proteine.flour","Humidity","Ash", "IDC","PMT","BEM"))]
data.rf.glut.DT.T7=data.glut.DT.T7[,-which(names(data.glut.DT.T7) %in% c("Protein.grain", "Proteine.flour","Humidity","Ash", "IDC","PMT","BEM"))]
data.rf.glut.DT.T7=na.omit(data.rf.glut.DT.T7)

# Data processing for Drought tolerant cultivar
rf.data.DT.T1<- cbind(otu.16S.DT.T1,otu.ITS.DT.T1,div.16S.DT.T1,div.ITS.DT.T1,biolog.DT.T1,qpcr.T1.DT,qual.DT.S1)
rf.data.DT.T1<-na.omit(rf.data.DT.T1)
rf.data.DT.T2<- cbind(otu.16S.DT.T2,otu.ITS.DT.T2,div.16S.DT.T2,div.ITS.DT.T2,biolog.DT.T2,qpcr.T2.DT,qual.DT.S2)
rf.data.DT.T2<-na.omit(rf.data.DT.T2)
rf.data.DT.T3<- cbind(otu.16S.DT.T3,otu.ITS.DT.T3,div.16S.DT.T3,div.ITS.DT.T3,biolog.DT.T3,qpcr.T3.DT,qual.DT.S3)
rf.data.DT.T3<-na.omit(rf.data.DT.T3)
rf.data.DT.T4<- cbind(otu.16S.DT.T4,otu.ITS.DT.T4,div.16S.DT.T4,div.ITS.DT.T4,biolog.DT.T4,qpcr.T4.DT,qual.DT.S4)
rf.data.DT.T4<-na.omit(rf.data.DT.T4)
rf.data.DT.T5<- cbind(otu.16S.DT.T5,otu.ITS.DT.T4,div.16S.DT.T5,div.ITS.DT.T5,biolog.DT.T5,qpcr.T5.DT,qual.DT.S5)
rf.data.DT.T5<-na.omit(rf.data.DT.T5)
rf.data.DT.T6<- cbind(otu.16S.DT.T6, otu.ITS.DT.T6,div.16S.DT.T6,div.ITS.DT.T6,biolog.DT.T6,qpcr.T6.DT,qual.DT.S6)
rf.data.DT.T6<-na.omit(rf.data.DT.T6)
rf.data.DT.T7<- cbind(otu.16S.DT.T7,otu.ITS.DT.T7,div.16S.DT.T7,div.ITS.DT.T7,biolog.DT.T7,qpcr.T7.DT,qual.DT.S7)
rf.data.DT.T7<-na.omit(rf.data.DT.T7)

rf.data.glut.DT.T1<-rf.data.DT.T1[,-which(names(rf.data.DT.T1) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.prot.DT.T1<-rf.data.DT.T1[,-which(names(rf.data.DT.T1) %in% c("Gluten", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.PMT.DT.T1<-rf.data.DT.T1[,-which(names(rf.data.DT.T1) %in% c("Protein.grain","Gluten", "Proteine.flour", "Humidity","Ash","IDC","BEM"))]
rf.data.BEM.DT.T1<-rf.data.DT.T1[,-which(names(rf.data.DT.T1) %in% c("Gluten","Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT"))]


rf.data.glut.DT.T2<-rf.data.DT.T2[,-which(names(rf.data.DT.T2) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.prot.DT.T2<-rf.data.DT.T2[,-which(names(rf.data.DT.T2) %in% c("Gluten", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.PMT.DT.T2<-rf.data.DT.T2[,-which(names(rf.data.DT.T2) %in% c("Protein.grain","Gluten", "Proteine.flour", "Humidity","Ash","IDC","BEM"))]
rf.data.BEM.DT.T2<-rf.data.DT.T2[,-which(names(rf.data.DT.T2) %in% c("Gluten","Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT"))]


rf.data.glut.DT.T3<-rf.data.DT.T3[,-which(names(rf.data.DT.T3) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.prot.DT.T3<-rf.data.DT.T3[,-which(names(rf.data.DT.T3) %in% c("Gluten", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.PMT.DT.T3<-rf.data.DT.T3[,-which(names(rf.data.DT.T3) %in% c("Protein.grain","Gluten", "Proteine.flour", "Humidity","Ash","IDC","BEM"))]
rf.data.BEM.DT.T3<-rf.data.DT.T3[,-which(names(rf.data.DT.T3) %in% c("Gluten","Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT"))]


rf.data.glut.DT.T4<-rf.data.DT.T4[,-which(names(rf.data.DT.T4) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.prot.DT.T4<-rf.data.DT.T4[,-which(names(rf.data.DT.T4) %in% c("Gluten", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.PMT.DT.T4<-rf.data.DT.T4[,-which(names(rf.data.DT.T4) %in% c("Protein.grain","Gluten", "Proteine.flour", "Humidity","Ash","IDC","BEM"))]
rf.data.BEM.DT.T4<-rf.data.DT.T4[,-which(names(rf.data.DT.T4) %in% c("Gluten","Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT"))]


rf.data.glut.DT.T5<-rf.data.DT.T5[,-which(names(rf.data.DT.T5) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.prot.DT.T5<-rf.data.DT.T5[,-which(names(rf.data.DT.T5) %in% c("Gluten", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.PMT.DT.T5<-rf.data.DT.T5[,-which(names(rf.data.DT.T5) %in% c("Protein.grain","Gluten", "Proteine.flour", "Humidity","Ash","IDC","BEM"))]
rf.data.BEM.DT.T5<-rf.data.DT.T5[,-which(names(rf.data.DT.T5) %in% c("Gluten","Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT"))]


rf.data.glut.DT.T6<-rf.data.DT.T6[,-which(names(rf.data.DT.T6) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.prot.DT.T6<-rf.data.DT.T6[,-which(names(rf.data.DT.T6) %in% c("Gluten", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.PMT.DT.T6<-rf.data.DT.T6[,-which(names(rf.data.DT.T6) %in% c("Protein.grain","Gluten", "Proteine.flour", "Humidity","Ash","IDC","BEM"))]
rf.data.BEM.DT.T6<-rf.data.DT.T6[,-which(names(rf.data.DT.T6) %in% c("Gluten","Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT"))]

rf.data.glut.DT.T7<-rf.data.DT.T7[,-which(names(rf.data.DT.T7) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.prot.DT.T7<-rf.data.DT.T7[,-which(names(rf.data.DT.T7) %in% c("Gluten", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.PMT.DT.T7<-rf.data.DT.T7[,-which(names(rf.data.DT.T7) %in% c("Protein.grain","Gluten", "Proteine.flour", "Humidity","Ash","IDC","BEM"))]
rf.data.BEM.DT.T7<-rf.data.DT.T7[,-which(names(rf.data.DT.T7) %in% c("Gluten","Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT"))]


## With drought  sensitive cultivar

rf.data.DS.T1<- cbind(otu.16S.DS.T1,otu.ITS.DS.T1,div.16S.DS.T1,div.ITS.DS.T1,biolog.DS.T1,qpcr.T1.DS,qual.DS.S1)
rf.data.DS.T1<-na.omit(rf.data.DS.T1)
rf.data.DS.T2<- cbind(otu.16S.DS.T2,otu.ITS.DS.T2,div.16S.DS.T2,div.ITS.DS.T2,biolog.DS.T2,qpcr.T2.DS,qual.DS.S2)
rf.data.DS.T2<-na.omit(rf.data.DS.T2)
rf.data.DS.T3<- cbind(otu.16S.DS.T3,otu.ITS.DS.T3,div.16S.DS.T3,div.ITS.DS.T3,biolog.DS.T3,qpcr.T3.DS,qual.DS.S3)
rf.data.DS.T3<-na.omit(rf.data.DS.T3)
rf.data.DS.T4<- cbind(otu.16S.DS.T4,otu.ITS.DS.T4,div.16S.DS.T4,div.ITS.DS.T4,biolog.DS.T4,qpcr.T4.DS,qual.DS.S4)
rf.data.DS.T4<-na.omit(rf.data.DS.T4)
rf.data.DS.T5<- cbind(otu.16S.DS.T5,otu.ITS.DS.T4,div.16S.DS.T5,div.ITS.DS.T5,biolog.DS.T5,qpcr.T5.DS,qual.DS.S5)
rf.data.DS.T5<-na.omit(rf.data.DS.T5)
rf.data.DS.T6<- cbind(otu.16S.DS.T6, otu.ITS.DS.T6,div.16S.DS.T6,div.ITS.DS.T6,biolog.DS.T6,qpcr.T6.DS,qual.DS.S6)
rf.data.DS.T6<-na.omit(rf.data.DS.T6)
rf.data.DS.T7<- cbind(otu.16S.DS.T7,otu.ITS.DS.T7,div.16S.DS.T7,div.ITS.DS.T7,biolog.DS.T7,qpcr.T7.DS,qual.DS.S7)
rf.data.DS.T7<-na.omit(rf.data.DS.T7)

rf.data.glut.DS.T1<-rf.data.DS.T1[,-which(names(rf.data.DS.T1) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.prot.DS.T1<-rf.data.DS.T1[,-which(names(rf.data.DS.T1) %in% c("Gluten", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.PMT.DS.T1<-rf.data.DS.T1[,-which(names(rf.data.DS.T1) %in% c("Protein.grain","Gluten", "Proteine.flour", "Humidity","Ash","IDC","BEM"))]
rf.data.BEM.DS.T1<-rf.data.DS.T1[,-which(names(rf.data.DS.T1) %in% c("Gluten","Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT"))]


rf.data.glut.DS.T2<-rf.data.DS.T2[,-which(names(rf.data.DS.T2) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.prot.DS.T2<-rf.data.DS.T2[,-which(names(rf.data.DS.T2) %in% c("Gluten", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.PMT.DS.T2<-rf.data.DS.T2[,-which(names(rf.data.DS.T2) %in% c("Protein.grain","Gluten", "Proteine.flour", "Humidity","Ash","IDC","BEM"))]
rf.data.BEM.DS.T2<-rf.data.DS.T2[,-which(names(rf.data.DS.T2) %in% c("Gluten","Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT"))]


rf.data.glut.DS.T3<-rf.data.DS.T3[,-which(names(rf.data.DS.T3) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.prot.DS.T3<-rf.data.DS.T3[,-which(names(rf.data.DS.T3) %in% c("Gluten", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.PMT.DS.T3<-rf.data.DS.T3[,-which(names(rf.data.DS.T3) %in% c("Protein.grain","Gluten", "Proteine.flour", "Humidity","Ash","IDC","BEM"))]
rf.data.BEM.DS.T3<-rf.data.DS.T3[,-which(names(rf.data.DS.T3) %in% c("Gluten","Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT"))]


rf.data.glut.DS.T4<-rf.data.DS.T4[,-which(names(rf.data.DS.T4) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.prot.DS.T4<-rf.data.DS.T4[,-which(names(rf.data.DS.T4) %in% c("Gluten", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.PMT.DS.T4<-rf.data.DS.T4[,-which(names(rf.data.DS.T4) %in% c("Protein.grain","Gluten", "Proteine.flour", "Humidity","Ash","IDC","BEM"))]
rf.data.BEM.DS.T4<-rf.data.DS.T4[,-which(names(rf.data.DS.T4) %in% c("Gluten","Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT"))]


rf.data.glut.DS.T5<-rf.data.DS.T5[,-which(names(rf.data.DS.T5) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.prot.DS.T5<-rf.data.DS.T5[,-which(names(rf.data.DS.T5) %in% c("Gluten", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.PMT.DS.T5<-rf.data.DS.T5[,-which(names(rf.data.DS.T5) %in% c("Protein.grain","Gluten", "Proteine.flour", "Humidity","Ash","IDC","BEM"))]
rf.data.BEM.DS.T5<-rf.data.DS.T5[,-which(names(rf.data.DS.T5) %in% c("Gluten","Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT"))]


rf.data.glut.DS.T6<-rf.data.DS.T6[,-which(names(rf.data.DS.T6) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.prot.DS.T6<-rf.data.DS.T6[,-which(names(rf.data.DS.T6) %in% c("Gluten", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.PMT.DS.T6<-rf.data.DS.T6[,-which(names(rf.data.DS.T6) %in% c("Protein.grain","Gluten", "Proteine.flour", "Humidity","Ash","IDC","BEM"))]
rf.data.BEM.DS.T6<-rf.data.DS.T6[,-which(names(rf.data.DS.T6) %in% c("Gluten","Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT"))]

rf.data.glut.DS.T7<-rf.data.DS.T7[,-which(names(rf.data.DS.T7) %in% c("Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.prot.DS.T7<-rf.data.DS.T7[,-which(names(rf.data.DS.T7) %in% c("Gluten", "Proteine.flour", "Humidity","Ash","IDC","PMT","BEM"))]
rf.data.PMT.DS.T7<-rf.data.DS.T7[,-which(names(rf.data.DS.T7) %in% c("Protein.grain","Gluten", "Proteine.flour", "Humidity","Ash","IDC","BEM"))]
rf.data.BEM.DS.T7<-rf.data.DS.T7[,-which(names(rf.data.DS.T7) %in% c("Gluten","Protein.grain", "Proteine.flour", "Humidity","Ash","IDC","PMT"))]


#Random Forest for T1

#Model for Protein Content-T2
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.glut.DT.T1.CV <- train(
  Gluten~ .,
  data = rf.data.glut.DT.T1,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.glut.DT.T1.CV


plot(model.glut.DT.T1.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.glut.DT.T1.final <- train(
  Gluten ~ .,
  data = rf.data.glut.DT.T1,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.glut.DT.T1.final

plot(model.glut.DT.T1.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.glut.DT.T1.final)
plot(rfImp,top = 20)


#Model for Protein Content-T1
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.prot.DT.T1.CV <- train(
  Protein.grain ~ .,
  data = rf.data.prot.DT.T1,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.prot.DT.T1.CV


plot(model.prot.DT.T1.CV)


#Hyper-parameter tuning
set.seed(1)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.prot.DT.T1.final <- train(
  Protein.grain~ .,
  data =rf.data.prot.DT.T1,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.prot.DT.T1.final

plot(model.prot.DT.T1.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.prot.DT.T1.final)
plot(rfImp,top = 25)



#Model for PMT Content-T1
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.PMT.DT.T1.CV <- train(
  PMT ~ .,
  data = rf.data.PMT.DT.T1,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.PMT.DT.T1.CV


plot(model.PMT.DT.T1.CV)



#Hyper-parameter tuning
set.seed(1)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.PMT.DT.T1.final <- train(
  PMT ~ .,
  data = rf.data.PMT.DT.T1,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.PMT.DT.T1.final

plot(model.PMT.DT.T1.final)

# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.PMT.DT.T1.final)
plot(rfImp,top = 10)






#Model for BEM Content-T1
#Cross-validation with 10 fold
# CV control function
set.seed(123)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.BEM.DT.T1.CV <- train(
  BEM ~ .,
  data = rf.data.BEM.DT.T1,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.BEM.DT.T1.CV


plot(model.BEM.DT.T1.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.BEM.DT.T1.final <- train(
  BEM ~ .,
  data = rf.data.BEM.DT.T1,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.BEM.DT.T1.final

plot(model.BEM.DT.T1.final)

# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.BEM.DT.T1.final)
plot(rfImp,top = 25)



############################################################################


#Random Forest for T2

#Model for Gluten Content-T2
#Cross-validation with 8 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.glut.DT.T2.CV <- train(
  Gluten~ .,
  data = rf.data.glut.DT.T2,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.glut.DT.T2.CV


plot(model.glut.DT.T2.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.glut.DT.T2.final <- train(
  Gluten ~ .,
  data = rf.data.glut.DT.T2,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.glut.DT.T2.final

plot(model.glut.DT.T2.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.glut.DT.T2.final)
plot(rfImp,top = 25)


#Model for Protein Content-T2
#Cross-validation with 10 fold
# CV control function
set.seed(123)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.prot.DT.T2.CV <- train(
  Protein.grain ~ .,
  data = rf.data.prot.DT.T2,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.prot.DT.T2.CV


plot(model.prot.DT.T2.CV)


#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.prot.DT.T2.final <- train(
  Protein.grain~ .,
  data =rf.data.prot.DT.T2,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.prot.DT.T2.final

plot(model.prot.DT.T2.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.prot.DT.T2.final)
plot(rfImp,top = 25)



#Model for PMT Content-T2
#Cross-validation with 8 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.PMT.DT.T2.CV <- train(
  PMT ~ .,
  data = rf.data.PMT.DT.T2,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.PMT.DT.T2.CV


plot(model.PMT.DT.T2.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.PMT.DT.T2.final <- train(
  PMT ~ .,
  data = rf.data.PMT.DT.T2,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.PMT.DT.T2.final

plot(model.PMT.DT.T2.final)

# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.PMT.DT.T2.final)
plot(rfImp,top = 25)




#Model for PMT Content-T2
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.BEM.DT.T2.CV <- train(
  BEM ~ .,
  data = rf.data.BEM.DT.T2,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.BEM.DT.T2.CV


plot(model.BEM.DT.T2.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.BEM.DT.T2.final <- train(
  BEM ~ .,
  data = rf.data.BEM.DT.T2,
  method = 'rf',
  importance = TRUE,
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.BEM.DT.T2.final

plot(model.BEM.DT.T2.final)

# Observed the top 10 important feature
library(mlbench) #for data set
library(caret)
library(tidyverse)
library(vimp)

rfImp <- varImp(model.BEM.DT.T2.final, SCALE=TRUE)
plot(rfImp,top = 25)


#imp.bem.DT.T1<-data.frame(model.BEM.DT.T2.final$finalModel$importance[,1]/model.BEM.DT.T2.final$finalModel$importanceSD)

imp.bem.DT.T2 <- data.frame(model.BEM.DT.T2.final$finalModel$importance)

imp.bem.DT.T2$Var.Names <-row.names(imp.glut.DT.T2)
#top_myvar <- head(arrange(imp.bem.DT.T1, desc(X.IncMSE)))


ggplot(top_myvar, aes(x=Var.Names, y=`IncNodePurity`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`IncNodePurity`), color="skyblue") +
  geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )


#Random Forest for T3

#Model for Gluten Content-T3
#Cross-validation with 8fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.glut.DT.T3.CV <- train(
  Gluten~ .,
  data = rf.data.glut.DT.T3,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.glut.DT.T3.CV


plot(model.glut.DT.T2.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.glut.DT.T3.final <- train(
  Gluten ~ .,
  data = rf.data.glut.DT.T3,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.glut.DT.T3.final

plot(model.glut.DT.T3.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.glut.DT.T3.final)
plot(rfImp,top = 20)


#Model for Protein Content-T3
#Cross-validation with 8 fold
# CV control function
set.seed(123)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(123)
model.prot.DT.T3.CV <- train(
  Protein.grain ~ .,
  data = rf.data.prot.DT.T3,
  method = 'rf',
  #preProcess = c("center", "scale"),
  trControl = ctrl
)
model.prot.DT.T3.CV


plot(model.prot.DT.T3.CV)


#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.prot.DT.T3.final <- train(
  Protein.grain~ .,
  data =rf.data.prot.DT.T3,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.prot.DT.T3.final

plot(model.prot.DT.T3.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.prot.DT.T3.final)
plot(rfImp,top = 25)



#Model for PMT Content-T2
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.PMT.DT.T3.CV <- train(
  PMT ~ .,
  data = rf.data.PMT.DT.T3,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.PMT.DT.T3.CV


plot(model.PMT.DT.T3.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.PMT.DT.T3.final <- train(
  PMT ~ .,
  data = rf.data.PMT.DT.T3,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.PMT.DT.T3.final

plot(model.PMT.DT.T3.final)

# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.PMT.DT.T3.final)
plot(rfImp,top = 25)




#Model for PMT Content-T3
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.BEM.DT.T3.CV <- train(
  BEM ~ .,
  data = rf.data.BEM.DT.T3,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.BEM.DT.T3.CV


plot(model.BEM.DT.T3.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.BEM.DT.T3.final <- train(
  BEM ~ .,
  data = rf.data.BEM.DT.T3,
  method = 'rf',
  importance = TRUE,
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.BEM.DT.T3.final

plot(model.BEM.DT.T3.final)


rfImp <- varImp(model.BEM.DT.T3.final)
plot(rfImp,top = 25)

#########################################################


#Random Forest for T4

#Model for Gluten Content-T3
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 10,
)

# Model test with CV
set.seed(1)
model.glut.DT.T4.CV <- train(
  Gluten~ .,
  data = rf.data.glut.DT.T4,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.glut.DT.T4.CV


plot(model.glut.DT.T4.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.glut.DT.T4.final <- train(
  Gluten ~ .,
  data = rf.data.glut.DT.T4,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.glut.DT.T4.final

plot(model.glut.DT.T4.final)

model.glut.DT.T4.final$finalModel$
# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.glut.DT.T4.final)
plot(rfImp,top = 10)

#Model for Protein Content-T4
#Cross-validation with 8 fold
# CV control function
set.seed(123)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(123)
model.prot.DT.T4.CV <- train(
  Protein.grain ~ .,
  data = rf.data.prot.DT.T4,
  method = 'rf',
  #preProcess = c("center", "scale"),
  trControl = ctrl
)
model.prot.DT.T4.CV


plot(model.prot.DT.T4.CV)


#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.prot.DT.T4.final <- train(
  Protein.grain~ .,
  data =rf.data.prot.DT.T4,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.prot.DT.T4.final

plot(model.prot.DT.T4.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.prot.DT.T4.final)
plot(rfImp,top = 25)



#Model for PMT Content-T2
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.PMT.DT.T4.CV <- train(
  PMT ~ .,
  data = rf.data.PMT.DT.T4,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.PMT.DT.T4.CV


plot(model.PMT.DT.T4.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.PMT.DT.T4.final <- train(
  PMT ~ .,
  data = rf.data.PMT.DT.T4,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.PMT.DT.T4.final

plot(model.PMT.DT.T4.final)

# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.PMT.DT.T4.final)
plot(rfImp,top = 25)




#Model for PMT Content-T3
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.BEM.DT.T4.CV <- train(
  BEM ~ .,
  data = rf.data.BEM.DT.T4,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.BEM.DT.T4.CV


plot(model.BEM.DT.T4.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.BEM.DT.T4.final <- train(
  BEM ~ .,
  data = rf.data.BEM.DT.T4,
  method = 'rf',
  importance = TRUE,
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.BEM.DT.T4.final

plot(model.BEM.DT.T4.final)


rfImp <- varImp(model.BEM.DT.T4.final)
plot(rfImp,top = 25)

###############################################


#Random Forest for T5

#Model for Gluten Content-T5
#Cross-validation with 10 fold
# CV control function
set.seed(123)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.glut.DT.T5.CV <- train(
  Gluten~ .,
  data = rf.data.glut.DT.T5,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.glut.DT.T5.CV


plot(model.glut.DT.T5.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.glut.DT.T5.final <- train(
  Gluten ~ .,
  data = rf.data.glut.DT.T5,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.glut.DT.T5.final

plot(model.glut.DT.T5.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.glut.DT.T5.final)
plot(rfImp,top = 20)


#Model for Protein Content-T3
#Cross-validation with 10 fold
# CV control function
set.seed(123)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(123)
model.prot.DT.T5.CV <- train(
  Protein.grain ~ .,
  data = rf.data.prot.DT.T5,
  method = 'rf',
  #preProcess = c("center", "scale"),
  trControl = ctrl
)
model.prot.DT.T5.CV


plot(model.prot.DT.T5.CV)


#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.prot.DT.T5.final <- train(
  Protein.grain~ .,
  data =rf.data.prot.DT.T5,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.prot.DT.T5.final

plot(model.prot.DT.T5.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.prot.DT.T5.final)
plot(rfImp,top = 25)



#Model for PMT Content-T2
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.PMT.DT.T5.CV <- train(
  PMT ~ .,
  data = rf.data.PMT.DT.T5,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.PMT.DT.T5.CV


plot(model.PMT.DT.T5.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.PMT.DT.T5.final <- train(
  PMT ~ .,
  data = rf.data.PMT.DT.T5,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.PMT.DT.T5.final

plot(model.PMT.DT.T5.final)

# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.PMT.DT.T5.final)
plot(rfImp,top = 25)




#Model for PMT Content-T3
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.BEM.DT.T5.CV <- train(
  BEM ~ .,
  data = rf.data.BEM.DT.T5,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.BEM.DT.T5.CV


plot(model.BEM.DT.T5.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.BEM.DT.T5.final <- train(
  BEM ~ .,
  data = rf.data.BEM.DT.T5,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.BEM.DT.T5.final

plot(model.BEM.DT.T5.final)


rfImp <- varImp(model.BEM.DT.T5.final)
plot(rfImp,top = 25)
#################################################


#Random Forest for T6

#Model for Gluten Content-T6
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.glut.DT.T6.CV <- train(
  Gluten~ .,
  data = rf.data.glut.DT.T6,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.glut.DT.T6.CV


plot(model.glut.DT.T6.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.glut.DT.T6.final <- train(
  Gluten ~ .,
  data = rf.data.glut.DT.T6,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.glut.DT.T6.final

plot(model.glut.DT.T6.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.glut.DT.T6.final)
plot(rfImp,top = 20)


#Model for Protein Content-T3
#Cross-validation with 10 fold
# CV control function
set.seed(123)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(123)
model.prot.DT.T6.CV <- train(
  Protein.grain ~ .,
  data = rf.data.prot.DT.T6,
  method = 'rf',
  #preProcess = c("center", "scale"),
  trControl = ctrl
)
model.prot.DT.T6.CV


plot(model.prot.DT.T6.CV)


#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.prot.DT.T6.final <- train(
  Protein.grain~ .,
  data =rf.data.prot.DT.T6,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.prot.DT.T6.final

plot(model.prot.DT.T6.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.prot.DT.T6.final)
plot(rfImp,top = 25)



#Model for PMT Content-T2
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.PMT.DT.T6.CV <- train(
  PMT ~ .,
  data = rf.data.PMT.DT.T6,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.PMT.DT.T6.CV


plot(model.PMT.DT.T6.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.PMT.DT.T6.final <- train(
  PMT ~ .,
  data = rf.data.PMT.DT.T6,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
)

model.PMT.DT.T6.final

plot(model.PMT.DT.T6.final)

# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.PMT.DT.T6.final)
plot(rfImp,top = 25)




#Model for PMT Content-T3
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.BEM.DT.T6.CV <- train(
  BEM ~ .,
  data = rf.data.BEM.DT.T6,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.BEM.DT.T6.CV


plot(model.BEM.DT.T6.CV)



#Hyper-parameter tuning
set.seed(1)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.BEM.DT.T6.final <- train(
  BEM ~ .,
  data = rf.data.BEM.DT.T6,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.BEM.DT.T6.final

plot(model.BEM.DT.T6.final)


rfImp <- varImp(model.BEM.DT.T6.final)
plot(rfImp,top = 25)




###################################################
#Random Forest for T7

#Model for Gluten Content-T7
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.glut.DT.T7.CV <- train(
  Gluten~ .,
  data = rf.data.glut.DT.T7,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.glut.DT.T7.CV


plot(model.glut.DT.T7.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.glut.DT.T7.final <- train(
  Gluten ~ .,
  data = rf.data.glut.DT.T7,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.glut.DT.T7.final

plot(model.glut.DT.T7.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.glut.DT.T7.final)
plot(rfImp,top = 25)


#Model for Protein Content-T3
#Cross-validation with 10 fold
# CV control function
set.seed(123)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(123)
model.prot.DT.T7.CV <- train(
  Protein.grain ~ .,
  data = rf.data.prot.DT.T7,
  method = 'rf',
  #preProcess = c("center", "scale"),
  trControl = ctrl
)
model.prot.DT.T7.CV


plot(model.prot.DT.T7.CV)


#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.prot.DT.T7.final <- train(
  Protein.grain~ .,
  data =rf.data.prot.DT.T7,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.prot.DT.T7.final

plot(model.prot.DT.T7.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.prot.DT.T7.final)
plot(rfImp,top = 25)



#Model for PMT Content-T2
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)




# Model test with CV
set.seed(1)
model.PMT.DT.T7.CV <- train(
  PMT ~ .,
  data = rf.data.PMT.DT.T7,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.PMT.DT.T7.CV


plot(model.PMT.DT.T7.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.PMT.DT.T7.final <- train(
  PMT ~ .,
  data = rf.data.PMT.DT.T7,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.PMT.DT.T7.final

plot(model.PMT.DT.T7.final)

# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.PMT.DT.T7.final)
plot(rfImp,top = 25)




#Model for PMT Content-T7
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.BEM.DT.T7.CV <- train(
  BEM ~ .,
  data = rf.data.BEM.DT.T7,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.BEM.DT.T6.CV


plot(model.BEM.DT.T7.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.BEM.DT.T7.final <- train(
  BEM ~ .,
  data = rf.data.BEM.DT.T7,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.BEM.DT.T7.final

plot(model.BEM.DT.T7.final)


rfImp <- varImp(model.BEM.DT.T7.final)
plot(rfImp,top = 25)

################################################################################

################################################################################

#Random Forest model fro drought sensitive cultivar

#Random Forest for T1

#Model for Protein Content-T2
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.glut.DS.T1.CV <- train(
  Gluten~ .,
  data = rf.data.glut.DS.T1,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.glut.DS.T1.CV


plot(model.glut.DS.T1.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.glut.DS.T1.final <- train(
  Gluten ~ .,
  data = rf.data.glut.DS.T1,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.glut.DS.T1.final

plot(model.glut.DS.T1.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.glut.DS.T1.final)
plot(rfImp,top = 20)


#Model for Protein Content-T1
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.prot.DS.T1.CV <- train(
  Protein.grain ~ .,
  data = rf.data.prot.DS.T1,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.prot.DS.T1.CV


plot(model.prot.DS.T1.CV)


#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.prot.DS.T1.final <- train(
  Protein.grain~ .,
  data =rf.data.prot.DS.T1,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.prot.DS.T1.final

plot(model.prot.DS.T1.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.prot.DS.T1.final)
plot(rfImp,top = 25)



#Model for PMT Content-T1
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.PMT.DS.T1.CV <- train(
  PMT ~ .,
  data = rf.data.PMT.DS.T1,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.PMT.DS.T1.CV


plot(model.PMT.DS.T1.CV)



#Hyper-parameter tuning
set.seed(1)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.PMT.DS.T1.final <- train(
  PMT ~ .,
  data = rf.data.PMT.DS.T1,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.PMT.DS.T1.final

plot(model.PMT.DS.T1.final)

# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.PMT.DS.T1.final)
plot(rfImp,top = 10)






#Model for BEM Content-T1
#Cross-validation with 10 fold
# CV control function
set.seed(123)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.BEM.DS.T1.CV <- train(
  BEM ~ .,
  data = rf.data.BEM.DS.T1,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.BEM.DS.T1.CV


plot(model.BEM.DS.T1.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.BEM.DS.T1.final <- train(
  BEM ~ .,
  data = rf.data.BEM.DS.T1,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.BEM.DS.T1.final

plot(model.BEM.DS.T1.final)

# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.BEM.DS.T1.final)
plot(rfImp,top = 25)

##############################################################
##############################################################


#Random Forest for T2

#Model for Gluten Content-T2
#Cross-validation with 8 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.glut.DS.T2.CV <- train(
  Gluten~ .,
  data = rf.data.glut.DS.T2,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.glut.DS.T2.CV


plot(model.glut.DS.T2.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.glut.DS.T2.final <- train(
  Gluten ~ .,
  data = rf.data.glut.DS.T2,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.glut.DS.T2.final

plot(model.glut.DS.T2.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.glut.DS.T2.final)
plot(rfImp,top = 25)


#Model for Protein Content-T2
#Cross-validation with 10 fold
# CV control function
set.seed(123)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.prot.DS.T2.CV <- train(
  Protein.grain ~ .,
  data = rf.data.prot.DS.T2,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.prot.DS.T2.CV


plot(model.prot.DS.T2.CV)


#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.prot.DS.T2.final <- train(
  Protein.grain~ .,
  data =rf.data.prot.DS.T2,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.prot.DS.T2.final

plot(model.prot.DS.T2.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.prot.DS.T2.final)
plot(rfImp,top = 25)



#Model for PMT Content-T2
#Cross-validation with 8 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.PMT.DS.T2.CV <- train(
  PMT ~ .,
  data = rf.data.PMT.DS.T2,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.PMT.DS.T2.CV


plot(model.PMT.DS.T2.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.PMT.DS.T2.final <- train(
  PMT ~ .,
  data = rf.data.PMT.DS.T2,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.PMT.DS.T2.final

plot(model.PMT.DS.T2.final)

# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.PMT.DS.T2.final)
plot(rfImp,top = 25)




#Model for PMT Content-T2
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.BEM.DS.T2.CV <- train(
  BEM ~ .,
  data = rf.data.BEM.DS.T2,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.BEM.DS.T2.CV


plot(model.BEM.DS.T2.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.BEM.DS.T2.final <- train(
  BEM ~ .,
  data = rf.data.BEM.DS.T2,
  method = 'rf',
  importance = TRUE,
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.BEM.DS.T2.final

plot(model.BEM.DS.T2.final)

# Observed the top 10 important feature
library(mlbench) #for data set
library(caret)
library(tidyverse)
library(vimp)

rfImp <- varImp(model.BEM.DS.T2.final, SCALE=TRUE)
plot(rfImp,top = 10)


#imp.bem.DT.T1<-data.frame(model.BEM.DT.T2.final$finalModel$importance[,1]/model.BEM.DT.T2.final$finalModel$importanceSD)
 
imp.bem.DS.T2 <- data.frame(model.BEM.DS.T2.final$finalModel$importance)

imp.bem.DS.T2$Var.Names <-row.names(imp.glut.DS.T2)
#top_myvar <- head(arrange(imp.bem.DT.T1, desc(X.IncMSE)))


ggplot(top_myvar, aes(x=Var.Names, y=`IncNodePurity`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`IncNodePurity`), color="skyblue") +
  geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )


######################################################################
######################################################################


#Random Forest for T3

#Model for Gluten Content-T3
#Cross-validation with 8fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(123)

model.glut.DS.T3.CV <- train(
  Gluten~ .,
  data = rf.data.glut.DS.T3,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.glut.DS.T3.CV


plot(model.glut.DS.T3.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.glut.DS.T3.final <- train(
  Gluten ~ .,
  data = rf.data.glut.DS.T3,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.glut.DS.T3.final

plot(model.glut.DS.T3.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.glut.DS.T3.final)
plot(rfImp,top = 20)


#Model for Protein Content-T3
#Cross-validation with 8 fold
# CV control function
set.seed(123)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(123)
model.prot.DS.T3.CV <- train(
  Protein.grain ~ .,
  data = rf.data.prot.DS.T3,
  method = 'rf',
  #preProcess = c("center", "scale"),
  trControl = ctrl
)
model.prot.DS.T3.CV


plot(model.prot.DS.T3.CV)


#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.prot.DS.T3.final <- train(
  Protein.grain~ .,
  data =rf.data.prot.DS.T3,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.prot.DS.T3.final

plot(model.prot.DS.T3.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.prot.DS.T3.final)
plot(rfImp,top = 10)



#Model for PMT Content-T2
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.PMT.DS.T3.CV <- train(
  PMT ~ .,
  data = rf.data.PMT.DS.T3,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.PMT.DS.T3.CV


plot(model.PMT.DS.T3.CV)




#Hyper-parameter tuning
set.seed(123)
memory.limit(size=56000)


tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.PMT.DS.T3.final <- train(
  PMT ~ .,
  data = rf.data.PMT.DS.T3,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.PMT.DS.T3.final

plot(model.PMT.DS.T3.final)

# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.PMT.DS.T3.final)
plot(rfImp,top = 25)




#Model for PMT Content-T3
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.BEM.DS.T3.CV <- train(
  BEM ~ .,
  data = rf.data.BEM.DS.T3,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.BEM.DS.T3.CV


plot(model.BEM.DS.T3.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.BEM.DS.T3.final <- train(
  BEM ~ .,
  data = rf.data.BEM.DS.T3,
  method = 'rf',
  importance = TRUE,
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.BEM.DS.T3.final

plot(model.BEM.DS.T3.final)


rfImp <- varImp(model.BEM.DS.T3.final)
plot(rfImp,top = 25)

############################################################################################



#Random Forest for T4

#Model for Gluten Content-T3
#Cross-validation with 10 fold
# CV control function
set.seed(123)
ctrl <- trainControl(
  method = "cv",
  number = 10,
)

# Model test with CV
set.seed(1)
model.glut.DS.T4.CV <- train(
  Gluten~ .,
  data = rf.data.glut.DS.T4,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.glut.DS.T4.CV


plot(model.glut.DS.T4.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.glut.DS.T4.final <- train(
  Gluten ~ .,
  data = rf.data.glut.DS.T4,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.glut.DS.T4.final

plot(model.glut.DS.T4.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.glut.DS.T4.final)
plot(rfImp,top = 25)


#Model for Protein Content-T4
#Cross-validation with 8 fold
# CV control function
set.seed(123)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(123)
model.prot.DS.T4.CV <- train(
  Protein.grain ~ .,
  data = rf.data.prot.DS.T4,
  method = 'rf',
  #preProcess = c("center", "scale"),
  trControl = ctrl
)
model.prot.DS.T4.CV


plot(model.prot.DS.T4.CV)


#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.prot.DS.T4.final <- train(
  Protein.grain~ .,
  data =rf.data.prot.DS.T4,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.prot.DS.T4.final

plot(model.prot.DS.T4.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.prot.DS.T4.final)
plot(rfImp,top = 25)



#Model for PMT Content-T2
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.PMT.DS.T4.CV <- train(
  PMT ~ .,
  data = rf.data.PMT.DS.T4,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.PMT.DS.T4.CV


plot(model.PMT.DS.T4.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.PMT.DS.T4.final <- train(
  PMT ~ .,
  data = rf.data.PMT.DS.T4,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.PMT.DS.T4.final

plot(model.PMT.DS.T4.final)

# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.PMT.DS.T4.final)
plot(rfImp,top = 25)




#Model for PMT Content-T3
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.BEM.DS.T4.CV <- train(
  BEM ~ .,
  data = rf.data.BEM.DS.T4,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.BEM.DS.T4.CV


plot(model.BEM.DS.T4.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.BEM.DS.T4.final <- train(
  BEM ~ .,
  data = rf.data.BEM.DS.T4,
  method = 'rf',
  importance = TRUE,
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.BEM.DS.T4.final

plot(model.BEM.DS.T4.final)


rfImp <- varImp(model.BEM.DS.T4.final)
plot(rfImp,top = 25)


##############################################################################



#Random Forest for T5

#Model for Gluten Content-T5
#Cross-validation with 10 fold
# CV control function
set.seed(123)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.glut.DS.T5.CV <- train(
  Gluten~ .,
  data = rf.data.glut.DS.T5,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.glut.DS.T5.CV


plot(model.glut.DS.T5.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.glut.DS.T5.final <- train(
  Gluten ~ .,
  data = rf.data.glut.DS.T5,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.glut.DS.T5.final

plot(model.glut.DS.T5.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.glut.DS.T5.final)
plot(rfImp,top = 20)


#Model for Protein Content-T3
#Cross-validation with 10 fold
# CV control function
set.seed(123)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(123)
model.prot.DS.T5.CV <- train(
  Protein.grain ~ .,
  data = rf.data.prot.DS.T5,
  method = 'rf',
  #preProcess = c("center", "scale"),
  trControl = ctrl
)
model.prot.DS.T5.CV


plot(model.prot.DS.T5.CV)


#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.prot.DS.T5.final <- train(
  Protein.grain~ .,
  data =rf.data.prot.DS.T5,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.prot.DS.T5.final

plot(model.prot.DS.T5.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.prot.DS.T5.final)
plot(rfImp,top = 25)



#Model for PMT Content-T2
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.PMT.DS.T5.CV <- train(
  PMT ~ .,
  data = rf.data.PMT.DS.T5,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.PMT.DS.T5.CV


plot(model.PMT.DS.T5.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.PMT.DS.T5.final <- train(
  PMT ~ .,
  data = rf.data.PMT.DS.T5,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.PMT.DS.T5.final

plot(model.PMT.DS.T5.final)

# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.PMT.DS.T5.final)
plot(rfImp,top = 25)




#Model for PMT Content-T3
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.BEM.DS.T5.CV <- train(
  BEM ~ .,
  data = rf.data.BEM.DS.T5,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.BEM.DS.T5.CV


plot(model.BEM.DS.T5.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.BEM.DS.T5.final <- train(
  BEM ~ .,
  data = rf.data.BEM.DS.T5,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.BEM.DS.T5.final

plot(model.BEM.DS.T5.final)


rfImp <- varImp(model.BEM.DS.T5.final)
plot(rfImp,top = 25)


##################################################################




#Random Forest for T6

#Model for Gluten Content-T6
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.glut.DS.T6.CV <- train(
  Gluten~ .,
  data = rf.data.glut.DS.T6,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.glut.DS.T6.CV


plot(model.glut.DS.T6.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.glut.DS.T6.final <- train(
  Gluten ~ .,
  data = rf.data.glut.DS.T6,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.glut.DS.T6.final

plot(model.glut.DS.T6.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.glut.DS.T6.final)
plot(rfImp,top = 25)


#Model for Protein Content-T3
#Cross-validation with 10 fold
# CV control function
set.seed(123)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(123)
model.prot.DS.T6.CV <- train(
  Protein.grain ~ .,
  data = rf.data.prot.DS.T6,
  method = 'rf',
  #preProcess = c("center", "scale"),
  trControl = ctrl
)
model.prot.DS.T6.CV


plot(model.prot.DS.T6.CV)


#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.prot.DS.T6.final <- train(
  Protein.grain~ .,
  data =rf.data.prot.DS.T6,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.prot.DS.T6.final

plot(model.prot.DS.T6.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.prot.DS.T6.final)
plot(rfImp,top = 25)



#Model for PMT Content-T2
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(123)
model.PMT.DS.T6.CV <- train(
  PMT ~ .,
  data = rf.data.PMT.DS.T6,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.PMT.DS.T6.CV


plot(model.PMT.DS.T6.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.PMT.DS.T6.final <- train(
  PMT ~ .,
  data = rf.data.PMT.DS.T6,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
)

model.PMT.DS.T6.final

plot(model.PMT.DS.T6.final)

# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.PMT.DS.T6.final)
plot(rfImp,top = 25)




#Model for PMT Content-T3
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.BEM.DS.T6.CV <- train(
  BEM ~ .,
  data = rf.data.BEM.DS.T6,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.BEM.DS.T6.CV


plot(model.BEM.DS.T6.CV)



#Hyper-parameter tuning
set.seed(1)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.BEM.DS.T6.final <- train(
  BEM ~ .,
  data = rf.data.BEM.DS.T6,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.BEM.DS.T6.final

plot(model.BEM.DS.T6.final)


rfImp <- varImp(model.BEM.DS.T6.final)
plot(rfImp,top = 25)




###################################################
#Random Forest for T7

#Model for Gluten Content-T7
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 5,
)

# Model test with CV
set.seed(1)
model.glut.DS.T7.CV <- train(
  Gluten~ .,
  data = rf.data.glut.DS.T7,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.glut.DS.T7.CV


plot(model.glut.DS.T7.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.glut.DS.T7.final <- train(
  Gluten ~ .,
  data = rf.data.glut.DS.T7,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.glut.DS.T7.final

plot(model.glut.DS.T7.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.glut.DS.T7.final)
plot(rfImp,top = 25)


#Model for Protein Content-T3
#Cross-validation with 10 fold
# CV control function
set.seed(123)
ctrl <- trainControl(
  method = "cv",
  number = 5,
)

# Model test with CV
set.seed(123)
model.prot.DS.T7.CV <- train(
  Protein.grain ~ .,
  data = rf.data.prot.DS.T7,
  method = 'rf',
  #preProcess = c("center", "scale"),
  trControl = ctrl
)
model.prot.DS.T7.CV


plot(model.prot.DS.T7.CV)


#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.prot.DS.T7.final <- train(
  Protein.grain~ .,
  data =rf.data.prot.DS.T7,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.prot.DS.T7.final

plot(model.prot.DS.T7.final)


# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.prot.DS.T7.final)
plot(rfImp,top = 25)



#Model for PMT Content-T2
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 5,
)




# Model test with CV
set.seed(1)
model.PMT.DS.T7.CV <- train(
  PMT ~ .,
  data = rf.data.PMT.DS.T7,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.PMT.DS.T7.CV


plot(model.PMT.DS.T7.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.PMT.DS.T7.final <- train(
  PMT ~ .,
  data = rf.data.PMT.DS.T7,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.PMT.DS.T7.final

plot(model.PMT.DS.T7.final)

# Observed the top 10 important feature
library(vimp)

rfImp <- varImp(model.PMT.DS.T7.final)
plot(rfImp,top = 25)




#Model for PMT Content-T3
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.BEM.DS.T7.CV <- train(
  BEM ~ .,
  data = rf.data.BEM.DS.T7,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.BEM.DS.T7.CV


plot(model.BEM.DS.T7.CV)

#Model for PMT Content-T7
#Cross-validation with 10 fold
# CV control function
set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 8,
)

# Model test with CV
set.seed(1)
model.BEM.DS.T7.CV <- train(
  BEM ~ .,
  data = rf.data.BEM.DS.T7,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl
)
model.BEM.DS.T6.CV


plot(model.BEM.DS.T7.CV)



#Hyper-parameter tuning
set.seed(123)

tuneGrid <- expand.grid(
  mtry = c(2:32)
)

model.BEM.DS.T7.final <- train(
  BEM ~ .,
  data = rf.data.BEM.DS.T7,
  method = 'rf',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid 
  #importance="permutation"
)
model.BEM.DS.T7.final

plot(model.BEM.DS.T7.final)


rfImp <- varImp(model.BEM.DS.T7.final)
plot(rfImp,top = 25)
