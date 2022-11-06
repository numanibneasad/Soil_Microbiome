#Author:Numan Ibne Asad
#2022/11/05
#Project:Predictive modeling

#Calculation of AIC and BIC  from the Lasso models
#AIC for BEM.DS.T1
tLL <- best_model.bem.DS.T1$nulldev - deviance(best_model.bem.DS.T1)
k <- best_model.bem.DS.T1$dim[2]
n <- best_model.bem.DS.T1$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL+2*k
AIC

BIC<-log(n)*k - tLL
BIC



#AIC for BEM.DS.T2
tLL <- best_model.bem.DS.T2$nulldev - deviance(best_model.bem.DS.T2)
k <- best_model.bem.DS.T2$dim[2]
n <- best_model.bem.DS.T2$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL+2*k
AIC

BIC<-log(n)*k - tLL
BIC

best_model.bem.DS.T2$dev.ratio

#AIC for BEM.DS.T3
tLL <- best_model.bem.DS.T3$nulldev - deviance(best_model.bem.DS.T3)
k <- best_model.bem.DS.T3$dim[2]
n <- best_model.bem.DS.T3$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL+2*k
AIC

BIC<-log(n)*k - tLL
BIC


#AIC for BEM.DS.T4
tLL <- best_model.bem.DS.T4$nulldev - deviance(best_model.bem.DS.T4)
k <- best_model.bem.DS.T4$dim[2]
n <- best_model.bem.DS.T4$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL+2*k
AIC

BIC<-log(n)*k - tLL
BIC


#AIC for BEM.DS.T5
tLL <- best_model.bem.DS.T5$nulldev - deviance(best_model.bem.DS.T5)
k <- best_model.bem.DS.T5$dim[2]
n <- best_model.bem.DS.T5$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL+2*k
AIC

BIC<-log(n)*k - tLL
BIC





#AIC for BEM.DS.T6
tLL <- best_model.bem.DS.T6$nulldev - deviance(best_model.bem.DS.T6)
k <- best_model.bem.DS.T6$dim[2]
n <- best_model.bem.DS.T6$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL+2*k
AIC

BIC<-log(n)*k - tLL
BIC

#AIC for BEM.DS.T7
tLL <- best_model.bem.DS.T7$nulldev - deviance(best_model.bem.DS.T7)
k <- best_model.bem.DS.T7$dim[2]
n <- best_model.bem.DS.T7$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL+2*k
AIC

BIC<-log(n)*k - tLL
BIC





#AIC for BEM.DT.T1
tLL <- best_model.bem.DT.T1$nulldev - deviance(best_model.bem.DT.T1)
k <- best_model.bem.DT.T1$dim[2]
n <- best_model.bem.DT.T1$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL+2*k
AIC

BIC<-log(n)*k - tLL
BIC

#AIC for BEM.DT.T2
tLL <- best_model.bem.DT.T2$nulldev - deviance(best_model.bem.DT.T2)
k <- best_model.bem.DT.T2$dim[2]
n <- best_model.bem.DT.T2$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL+2*k
AIC

BIC<-log(n)*k - tLL
BIC




#AIC for BEM.DT.T3
tLL <- best_model.bem.DT.T3$nulldev - deviance(best_model.bem.DT.T3)
k <- best_model.bem.DT.T3$dim[2]
n <- best_model.bem.DT.T3$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL+2*k
AIC

BIC<-log(n)*k - tLL
BIC



#AIC for BEM.DT.T5
tLL <- best_model.bem.DT.T5$nulldev - deviance(best_model.bem.DT.T5)
k <- best_model.bem.DT.T5$dim[2]
n <- best_model.bem.DT.T5$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL+2*k
AIC

BIC<-log(n)*k - tLL
BIC


#AIC for BEM.DT.T6
tLL <- best_model.bem.DT.T6$nulldev - deviance(best_model.bem.DT.T6)
k <- best_model.bem.DT.T6$dim[2]
n <- best_model.bem.DT.T6$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL+2*k
AIC

BIC<-log(n)*k - tLL
BIC


#AIC for BEM.DT.T7
tLL <- best_model.bem.DT.T7$nulldev - deviance(best_model.bem.DT.T7)
k <- best_model.bem.DT.T7$dim[2]
n <- best_model.bem.DT.T7$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL+2*k
AIC

BIC<-log(n)*k - tLL
BIC


##AIC for Gluten.DT.T1

tLL <- best_model.glut.DT.T1$nulldev - deviance(best_model.glut.DT.T1)
k <- best_model.glut.DT.T1$dim[2]
n <- best_model.glut.DT.T1$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL+2*k
AIC

BIC<-log(n)*k - tLL
BIC 


##AIC for Gluten.DT.T2

tLL.DT.T2 <- best_model.glut.DT.T2$nulldev - deviance(best_model.glut.DT.T2)
k <- best_model.glut.DT.T2$dim[2]
n <- best_model.glut.DT.T2$nobs
AICc <- -tLL.DT.T2+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.DT.T2+2*k
AIC

BIC<-log(n)*k - tLL.DT.T2
BIC 


##AIC for Gluten.DT.T3

tLL <- best_model.glut.DT.T3$nulldev - deviance(best_model.glut.DT.T3)
k <- best_model.glut.DT.T3$dim[2]
n <- best_model.glut.DT.T3$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL+2*k
AIC

BIC<-log(n)*k - tLL
BIC 

##AIC for Gluten.DT.T4

tLL <- best_model.glut.DT.T4$nulldev - deviance(best_model.glut.DT.T4)
k <- best_model.glut.DT.T4$dim[2]
k <- best_model.glut.DT.T4$df
n <- best_model.glut.DT.T4$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL+2*k
AIC

BIC<-log(n)*k - tLL
BIC 

##AIC for Gluten.DT.T5

tLL <- best_model.glut.DT.T5$nulldev - deviance(best_model.glut.DT.T5)
k <- best_model.glut.DT.T5$dim[2]
k <- best_model.glut.DT.T5$df
n <- best_model.glut.DT.T5$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL+2*k
AIC

BIC<-log(n)*k - tLL
BIC 





##AIC for Gluten.DT.T6

tLL.T6 <- best_model.glut.DT.T6$nulldev - deviance(best_model.glut.DT.T6)
k <- best_model.glut.DT.T6$dim[2]
n <- best_model.glut.DT.T6$nobs
AICc <- -tLL.T6+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.T6+2*k
AIC

BIC<-log(n)*k - tLL.T6
BIC 


##AIC for Gluten.DT.T7

tLL.T7 <- best_model.glut.DT.T7$nulldev - deviance(best_model.glut.DT.T7)
k <- best_model.glut.DT.T7$dim[2]
n <- best_model.glut.DT.T7$nobs
AICc <- -tLL.T6+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.T7+2*k
AIC

BIC<-log(n)*k - tLL.T7
BIC 


##AIC for Gluten.DS.T1

tLL <- best_model.glut.DS.T1$nulldev - deviance(best_model.glut.DS.T1)
k <- best_model.glut.DS.T1$dim[2]
n <- best_model.glut.DS.T1$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL+2*k
AIC

BIC<-log(n)*k - tLL
BIC 


##AIC for Gluten.DS.T2

tLL.T2 <- best_model.glut.DS.T2$nulldev - deviance(best_model.glut.DS.T2)
k <- best_model.glut.DS.T2$dim[2]
n <- best_model.glut.DS.T2$nobs
AICc <- -tLL.T2+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.T2+2*k
AIC

BIC<-log(n)*k - tLL.T2
BIC 


##AIC for Gluten.DS.T3

tLL <- best_model.glut.DS.T3$nulldev - deviance(best_model.glut.DS.T3)
k <- best_model.glut.DS.T3$dim[2]
n <- best_model.glut.DS.T3$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL+2*k
AIC

BIC<-log(n)*k - tLL
BIC 

##AIC for Gluten.DS.T4

tLL.T4 <- best_model.glut.DS.T4$nulldev - deviance(best_model.glut.DS.T4)
k <- best_model.glut.DS.T4$dim[2]
n <- best_model.glut.DS.T4$nobs
AICc <- -tLL.T4+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.T4+2*k
AIC

BIC<-log(n)*k - tLL.T4
BIC 

##AIC for Gluten.DS.T5

tLL.T5 <- best_model.glut.DS.T5$nulldev - deviance(best_model.glut.DS.T5)
k <- best_model.glut.DS.T5$dim[2]
n <- best_model.glut.DS.T5$nobs
AICc <- -tLL.T5+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.T5+2*k
AIC

BIC<-log(n)*k - tLL.T5
BIC 



##AIC for Gluten.DS.T6

tLL.T6 <- best_model.glut.DS.T6$nulldev - deviance(best_model.glut.DS.T6)
k <- best_model.glut.DS.T6$dim[2]
n <- best_model.glut.DS.T6$nobs
AICc <- -tLL.T6+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.T6+2*k
AIC

BIC<-log(n)*k - tLL.T6
BIC 


##AIC for Gluten.DS.T7

tLL.T7 <- best_model.glut.DS.T7$nulldev - deviance(best_model.glut.DS.T7)
k <- best_model.glut.DS.T7$dim[2]
n <- best_model.glut.DS.T7$nobs
AICc <- -tLL.T7+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.T7+2*k
AIC

BIC<-log(n)*k - tLL.T7
BIC 





##AIC for protein.DS.T1

tLL.pro <-best_model.prot.DT.T1$nulldev - deviance(best_model.prot.DT.T1)
k <- best_model.prot.DT.T1$dim[2]
n <- best_model.prot.DT.T1$nobs
AICc <- -tLL.pro+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.pro+2*k
AIC

BIC<-log(n)*k - tLL.pro
BIC 


##AIC for protein.DT.T2
tLL.pro <-best_model.prot.DT.T2$nulldev - deviance(best_model.prot.DT.T2)
k <- best_model.prot.DT.T2$dim[2]
n <- best_model.prot.DT.T2$nobs
AICc <- -tLL.pro+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.pro+2*k
AIC

BIC<-log(n)*k - tLL.pro
BIC 



##AIC for protein.DT.T3

tLL.pro <-best_model.prot.DT.T3$nulldev - deviance(best_model.prot.DT.T3)
k <- best_model.prot.DT.T3$dim[2]
n <- best_model.prot.DT.T3$nobs
AICc <- -tLL.pro+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.pro+2*k
AIC

BIC<-log(n)*k - tLL.pro
BIC 


##AIC for protein.DS.T4

tLL.pro <-best_model.prot.DT.T4$nulldev - deviance(best_model.prot.DT.T4)
k <- best_model.prot.DT.T4$dim[2]
n <- best_model.prot.DT.T4$nobs
AICc <- -tLL.pro+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.pro+2*k
AIC

BIC<-log(n)*k - tLL.pro
BIC 
#AIC for protein.DS.T5

tLL.pro <-best_model.prot.DT.T5$nulldev - deviance(best_model.prot.DT.T5)
k <- best_model.prot.DT.T5$dim[2]
n <- best_model.prot.DT.T5$nobs
AICc <- -tLL.pro+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.pro+2*k
AIC

BIC<-log(n)*k - tLL.pro
BIC 



##AIC for Protein.DS.T6

tLL.pro.t6 <- best_model.prot.DT.T6 $nulldev - deviance(best_model.prot.DT.T6 )
k <- best_model.prot.DT.T6$dim[2]
n <- best_model.prot.DT.T6$nobs
AICc <- -tLL.pro.t6+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.pro.t6+2*k
AIC

BIC<-log(n)*k - tLL.pro.t6
BIC 


##AIC for Protein.DS.T7

tLL.T7 <- best_model.prot.DT.T7$nulldev - deviance(best_model.prot.DT.T7)
k <- best_model.prot.DS.T7$dim[2]
n <- best_model.prot.DS.T7$nobs
AICc <- -tLL.T7+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.T7+2*k
AIC

BIC<-log(n)*k - tLL.T7
BIC 


##AIC for protein.DS.T1

tLL.pro <-best_model.prot.DS.T1$nulldev - deviance(best_model.prot.DS.T1)
k <- best_model.prot.DS.T1$dim[2]
n <- best_model.prot.DS.T1$nobs
AICc <- -tLL.pro+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.pro+2*k
AIC

BIC<-log(n)*k - tLL.pro
BIC 


#AIC for protein.DS.T5

tLL.pro <-best_model.prot.DS.T5$nulldev - deviance(best_model.prot.DS.T5)
k <- best_model.prot.DS.T5$dim[2]
n <- best_model.prot.DS.T5$nobs
AICc <- -tLL.pro+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.pro+2*k
AIC

BIC<-log(n)*k - tLL.pro
BIC 

#AIC for protein.DS.T6

tLL.pro <-best_model.prot.DS.T6$nulldev - deviance(best_model.prot.DS.T6)
k <- best_model.prot.DS.T6$dim[2]
n <- best_model.prot.DS.T6$nobs
AICc <- -tLL.pro+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.pro+2*k
AIC

BIC<-log(n)*k - tLL.pro
BIC 


#AIC for protein.DS.T7

tLL.pro <-best_model.prot.DS.T7$nulldev - deviance(best_model.prot.DS.T7)
k <- best_model.prot.DS.T7$dim[2]
n <- best_model.prot.DS.T7$nobs
AICc <- -tLL.pro+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.pro+2*k
AIC

BIC<-log(n)*k - tLL.pro
BIC 



##AIC for PMT.DT.T1

tLL.pmt <-best_model.pmt.DT.T1$nulldev - deviance(best_model.pmt.DT.T1)
k <- best_model.pmt.DT.T1$dim[2]
n <- best_model.pmt.DT.T1$nobs
AICc <- -tLL.pmt+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.pmt+2*k
AIC

BIC<-log(n)*k - tLL.pmt
BIC 


##AIC for protein.DS.T2

tLL.pmt <-best_model.pmt.DT.T2$nulldev - deviance(best_model.pmt.DT.T2)
k <- best_model.pmt.DT.T2$dim[2]
n <- best_model.pmt.DT.T2$nobs
AICc <- -tLL.pmt+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.pmt+2*k
AIC

BIC<-log(n)*k - tLL.pmt
BIC 



##AIC for protein.DS.T4

tLL.pmt <-best_model.pmt.DT.T4$nulldev - deviance(best_model.pmt.DT.T4)
k <- best_model.pmt.DT.T4$dim[2]
n <- best_model.pmt.DT.T4$nobs
AICc <- -tLL.pmt+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.pmt+2*k
AIC

BIC<-log(n)*k - tLL.pmt
BIC 

##AIC for protein.DS.T7

tLL.pmt <-best_model.pmt.DT.T7$nulldev - deviance(best_model.pmt.DT.T7)
k <- best_model.pmt.DT.T7$dim[2]
n <- best_model.pmt.DT.T7$nobs
AICc <- -tLL.pmt+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.pmt+2*k
AIC

BIC<-log(n)*k - tLL.pmt
BIC 




##AIC for PMT.DS.T1

tLL.pmt <-best_model.pmt.DS.T1$nulldev - deviance(best_model.pmt.DS.T1)
k <- best_model.pmt.DS.T1$dim[2]
n <- best_model.pmt.DS.T1$nobs
AICc <- -tLL.pmt+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.pmt+2*k
AIC

BIC<-log(n)*k - tLL.pmt
BIC 



##AIC for protein.DS.T2

tLL.pmt <-best_model.pmt.DS.T2$nulldev - deviance(best_model.pmt.DS.T2)
k <- best_model.pmt.DS.T2$dim[2]
n <- best_model.pmt.DS.T2$nobs
AICc <- -tLL.pmt+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.pmt+2*k
AIC

BIC<-log(n)*k - tLL.pmt
BIC 

##AIC for protein.DS.T3

tLL.pmt <-best_model.pmt.DS.T3$nulldev - deviance(best_model.pmt.DS.T3)
k <- best_model.pmt.DS.T3$dim[2]
n <- best_model.pmt.DS.T3$nobs
AICc <- -tLL.pmt+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.pmt+2*k
AIC

BIC<-log(n)*k - tLL.pmt
BIC 

##AIC for protein.DS.T5

tLL.pmt <-best_model.pmt.DS.T5$nulldev - deviance(best_model.pmt.DS.T5)
k <- best_model.pmt.DS.T5$dim[2]
n <- best_model.pmt.DS.T5$nobs
AICc <- -tLL.pmt+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.pmt+2*k
AIC

BIC<-log(n)*k - tLL.pmt
BIC 



##AIC for protein.DS.T6

tLL.pmt <-best_model.pmt.DS.T6$nulldev - deviance(best_model.pmt.DS.T6)
k <- best_model.pmt.DS.T6$dim[2]
n <- best_model.pmt.DS.T6$nobs
AICc <- -tLL.pmt+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.pmt+2*k
AIC

BIC<-log(n)*k - tLL.pmt
BIC 

##AIC for protein.DS.T7

tLL.pmt <-best_model.pmt.DS.T7$nulldev - deviance(best_model.pmt.DS.T7)
k <- best_model.pmt.DS.T7$dim[2]
n <- best_model.pmt.DS.T7$nobs
AICc <- -tLL.pmt+2*k+2*k*(k+1)/(n-k-1)
AICc
AIC<- -tLL.pmt+2*k
AIC

BIC<-log(n)*k - tLL.pmt
BIC 


