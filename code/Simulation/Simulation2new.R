#### Simulation 2 #####

# This simulation evaluate the performance of the method without considering
# measurement error in covariates. Different from the Simulation 1, this simulation
# evaluate the case using GEE method rather than likelihood approach.


## Update: Jan 15, 2020
# 1. the Theta is the precision matrix rather than covariance matrix
# 2. make the parameter beta fixed in each simulation


#### 1. Set up ####
### 1.1 Global parameters  ####
set.seed(2019)
seed_i <- sample(1000000,1000)
ProjectName <- paste0("Simulation2")


### 1.2 R packages ####
library(doParallel)
library(scales)
library(huge)
library(MASS)
library(glmnet)
library(GeneErrorMis)
library(ggplot2)
library(nleqslv)
library(xtable)
library(tidyr)
library(tidyverse)


### 1.3 Global Functions ####
expit <- function(x){
  value <- exp(x)/(1+exp(x))
  ifelse(is.na(value),1,value) 
}


remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}


getAUC <- function(fp,tp){
  ord.fp = order(fp)
  tmp1 = fp[ord.fp]
  tmp2 = tp[ord.fp]
  AUC = sum(diff(tmp1) * (tmp2[-1] + tmp2[-length(tmp2)]))/2
  return(AUC)
}

likelihood_GQ<-function(Theta, weights, nodes, Y1, Y2, Covariates, R){
  nbeta <- dim(Covariates)[2]
  rett<-logLik_GQapprox_noerror(weights, nodes, Y1, Y2,Covariates, R,
                  Theta[1:nbeta],Theta[(nbeta+1):(2*nbeta)],Theta[(2*nbeta)+1],Theta[2*nbeta+2])
  return(rett)
}

score_GQ<-function(Theta, weights, nodes, Y1, Y2, Covariates, R){
  nbeta <- dim(Covariates)[2]
  score<-IL_score_noerror(weights, nodes, Y1, Y2, Covariates, R,  
                   Theta[1:nbeta],Theta[(nbeta+1):(2*nbeta)],Theta[2*nbeta+1],Theta[2*nbeta+2])
  return(score)
}

infomat_GQ<-function(Theta, weights, nodes, Y1, Y2, Covariates, R){
  nbeta <- dim(Covariates)[2]
  score<-IL_infomat_noerror(weights, nodes, Y1, Y2, Covariates, R,  
                            Theta[1:nbeta],Theta[(nbeta+1):(2*nbeta)],Theta[2*nbeta+1],Theta[2*nbeta+2])
  return(score)
}

GEE_UI <- function(Theta, Y1star, Y2star, Covariates){
  # cat(theta, " \n")
  nbeta <- dim(Covariates)[2]
  return(GEE_UfuncIns(Y1star, Y2star, DesignMatrix1=as.matrix(Covariates), 
                      DesignMatrix2=as.matrix(Covariates), 
                      CovMis1 = matrix(rep(0,dim(Covariates)[1]*2),ncol=2), 
                      CovMis2 = as.matrix(rep(1,dim(Covariates)[1])),
                      beta1 = Theta[1:nbeta], 
                      beta2 = Theta[(nbeta+1):(2*nbeta)], 
                      sigma = Theta[2*nbeta+1], xi = Theta[2*nbeta+2], 
                      gamma1=1, gamma=c(0,0), alpha1=-Inf, alpha0=-Inf,
                      sigma_e=0))
}

GEE_SIGMA <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, 
                      gamma1, gamma, alpha1, alpha0, sigma_e){
  nbeta <- dim(DesignMatrix1)[2]
  return(GEE_SIGMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                      theta[1:nbeta], theta[(nbeta+1):(2*nbeta)], 
                      sigma = theta[2*nbeta+1], xi = theta[2*nbeta+2], 
                      gamma1, gamma, alpha1, alpha0, sigma_e))
}

GEE_GAMMA <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2){
  nbeta <- dim(DesignMatrix1)[2]
  GAMMA <- GEE_GAMMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, beta1=theta[1:nbeta], 
                        beta2=theta[(nbeta+1):(2*nbeta)],
                        sigma = theta[2*nbeta+1], xi = theta[2*nbeta+2])
  return(GAMMA)
}

GEE_GAMMA.inv <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2){
  nbeta <- dim(DesignMatrix1)[2]
  GAMMA <- GEE_GAMMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, beta1=theta[1:nbeta], 
                        beta2=theta[(nbeta+1):(2*nbeta)],
                        sigma = theta[2*nbeta+1], xi = theta[2*nbeta+2])
  GAMMA.inv <- solve(GAMMA)
  return(GAMMA.inv)
}

GEE_cov <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                    gamma1, gamma, alpha1, alpha0, sigma_e){
  GAMMA.inv <- GEE_GAMMA.inv(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
  SIGMA <- GEE_SIGMA(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                     gamma1, gamma, alpha1, alpha0, sigma_e)
  covmatrix <- GAMMA.inv %*% SIGMA %*% t(as.matrix(GAMMA.inv))
  return(covmatrix)
}



ROC_curve_GNMM <- function (path, theta) 
{ ROC = list()
  theta = as.matrix(theta)
  d = ncol(theta)
  diag(theta) <- 0
  pos.total = sum(theta != 0)
  neg.total = d * (d - 1) - pos.total
  
  ROC$tp = rep(0, length(path))
  ROC$fp = rep(0, length(path))
  ROC$F1 = rep(0, length(path))
  for (r in 1:length(path)) {
    tmp = as.matrix(path[[r]])
    tp.all = (theta != 0) * (tmp != 0)
    diag(tp.all) = 0
    ROC$tp[r] <- sum(tp.all != 0)/pos.total
    fp.all = (theta == 0) * (tmp != 0)
    diag(fp.all) = 0
    ROC$fp[r] <- sum(fp.all != 0)/neg.total
    fn = 1 - ROC$tp[r]
    precision = ROC$tp[r]/(ROC$tp[r] + ROC$fp[r])
    recall = ROC$tp[r]/(ROC$tp[r] + fn)
    ROC$F1[r] = 2 * precision * recall/(precision + recall)
    if (is.na(ROC$F1[r])) 
      ROC$F1[r] = 0
  }
  # rm(precision, recall, tp.all, fp.all, path, theta, fn)
  ord.fp = order(ROC$fp)
  tmp1 = ROC$fp[ord.fp]
  tmp2 = ROC$tp[ord.fp]
  ROC$AUC = sum(diff(tmp1) * (tmp2[-1] + tmp2[-length(tmp2)]))/2
  return(list(fp = tmp1, tp = tmp2, AUC = ROC$AUC, F1= ROC$F1))
}

# LASSO <- LASSO1$beta[7:dim(cordcomp)[1],]
# theta <- cordcomp[7:dim(cordcomp)[1],3]

ROC_curve_LASSO <- function (LASSO, theta) 
{ ROC = list()
  pos.total = sum(theta != 0)
  neg.total = sum(theta == 0)
  npath <- dim(LASSO)[2]
  ROC$tp = rep(0, npath)
  ROC$fp = rep(0, npath)
  ROC$F1 = rep(0, npath)
  for (r in 1:npath) {
     tmp = LASSO[,r]
     tp.all = (theta != 0) * (tmp != 0)
     ROC$tp[r] <- sum(tp.all != 0)/pos.total
     fp.all = (theta == 0) * (tmp != 0)
     ROC$fp[r] <- sum(fp.all != 0)/neg.total
     fn = 1 - ROC$tp[r]
     precision = ROC$tp[r]/(ROC$tp[r] + ROC$fp[r])
     recall = ROC$tp[r]/(ROC$tp[r] + fn)
     ROC$F1[r] = 2 * precision * recall/(precision + recall)
     if (is.na(ROC$F1[r])) 
        ROC$F1[r] = 0
  }
  # rm(precision, recall, tp.all, fp.all, path, theta, fn)
  ord.fp = order(ROC$fp)
  tmp1 = ROC$fp[ord.fp]
  tmp2 = ROC$tp[ord.fp]
  ROC$AUC = sum(diff(tmp1) * (tmp2[-1] + tmp2[-length(tmp2)]))/2
  return(list(fp = tmp1, tp = tmp2, AUC = ROC$AUC, F1= ROC$F1))
}

# graphhub <- huge.generator(n=1000,d=6,graph ="hub",g=2,vis=T)
# graphcluster <- huge.generator(n=1000,d=6,graph ="cluster",g=2,vis=T)
# graphscalefree <- huge.generator(n=1000,d=6,graph ="scale-free",vis=T)
# graphhub$theta
# graphcluster$theta
# graphscalefree$theta
thetas <- 0.2
betas <- 0.5


#### 2. Implementation Function #####

### For dubuging puprose (comment when not debugging)

i <- 1
nsample <-300
graphtype <- "hub"

### Generate true beta:
## true parameters (assume no intercept)
nbeta <- 15+6
#(main para)   (inter.) 


### Function:

SIM2_main <- function(i,thetas, betas, nsample, graphtype){

#### 1. Set up ####
### 1.1 Global parameters  ####
  set.seed(2019)
  seed_i <- sample(1000000,1000)
  set.seed(seed_i[i])
  
  
### 1.2 R packages ####
  library(parallel)
  library(scales)
  library(huge)
  library(MASS)
  library(glmnet)
  library(GeneErrorMis)
  library(nleqslv)
  library(xtable)
  

#### 2.1 Data Generation ####
### 2.1.1 Specification of correlation matrix ####

## hub plot
Theta_hab <- thetas * (2*rbinom(5,1,0.5)-1)

THETA_hub <-matrix(c(1, Theta_hab[1], Theta_hab[2], Theta_hab[5], 0, 0,
              Theta_hab[1], 1, 0, 0, 0, 0,
              Theta_hab[2], 0, 1, 0, 0, 0,
              Theta_hab[5], 0, 0, 1, Theta_hab[3], Theta_hab[4],
              0, 0, 0, Theta_hab[3], 1, 0,
              0, 0, 0, Theta_hab[4], 0, 1),ncol=6)

## scale-free plot
Theta_scfr <- thetas * (2*rbinom(5,1,0.5)-1)

THETA_scfr <-matrix(c(1, Theta_scfr[1], Theta_scfr[2], Theta_scfr[3], 0, Theta_scfr[4],
                      Theta_scfr[1], 1, 0, 0, 0, 0,
                      Theta_scfr[2], 0, 1, 0, 0, 0,
                      Theta_scfr[3], 0, 0, 1, 0, 0,
                      0, 0, 0, 0, 1, Theta_scfr[5],
                      Theta_scfr[4], 0, 0, 0, Theta_scfr[5], 1),ncol=6)

## block plot
Theta_bloc <- thetas * (2*rbinom(3,1,0.5)-1)

THETA_bloc <-matrix(c(1, Theta_bloc[1], Theta_bloc[2], 0, 0, 0,
                      Theta_bloc[1], 1, Theta_bloc[3], 0, 0, 0,
                      Theta_bloc[2], Theta_bloc[3], 1, 0, 0, 0,
                      0, 0, 0, 1, 0, 0,
                      0, 0, 0, 0, 1, 0,
                      0, 0, 0, 0, 0, 1),ncol=6)

if (graphtype=="hub") {
  THETAgraph <- THETA_hub
 }else {
   if (graphtype=="bloc") {
     THETAgraph <- THETA_bloc
    }else{
      THETAgraph <- THETA_scfr
   }
}



Xcov <- mvrnorm(n=nsample, mu=rep(0,6), Sigma=solve(THETAgraph)) 
#### 2.2 Step 1: variable selection ####
### 2.2.1 Parameters

### 2.2.2 Selection Algorithm
 # use invisible() to suppress the function message
invisible(capture.output(varsel <- huge(Xcov, lambda = NULL, nlambda = 30, lambda.min.ratio = NULL, method = "mb",
                                        scr = F, scr.num = NULL, sym = "or", verbose = TRUE,cov.output =T)))

panelty <-  seq(from = max(varsel$lambda)*3, to = 0, length.out = 30)
panelty2 <-  seq(from = max(varsel$lambda)/3, to = 0, length.out = 30)

#### 2.3 Step 2: regression analysis ####
## 2.3.1 Data Generation ####

EdgeTrue <- NULL
EdgeHash <- rep(T,6)

for (i in 1:5){
  for (j in (i+1):6){
    if (THETAgraph[i,j]!=0){
      EdgeTrue <- rbind(EdgeTrue,c(i,j))
      EdgeHash <- c(EdgeHash,T)
    } else { EdgeHash <- c(EdgeHash,F)}
  }
}

DesMatrix <- Xcov

for (i in 1:dim(EdgeTrue)[1]){
  DesMatrix <- cbind(DesMatrix,Xcov[,EdgeTrue[i,1]]*Xcov[,EdgeTrue[i,2]])
}


cordinates_true <- matrix(c(rep(0,6),1:6),ncol=2)
cordinates_true <- rbind(cordinates_true, EdgeTrue)


## true parameters (assume no intercept)
nbeta <- dim(cordinates_true)[1]
#(main para)   (inter.) 
set.seed(2019)
beta1 <- rnorm(nbeta, betas, 0.3) * (2*rbinom(nbeta,1,0.5)-1)
beta2 <- rnorm(nbeta, betas, 0.3) * (2*rbinom(nbeta,1,0.5)-1)

sigma <- 1

cordtrue <- data.frame(cordinates_true, beta1 = beta1, beta2 = beta2)

### Generate the true data sets

mu1 <- DesMatrix %*% t(t(beta1)) 
mu2 <- DesMatrix %*% t(t(beta2))


## Response
epsilon <- rnorm(nsample,0,1)
U <- runif(nsample,0,1)
mu2expit <- expit(mu2)

Y1 <- mu1 +  epsilon
Y2 <- ifelse(U < mu2expit,1,0)

## measurement error and misclassification
e <- rnorm(nsample,0,sigma_e)
U2 <- runif(nsample,0,1)

Y1star <- Y1 + gammas * Y2 + e
Y2star <- ifelse(U2>expit(alphas),Y2,1-Y2)

# ## Naive model
# naive.model1 <- lm(Y1star ~ X + W)
# true.model1 <- lm(Y1 ~ X + W)
# naive.model2 <- glm(Y2star ~ X + W, family = binomial(link = logit))
# true.model2 <- glm(Y2star ~ X + W, family = binomial(link = logit))

CovMis1 <- cbind(rep(0,length(Y1star)),rep(1,length(Y1star)))
CovMis2 <- c(rep(1,length(Y1star)))

CovMis1 <- as.matrix(CovMis1)
CovMis2 <- as.matrix(CovMis2)




tryCatch({
## 2.3.2 Algorithm Implementation -- Proposed Approach ####

## Create the mismeasured data and the validation data
## We start with a choice from the previous step ``varsel$path[[15]]"
graph <- varsel$path[[15]]
cordinates <- matrix(c(rep(0,6),1:6),ncol=2)

for (i in 1:5){
  for (j in (i+1):6){
    if (graph[i,j]==1){
      cordinates <- rbind(cordinates,c(i,j))
    }
  }
}

DesMatrixhat <- Xcov

for (i in 1:dim(cordinates)[1]){
  DesMatrixhat <- cbind(DesMatrixhat,Xcov[,cordinates[i,1]]*Xcov[,cordinates[i,2]])
}

theta0 <- c(rep(0,(dim(cordinates)[1])*2),1,0)
# theta1 <- c(beta1,beta2,1,1)

## 2.3.2.1 Simulation - point estimation ####

NR <- nleqslv(theta0, GEE_UI, Y1star=Y1, Y2star=Y2, Covariates=DesMatrixhat,
              jacobian=T, control=list(maxit=2000))

betahat <- ifelse(abs(NR$x)<10,NR$x,NA)



## 2.3.2.2 Results Rearrangement ####

thetaresults <- betahat
betaresults <- data.frame(cordinates, beta1 = thetaresults[1:dim(cordinates)[1]], 
                          beta2 = thetaresults[(dim(cordinates)[1]+1):(2*dim(cordinates)[1])])


final_GNNM <- merge(cordtrue, betaresults,by = c("X1","X2"), all.x = T, all.y = F)

## 2.3.2.3 Variance Estimationo ####

betaI <- c(final_GNNM$beta1.y,final_GNNM$beta2.y,betahat[(length(betahat)-1):length(betahat)])
betaI <- ifelse(is.na(betaI),0,betaI)


if (!any(is.na(betahat))) {
  cov <- GEE_cov(betaI,Y1star = Y1, Y2star = Y2, 
                 DesignMatrix1=as.matrix(DesMatrix),
                 DesignMatrix2 = as.matrix(DesMatrix), 
                 CovMis1 = matrix(rep(0,dim(DesMatrix)[1]*2),ncol=2), 
                 CovMis2 = as.matrix(rep(1,dim(DesMatrix)[1])),
                 gamma1 = 1, gamma=c(0,0), alpha1= -Inf, alpha0= -Inf, sigma_e = 0)
  betaIsd <- sqrt(diag(cov))} else {
    betaIsd <- rep(NA,length(betaI))
  }



## 2.3.3 LASSO ####
### 2.3.3.1 First Component ####

cordcomp <- matrix(c(rep(0,6),1:6,rep(1,6)),ncol=3)

for (i in 1:5){
  for (j in (i+1):6){
    cordcomp <- rbind(cordcomp,c(i,j,(THETAgraph[i,j]!=0)*1))
  }
}


DesMatrixComp <- Xcov

for (i in (dim(Xcov)[2]+1):(dim(cordcomp)[1])){
  DesMatrixComp <- cbind(DesMatrixComp,Xcov[,cordcomp[i,1]]*Xcov[,cordcomp[i,2]])
}


LASSO1 <- glmnet(x = DesMatrixComp, y = Y1, family="gaussian",  alpha = 1, lambda = panelty,
       intercept=FALSE) 

cv_LASSO1 <- cv.glmnet(x = DesMatrixComp, y = Y1, lambda = panelty, nfolds=10,intercept=FALSE)
diff1 <- abs(panelty - cv_LASSO1$lambda.1se)
LASSO1best <- which(diff1==min(diff1))

### 2.3.3.2 Second Component ####
LASSO2 <- glmnet(x = DesMatrixComp, y = Y2, family="binomial",  alpha = 1, lambda = panelty2,
                 intercept=FALSE) 

cv_LASSO2 <- cv.glmnet(x = DesMatrixComp, y = Y2, family = "binomial", type.measure = "class",
                       lambda = panelty2, nfolds=10,intercept=FALSE)
diff2 <- abs(panelty2 - cv_LASSO2$lambda.1se)
LASSO2best <- which(diff2==min(diff2))



LASSOresults <- data.frame(cordcomp, beta1 = LASSO1$beta[,LASSO1best] , beta2 = LASSO2$beta[,LASSO2best])
final_LASSO <- merge(cordtrue,LASSOresults,by = c("X1","X2"),all.x = T, all.y = F)


#### 2.3.4 Comparing the results of the first stage ####
## 2.3.4.1 ROC curve for GNMM ####
ROCGNMM <- ROC_curve_GNMM( varsel$path,  THETAgraph)


ROCGNMMLASSO1 <- ROC_curve_LASSO (LASSO1$beta[7:dim(cordcomp)[1],],cordcomp[7:dim(cordcomp)[1],3])
ROCGNMMLASSO2 <- ROC_curve_LASSO (LASSO2$beta[7:dim(cordcomp)[1],],cordcomp[7:dim(cordcomp)[1],3])




#### 2.3.5 Comparing the results of the second stage ####
Norm1GNNM <- (colSums(abs(final_GNNM[,5:6]-final_GNNM[,3:4]),na.rm=T))
Norm1LASSO <- (colSums(abs(final_LASSO[,6:7]-final_LASSO[,3:4]),na.rm=T))

Norm2GNNM <- sqrt(colSums((final_GNNM[,5:6]-final_GNNM[,3:4])^2,na.rm=T))
Norm2LASSO <- sqrt(colSums((final_LASSO[,6:7]-final_LASSO[,3:4])^2,na.rm=T))


return(list( betaI, betaIsd, final_GNNM=final_GNNM, final_LASSO = final_LASSO, ROCGNMM = ROCGNMM, ROCGNMMLASSO1 = ROCGNMMLASSO1, ROCGNMMLASSO2 = ROCGNMMLASSO2,
    Norm1GNNM=Norm1GNNM, Norm1LASSO=Norm1LASSO,
    Norm2GNNM=Norm2GNNM, Norm2LASSO=Norm2LASSO ))

}, error = function(e)   return(NULL) )
} 




#### 3. Implementation ####
### 3.1 Example ####
### block plot ####
SIM2_main(i = 1, thetas = 0.2, betas = 0.4, nsample = 500, graphtype="bloc")

### hub plot ####
SIM2_main(i = 1, thetas = 0.2, betas = 0.4, nsample = 500, graphtype="hub")

### scale-free plot ####
SIM2_main(i = 1, thetas = 0.2, betas = 0.4, nsample = 500, graphtype="scfr")

### 3.2 parallel computing ####
cl = makeCluster( 40 ,outfile=paste0(ProjectName,".txt")) 

clusterExport(cl=cl, varlist=c("seed_i", "expit","likelihood_GQ",
                               "score_GQ","infomat_GQ","ROC_curve_GNMM","ROC_curve_LASSO",
                               "GEE_UI","GEE_SIGMA","GEE_GAMMA","GEE_GAMMA.inv","GEE_cov"))

registerDoParallel(cl)



#### 4.  Simulation Studies ####

### 4.1 Study 1: Sample Size ####
SIM1  <- foreach(nsample=c(200,500)) %:% foreach(grty=c("bloc","hub","scfr")) %:% foreach(var1=1:1000) %dopar% {
  SIM2_main(i = var1, thetas = 0.2, betas = 0.4, nsample = nsample, graphtype = grty)
}

save(SIM1,file="output/SIMGEE1.RData")

SIM1[[1]][[1]][[1]]

ResTable <- NULL

graphtype <- c("block","hub","scale-free")
data.ROC <- NULL
nsample <- c(200,500)
GNMMAUClist <- NULL
LASSO1AUClist <- NULL
LASSO2AUClist <- NULL


for (k in 1:2){
  for (i in 1:3){
    GNMMtptab <- NULL
    GNMMfptab <- NULL
    GNMMAUCvec <- NULL
    LASSO1tptab <- NULL
    LASSO1fptab <- NULL
    LASSO1AUCvec <- NULL
    LASSO2tptab <- NULL
    LASSO2fptab <- NULL
    LASSO2AUCvec <- NULL
    for (j in 1:1000){
      GNMMtptab <- rbind(GNMMtptab,SIM1[[k]][[i]][[j]]$ROCGNMM$tp)
      GNMMfptab <- rbind(GNMMfptab,SIM1[[k]][[i]][[j]]$ROCGNMM$fp)
      GNMMAUCvec <- c(GNMMAUCvec,SIM1[[k]][[i]][[j]]$ROCGNMM$AUC)

      LASSO1tptab <- rbind(LASSO1tptab,SIM1[[k]][[i]][[j]]$ROCGNMMLASSO1$tp)
      LASSO1fptab <- rbind(LASSO1fptab,SIM1[[k]][[i]][[j]]$ROCGNMMLASSO1$fp)
      LASSO1AUCvec <- c(LASSO1AUCvec,SIM1[[k]][[i]][[j]]$ROCGNMMLASSO1$AUC)

      LASSO2tptab <- rbind(LASSO2tptab,SIM1[[k]][[i]][[j]]$ROCGNMMLASSO2$tp)
      LASSO2fptab <- rbind(LASSO2fptab,SIM1[[k]][[i]][[j]]$ROCGNMMLASSO2$fp)
      LASSO2AUCvec <- c(LASSO2AUCvec,SIM1[[k]][[i]][[j]]$ROCGNMMLASSO2$AUC)
    }
    GNMMtp <- c(0,colMeans(GNMMtptab,na.rm=T),1)
    GNMMfp <- c(0,colMeans(GNMMfptab,na.rm=T),1)
    GNMMAUC <- getAUC(GNMMfp,GNMMtp)

    LASSO1tp <- c(0,colMeans(LASSO1tptab,na.rm=T),1)
    LASSO1fp <- c(0,colMeans(LASSO1fptab,na.rm=T),1)
    LASSO1AUC <- getAUC(LASSO1fp,LASSO1tp)

    LASSO2tp <- c(0,colMeans(LASSO2tptab,na.rm=T),1)
    LASSO2fp <- c(0,colMeans(LASSO2fptab,na.rm=T),1)
    LASSO2AUC <- getAUC(LASSO2fp,LASSO2tp)

    method <- c(rep("GNMM",length(GNMMtp)),rep("LASSO1",length(LASSO1tp)),rep("LASSO2",length(LASSO2tp)))
    tps <- c(GNMMtp,LASSO1tp,LASSO2tp)
    fps <- c(GNMMfp,LASSO1fp,LASSO2fp)

    data.ROC.this <- data.frame(method=method, tps=tps, fps=fps, graphtype = graphtype[i], nsample=nsample[k])
    data.ROC <- rbind(data.ROC,data.ROC.this)

    GNMMAUClist <- c(GNMMAUClist, GNMMAUC)
    LASSO1AUClist <- c(LASSO1AUClist, LASSO1AUC)
    LASSO2AUClist <- c(LASSO2AUClist, LASSO2AUC)
  }
}


dat_text0 <- data.frame(
  label = unlist(lapply(GNMMAUClist, FUN = function(x){
    paste0("AUC")
  })),
  nsample   = rep(c(200,500), each=3),
  graphtype = rep(c("block","hub","scale-free"),times =2 ),
  fps = 0.75,
  tps = 0.35,
  stringsAsFactors = FALSE
)

dat_text1 <- data.frame(
  label1 = unlist(lapply(GNMMAUClist, FUN = function(x){
    paste0("GNM    ",sprintf("%.3f", round(x,digits =3)))
    # paste0("GNMM")
  })),
  nsample   = rep(c(200,500), each=3),
  graphtype = rep(c("block","hub","scale-free"),times =2 ),
  fps = 0.72,
  tps = 0.25,
  stringsAsFactors = FALSE
)

dat_text2 <- data.frame(
  label1 = unlist(lapply(LASSO1AUClist, FUN = function(x){
    paste0("LASSO1  ",sprintf("%.3f", round(x,digits =3)))
    # paste0("GNMM")
  })),
  nsample   = rep(c(200,500), each=3),
  graphtype = rep(c("block","hub","scale-free"),times =2 ),
  fps = 0.7,
  tps = 0.15,
  stringsAsFactors = FALSE
)

dat_text3 <- data.frame(
  label1 = unlist(lapply(LASSO2AUClist, FUN = function(x){
    paste0("LASSO2  ",sprintf("%.3f", round(x,digits =3)))
    # paste0("GNMM")
  })),
  nsample   = rep(c(200,500), each=3),
  graphtype = rep(c("block","hub","scale-free"),times =2 ),
  fps = 0.7,
  tps = 0.05,
  stringsAsFactors = FALSE
)

pdf(file = "output/fig1sample.pdf",height = 5.6, width = 10)
p1 <- ggplot(data.ROC, aes(x=fps, y=tps, group=method)) +
  geom_line(aes(linetype=method,color=method),size=1.1)+
  theme(text=element_text(size=13, family="mono"))+
  scale_color_brewer(palette="Dark2")+
  geom_point(aes(shape=method,color=method))+
  facet_grid(nsample~graphtype,switch ="y") +
  labs(x = "1 - Specificity", y="Sensitivity")


p1 +
  geom_text(
    family = "mono", fontface="bold", size =5,
    data    = dat_text0,
    mapping = aes(x = fps, y = tps, label = label),
    inherit.aes=F
 )+
  geom_text(
    family = "mono", size =3,
    data    = dat_text1,
    mapping = aes(x = fps, y = tps, label = label1),
    inherit.aes=F
  ) + geom_text(
    family = "mono", size =3,
    data    = dat_text2,
    mapping = aes(x = fps, y = tps, label = label1),
    inherit.aes=F
  ) + geom_text(
    family = "mono", size =3,
    data    = dat_text3,
    mapping = aes(x = fps, y = tps, label = label1),
    inherit.aes=F
  )

dev.off()


#### 4.1.2 Step 2 Evaluation ####

Results <- NULL

samplesize_range <- c(200,500)
graph_range <- c("block","hub","scale-free")

for (j in 1:3){
    for (k in 1:2) {
    results <- SIM1[[k]][[j]]
    
    Norm1GNNM <- NULL
    Norm1LASSO <- NULL
    Norm2GNNM <- NULL
    Norm2LASSO <- NULL
    
    for (i in 1:1000){
      Norm1GNNMi <- results[[i]]$Norm1GNNM
      Norm1LASSOi <- results[[i]]$Norm1LASSO
      Norm2GNNMi  <- results[[i]]$Norm2GNNM
      Norm2LASSOi <- results[[i]]$Norm2LASSO
      
      Norm1GNNM <- rbind(Norm1GNNM,Norm1GNNMi)
      Norm1LASSO <- rbind(Norm1LASSO,Norm1LASSOi)
      Norm2GNNM <- rbind(Norm2GNNM,Norm2GNNMi)
      Norm2LASSO <- rbind(Norm2LASSO,Norm2LASSOi)
    }
    
    Norm1GNNMv <- colMeans(na.omit(Norm1GNNM),na.rm = T)
    Norm1LASSOv <- colMeans(na.omit(Norm1LASSO),na.rm = T)
    Norm2GNNMv <- colMeans(na.omit(Norm2GNNM),na.rm = T)
    Norm2LASSOv <- colMeans(na.omit(Norm2LASSO),na.rm = T)
    
    Results0 <- data.frame(samplesize = samplesize_range[k],
                           graph = graph_range[j],
                           type= c("beta1","beta2"),
                           Norm1GNNM=round(Norm1GNNMv,3),Norm1LASSO=round(Norm1LASSOv,3),Norm2GNNM=round(Norm2GNNMv,3),Norm2LASSO=round(Norm2LASSOv,3))
    
    Results <- rbind(Results,Results0)
  }
}

Results_Long <- gather(Results, key = Method, value = Measurements, Norm1GNNM:Norm2LASSO)
Results_Long <- cbind(Results_Long, Norm = rep(c("Norm1","Norm2"), each = dim(Results_Long)[1]/2))

Results_rar <- Results_Long  %>% arrange((graph)) %>% arrange((samplesize)) %>% arrange((type)) 

# Results_c <- data.frame(Results_rar[Results_rar$Norm=="Norm1",], Norm2 = Results_rar[Results_rar$Norm=="Norm2","Measurements"])
Results_rar$Method <- ifelse (Results_rar$Method %in% c("Norm1GNNM","Norm2GNNM"),"GNNM","LASSO")


Results_Wide <- spread(Results_rar, key = Norm, value = Measurements)
Results_Wide2 <- cbind(Results_Wide[Results_Wide$type=="beta1",],
                       Results_Wide[Results_Wide$type=="beta2",c("Norm1","Norm2")])
Results_final <- Results_Wide2[,names(Results_Wide2)!="type"]

xtable(Results_final)



# Results_output <- cbind(Results[Results$type=="beta1",],Results[Results$type=="beta2",4:7])

# ### 4.1 Study 1: Sample Size ####
# SIM2  <- foreach(nsample=c(200,500)) %:% foreach(grty=c("bloc","hub","scfr")) %:% foreach(var1=1:1000) %dopar% {
#   SIM2_main(i = var1, thetas = 0.2, betas = 0.4, nsample = nsample, graphtype = grty)
# }
# 
# save(SIM2,file="output/SIMGEE2.RData")
# 
# 
# ### Summarize the performance of the second step
# Results <- NULL
# 
# samplesize_range <- c(200,1000)
# 
# granum <- 2
# for (k in 1:2) {
#   results <- SIM2[[k]][[granum]]
# 
#   betas <- NULL
#   betabias <- NULL
#   sds <- NULL
#   CIs <- NULL
#   for (i in 1:1000){
#     betaItrue <- c(results[[i]]$final_GNNM$beta1.x,results[[i]]$final_GNNM$beta2.x,1,0)
#     betahat0 <- results[[i]][[1]]
#     sd0 <- results[[i]][[2]]
#     # sd0 <- ifelse(abs(sd0)<100,sd0,NA)
#     betas <- rbind(betas,betahat0)
#     betabias <- rbind(betabias, betahat0-betaItrue)
#     sds <- rbind(sds, sd0)
# 
#     CILB <- betahat0 - 1.96 *(sd0)
#     CIUB <- betahat0 + 1.96 *(sd0)
#     CIs <- rbind(CIs,ifelse((betaItrue<as.vector(CIUB)) & (betaItrue>as.vector(CILB)),1,0))
#   }
# 
#   bias1 <- colMeans(na.omit(betabias),na.rm = T)
# 
#   sd_emp <- apply(betas,MARGIN = 2, FUN = sd, na.rm = T)
# 
#   # sd_emp <- apply(na.omit(betas),MARGIN = 2, FUN = function(x){
#   #   x.noout <-  remove_outliers(x)
#   #   return( sd(x.noout,na.rm = T ) )
#   # })
# 
#   sd_mod <- colMeans(na.omit(sds),na.rm = T)
# 
#   # sd_mod <- sd_mod <- apply(na.omit(sds),MARGIN = 2, FUN = function(x){
#   #   x.noout <-  remove_outliers(x)
#   #   return( mean(x.noout,na.rm = T ))
#   # })
#   CIrate <- colMeans(na.omit(CIs),na.rm = T)
# 
#   Results0 <- data.frame(samplesize = samplesize_range[k],
#                          biasprop=round(bias1,3),propose_esd=round(sd_emp,3),sdpropose=round(sd_mod,3),CI_propose=percent(round(CIrate,3)))
# 
#   Results <- rbind(Results,Results0)
# }

  


### 4.2 Study 2: Sample Size ####
# SIM2  <- foreach(thetas=c(0.2,0.4)) %:% foreach(betas=c(0.2,0.5,0.9)) %:% foreach(var1=1:1000) %dopar% {
#   SIM2_main(i = var1, thetas = thetas, betas = betas, nsample = 200, graphtype = "scale-free")
# }
# 
# save(SIM2,file="output/SIMGEE2.RData")
# 
# SIM2[[1]][[1]][[1]]
# 
# ResTable <- NULL
# 
# thetas=c(0.2,0.4)
# betas=c(0.2,0.5,0.9)
# 
# data.ROC <- NULL
# nsample <- c(200,500)
# GNMMAUClist <- NULL
# LASSO1AUClist <- NULL
# LASSO2AUClist <- NULL
# 
# 
# for (k in 1:2){
#   for (i in 1:3){
#     GNMMtptab <- NULL
#     GNMMfptab <- NULL
#     GNMMAUCvec <- NULL
#     LASSO1tptab <- NULL
#     LASSO1fptab <- NULL
#     LASSO1AUCvec <- NULL
#     LASSO2tptab <- NULL
#     LASSO2fptab <- NULL
#     LASSO2AUCvec <- NULL
#     for (j in 1:1000){
#       GNMMtptab <- rbind(GNMMtptab,SIM2[[k]][[i]][[j]]$ROCGNMM$tp)
#       GNMMfptab <- rbind(GNMMfptab,SIM2[[k]][[i]][[j]]$ROCGNMM$fp)
#       GNMMAUCvec <- c(GNMMAUCvec,SIM2[[k]][[i]][[j]]$ROCGNMM$AUC)
#       
#       LASSO1tptab <- rbind(LASSO1tptab,SIM2[[k]][[i]][[j]]$ROCGNMMLASSO1$tp)
#       LASSO1fptab <- rbind(LASSO1fptab,SIM2[[k]][[i]][[j]]$ROCGNMMLASSO1$fp)
#       LASSO1AUCvec <- c(LASSO1AUCvec,SIM2[[k]][[i]][[j]]$ROCGNMMLASSO1$AUC)
#       
#       LASSO2tptab <- rbind(LASSO2tptab,SIM2[[k]][[i]][[j]]$ROCGNMMLASSO2$tp)
#       LASSO2fptab <- rbind(LASSO2fptab,SIM2[[k]][[i]][[j]]$ROCGNMMLASSO2$fp)
#       LASSO2AUCvec <- c(LASSO2AUCvec,SIM2[[k]][[i]][[j]]$ROCGNMMLASSO2$AUC)
#     }
#     GNMMtp <- c(0,colMeans(GNMMtptab,na.rm=T),1)
#     GNMMfp <- c(0,colMeans(GNMMfptab,na.rm=T),1)
#     GNMMAUC <- getAUC(GNMMfp,GNMMtp)
#     
#     LASSO1tp <- c(0,colMeans(LASSO1tptab,na.rm=T),1)
#     LASSO1fp <- c(0,colMeans(LASSO1fptab,na.rm=T),1)
#     LASSO1AUC <- getAUC(LASSO1fp,LASSO1tp)
#     
#     LASSO2tp <- c(0,colMeans(LASSO2tptab,na.rm=T),1)
#     LASSO2fp <- c(0,colMeans(LASSO2fptab,na.rm=T),1)
#     LASSO2AUC <- getAUC(LASSO2fp,LASSO2tp)
#     
#     method <- c(rep("GNMM",length(GNMMtp)),rep("LASSO1",length(LASSO1tp)),rep("LASSO2",length(LASSO2tp)))
#     tps <- c(GNMMtp,LASSO1tp,LASSO2tp)
#     fps <- c(GNMMfp,LASSO1fp,LASSO2fp)
#     
#     data.ROC.this <- data.frame(method=method, tps=tps, fps=fps, betas = betas[i],thetas=thetas[k])
#     data.ROC <- rbind(data.ROC,data.ROC.this)
#     
#     GNMMAUClist <- c(GNMMAUClist, GNMMAUC)
#     LASSO1AUClist <- c(LASSO1AUClist, LASSO1AUC)
#     LASSO2AUClist <- c(LASSO2AUClist, LASSO2AUC)
#   }
# }
# 
# 
# dat_text0 <- data.frame(
#   label = unlist(lapply(GNMMAUClist, FUN = function(x){
#     paste0("AUC")
#   })),
#   thetas   = rep(thetas, each=3),
#   betas = rep(betas,times =2 ),
#   fps = 0.75,
#   tps = 0.35,
#   stringsAsFactors = FALSE
# )
# 
# dat_text1 <- data.frame(
#   label1 = unlist(lapply(GNMMAUClist, FUN = function(x){
#     paste0("GNMM   ",sprintf("%.3f", round(x,digits =3)))
#     # paste0("GNMM")
#   })),
#   thetas   = rep(thetas, each=3),
#   betas = rep(betas,times =2 ),
#   fps = 0.72,
#   tps = 0.25,
#   stringsAsFactors = FALSE
# )
# 
# dat_text2 <- data.frame(
#   label1 = unlist(lapply(LASSO1AUClist, FUN = function(x){
#     paste0("LASSO1  ",sprintf("%.3f", round(x,digits =3)))
#     # paste0("GNMM")
#   })),
#   thetas   = rep(thetas, each=3),
#   betas = rep(betas,times =2 ),
#   fps = 0.7,
#   tps = 0.15,
#   stringsAsFactors = FALSE
# )
# 
# dat_text3 <- data.frame(
#   label1 = unlist(lapply(LASSO2AUClist, FUN = function(x){
#     paste0("LASSO2  ",sprintf("%.3f", round(x,digits =3)))
#     # paste0("GNMM")
#   })),
#   thetas   = rep(thetas, each=3),
#   betas = rep(betas,times =2 ),
#   fps = 0.7,
#   tps = 0.05,
#   stringsAsFactors = FALSE
# )
# 
# pdf(file = "output/fig2strength.pdf",height = 5.6, width = 10)
# p1 <- ggplot(data.ROC, aes(x=fps, y=tps, group=method)) +
#   geom_line(aes(linetype=method,color=method),size=1.1)+
#   theme(text=element_text(size=13, family="mono"))+
#   scale_color_brewer(palette="Dark2")+
#   geom_point(aes(shape=method,color=method))+
#   facet_grid(thetas~betas,switch ="y") +
#   labs(x = "1 - Specificity", y="Sensitivity")
#  
# 
# p1 +
#   geom_text(
#     family = "mono", fontface="bold", size =5,
#     data    = dat_text0,
#     mapping = aes(x = fps, y = tps, label = label),
#     inherit.aes=F
#   )+
#   geom_text(
#     family = "mono", size =3,
#     data    = dat_text1,
#     mapping = aes(x = fps, y = tps, label = label1),
#     inherit.aes=F
#   ) + geom_text(
#     family = "mono", size =3,
#     data    = dat_text2,
#     mapping = aes(x = fps, y = tps, label = label1),
#     inherit.aes=F
#   ) + geom_text(
#     family = "mono", size =3,
#     data    = dat_text3,
#     mapping = aes(x = fps, y = tps, label = label1),
#     inherit.aes=F
#   )
# 
# 
# dev.off() 