#### Simulation 1 #####

# This simulation evaluate the performance of the method without considering
# measurement error in covariates.



#### 1. Set up ####
### 1.1 Global parameters  ####
set.seed(2019)
seed_i <- sample(1000000,1000)



### 1.2 R packages ####
library(parallel)
library(scales)
library(huge)
library(MASS)
library(glmnet)
library(GeneErrorMis)
library(ggplot2)
library(pROC)


### 1.3 Global Functions ####
expit <- function(x){
  value <- exp(x)/(1+exp(x))
  ifelse(is.na(value),1,value) 
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
SIM1_main <- function(i,thetas, betas, nsample, THETA){

#### 2. Data Generation ####
### 2.1 Specification of correlation matrix ####

## hub plot
Theta_hab <- thetas * (2*rbinom(4,1,0.5)-1)

THETA_hub <-matrix(c(1, Theta_hab[1], Theta_hab[2], 0, 0, 0,
              Theta_hab[1], 1, 0, 0, 0, 0,
              Theta_hab[2], 0, 1, 0, 0, 0,
              0, 0, 0, 1, Theta_hab[3], Theta_hab[4],
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




nsample <- 500
Xcov <- mvrnorm(n=nsample,mu=rep(0,6),Sigma=THETA_bloc) 
#### 2. Step 1: variable selection ####
### 2.1 Parameters

### 2.2 Selection Algorithm

varsel <- huge(Xcov, lambda = NULL, nlambda = 30, lambda.min.ratio = NULL, method = "mb",
        scr = F, scr.num = NULL, sym = "or", verbose = TRUE,cov.output =T)

panelty <-  seq(from = max(varsel$lambda)*3, to = min(varsel$lambda), length.out = 30)

# varsel_opt <- huge.select(varsel, criterion = "stars", stars.thresh = 0.05, 
#                                         stars.subsample.ratio = NULL, rep.num = 50, verbose = TRUE)
# varsel_opt <- huge.select(varsel, criterion = "ric", stars.thresh = 0.1,
#             stars.subsample.ratio = NULL, rep.num = 100, verbose = TRUE)
# varsel_opt <- huge.select(varsel, criterion = "stars", stars.thresh = 0.5,
#                           stars.subsample.ratio = NULL, rep.num = 20, verbose = TRUE)


# huge.roc(varsel_opt$path,L$theta)



#### 3. Step 2: regression analysis ####
## 3.1 Data Generation ####
Edge_bloc <- rbind(c(1,2),c(1,3),c(2,3))

DesMatrix <- Xcov

for (i in 1:dim(Edge_bloc)[1]){
  DesMatrix <- cbind(DesMatrix,Xcov[,Edge_bloc[i,1]]*Xcov[,Edge_bloc[i,2]])
}


FF <- c(rep(1,250),rep(2,250))


cordinates_true <- matrix(c(rep(0,6),1:6),ncol=2)
cordinates_true <- rbind(cordinates_true, Edge_bloc)

## true parameters (assume no intercept)
nbeta <- dim(cordinates_true)[1]
#(main para)   (inter.) 
beta1 <- rnorm(nbeta,0.4,0.3) * (2*rbinom(nbeta,1,0.5)-1)
beta2 <- rnorm(nbeta,0.4,0.3) * (2*rbinom(nbeta,1,0.5)-1)

sigma <- 1
sigma_c <- 1

EF <- rnorm(nsample,0,sigma_c)

cordtrue <- data.frame(cordinates_true, beta1 = beta1, beta2 = beta2)

### Generate the true data sets

mu1 <- DesMatrix %*% t(t(beta1)) + FF * EF
mu2 <- DesMatrix %*% t(t(beta2)) + FF * EF


## Response
epsilon <- rnorm(nsample,0,1)
U <- runif(nsample,0,1)
mu2expit <- expit(mu2)

Y1 <- mu1 +  epsilon
Y2 <- ifelse(U < mu2expit,1,0)



## 3.2 Algorithm Implementation -- Proposed Approach ####

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

theta0 <- c(rep(0,(dim(cordinates)[1])*2),1,1)
theta1 <- c(beta1,beta2,1,1)

## 3.2.1 Hermite Quadracture Weights ####
library(statmod)
quadrature <- gauss.quad(20,kind="hermite")
nodes <- quadrature$nodes
weights <- quadrature$weights

## 3.2.2 Simulation ####
R <- c(rep(1,500),rep(4,500))

Results_IL_ini1<-tryCatch({optim(theta0,likelihood_GQ, gr = score_GQ,
                                 method="BFGS",control=list(fnscale=-1),hessian = TRUE,
                                 weights = weights, nodes = nodes, Y1 = Y1, R = R,
                                 Y2 = Y2, Covariates = DesMatrixhat
)}, error = function(e)
{return(list(par=NULL,value=NULL,vcovinitial1=NULL))}
)

if (!is.null(Results_IL_ini1$par)){
  infomatini1<-infomat_GQ(Results_IL_ini1$par, weights = weights, nodes = nodes, Y1 = Y1,
                                Y2 = Y2, Covariates = DesMatrixhat, R = R)
  vcovinitial1 <- solve(infomatini1,tol=1e-100)
  # beta <- Results_IL_ini1$par[1:6]
  # sd <- sqrt(diag(vcovinitial1[1:6,1:6]))
  # testsig<-ifelse(abs(beta/sd)>1.96,1,0)
} else{
  infomatini1<- NULL
  # vcovinitial1 <- Results_IL_ini1$vcovinitial1
}


## 3.2.3 Results Rearrangement ####

thetaresults <- Results_IL_ini1$par
betaresults <- data.frame(cordinates, beta1 = thetaresults[1:dim(cordinates)[1]], 
                          beta2 = thetaresults[(dim(cordinates)[1]+1):(2*dim(cordinates)[1])])


final_GNNM <- merge(cordtrue,betaresults,by = c("X1","X2"),all.x = T, all.y = F)


## 3.3 LASSO ####
### 3.3.1 First Component ####

cordcomp <- matrix(c(rep(0,6),1:6,rep(1,6)),ncol=3)

for (i in 1:5){
  for (j in (i+1):6){
    cordcomp <- rbind(cordcomp,c(i,j,(THETA_bloc[i,j]!=0)*1))
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

### 3.3.2 Second Component ####
LASSO2 <- glmnet(x = DesMatrixComp, y = Y2, family="binomial",  alpha = 1, lambda = panelty,
                 intercept=FALSE) 

cv_LASSO2 <- cv.glmnet(x = DesMatrixComp, y = Y2, family = "binomial", type.measure = "class",
                       lambda = panelty, nfolds=10,intercept=FALSE)
diff2 <- abs(panelty - cv_LASSO2$lambda.1se)
LASSO2best <- which(diff2==min(diff2))



LASSOresults <- data.frame(cordcomp, beta1 = LASSO1$beta[,LASSO1best] , beta2 = LASSO2$beta[,LASSO2best])
final_LASSO <- merge(cordtrue,LASSOresults,by = c("X1","X2"),all.x = T, all.y = F)


#### 3.4 Comparing the results of the first stage ####
## 3.4.1 ROC curve for GNMM ####
ROCGNMM <- ROC_curve_GNMM( varsel$path,  THETA_bloc)

LASSO <- LASSO1$beta[7:dim(cordcomp)[1],]
theta <- cordcomp[7:dim(cordcomp)[1],3]

          
ROCGNMMLASSO1 <- ROC_curve_LASSO (LASSO1$beta[7:dim(cordcomp)[1],],cordcomp[7:dim(cordcomp)[1],3])
ROCGNMMLASSO2 <- ROC_curve_LASSO (LASSO2$beta[7:dim(cordcomp)[1],],cordcomp[7:dim(cordcomp)[1],3])




#### 3.5 Comparing the results of the second stage ####
Norm1GNNM <- colSums(abs(final_GNNM[,5:6]-final_GNNM[,3:4]))
Norm1LASSO <- colSums(abs(final_LASSO[,6:7]-final_LASSO[,3:4]))


return(list( ROCGNMM = ROCGNMM, ROCGNMMLASSO1 = ROCGNMMLASSO1,ROCGNMMLASSO2 = ROCGNMMLASSO2,
    Norm1GNNM=Norm1GNNM, Norm1LASSO=Norm1LASSO))
} 
