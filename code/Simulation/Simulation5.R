#### Simulation 5 #####

# This simulation evaluate the performance of the method WITH 
# measurement error in responses. Different from the Simulation 1, this simulation
# evaluate the case using GEE method rather than likelihood approach.
# And different from simulation 3, this study consider the case with 
# EXTERNAL VALIDATION data available

## Update Mar 17, 2021
# Change covariate in the misclassification process to be part of covariate

## Update Feb 7th 
#  Add covariates in the misclassification process


#### 1. Set up ####
### 1.1 Global parameters  ####
set.seed(2020)
seed_i <- sample(1000000,1000)
ProjectName <- paste0("Simulation5")


ncore <- 40

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

GEE_UI <- function(Theta, Y1star, Y2star, Covariates){
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

GEE_UI_ErrMis <- function(Theta, Y1star, Y2star, Covariates, CovMis1, CovMis2,
                   gamma1, gamma, alpha1, alpha0, sigma_e){
  nbeta <- dim(Covariates)[2]
  return(GEE_UfuncIns(Y1star, Y2star, DesignMatrix1=as.matrix(Covariates), 
                      DesignMatrix2=as.matrix(Covariates),  CovMis1, CovMis2,
                      beta1 = Theta[1:nbeta], 
                      beta2 = Theta[(nbeta+1):(2*nbeta)], 
                      sigma = Theta[2*nbeta+1], xi = Theta[2*nbeta+2], 
                      gamma1, gamma, alpha1=alpha1, alpha0=alpha0, sigma_e))
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

## 1.3.2 Functions with Internal Validation  ####
GEE_UI_EV <- function(theta, ncov, data.mismeasure,
                      gamma1, gamma, alpha1, alpha0, sigma_e){
  return(GEE_UfuncIns(Y1star=data.mismeasure$Y1star,
                      Y2star=data.mismeasure$Y2star,
                      DesignMatrix1 = as.matrix(data.mismeasure[,3:(2+ncov)]),
                      DesignMatrix2 = as.matrix(data.mismeasure[,3:(2+ncov)]), 
                      CovMis1 = as.matrix(data.mismeasure[,(3+ncov):(4+ncov)]),
                      CovMis2 = as.matrix(data.mismeasure[,(5+ncov):(6+ncov)]),
                      beta1=theta[1:ncov], beta2=theta[(ncov+1):(ncov*2)], sigma = theta[ncov*2+1], xi = theta[ncov*2+2], 
                      gamma1=gamma1, gamma=gamma, alpha1=alpha1, alpha0=alpha0, sigma_e=sigma_e))
}

GEE_GAMMA_EV0 <- function(theta, ncov, Y1star, Y2star, Y1, Y2, CovMis1, CovMis2, 
                          gamma1, gamma, alpha1, alpha0, sigma_e,
                          fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  return(GEE_GAMMAInsEV0(Y1star, Y2star, Y1, Y2,
                         CovMis1, CovMis2, ncov1=ncov, ncov2=ncov,
                         gamma1, gamma, alpha1, alpha0, sigma_e, 
                         fixgamma1=fixgamma1, fixgamma=fixgamma, fixsigma_e=fixsigma_e, fixalpha1=fixalpha1, fixalpha0=fixalpha0)
  )
}


GEE_GAMMA_EVI <- function(theta, ncov, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, 
                          gamma1, gamma, alpha1, alpha0, sigma_e,
                          fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  return(GEE_GAMMAInsEVI(Y1star, Y2star, DesignMatrix1, DesignMatrix2,
                         CovMis1, CovMis2,
                         beta1=theta[1:ncov], beta2=theta[(ncov+1):(ncov*2)], xi=theta[ncov*2+2], sigma=theta[ncov*2+1],
                         gamma1, gamma, alpha1, alpha0, sigma_e,
                         fixgamma1=fixgamma1, fixgamma=fixgamma, fixsigma_e=fixsigma_e, fixalpha1=fixalpha1, fixalpha0=fixalpha0)
  )
}

GEE_SIGMA_EV0 <- function(theta, ncov, Y1star, Y2star, Y1, Y2, CovMis1, CovMis2, 
                          gamma1, gamma, alpha1, alpha0, sigma_e,
                          fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  return(GEE_SIGMAInsEV0(Y1star, Y2star, Y1, Y2,
                         CovMis1, CovMis2, ncov1 =ncov, ncov2=ncov,
                         gamma1, gamma, alpha1, alpha0, sigma_e,
                         fixgamma1=fixgamma1, fixgamma=fixgamma, fixsigma_e=fixsigma_e, fixalpha1=fixalpha1, fixalpha0=fixalpha0)
  )
}

GEE_SIGMA_EVI <- function(theta, ncov, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, 
                          gamma1, gamma, alpha1, alpha0, sigma_e,
                          fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  return(GEE_SIGMAInsEVI(Y1star, Y2star, DesignMatrix1, DesignMatrix2,
                         CovMis1, CovMis2,
                         beta1=theta[1:ncov], beta2=theta[(ncov+1):(ncov*2)], xi=theta[ncov*2+2], sigma=theta[ncov*2+1],
                         gamma1, gamma, alpha1, alpha0, sigma_e,
                         fixgamma1=fixgamma1, fixgamma=fixgamma, fixsigma_e=fixsigma_e, fixalpha1=fixalpha1, fixalpha0=fixalpha0)
  )
}

GEE_covEV <- function(theta, ncov, data.validation, data.mismeasure,
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  nvalidation <- dim(data.validation)[1]
  nsample <- dim(data.mismeasure)[1] + nvalidation
  
  M0 <- GEE_GAMMA_EV0(theta, ncov, 
                      Y1star=data.validation$Y1star, 
                      Y2star=data.validation$Y2star, 
                      Y1 = data.validation$Y1, 
                      Y2 = data.validation$Y2, 
                      CovMis1 = as.matrix(data.validation[,5:6]), 
                      CovMis2 = as.matrix(data.validation[,7:8]), 
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0)
  
  M1 <- GEE_GAMMA_EVI(theta, ncov, 
                      Y1star=data.mismeasure$Y1star, 
                      Y2star=data.mismeasure$Y2star, 
                      DesignMatrix1 = as.matrix(data.mismeasure[,3:(2+ncov)]),
                      DesignMatrix2 = as.matrix(data.mismeasure[,3:(2+ncov)]), 
                      CovMis1 = as.matrix(data.mismeasure[,(3+ncov):(4+ncov)]), 
                      CovMis2 = as.matrix(data.mismeasure[,(5+ncov):(6+ncov)]), 
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0)
  
  GAMMA_EV <- M1 + M0
  
  B0 <- GEE_SIGMA_EV0(theta, ncov, 
                      Y1star=data.validation$Y1star, 
                      Y2star=data.validation$Y2star, 
                      Y1 = data.validation$Y1, 
                      Y2 = data.validation$Y2,
                      CovMis1 = as.matrix(data.validation[,5:6]), 
                      CovMis2 = as.matrix(data.validation[,7:8]), 
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0)
  
  B1 <- GEE_SIGMA_EVI(theta, ncov, 
                      Y1star=data.mismeasure$Y1star, 
                      Y2star=data.mismeasure$Y2star, 
                      DesignMatrix1 = as.matrix(data.mismeasure[,3:(2+ncov)]),
                      DesignMatrix2 = as.matrix(data.mismeasure[,3:(2+ncov)]), 
                      CovMis1 = as.matrix(data.mismeasure[,(3+ncov):(4+ncov)]), 
                      CovMis2 = as.matrix(data.mismeasure[,(5+ncov):(6+ncov)]), 
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0)
  
  SIGMA_EV <- B1 + B0
  
  GAMMA.inv <- solve(GAMMA_EV,tol=1e-200)
  covmatrix <- GAMMA.inv %*% SIGMA_EV %*% t(as.matrix(GAMMA.inv))
  return(covmatrix)
}





## 1.3.3 Functions with ROC curve  ####

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
  
  ord.fp = order(ROC$fp)
  tmp1 = ROC$fp[ord.fp]
  tmp2 = ROC$tp[ord.fp]
  ROC$AUC = sum(diff(tmp1) * (tmp2[-1] + tmp2[-length(tmp2)]))/2
  return(list(fp = tmp1, tp = tmp2, AUC = ROC$AUC, F1= ROC$F1))
}


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

  ord.fp = order(ROC$fp)
  tmp1 = ROC$fp[ord.fp]
  tmp2 = ROC$tp[ord.fp]
  ROC$AUC = sum(diff(tmp1) * (tmp2[-1] + tmp2[-length(tmp2)]))/2
  return(list(fp = tmp1, tp = tmp2, AUC = ROC$AUC, F1= ROC$F1))
}

thetas <- 0.2
betas <- 0.5


#### 2. Implementation Function #####

### For dubuging puprose (comment when not debugging)

i <- 1
nsample <-1500
nvalidation <- 500
graphtype <- "scale-free"
sigma_e <- 0.1
alphas <-  c(-3,1)
gammas <- 0.5

### Generate true beta:
## true parameters (assume no intercept)
nbeta <- 15+6
#(main para)   (inter.) 
set.seed(2019)
beta1pool <- runif(nbeta, 0.1, 0.7) * (2*rbinom(nbeta,1,0.5)-1)
beta2pool <- runif(nbeta, 0.1, 0.7) * (2*rbinom(nbeta,1,0.5)-1)


### Function:


SIM5_main <- function(i,thetas, betas, nsample, nvalidation, graphtype, sigma_e, alphas, gammas){
  
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
  
  
  Xcov <- mvrnorm(n=nsample, mu=rep(0,6), Sigma=THETAgraph) 
  
  #### 2.2 Step 1: variable selection ####
  ### 2.2.1 Parameters
  
  ### 2.2.2 Selection Algorithm
  # use invisible() to suppress the function message
  invisible(capture.output(varsel <- huge(Xcov, lambda = NULL, nlambda = 30, lambda.min.ratio = NULL, method = "mb",
                                          scr = F, scr.num = NULL, sym = "or", verbose = TRUE, cov.output =T)))
  
  panelty <-  seq(from = max(varsel$lambda)*3, to = 0, length.out = 30)
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
  
  
  sigma <- 1
  
  
  
  ### Generate the true data sets
  
  beta1 <- beta1pool[EdgeHash]
  beta2 <- beta2pool[EdgeHash]
  
  cordtrue <- data.frame(cordinates_true, beta1 = beta1, beta2 = beta2)
  
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
  
  Z <- Xcov[,1]
  CovMis2 <- data.frame(intercept= rep(1,nsample), Z=Z)
  CovMis2 <- as.matrix(CovMis2)
  
  
  # alphas <- c(-4,-1)   #1%
  # alphas <- c(-3,0)    #5%
  # alphas <- c(-2.8,-1)    #10%
  alphastcov <- as.matrix(CovMis2) %*% t(t(alphas))
  # mean(expit(alphastcov))
  
  
  Y1star <- Y1 + gammas * Y2 + e
  
  Y2star <- ifelse(U2>expit(alphastcov),Y2,1-Y2)
  
  CovMis1 <- cbind(rep(0,length(Y1star)), rep(1,length(Y1star)))
  
  CovMis1 <- as.matrix(CovMis1)
  
  
  
  tryCatch({
    ## 2.3.2 Algorithm Implementation -- Naive Approach ####
    
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
    
    
    ## Create the mismeasured data and the validation data
    
    data.mismeasure <- data.frame(Y1star=Y1star[1:(nsample - nvalidation)],Y2star=Y2star[1:(nsample - nvalidation)], DesMatrixhat[1:(nsample - nvalidation),],CovMis1[1:(nsample - nvalidation),],CovMis2[1:(nsample - nvalidation),]) 
    data.validation <- data.frame(Y1=Y1[(nsample - nvalidation+1):nsample],Y2=Y2[(nsample - nvalidation+1):nsample], Y1star=Y1star[(nsample - nvalidation+1):nsample],Y2star=Y2star[(nsample - nvalidation+1):nsample],
                                  CovMis1[(nsample - nvalidation+1):nsample,],CovMis2[(nsample - nvalidation+1):nsample,]) 
    
    ## 2.3.2.1 Simulation - point estimation ####
    
    NR_Naive <- nleqslv(theta0, GEE_UI, Y1star=Y1star, Y2star=Y2star, Covariates=DesMatrixhat,
                        jacobian=T, control=list(maxit=2000))
    
    betahat_Naive <- ifelse(abs(NR_Naive$x)<10,NR_Naive$x,NA)
    
    
    
    ## 2.3.2.2 Results Rearrangement ####
    
    thetaresults <- betahat_Naive
    
    betaresults <- data.frame(cordinates, beta1 = thetaresults[1:dim(cordinates)[1]], 
                              beta2 = thetaresults[(dim(cordinates)[1]+1):(2*dim(cordinates)[1])])
    
    
    naive_GNNM <- merge(cordtrue, betaresults,by = c("X1","X2"), all.x = T, all.y = F)
    
    ## 2.3.2.3 Variance Estimation ####
    
    betaI_naive <- c(naive_GNNM$beta1.y,naive_GNNM$beta2.y,
                     betahat_Naive[(length(betahat_Naive)-1):length(betahat_Naive)])
    betaI_naive <- ifelse(is.na(betaI_naive),0,betaI_naive)
    
    
    if (!any(is.na(betahat_Naive))) {
      cov_naive <- GEE_cov(betaI_naive,Y1star = Y1star, Y2star = Y2star, 
                           DesignMatrix1=as.matrix(DesMatrix),
                           DesignMatrix2 = as.matrix(DesMatrix), 
                           CovMis1 = matrix(rep(0,dim(DesMatrix)[1]*2),ncol=2), 
                           CovMis2 = as.matrix(rep(1,dim(DesMatrix)[1])),
                           gamma1 = 1, gamma=c(0,0), alpha1= -Inf, alpha0= -Inf, sigma_e = 0)
      betaIsd_naive <- sqrt(diag(cov_naive))} else {
        betaIsd_naive <- rep(NA,length(betaI_naive))
      }
    s
    
    ## 2.3.3 Algorithm Implementation -- Proposed Approach ####s
    
    ## 2.3.3.1 Measurement Error and Misclassification Parameters
    model.measure <- lm(Y1star ~ -1 + offset(Y1) + Y2,data = data.validation) 
    model.class1 <- glm((1-Y2star) ~ Z, data = data.validation[data.validation$Y2==1,],family = binomial(link="logit")) 
    model.class0 <- glm(Y2star ~ Z, data = data.validation[data.validation$Y2==0,],family = binomial(link="logit")) 
    
    gamma2 <- model.measure$coefficients
    sigma_e <- sigma(model.measure)
    alpha1 <- model.class1$coefficients
    alpha0 <- model.class0$coefficients
    
    
    theta0 <- c(rep(0,(dim(cordinates)[1])*2),1,0)
    
    ## 2.4.2.1 Simulation - point estimation ####
    
    NR <- nleqslv(theta0, GEE_UI_EV, jacobian=T, control=list(maxit=2000),
                  data.mismeasure = data.mismeasure, 
                  ncov= dim(cordinates)[1], 
                  gamma1 = 1, gamma=c(0,gamma2), alpha1= alpha1, alpha0= alpha0, sigma_e = sigma_e)
    
    betahat <- ifelse(abs(NR$x)<10,NR$x,NA)
    
    ## 2.3.2.2 Results Rearrangement ####
    
    thetaresults <- betahat
    
    betaresults <- data.frame(cordinates, beta1 = thetaresults[1:dim(cordinates)[1]], 
                              beta2 = thetaresults[(dim(cordinates)[1]+1):(2*dim(cordinates)[1])])
    
    
    proposed_GNNM <- merge(cordtrue, betaresults,by = c("X1","X2"), all.x = T, all.y = F)
    
    ## 2.3.2.3 Variance Estimation ####
    
    betaI_proposed <- c(proposed_GNNM$beta1.y,proposed_GNNM$beta2.y,betahat[(length(betahat)-1):length(betahat)])
    
    betaI_proposed <- ifelse(is.na(betaI_proposed),0,betaI_proposed)
    
    ## Create the mismeasured data and the validation data under true dataset
    
    dataE.mismeasure <- data.frame(Y1star=Y1star[1:(nsample - nvalidation)],Y2star=Y2star[1:(nsample - nvalidation)], DesMatrix[1:(nsample - nvalidation),],CovMis1[1:(nsample - nvalidation),],CovMis2[1:(nsample - nvalidation),]) 
    
    dataE.validation <- data.frame(Y1=Y1[(nsample - nvalidation+1):nsample],Y2=Y2[(nsample - nvalidation+1):nsample], Y1star=Y1star[(nsample - nvalidation+1):nsample],Y2star=Y2star[(nsample - nvalidation+1):nsample],
                                   CovMis1[(nsample - nvalidation+1):nsample,],CovMis2[(nsample - nvalidation+1):nsample,]) 
    
    
    
    ### variance estimation with validation data
    if (!any(is.na(betahat))) {
      cov <- GEE_covEV (betaI_proposed, ncov= length(betaI_proposed)/2-1,  dataE.validation, dataE.mismeasure, 
                        gamma1=1, gamma = c(0,gamma2), alpha1, alpha0, sigma_e,
                        fixgamma1=1, fixgamma=c(1,0), fixsigma_e=0, fixalpha1=c(0,0), fixalpha0=c(0,0))
      betaIsd_proposed <- sqrt(diag(cov))} else {
        betaIsd_proposed <- rep(NA,length(betahat))
      }
    
    
    betatrue <- c(cordtrue$beta1, cordtrue$beta2, 1, 0)
    
    return(list( betaI_naive=betaI_naive, 
                 betaIsd_naive=betaIsd_naive, 
                 naive_GNNM=naive_GNNM, 
                 betaI_proposed=betaI_proposed, 
                 betaIsd_proposed=betaIsd_proposed, 
                 proposed_GNNM=proposed_GNNM,
                 betatrue=betatrue))
    
  }, error = function(e) return(NULL))
} 




#### 3. Implementation ####
### 3.1 Example ####
### block plot ####
SIM5_main(i = 1, thetas = 0.2, betas = 0.4, nsample = 1500, nvalidation=500,
          graphtype="bloc", sigma_e = 0.1, alphas = c(-3,1), gammas = 0.5)

### hub plot ####
SIM5_main(i = 1, thetas = 0.2, betas = 0.4, nsample = 1500, nvalidation=500, 
          graphtype="hub", sigma_e = 0.1, alphas = c(-3,1), gammas = 0.5)

### scale-free plot ####
SIM5_main(i = 1, thetas = 0.2, betas = 0.4, nsample = 1500, nvalidation=500, 
          graphtype="scfr", sigma_e = 0.1, alphas = c(-3,1), gammas = 0.5)

### 3.2 parallel computing ####
cl = makeCluster( ncore ,outfile=paste0(ProjectName,".txt"))

clusterExport(cl=cl, varlist=c("seed_i", "expit","likelihood_GQ", "GEE_UI_ErrMis",
                               "score_GQ","infomat_GQ","ROC_curve_GNMM","ROC_curve_LASSO",
                               "GEE_UI","GEE_SIGMA","GEE_GAMMA","GEE_GAMMA.inv","GEE_cov",
                               "GEE_UI_EV","GEE_GAMMA_EV0","GEE_GAMMA_EVI","GEE_SIGMA_EV0",
                               "GEE_SIGMA_EVI", "GEE_covEV"))

registerDoParallel(cl)


#### 4.  Simulation Studies ####
# alphas <- c(-5,-1)   #1%
# alphas <- c(-3,0)    #5%
# alphas <- c(c(-2.8,0.5))    #10%
alphasset <- c(-5,-1,-3,0,-2.8,0.5)
alphassetM <-matrix(alphasset,nrow=2)

### 4.1 Study 1: Sample Size ####
SIM5  <- foreach(sige=c(0.2,0.7)) %:% foreach(alas=1:3) %:% foreach(grty=c("bloc","hub","scfr")) %:% foreach(var1=1:1000) %dopar% {
  SIM5_main(i = var1, thetas = 0.2, betas = 0.4, nsample = 1500, nvalidation=500,
            graphtype = grty,sigma_e = sige, alphas = alphassetM[,alas], gammas = 0.5)
}

load(file="output/SIMGEE5.RData")

# Example
SIM5[[2]][[3]][[3]][[1]]

### Summarize the performance of the second step
Results <- NULL

samplesize_range <- c(1000)
sigmae_range <- c(0.2,0.7)
alpha_range <- c(-4.595,-2.944,-2.197)
graph_range <- c("bloc","hub","scfr")

for (i in 1:2) {
  for (j in 1:3){
    for (k in 1:3){
      results <- SIM5[[i]][[j]][[k]]

      betas_naive <- NULL
      betabias_naive <- NULL
      sds_naive <- NULL
      CIs_naive <- NULL

      betas <- NULL
      betabias <- NULL
      sds <- NULL
      CIs <- NULL
      for (num in 1:1000){
        if (is.null(results[[num]])) {next}
        betaItrue <- results[[num]]$betatrue

        betahat0_naive <- results[[num]][[1]]
        sd0_naive <- results[[num]][[2]]
        betas_naive <- rbind(betas_naive,betahat0_naive)
        betabias_naive <- rbind(betabias_naive, (betahat0_naive-betaItrue))
        sds_naive <- rbind(sds_naive, sd0_naive)

        CILB_naive <- betahat0_naive - 1.96 *(sd0_naive)
        CIUB_naive <- betahat0_naive + 1.96 *(sd0_naive)
        CIs_naive <- rbind(CIs_naive,ifelse((betaItrue<as.vector(CIUB_naive)) & (betaItrue>as.vector(CILB_naive)),1,0))


        betahat0_props <- results[[num]][[4]]
        sd0_props <- results[[num]][[5]]
        sd0_props <- ifelse(abs(sd0_props)<3,sd0_props,NA)
        betas <- rbind(betas,betahat0_props)
        betabias <- rbind(betabias, (betahat0_props-betaItrue))
        sds <- rbind(sds, sd0_props)

        CILB_props <- betahat0_props - 1.96 *(sd0_props)
        CIUB_props <- betahat0_props + 1.96 *(sd0_props)
        CIs <- rbind(CIs,ifelse((betaItrue<as.vector(CIUB_props)) & (betaItrue>as.vector(CILB_props)),1,0))
      }

      biasnaive <- colMeans(na.omit(betabias_naive),na.rm = T)


      sd_empnaive <- apply(na.omit(betas_naive),MARGIN = 2, FUN = function(x){
        x.noout <-  remove_outliers(x)
        return( sd(x.noout,na.rm = T ) )
      })

      sd_modnaive <- colMeans(na.omit(sds_naive),na.rm = T)

      CIrate_naive <- colMeans(na.omit(CIs_naive),na.rm = T)

      bias1 <- colMeans(na.omit(betabias),na.rm = T)

      sd_emp <- apply(na.omit(betas),MARGIN = 2, FUN = function(x){
        x.noout <-  remove_outliers(x)
        return( sd(x.noout,na.rm = T ) )
      })

      sd_mod <- colMeans(na.omit(sds),na.rm = T)

      
      CIrate <- colMeans(na.omit(CIs),na.rm = T)
      
      ### Only focus on the report for beta
      biasnaive <- biasnaive[1:(length(biasnaive)-2)]
      sd_modnaive <- sd_modnaive[1:(length(sd_modnaive)-2)]
      sd_empnaive <- sd_empnaive[1:(length(sd_empnaive)-2)]
      CIrate_naive <- CIrate_naive[1:(length(CIrate_naive)-2)]
      
      bias1 <- bias1[1:(length(bias1)-2)]
      sd_mod <- sd_mod[1:(length(sd_mod)-8)]
      sd_emp <- sd_emp[1:(length(sd_emp)-2)]
      CIrate <- CIrate[1:(length(CIrate)-8)]
      

      Results0 <- data.frame(sigmae = sigmae_range[i],
                             alpha = alpha_range[j],
                             graph = graph_range[k],
                             biasnaive=round(abs(biasnaive),3),sdempnaive=round(sd_empnaive,3),sdmodnaive=round(sd_modnaive,3),CI_naive=(round(CIrate_naive,3)),
                             biasprop=round(abs(bias1),3),sdempprop=round(sd_emp,3),sdmodprop=round(sd_mod,3),CI_prop=(round(CIrate,3)))
      Results0 <- cbind(rep(c("beta1","beta2"),each=dim(Results0)[1]/2),Results0)
      colnames(Results0)[1] <- "type"
      Results <- rbind(Results,Results0)
    }
  }
}

results_agg <- aggregate(Results[,5:12], by=list(Results$type,Results$sigmae,Results$alpha,Results$graph), FUN=mean,na.rm=T)

results_agg <- results_agg[order(results_agg$Group.1),]

colnames(results_agg) <- c("Type","Sigma_e","alpha","graph",
                           "biasnaive", "sdempnaive", "sdmodnaive", "CI_naive",
                           "biasprop", "sdempprop", "sdmodprop", "CI_prop")

results_agg2 <- results_agg

results_agg2[,1:4] <- results_agg2[,c(1,4,3,2)]

colnames(results_agg2) <- c("Type","Sigma_e","alpha","graph",
                            "biasnaive", "sdempnaive", "sdmodnaive", "CI_naive",
                            "biasprop", "sdempprop", "sdmodprop", "CI_prop")

results_agg2[,8]<-percent(results_agg2[,8], accuracy = .1)
results_agg2[,12]<-percent(results_agg2[,12], accuracy = .1)

xtable(results_agg2,digits = 3)
save(results_agg2,file="results/Table2/Table2_EXV.RData")

