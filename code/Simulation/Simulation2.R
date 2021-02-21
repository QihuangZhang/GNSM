#### Simulation 2 #####

# This simulation evaluate the performance of the method without considering
# measurement error in covariates. Different from the Simulation 1, this simulation
# evaluate the case using GEE method rather than likelihood approach.

## Update: Feb 8, 2021
# 1. the Theta is the precision matrix rather than covariance matrix
# 2. make the parameter beta fixed in each simulation


## Update: Jan 15, 2020
# 1. the Theta is the precision matrix rather than covariance matrix
# 2. make the parameter beta fixed in each simulation


## Update: Feb 5, 2020
# 1. Remove the LASSO comparison and compare the method by sample size


#### 1. Set up ####
### 1.1 Global parameters  ####
set.seed(2019)
seed_i <- sample(1000000,1000)
ProjectName <- paste0("Simulation2")
ncore <- 40

## Parameter for parallel purpose
TargetRegion <- c(1,1000) ## 1-200 201-400 401-600 601-800 800-1000
scanregion<- paste0(TargetRegion[1],"-",TargetRegion[2])
ProjectName <- paste0(ProjectName, scanregion)


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
library(patchwork)


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


ROC_curve_RetNet <- function (TP, FP, theta) 
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

SIM2_main <- function(i,thetas, nsample, graphtype){
  
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
  Theta_hab <- runif(5,0, 2*thetas) * (2*rbinom(5,1,0.5)-1)
  
  THETA_hub <-matrix(c(1, Theta_hab[1], Theta_hab[2], Theta_hab[5], 0, 0,
                       Theta_hab[1], 1, 0, 0, 0, 0,
                       Theta_hab[2], 0, 1, 0, 0, 0,
                       Theta_hab[5], 0, 0, 1, Theta_hab[3], Theta_hab[4],
                       0, 0, 0, Theta_hab[3], 1, 0,
                       0, 0, 0, Theta_hab[4], 0, 1),ncol=6)
  
  ## scale-free plot
  Theta_scfr <- runif(5,0, 2*thetas) * (2*rbinom(5,1,0.5)-1)
  
  THETA_scfr <-matrix(c(1, Theta_scfr[1], Theta_scfr[2], Theta_scfr[3], 0, Theta_scfr[4],
                        Theta_scfr[1], 1, 0, 0, 0, 0,
                        Theta_scfr[2], 0, 1, 0, 0, 0,
                        Theta_scfr[3], 0, 0, 1, 0, 0,
                        0, 0, 0, 0, 1, Theta_scfr[5],
                        Theta_scfr[4], 0, 0, 0, Theta_scfr[5], 1),ncol=6)
  
  ## block plot
  Theta_bloc <- runif(3,0, 2*thetas) * (2*rbinom(3,1,0.5)-1)
  
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
  cutpointrho <-  seq(from = 0, to = 1, length.out = 30)
  cutpointrho <- cutpointrho^2
  
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
  # beta1 <- rnorm(nbeta, betas, 0.3) * (2*rbinom(nbeta,1,0.5)-1)
  # beta2 <- rnorm(nbeta, betas, 0.3) * (2*rbinom(nbeta,1,0.5)-1)
  beta1 <- runif(nbeta, 0.1, 0.7) * (2*rbinom(nbeta,1,0.5)-1)
  beta2 <- runif(nbeta, 0.1, 0.7) * (2*rbinom(nbeta,1,0.5)-1)
  
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
    
    
    # if (!any(is.na(betahat))) {
    #   cov <- GEE_cov(betaI,Y1star = Y1, Y2star = Y2, 
    #                  DesignMatrix1=as.matrix(DesMatrix),
    #                  DesignMatrix2 = as.matrix(DesMatrix), 
    #                  CovMis1 = matrix(rep(0,dim(DesMatrix)[1]*2),ncol=2), 
    #                  CovMis2 = as.matrix(rep(1,dim(DesMatrix)[1])),
    #                  gamma1 = 1, gamma=c(0,0), alpha1= -Inf, alpha0= -Inf, sigma_e = 0)
    #   betaIsd <- sqrt(diag(cov))} else {
    #     betaIsd <- rep(NA,length(betaI))
    #   }
    
    
    
    ## 2.3.3 RelNet ####
    ### 2.3.3.1 First Component ####
    
    Cov_Matrix <- cov(Xcov)
    
    results_RelNet <- lapply(cutpointrho, FUN= function(rho) {
      keep <- abs(Cov_Matrix) > rho
      
      ## Get True Positive and True Negative for ploting ROC
      
      TP <- sum((THETAgraph!=0) * keep )
      FP <- sum((THETAgraph==0) * keep )
      
      TN <- sum((THETAgraph==0) * (1-keep) )
      FN <- sum((THETAgraph!=0) * (1-keep) )
       
      return(list(RelNetTP=TP/(TP+FN), RelNetFP=FP/(FP+TN)))
    })

    
    keep <- abs(Cov_Matrix) > 0.5
    
    cordcomp <- matrix(c(rep(0,6),1:6),ncol=2)
    
    for (i in 1:5){
      for (j in (i+1):6){
        if (keep[i,j]){
          cordcomp <- rbind(cordcomp,c(i,j))
        }
      }
    }
    
    DesMatrixComp <- Xcov
    
    if ((dim(Xcov)[2]+1) < (dim(cordcomp)[1])) {
      for (i in (dim(Xcov)[2]+1):(dim(cordcomp)[1])){
        DesMatrixComp <- cbind(DesMatrixComp,Xcov[,cordcomp[i,1]]*Xcov[,cordcomp[i,2]])
      }
    }
    
    
    theta0_RelNet <- c(rep(0,(dim(cordcomp)[1])*2),1,0)
    # theta1 <- c(beta1,beta2,1,1)
    
    ## 2.3.2.1 Simulation - point estimation ####
    
    NR_RelNet <- nleqslv(theta0_RelNet, GEE_UI, Y1star=Y1, Y2star=Y2, Covariates=DesMatrixComp,
                  jacobian=T, control=list(maxit=2000))
    
    betahat_RelNet <- ifelse(abs(NR_RelNet$x)<10,NR_RelNet$x,NA)
    
    
    
    ## 2.3.2.2 Results Rearrangement ####
    
    thetaresults_RelNet <- betahat_RelNet
    betaresults_RelNet <- data.frame(cordcomp, beta1 = thetaresults_RelNet[1:dim(cordcomp)[1]], 
                              beta2 = thetaresults_RelNet[(dim(cordcomp)[1]+1):(2*dim(cordcomp)[1])])
    
    
    final_RelNet <- merge(cordtrue, betaresults_RelNet,by = c("X1","X2"), all.x = T, all.y = F)
    
    ## 2.3.2.3 Variance Estimationo ####
    
    betaI_RelNet <- c(final_RelNet$beta1.y,final_RelNet$beta2.y,betahat_RelNet[(length(betahat_RelNet)-1):length(betahat_RelNet)])
    betaI_RelNet <- ifelse(is.na(betaI),0,betaI)
    
    
    # if (!any(is.na(betaI_RelNet))) {
    #   cov_RelNet <- GEE_cov(betaI_RelNet,Y1star = Y1, Y2star = Y2, 
    #                  DesignMatrix1=as.matrix(DesMatrixComp),
    #                  DesignMatrix2 = as.matrix(DesMatrixComp), 
    #                  CovMis1 = matrix(rep(0,dim(DesMatrixComp)[1]*2),ncol=2), 
    #                  CovMis2 = as.matrix(rep(1,dim(DesMatrixComp)[1])),
    #                  gamma1 = 1, gamma=c(0,0), alpha1= -Inf, alpha0= -Inf, sigma_e = 0)
    #   betaIsd_RelNet <- sqrt(diag(cov_RelNet))} else {
    #     betaIsd_RelNet <- rep(NA,length(theta0_RelNet))
    #   }
    
     
    
    #### 2.3.4 Comparing the results of the first stage ####
    ## 2.3.4.1 ROC curve for GNMM ####
    ROCGNMM <- ROC_curve_GNMM( varsel$path,  THETAgraph)
    
    RelNetTP <- unlist(lapply(results_RelNet, `[[`, 1))
    RelNetFP <- unlist(lapply(results_RelNet, `[[`, 2))
    
    
    #### 2.3.5 Comparing the results of the second stage ####
    Norm1GNNM <- (colSums(abs(final_GNNM[,5:6]-final_GNNM[,3:4]),na.rm=T))
    Norm1RelNet <- (colSums(abs(final_RelNet[,5:6]-final_RelNet[,3:4]),na.rm=T))
    
    Norm2GNNM <- sqrt(colSums((final_GNNM[,5:6]-final_GNNM[,3:4])^2,na.rm=T))
    Norm2RelNet <- sqrt(colSums((final_RelNet[,5:6]-final_RelNet[,3:4])^2,na.rm=T))
    
    
    return(list( betaI, #betaIsd, 
                 final_GNNM=final_GNNM,   
                 Norm1GNNM=Norm1GNNM, Norm2GNNM=Norm2GNNM,
                 final_RelNet = final_RelNet, 
                 Norm1RelNet=Norm1RelNet, Norm2RelNet=Norm2RelNet,
                 ROCGNMM = ROCGNMM,
                 RelNetTP = RelNetTP,
                 RelNetFP = RelNetFP
                 #,path= varsel$path
                 ))
    
  }, error = function(e)   return(NULL) )
} 



SIM2_high <- function(i,thetas, nsample){
  
  #### 2.4 Set up ####
  ### 2.4.1 Global parameters  ####
  set.seed(2019)
  seed_i <- sample(1000000,1000)
  set.seed(seed_i[i])
  ndim <- 500 ## p=1000 by combining both outcome and the interaction terms
  
  
  ### 2.4.2 R packages ####
  library(parallel)
  library(scales)
  library(huge)
  library(MASS)
  library(glmnet)
  library(GeneErrorMis)
  library(nleqslv)
  library(xtable)
  library(Matrix)
  
  tryCatch({
  ## hub plot
  Theta_g <- runif(ndim^2,0, 4*thetas) * (rbinom(ndim^2,1,0.05)-rbinom(ndim^2,1,0.05))

  THETAgraph_high <- matrix(Theta_g, ncol = ndim)
  Higher.Trian <-  THETAgraph_high

  ## scale-free plot
  # Theta_scfr <- runif(5,0, 2*thetas) * (2*rbinom(5,1,0.5)-1)
  # 
  # THETA_scfr <-matrix(c(1, Theta_scfr[1], Theta_scfr[2], Theta_scfr[3], 0, Theta_scfr[4], 
  #                       Theta_scfr[1], 1, 0, 0, 0, 0,
  #                       Theta_scfr[2], 0, 1, 0, 0, 0,
  #                       Theta_scfr[3], 0, 0, 1, 0, 0,
  #                       0, 0, 0, 0, 1, Theta_scfr[5],
  #                       Theta_scfr[4], 0, 0, 0, Theta_scfr[5], 1),ncol=6)
  
  # if (graphtype=="hub") {
  #   THETAgraph <- THETA_hub
  # }else {
  #   if (graphtype=="bloc") {
  #     THETAgraph <- THETA_bloc
  #   }else{
  #     THETAgraph <- THETA_scfr
  #   }
  # }
  # 
  
  # THETAgraph_high1 <- cbind(THETA_hub,matrix(0,nrow=6,ncol=ndim-6))
  # THETAgraph_high2 <- cbind(matrix(0,nrow=6,ncol=6),THETA_scfr,matrix(0,nrow=6,ncol=ndim-12))
  # THETAgraph_high <- rbind (THETAgraph_high1, THETAgraph_high2, matrix(0,nrow=ndim-12,ncol=ndim))  
  
  Higher.Trian[lower.tri(THETAgraph_high, diag = F)] <- 0
  diag(Higher.Trian) <- rep(1,ndim)
  #  <- t(Higher.Trian)
  
  THETAgraph_high<- Higher.Trian %*%  t(Higher.Trian)

  # diag(THETAgraph_high) <- rep(2,ndim)
  # THETAgraph_high <- THETAgraph_high/2
  
  Xcov <- mvrnorm(n=nsample, mu=rep(0,ndim), Sigma=solve(THETAgraph_high)) 
  #### 2.4.3.2 Step 1: variable selection ####
  ### Parameters
  
  ### Selection Algorithm
  # use invisible() to suppress the function message
  lambdas <- seq(from = 0.84, to = 0, length.out = 30)
  
  invisible(capture.output(varsel <- huge(Xcov, lambda = lambdas^4, nlambda = 30, lambda.min.ratio = NULL, method = "mb",
                                          scr = F, scr.num = NULL, sym = "or", verbose = TRUE,cov.output =T)))
  
  panelty <-  seq(from = max(varsel$lambda)*30, to = 0, length.out = 30)
  cutpointrho <-  (seq(from = 0, to = 1, length.out = 30))
  cutpointrho <- cutpointrho^2
  
  #### 2.4.3.3 Step 2: regression analysis ####
  ##  Data Generation ####
  
  ### Generate the true data sets
  
  # mu1 <- DesMatrix %*% t(t(beta1)) 
  # mu2 <- DesMatrix %*% t(t(beta2))
  # 
  # 
  # ## Response
  # epsilon <- rnorm(nsample,0,1)
  # U <- runif(nsample,0,1)
  # mu2expit <- expit(mu2)
  # 
  # Y1 <- mu1 +  epsilon
  # Y2 <- ifelse(U < mu2expit,1,0)
  # 
    ## 2.4.3.4 Algorithm Implementation -- Proposed Approach ####
    
    ## Create the mismeasured data and the validation data
    ## We start with a choice from the previous step ``varsel$path[[15]]"
    # graph <- varsel$path[[15]]
    # cordinates <- matrix(c(rep(0,ndim),1:ndim),ncol=2)
    # 
    # for (i in 1:(ndim-1)){
    #   for (j in (i+1):ndim){
    #     if (graph[i,j]==1){
    #       cordinates <- rbind(cordinates,c(i,j))
    #     }
    #   }
    # }
    # 
    # DesMatrixhat <- Xcov
    # 
    # for (i in 1:dim(cordinates)[1]){
    #   DesMatrixhat <- cbind(DesMatrixhat,Xcov[,cordinates[i,1]]*Xcov[,cordinates[i,2]])
    # }
    # 
    
    # ## Naive models
    # model1 <- lm(Y1~-1+DesMatrixhat)
    # model2 <- glm(Y2~-1+DesMatrixhat, family = binomial(link="logit"))
    # 
    # coef1 <- coef(model1) + rnorm(length(coef(model1)),0,0.1)
    # coef2 <- coef(model2) + rnorm(length(coef(model2)),0,1)
    # # theta0 <- c(rep(0,(dim(cordinates)[1])*2),1,0)
    # # theta0 <- c(coef1,rep(0,dim(cordinates)[1]),1,0)
    # theta0 <- c(coef1,coef2/10,1,0)
    # # theta0 <- c(coef1,beta2,1,0)
    
    # ## 2.4.3.5 Simulation - point estimation ####
    # 
    # NR <- nleqslv(theta0, GEE_UI, Y1star=Y1, Y2star=Y2, Covariates=DesMatrixhat,
    #               jacobian=T, control=list(maxit=50))
    # 
    # betahat <- ifelse(abs(NR$x)<10,NR$x,NA)
    # 
    # 
    # 
    # ## 2.4.3.6 Results Rearrangement ####
    # 
    # thetaresults <- betahat
    # betaresults <- data.frame(cordinates, beta1 = thetaresults[1:dim(cordinates)[1]], 
    #                           beta2 = thetaresults[(dim(cordinates)[1]+1):(2*dim(cordinates)[1])])
    # 
    # 
    # final_GNNM <- merge(cordtrue, betaresults,by = c("X1","X2"), all.x = T, all.y = F)
    # 
    # ## 2.4.3.7 Variance Estimation ####
    # cat("Variance Estimation - GNSM /n")
    # 
    # betaI <- c(final_GNNM$beta1.y,final_GNNM$beta2.y,betahat[(length(betahat)-1):length(betahat)])
    # betaI <- ifelse(is.na(betaI),0,betaI)
    # 
    # 
    # if (!any(is.na(betahat))) {
    #   cov <- GEE_cov(betaI,Y1star = Y1, Y2star = Y2, 
    #                  DesignMatrix1=as.matrix(DesMatrix),
    #                  DesignMatrix2 = as.matrix(DesMatrix), 
    #                  CovMis1 = matrix(rep(0,dim(DesMatrix)[1]*2),ncol=2), 
    #                  CovMis2 = as.matrix(rep(1,dim(DesMatrix)[1])),
    #                  gamma1 = 1, gamma=c(0,0), alpha1= -Inf, alpha0= -Inf, sigma_e = 0)
    #   betaIsd <- sqrt(diag(cov))} else {
    #     betaIsd <- rep(NA,length(betaI))
    #   }
    
    
    
    ## 2.4.4 RelNet ####
    cat("RelNet Method \n")
    ### 2.4.4.1 First Component ####
    
    Cov_Matrix <- cor(Xcov)
    
    results_RelNet <- lapply(cutpointrho, FUN= function(rho) {
      keep <- abs(Cov_Matrix) > rho
      
      ## Get True Positive and True Negative for ploting ROC
      
      TP <- sum((THETAgraph_high!=0) * keep )
      FP <- sum((THETAgraph_high==0) * keep )
      
      TN <- sum((THETAgraph_high==0) * (1-keep) )
      FN <- sum((THETAgraph_high!=0) * (1-keep) )
      
      return(list(RelNetTP=TP/(TP+FN), RelNetFP=FP/(FP+TN)))
    })
    
    
    #### 2.4.5 Comparing the results of the first stage ####
    ## 2.4.5.1 ROC curve for GNMM ####
    ROCGNMM <- ROC_curve_GNMM( varsel$path,  THETAgraph_high)
    
    RelNetTP <- unlist(lapply(results_RelNet, `[[`, 1))
    RelNetFP <- unlist(lapply(results_RelNet, `[[`, 2))
    
    
    #### 2.4.5.2 Comparing the results of the second stage ####
    # Norm1GNNM <- (colSums(abs(final_GNNM[,5:6]-final_GNNM[,3:4]),na.rm=T))
    # Norm1RelNet <- (colSums(abs(final_RelNet[,5:6]-final_RelNet[,3:4]),na.rm=T))
    # 
    # Norm2GNNM <- sqrt(colSums((final_GNNM[,5:6]-final_GNNM[,3:4])^2,na.rm=T))
    # Norm2RelNet <- sqrt(colSums((final_RelNet[,5:6]-final_RelNet[,3:4])^2,na.rm=T))
    
    
    return(list( ROCGNMM = ROCGNMM, 
                 RelNetTP = RelNetTP,
                 RelNetFP = RelNetFP))
    
   }, error = function(e)   return(NULL) )
} 



#### 3. Implementation ####
### 3.1 Example ####
### block plot ####
# SIM2_main(i = 1, thetas = 0.3, nsample = 200, graphtype="bloc")
# 
# ### hub plot ####
# SIM2_main(i = 1, thetas = 0.3, nsample = 250, graphtype="hub")
# 
# ### scale-free plot ####
# SIM2_main(i = 1, thetas = 0.3, nsample = 250, graphtype="scfr")

### scale-free high-dimension ####
# start <- proc.time()
# SIM2_high(i = 1, thetas = 0.3, nsample = 4000)
# end <- proc.time() - start
# end[3]

### 3.2 parallel computing ####
cl = makeCluster( ncore ,outfile=paste0("computelog/",ProjectName,".txt")) 

clusterExport(cl=cl, varlist=c("seed_i", "expit","likelihood_GQ",
                               "score_GQ","infomat_GQ","ROC_curve_GNMM","ROC_curve_LASSO",
                               "GEE_UI","GEE_SIGMA","GEE_GAMMA","GEE_GAMMA.inv","GEE_cov"))

registerDoParallel(cl)



#### 4.  Simulation Studies ####

# ### 4.1 Study 1: Low Dimension ####
# c("bloc","hub","scfr")
SIM1  <- foreach(nsample=c(50, 250, 500)) %:% foreach(grty=c("bloc","hub","scfr")) %:% foreach(var1=1:1000) %dopar% {
  tryCatch({
    SIM2_main(i = var1, thetas = 0.35, nsample = nsample, graphtype = grty)
  }, error = function(e)   return(NULL) )
}

save(SIM1,file="output/SIMGEE1_2.RData")

load(file="output/SIMGEE1_2.RData")

SIM1[[1]][[2]][[1]]

ResTable <- NULL


graphtype <- c("bloc","hub","scfr")
data.ROC <- NULL
nsample <- c(50, 250, 500)
GNMMAUClist <- NULL
RelNetAUClist <- NULL


for (k in 1:3){
  for (i in 1:3){
    GNMMtptab <- NULL
    GNMMfptab <- NULL
    GNMMAUCvec <- NULL
    
    for (j in 1:1000){
      GNMMtptab <- rbind(GNMMtptab,SIM1[[k]][[i]][[j]]$ROCGNMM$tp)
      GNMMfptab <- rbind(GNMMfptab,SIM1[[k]][[i]][[j]]$ROCGNMM$fp)
      GNMMAUCvec <- c(GNMMAUCvec,SIM1[[k]][[i]][[j]]$ROCGNMM$AUC)
    }
    GNMMtp <- c(0,colMeans(GNMMtptab,na.rm=T),1)
    GNMMfp <- c(0,colMeans(GNMMfptab,na.rm=T),1)
    GNMMAUC <- getAUC(GNMMfp,GNMMtp)
    GNMMAUClist <- c(GNMMAUClist, GNMMAUC)
    
    # method <- c(rep("GNMM",length(GNMMtp)),rep("LASSO1",length(LASSO1tp)),rep("LASSO2",length(LASSO2tp)))
    tps <- GNMMtp
    fps <- GNMMfp
    
    data.ROC.this <- data.frame(tps=tps, fps=fps, graphtype = graphtype[i], nsample=nsample[k], method = "GNSM")
    
    data.ROC <- rbind(data.ROC,data.ROC.this)
    
    RelNettptab <- NULL
    RelNetfptab <- NULL
    
    for (j in 1:1000){
      RelNettptab <- rbind(RelNettptab,SIM1[[k]][[i]][[j]]$RelNetTP)
      RelNetfptab <- rbind(RelNetfptab,SIM1[[k]][[i]][[j]]$RelNetFP)
    }
    RelNettptab <- c(0,colMeans(RelNettptab,na.rm=T),1)
    RelNetfptab <- c(0,colMeans(RelNetfptab,na.rm=T),1)
    RelNetAUC <- getAUC(RelNetfptab,RelNettptab)
    
    # method <- c(rep("GNMM",length(GNMMtp)),rep("LASSO1",length(LASSO1tp)),rep("LASSO2",length(LASSO2tp)))
    # tps <- GNMMtp
    # fps <- GNMMfp
    
    data.ROC.this <- data.frame(tps=RelNettptab, fps=RelNetfptab, graphtype = graphtype[i], 
                                nsample=nsample[k], method = "RelNet")
    
    data.ROC <- rbind(data.ROC,data.ROC.this)
    
    RelNetAUClist <- c(RelNetAUClist, RelNetAUC)
  }
}


### 4.2 Study 2: High Dimension ####


SIM1_High  <- foreach(nsample=c(4000)) %:% foreach(grty=c("scfr")) %:% foreach(var1=1:1000) %dopar% {
  start <- proc.time()
  results <-  SIM2_high(i = var1, thetas = 0.3, nsample = 5000)
  end <- proc.time() - start
  cat(var1, ":", end[3]," \n")
  return(results)
}

save(SIM1_High,file=paste0("output/SIM1_High.RData"))

load(file="output/SIM1_High.RData")



SIM1_High[[1]][[1]][[1]]

ResTable <- NULL


graphtype <- c("scfr")
# data.ROC <- NULL
nsample <- c(4000)
GNMMAUClist <- NULL
RelNetAUClist <- NULL
data.ROC.high <- NULL


for (k in 1:1){
  for (i in 1:1){
    GNMMtptab <- NULL
    GNMMfptab <- NULL
    GNMMAUCvec <- NULL
    
    for (j in 1:1000){
      GNMMtptab <- rbind(GNMMtptab,SIM1_High[[k]][[i]][[j]]$ROCGNMM$tp)
      GNMMfptab <- rbind(GNMMfptab,SIM1_High[[k]][[i]][[j]]$ROCGNMM$fp)
      GNMMAUCvec <- c(GNMMAUCvec,SIM1_High[[k]][[i]][[j]]$ROCGNMM$AUC)
    }
    GNMMtp <- c(0,colMeans(GNMMtptab,na.rm=T),1)
    GNMMfp <- c(0,colMeans(GNMMfptab,na.rm=T),1)
    GNMMAUC <- getAUC(GNMMfp,GNMMtp)
    GNMMAUClist <- c(GNMMAUClist, GNMMAUC)
    
    # method <- c(rep("GNMM",length(GNMMtp)),rep("LASSO1",length(LASSO1tp)),rep("LASSO2",length(LASSO2tp)))
    tps <- GNMMtp
    fps <- GNMMfp
    
    data.ROC.this <- data.frame(tps=tps, fps=fps, graphtype = graphtype[i], nsample=nsample[k], method = "GNSM")
    
    data.ROC.high <- rbind(data.ROC.high,data.ROC.this)
    
    RelNettptab <- NULL
    RelNetfptab <- NULL
    
    for (j in 1:1000){
      RelNettptab <- rbind(RelNettptab,SIM1_High[[k]][[i]][[j]]$RelNetTP)
      RelNetfptab <- rbind(RelNetfptab,SIM1_High[[k]][[i]][[j]]$RelNetFP)
    }
    RelNettptab <- c(0,colMeans(RelNettptab,na.rm=T),1)
    RelNetfptab <- c(0,colMeans(RelNetfptab,na.rm=T),1)
    RelNetAUC <- getAUC(RelNetfptab,RelNettptab)
    
    # method <- c(rep("GNMM",length(GNMMtp)),rep("LASSO1",length(LASSO1tp)),rep("LASSO2",length(LASSO2tp)))
    # tps <- GNMMtp
    # fps <- GNMMfp
    
    data.ROC.this <- data.frame(tps=RelNettptab, fps=RelNetfptab, graphtype = graphtype[i], 
                                nsample=nsample[k], method = "RelNet")
    
    data.ROC.high <- rbind(data.ROC.high,data.ROC.this)
    
    RelNetAUClist <- c(RelNetAUClist, RelNetAUC)
  }
}





######## Present the AUG values in the curve #####
### remove this part as suggested by Grace

# dat_text0 <- data.frame(
#   label = unlist(lapply(GNMMAUClist, FUN = function(x){
#     paste0("AUC")
#   })),
#   nsample   = rep(c(200,1000), each=3),
#   graphtype = rep(c("block","hub","scale-free"),times =2 ),
#   fps = 0.75,
#   tps = 0.35,
#   stringsAsFactors = FALSE
# )
#
# dat_text1 <- data.frame(
#   label1 = unlist(lapply(GNMMAUClist, FUN = function(x){
#     paste0("GNM    ",sprintf("%.3f", round(x,digits =3)))
#     # paste0("GNMM")
#   })),
#   nsample   = rep(c(200,1000), each=3),
#   graphtype = rep(c("block","hub","scale-free"),times =2 ),
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
#   nsample   = rep(c(200,1000), each=3),
#   graphtype = rep(c("block","hub","scale-free"),times =2 ),
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
#   nsample   = rep(c(200,1000), each=3),
#   graphtype = rep(c("block","hub","scale-free"),times =2 ),
#   fps = 0.7,
#   tps = 0.05,
#   stringsAsFactors = FALSE
# )

pdf(file = "output/fig1_varselect.pdf",height = 14, width = 18)


datp1 <- ggplot(data.ROC %>%
               filter(method == "GNSM" & graphtype  == "hub") %>%
               mutate(nsample = factor(nsample))
               , aes(x=fps, y=tps, group=nsample)) +
  geom_line(aes(linetype=nsample,color=nsample),size=1.5)+
  scale_color_manual(values = c("#7A989A", "#CF9546", "#C67052"))+
  # scale_color_brewer(type = 'div', palette = 'Blues', direction = 1) +
  geom_point(aes(shape=nsample,color=nsample))+
  labs(x = "False Negative Rate", y="True Positive Rate", shape = "Sample Size", linetype = "Sample Size", color = "Sample Size") +
  labs(title="(a)") +
  theme_bw() +
  theme(text=element_text(size=12, family="mono"), axis.text = element_text(size = 14, family="mono"), panel.spacing = unit(1, "lines"),
        legend.position="bottom") 
# c("bloc","hub","scfr")

datp2 <- ggplot(data.ROC %>%
                  filter(nsample == 250 & method == "GNSM") %>%
                  mutate(graphtype = factor (graphtype, levels = c("bloc","hub","scfr"), labels = c("block","hub","scale-free")))
                , aes(x= fps, y= tps, group= graphtype)) +
  geom_line(aes(linetype=graphtype,color=graphtype),size=1.5)+
  scale_color_manual(values = c("#7A989A", "#C1AE8D", "#C67052")) +
  # scale_color_brewer(type = 'div', palette = 'Blues', direction = 1) +
  geom_point(aes(shape=graphtype,color=graphtype))+
  labs(x = "False Negative Rate", y="True Positive Rate", shape = "Graph Type", linetype = "Graph Type", color = "Graph Type") +
  labs(title="(c)") +
  theme_bw() +
  theme(text=element_text(size=12, family="mono"), axis.text = element_text(size = 14, family="mono"), panel.spacing = unit(1, "lines"),
        legend.position="bottom") 

datp3 <- ggplot(data.ROC %>%
                  filter(nsample == 250 & graphtype  == "hub"), aes(x=fps, y=tps, group=method)) +
  geom_line(aes(linetype=method,color=method),size=1.5)+
  scale_color_manual(values = c("#7A989A", "#CF9546", "#C67052"))+
  # scale_color_brewer(type = 'div', palette = 'Blues', direction = 1) +
  geom_point(aes(shape=method,color=method))+
  labs(x = "False Negative Rate", y="True Positive Rate", shape = "Method", linetype = "Method", color = "Method") +
  labs(title="(b)") +
  theme_bw() +
  theme(text=element_text(size=12, family="mono"), axis.text = element_text(size = 14, family="mono"), panel.spacing = unit(1, "lines"),
        legend.position="bottom") 

datp4 <- ggplot(data.ROC.high, aes(x=fps, y=tps, group=method)) +
  geom_line(aes(linetype=method,color=method),size=1.5)+
  scale_color_manual(values = c("#7A989A", "#CF9546", "#C67052"))+
  # scale_color_brewer(type = 'div', palette = 'Blues', direction = 1) +
  geom_point(aes(shape=method,color=method))+
  labs(x = "False Negative Rate", y="True Positive Rate", shape = "Method", linetype = "Method", color = "Method") +
  labs(title="(d)") +
  theme_bw() +
  theme(text=element_text(size=12, family="mono"), axis.text = element_text(size = 14, family="mono"), panel.spacing = unit(1, "lines"),
        legend.position="bottom") 


(datp1 + datp3) / ( datp2 + datp4 )

# p1 +
#   geom_text(
#     family = "mono", fontface="bold", size =5,
#     data    = dat_text0,
#     mapping = aes(x = fps, y = tps, label = label),
#     inherit.aes=F
#  )+
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
dev.off()


#### 4.1.2 Step 2 Evaluation ####

Results <- NULL

samplesize_range <- c(50, 250,500)
graph_range <- c("block","hub","scale-free")

for (j in 1:3){
    for (k in 1:3) {
    results <- SIM1[[k]][[j]]

    Norm1GNNM <- NULL
    Norm1RelNet <- NULL
    Norm2GNNM <- NULL
    Norm2RelNet <- NULL

    for (i in 1:1000){
      Norm1GNNMi <- results[[i]]$Norm1GNNM
      Norm1RelNeti <- results[[i]]$Norm1RelNet
      Norm2GNNMi  <- results[[i]]$Norm2GNNM
      Norm2RelNeti <- results[[i]]$Norm2RelNet

      Norm1GNNM <- rbind(Norm1GNNM,Norm1GNNMi)
      Norm1RelNet <- rbind(Norm1RelNet,Norm1RelNeti)
      Norm2GNNM <- rbind(Norm2GNNM,Norm2GNNMi)
      Norm2RelNet <- rbind(Norm2RelNet,Norm2RelNeti)
    }

    Norm1GNNMv <- colMeans(na.omit(Norm1GNNM),na.rm = T)
    Norm1RelNetv <- colMeans(na.omit(Norm1RelNet),na.rm = T)
    Norm2GNNMv <- colMeans(na.omit(Norm2GNNM),na.rm = T)
    Norm2RelNetv <- colMeans(na.omit(Norm2RelNet),na.rm = T)

    Results0 <- data.frame(samplesize = samplesize_range[k],
                           graph = graph_range[j],
                           type= c("beta1","beta2"),
                           Norm1GNNM=round(Norm1GNNMv,3),Norm2GNNM=round(Norm2GNNMv,3),
                           Norm1RelNet=round(Norm1RelNetv,3),Norm2RelNet=round(Norm2RelNetv,3))

    Results <- rbind(Results,Results0)
  }
}

# Results_Long <- gather(Results, key = type, value = Measurements, Norm1GNNM:Norm2GNNM)
# Results_Long <- cbind(Results_Long, Norm = rep(c("Norm1","Norm2"), each = dim(Results_Long)[1]/2))
#
# Results_rar <- Results_Long  %>% arrange((graph)) %>% arrange((samplesize)) %>% arrange((type))
#
# # Results_c <- data.frame(Results_rar[Results_rar$Norm=="Norm1",], Norm2 = Results_rar[Results_rar$Norm=="Norm2","Measurements"])
# Results_rar$Method <- ifelse (Results_rar$Method %in% c("Norm1GNNM","Norm2GNNM"),"GNNM","LASSO")


Results_Wide1 <- spread(Results[,c("samplesize","graph", "type", "Norm1GNNM")], key = type , value = Norm1GNNM )
Results_Wide2 <- spread(Results[,c("samplesize","graph", "type", "Norm2GNNM")], key = type , value = Norm2GNNM )
Results_Wide3 <- spread(Results[,c("samplesize","graph", "type", "Norm1RelNet")], key = type , value = Norm1RelNet )
Results_Wide4 <- spread(Results[,c("samplesize","graph", "type", "Norm2RelNet")], key = type , value = Norm2RelNet )



Results_final <- cbind(Results_Wide1,Results_Wide2[,3:4],Results_Wide3[,3:4],Results_Wide4[,3:4])

xtable(Results_final, digits = 3)

