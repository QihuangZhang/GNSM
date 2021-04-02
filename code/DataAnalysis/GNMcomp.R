
### Project GEEmix

### This is the project with internal validation data

#### 0 Project Information and Configuration ####
ProjectName<- paste0("DataAnalysis")
cat(paste0("Note: This is for the data analysis of SNP 109 \n")) 
set.seed(2019)

## 0.1 Set file path ####
# setwd("N:/Work/waterloo/2019/GNMM")
# WD=paste0(getwd(),"/results/EMsensi/")


## Finished part -- no need to run again
# genotype<-read.table(file="data/geno.txt",header=T)
# save(genotype,file="genotype.RData")

### Preperation of the Genotype Data

load("~/Work/waterloo/2018/mixResponse/file/genotype.RData")
# 
SNPlist <- c("id","discard","rs45690064","rs27338905","rs32962338","rs33583459","rs224051056",
             "rs33217671","rs38916331","rs47869247","rs217439518","rs29477109",
             "rs252503010","rs265727287","rs246035173","rs231489766","rs46826545",
             "rs51809856","rs6279141","rs30535702","rs30201629","rs30549753")
# 
genotypes <- genotype[,colnames(genotype) %in% SNPlist]
# 
# save(genotypes,file="~/Work/waterloo/2019/GNMM/file/genotypes.RData")
  
load("~/Work/waterloo/2019/GNMM/file/genotypes.RData")
phenotype<-read.csv(file="~/Work/waterloo/2018/mixResponse/data/pheno.csv")

### Note:!! if use windows, change ~/ back to P:/ or N:/


## 0.2 Prepare package of GWAS toolbox ####
# library(maxLik)
library(MASS)
library(GeneErrorMis)
library(parallel)
library(nleqslv)
library(huge)
library(igraph)
# source("code/CGRM.R")

## 0.3 Prepare the datasets ####

## Predict the TAstar
# modelTA <- lm(TA~gastroc+EDL+plantaris, data=phenotype)
# phenotype$TAstar <- predict(modelTA, newdata=phenotype)
# 
BMD90 <- quantile(phenotype$BMD,na.rm=T, probs = 0.9)
phenotype$BCstar <- ifelse(phenotype$BMD<BMD90,0,1)

## Association Data Frame
data_GWA_Pre <- data.frame(id = phenotype$id,Y1 = phenotype$tibia,Y2 = phenotype$BCstar,
                           discard = ifelse(phenotype$discard=="no",1,0),
                           mixup = ifelse(phenotype$mixup=="no",0,1), 
                           BMD=phenotype$BMD, Y1star=phenotype$tibia, Y2star=phenotype$abnormalbone)


data_GWAS_unsorted <- data_GWA_Pre[data_GWA_Pre$id %in% genotype$id,]
data_GWAS <- data_GWAS_unsorted[order(data_GWAS_unsorted$id),]

# Response
Y1star <- data_GWAS$Y1star
Y2star <- data_GWAS$Y2star
# Y1 <- data_GWAS$Y1
# Y2 <- data_GWAS$Y2

# ## Covariates for main model
# Covariates<-data_GWAS[,c("BW","BMD","mixup","batch")]
# Covariates<-as.matrix(Covariates)



### Deciding covariates dependence structure

### Standardize the candidate genotypes 

scaled.gene <- scale(genotypes[,3:dim(genotypes)[2]])


varsel <- huge(scaled.gene, lambda = NULL, nlambda = 30, lambda.min.ratio = NULL, method = "mb",
               scr = F, scr.num = NULL, sym = "or", verbose = TRUE, cov.output =T)
out.select = huge.select(varsel, criterion = "stars", stars.thresh = 0.05,rep.num=30)


THETAgraph <- varsel$path[[out.select$opt.index]]

EdgeTrue <- NULL  
EdgeHash <- rep(T,20)

for (i in 1:19){
  for (j in (i+1):20){
    if (THETAgraph[i,j]!=0){
      EdgeTrue <- rbind(EdgeTrue,c(i,j))
      EdgeHash <- c(EdgeHash,T)
    } else { EdgeHash <- c(EdgeHash,F)}
  }
}

cat(EdgeTrue)

DesMatrix <- scale(genotypes[,3:dim(genotypes)[2]])


for (i in 1:dim(EdgeTrue)[1]){
  DesMatrix <- cbind(DesMatrix,DesMatrix[,EdgeTrue[i,1]]*DesMatrix[,EdgeTrue[i,2]])
}

# DesMatrix <- scale(DesMatrix)


#### plot the selected graph
plot.selectgraph = function(x, ...){
  if(x$cov.input){
    cat("Model selection is not available when using the covariance matrix as input.")
  }
  if(!x$cov.input)
  {
    par(mfrow=c(1,2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
    
    g = graph.adjacency(as.matrix(x$refit), mode="undirected", diag=FALSE)
    layout.grid = layout.fruchterman.reingold(g)
    
    plot(g, layout=layout.grid, edge.color='gray50',vertex.color="white", vertex.size=10, vertex.label=SNPlist[-1*1:2], vertex.label.dist=0.5,vertex.label.degree =pi)
    plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameter", ylab = "Sparsity Level", type = "l",xlim = rev(range(x$lambda)), main = "Solution path sparsity levels")
    lines(x$opt.lambda,x$opt.sparsity,type = "p")
  }
}


pdf(file = "output/figdatacov.pdf",height = 5.6, width = 10)
set.seed(2017)
plot.selectgraph(out.select)
dev.off()


## Mark the complete sample
comp<-data.frame(Y1star,Y2star,DesMatrix)
index.notna<-complete.cases(comp)
index.notna <- ifelse(data_GWAS$discard == 1, index.notna , FALSE)

Y1star<-Y1star[index.notna]
Y2star<-Y2star[index.notna]
# Y1<-Y1[index.notna]
# Y2<-Y2[index.notna]
Covariates_all<-comp[index.notna,]
nsample<-dim(Covariates_all)[1]


CovMis1 <- as.matrix(data.frame(intercept=rep(0,nsample),cov=rep(0,nsample)))
  
CovMis2 <- as.matrix(data.frame(intercept=rep(1,nsample)))




## 0.4 Function set up ####

## 0.4.1 Functions of Original Version ####
GEE_UI <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, alpha1, alpha0, sigma_e){
  # cat(theta, " \n")
  return(GEE_UfuncIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                      theta[1:4], theta[5:8], sigma = theta[9], xi = theta[10], 
                      gamma1 = 1, gamma=c(0,0),  alpha1, alpha0, sigma_e
  ))
}

GEE_SIGMA <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, gamma, alpha1, alpha0,sigma_e){
  nbeta <- dim(DesignMatrix1)[2]
  return(GEE_SIGMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                      theta[1:nbeta], theta[(nbeta+1):(2*nbeta)], sigma = theta[2*nbeta+1], xi = theta[2*nbeta+2], 
                      gamma1 = 1, gamma, alpha1, alpha0, sigma_e))
}

GEE_GAMMA <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2){
  GAMMA <- GEE_GAMMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, 
                        beta1=theta[1:nbeta], beta2=theta[(nbeta+1):(2*nbeta)], sigma = theta[2*nbeta+1], xi = theta[2*nbeta+2])
  return(GAMMA)
}

GEE_GAMMA.inv <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2){
  nbeta <- dim(DesignMatrix1)[2]
  GAMMA <- GEE_GAMMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, 
                        beta1=theta[1:nbeta], beta2=theta[(nbeta+1):(2*nbeta)], sigma = theta[2*nbeta+1], xi = theta[2*nbeta+2])
  GAMMA.inv <- solve(GAMMA,tol=1e-200)
  return(GAMMA.inv)
}

GEE_cov <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, gamma, alpha1, alpha0, sigma_e){
  GAMMA.inv <- GEE_GAMMA.inv(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
  SIGMA <- GEE_SIGMA(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, gamma, alpha1, alpha0, sigma_e)
  covmatrix <- GAMMA.inv %*% SIGMA %*% t(as.matrix(GAMMA.inv))
  return(covmatrix)
}

GEE_UI_ErrMis <- function(Theta, Y1star, Y2star, Covariates, CovMis1, CovMis2,
                          gamma1, gamma, alpha1, alpha0, sigma_e){
  # cat(theta, " \n")
  nbeta <- dim(Covariates)[2]
  return(GEE_UfuncIns(Y1star, Y2star, DesignMatrix1=as.matrix(Covariates), 
                      DesignMatrix2=as.matrix(Covariates),  CovMis1, CovMis2,
                      beta1 = Theta[1:nbeta], 
                      beta2 = Theta[(nbeta+1):(2*nbeta)], 
                      sigma = Theta[2*nbeta+1], xi = Theta[2*nbeta+2], 
                      gamma1, gamma, alpha1=alpha1, alpha0=alpha0, sigma_e))
}

#### 1. Data Preparation ####

Covariates1 <- Covariates2 <- as.matrix(cbind(Covariates_all[,! names(Covariates_all) %in% c("Y1star","Y2star","BW")])) 

CovMis1 <- as.matrix(data.frame(intercept=rep(1,nsample),cov=rep(1,nsample)))

CovMis2 <- as.matrix(data.frame(intercept=rep(1,nsample)))

# indexVal <- ifelse(is.na(Covariates_all[,"BMD"])|(Covariates_all[,"mixup"]==1)|(Covariates_all[,"batch"]==1),F,T)

# mismeasuredata <- data.frame(Y1star=Y1star[!indexVal], Y2star=Y2star[!indexVal],
#                              Covariates1[!indexVal,],CovMis1[!indexVal,],CovMis2[!indexVal,])
# validationdata <- data.frame(Y1star=Y1star[indexVal], Y2star=Y2star[indexVal],
#                              Y1=Y1[indexVal], Y2=Y2[indexVal],
#                              Covariates1[indexVal,],CovMis1[indexVal,],CovMis2[indexVal,])

completedata <- as.matrix(cbind(Covariates_all[,! names(Covariates_all) %in% c("BW")])) 

#### 2. Analysis  of 109 ####


## 2.1 Naive Analysis of 109 ####


# beta_Y1_0 <- mean(mismeasuredata$Y1star)
# beta_Y2_0 <- log(mean(mismeasuredata$Y2star)/(1-mean(mismeasuredata$Y2star)))



## 2.2 Estimating Procedure ####

# naive.model1 <- lm(Y1star ~ 1)
# naive.model2 <- glm(Y2star ~ 1)

## 2.2.1 Naive model ####
naive.model1 <- lm(Y1star ~ -1+ ., data = data.frame(completedata[,-2]))
naive.model2 <- glm(Y2star ~ -1+ ., family = binomial(link = logit), data = data.frame(completedata[,-1]))

summary(naive.model1)
summary(naive.model2)


initial2 <- c(naive.model1$coefficients,naive.model2$coefficients,0.03,0)
cat("Initial value:", initial2, "\n")


# naive.model3 <- lm(Y1 ~ X905 + PC1 + PC2, data = validationdata)
# naive.model4 <- glm(Y2 ~ X905 + PC1 + PC2, family = binomial(link = logit), data = validationdata)
# 
# intial3 <- c(naive.model3$coefficients,naive.model4$coefficients,1,0)
# cat("Initial value:", intial3, "\n")


## 2.2.2 the proposed approach ####

## Measurement Error and Misclassification Parameters
# model.measure <- lm(Y1star ~ Y2 + offset(Y1),data = validationdata) 
# model.class1 <- glm((1-Y2star) ~ BW, data = validationdata[validationdata$Y2==1,],family = binomial(link="logit")) 
# model.class0 <- glm(Y2star ~ BW, data = validationdata[validationdata$Y2==0,],family = binomial(link="logit")) 
# 
# gamma0 <- model.measure$coefficients[1]
# gamma2 <- model.measure$coefficients[2] 
# sigma_e <- sigma(model.measure)
# alpha1 <- model.class1$coefficients
# alpha0 <- model.class0$coefficients


GEE_UI_ErrMis <- function(Theta, Y1star, Y2star, Covariates, CovMis1, CovMis2,
                          gamma1, gamma, alpha1, alpha0, sigma_e){
  # cat(theta, " \n")
  nbeta <- dim(Covariates)[2]
  return(GEE_UfuncIns(Y1star, Y2star, DesignMatrix1=as.matrix(Covariates), 
                      DesignMatrix2=as.matrix(Covariates),  CovMis1, CovMis2,
                      beta1 = Theta[1:nbeta], 
                      beta2 = Theta[(nbeta+1):(2*nbeta)], 
                      sigma = Theta[2*nbeta+1], xi = Theta[2*nbeta+2], 
                      gamma1, gamma, alpha1=alpha1, alpha0=alpha0, sigma_e))
}


SensiAnalysis <- function(alphasensi,sigma_esensi){
  
  NR <- nleqslv(initial2, GEE_UI_ErrMis, Y1star=Y1star, Y2star=Y2star, Covariates=Covariates1,
                jacobian=T, control=list(maxit=10000),
                CovMis1=CovMis1, CovMis2=CovMis2,
                gamma1 = 1, gamma = c(0,0), alpha1=alphasensi, alpha0=alphasensi, sigma_e= sigma_esensi)
  
  betahat <- ifelse(abs(NR$x)<100,NR$x,NA)
  
  
  ## 2.3.2.2 Results Rearrangement ####
  
  # thetaresults <- betahat
  # betaresults <- data.frame(cordinates, beta1 = thetaresults[1:dim(cordinates)[1]], 
  #                           beta2 = thetaresults[(dim(cordinates)[1]+1):(2*dim(cordinates)[1])])
  # 
  # 
  # proposed_GNNM <- merge(cordtrue, betaresults,by = c("X1","X2"), all.x = T, all.y = F)
  # 
  # ## 2.3.2.3 Variance Estimation ####
  # 
  # betaI_proposed <- c(proposed_GNNM$beta1.y,proposed_GNNM$beta2.y,betahat[(length(betahat)-1):length(betahat)])
  # betaI_proposed <- ifelse(is.na(betaI_proposed),0,betaI_proposed)
  
  
  if (!any(is.na(betahat))) {
    cov <- GEE_cov(betahat,Y1star = Y1star, Y2star = Y2star, 
                   DesignMatrix1 = as.matrix(Covariates1),
                   DesignMatrix2 = as.matrix(Covariates1), 
                   CovMis1 = matrix(rep(0,dim(Covariates1)[1]*2),ncol=2), 
                   CovMis2 = as.matrix(rep(1,dim(Covariates1)[1])),
                   gamma=c(0,0), alpha1= alphasensi, alpha0= alphasensi, sigma_e = sigma_esensi)
    betaIsd_proposed <- sqrt(diag(cov))} else {
      betaIsd_proposed <- rep(NA,length(betahat))
    }
  
  
  # betahat <- c(betahat,gamma0,gamma2,sigma_e,alpha1,alpha0)
  
  Zvalue <- betahat/betaIsd_proposed
  pvalue <- 2*(pnorm(-abs(Zvalue)))
  pscale <- -log(2*(pnorm(-abs(Zvalue))), base = 10)
  # pscale <- -log(2*(1-pnorm(abs(Zvalue))))
  
  Table_numeric <- data.frame(SNPnames = c(rep(c(SNPlist[-1*1:2],"v1","v2","v3","c4"),2),"sigma","phi"),
                              propobeta=betahat,
                              proposd = betaIsd_proposed,
                              propoZ = Zvalue,
                              propop = pvalue,
                              pscale = pscale)

}

library(xtable)


Tablelow0.77 <- SensiAnalysis(-2.944,0.77) 

xtable(Tablelow0.77,digits = 3)

Tablemedian0.77 <- SensiAnalysis(-2.197,0.77) 

xtable(Tablemedian0.77,digits = 3)

Tablehigh0.77 <- SensiAnalysis( -1.386,0.77) 

xtable(Tablehigh0.77,digits = 3)

