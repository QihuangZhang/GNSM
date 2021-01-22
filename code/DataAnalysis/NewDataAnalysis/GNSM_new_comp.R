###
### Project GNM
## As advised by referee, we conduct a two-stage GWAS, where in the first stage, we filter the SNPs with the highest association.
## A sensitivity study is conducted to evaluate the performance of the method.

#### 0 Project Information and Configuration ####
ProjectName<- paste0("DataAnalysis")
set.seed(2019)

## 0.1 Set file path ####
# setwd("N:/Work/waterloo/2019/GNMM")
# WD=paste0(getwd(),"/results/EMsensi/")


## Finished part -- no need to run again
# genotype<-read.table(file="data/geno.txt",header=T)
# save(genotype,file="genotype.RData")

### Preperation of the Genotype Data

load("~/Work/2018/mixResponse(conclude)/file/genotype.RData")
# 

# 
# save(genotypes,file="~/Work/waterloo/2019/GNMM/file/genotypes.RData")

# load("~/Work/waterloo/2019/GNMM/file/genotypes.RData")
phenotype<-read.csv(file="~/Work/2018/mixResponse(conclude)/data/pheno.csv")

### Note:!! if use windows, change ~/ back to P:/ or N:/


## 0.2 Prepare package of GWAS toolbox ####
# library(maxLik)
library(MASS)
library(GeneErrorMis)
library(parallel)
library(nleqslv)
library(huge)
library(igraph)
library(dplyr)
library(geomnet)
library(patchwork)
library(stringr)
library(ggraph)
library(xtable)



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
  
  
  Zvalue <- betahat/betaIsd_proposed
  pvalue <- 2*(pnorm(-abs(Zvalue)))
  pscale <- -log(2*(pnorm(-abs(Zvalue))), base = 10)
  # pscale <- -log(2*(1-pnorm(abs(Zvalue))))
  
  Table_numeric <- data.frame(SNPnames = c(SNPnames, "sigma", "phi"),
                              propobeta=betahat,
                              proposd = betaIsd_proposed,
                              propoZ = Zvalue,
                              propop = pvalue,
                              pscale = pscale)
}

integrate <- function(small,median,large){
  Tablelist <- list(small,median,large)
  FinalTable <- NULL
  
  FinalTablelist <- lapply(1:3, FUN=function(i){
    TableTemp <- data.frame(SNPnames=Tablelist[[i]][,1], round(Tablelist[[i]][,2:dim(Tablelist[[i]])[2]],3))
    TableTemp <- TableTemp[-1*c(dim(Tablelist[[i]])[1]-1,dim(Tablelist[[i]])[1]),]
    paste(TableTemp[,2:3],sep="")
    TableTemp$parabias <-apply(TableTemp[,2:3], MARGIN = 1, FUN = function(x){
      return(paste0(sprintf("%.3f",x[1])," (", sprintf("%.3f",x[2]),")"))
    })
    
    TableTemp <- TableTemp[,c("SNPnames","parabias","propop")]
  })
  
  FinalTable <- do.call(cbind, FinalTablelist)
  FinalTable <- FinalTable[,-c(4,7)]
  
  FinalTable2 <- cbind(FinalTable[1:(dim(FinalTable)[1]/2),],
                       FinalTable[(dim(FinalTable)[1]/2+1):dim(FinalTable)[1],])
  return(FinalTable2[,-8])
}


integrate2 <- function(small,median,large, titles){
  
  Tablemains1 <- small %>%
    filter(!SNPnames %in% c("sigma", "phi")) %>%
    mutate(method = titles[1]) 
  Tablemains1 <- Tablemains1 %>% 
    mutate(type = rep(c("continuous", "discrete"),each = dim(Tablemains1)[1]/2))
  
  Tablemains2 <- median %>%
    filter(!SNPnames %in% c("sigma", "phi")) %>%
    mutate(method = titles[2]) 
  Tablemains2 <- Tablemains2 %>% 
    mutate(type = rep(c("continuous", "discrete"),each = dim(Tablemains2)[1]/2))
  
  Tablemains3 <- large %>%
    filter(!SNPnames %in% c("sigma", "phi")) %>%
    mutate(method = titles[3]) 
  Tablemains3 <- Tablemains3 %>% 
    mutate(type = rep(c("continuous", "discrete"),each = dim(Tablemains3)[1]/2))
  
  outputtable <- rbind(Tablemains1, Tablemains2, Tablemains3)
  
  outputtable$method <- factor(outputtable$method, levels = titles)
  
  return(outputtable)
}

plotintegrate2 <- function(Atable) {
  TablefirstPara <- Atable %>% 
    mutate(HLPara = ifelse(abs(propobeta)>1.5,"Highlight","Not Highlight")) %>%
    mutate(HLParalabel = ifelse(HLPara=="Highlight",SNPnames,"")) %>%
    mutate(HLParalabel = gsub("[.]x[.]", " x ", HLParalabel)) %>%
    mutate(HLParalabel = gsub("[.]", "-", HLParalabel)) 
  
  p1 <- ggplot(TablefirstPara, aes(x=type,y=propobeta)) +
    geom_point(aes(color = HLPara), size =2) +
    facet_wrap(.~method) + 
    geom_text_repel(aes(label = HLParalabel), color="#C1AE8D", fontface=2, alpha = 0.9) +
    scale_color_manual(values = c("#C67052", "#7A989A")) +
    theme_bw() + 
    theme(text=element_text(size=12, family="mono"), axis.text = element_text(size = 14, family="mono"),
          legend.position = "none") +
    labs(x = "Data Type", y = "Parameter Estimate") +
    labs(title="(a)") +
    theme(strip.background =element_rect(fill="#4F534A",color="#4F534A"))+ # #535b44
    theme(strip.text = element_text(colour = 'white')) +
    theme(panel.border = element_rect(colour = "#4F534A")) 
  
  TablefirstPvalue <- Atable %>% 
    mutate(HLPvalue = ifelse(abs(pscale)>3.6,"Highlight","Not Highlight")) %>%
    mutate(HLPvaluelabel = ifelse(HLPvalue=="Highlight",SNPnames,""))  %>%
    mutate(HLPvaluelabel = gsub("[.]x[.]", " x ", HLPvaluelabel)) %>%
    mutate(HLPvaluelabel = gsub("[.]", "-", HLPvaluelabel)) 
  
  p2 <- ggplot(TablefirstPvalue, aes(x=type,y=pscale)) +
    geom_point(aes(color = HLPvalue), size =2) +
    facet_wrap(.~method) + 
    geom_text_repel(aes(label = HLPvaluelabel), color="#C1AE8D", fontface=2, alpha = 0.9) +
    scale_color_manual(values = c("#C67052", "#7A989A")) +
    theme_bw() + 
    theme(text=element_text(size=12, family="mono"), axis.text = element_text(size = 14, family="mono"),
          legend.position = "none") +
    labs(x = "Data Type", y = bquote("-"~log[10]~"P-value")) +
    labs(title="(b)") +
    theme(strip.background =element_rect(fill="#4F534A",color="#4F534A"))+ # #535b44
    theme(strip.text = element_text(colour = 'white')) +
    theme(panel.border = element_rect(colour = "#4F534A"))
  
  print(p1 / p2 +  plot_layout(heights = c(2, 1)))
   
}

EdgeBundling <- function(Tablefirstat, outcomearg, methodarg) {
  
  ## The function of preparing the data set for edge bundling plot
  
  connect <- Tablefirstat %>% 
    filter(method==methodarg) %>% 
    filter(type==outcomearg) %>% 
    filter(!SNPname2=="") %>% 
    mutate(from=SNPname1) %>% 
    mutate(to=SNPname2)%>%
    select(-c("SNPname1", "SNPname2"))
  
  vertices <- Tablefirstat %>% 
    filter(method==methodarg) %>% 
    filter(type==outcomearg) %>% 
    filter(SNPname2=="") 
  
  d1 <- data.frame(from = "origin", to = "group1")
  d2 <- data.frame(from = "group1", to = vertices$SNPname1)
  myedges <- rbind(d1, d2)
  
  
  vertices[nrow(vertices)+1:2,] <- NA
  vertices[nrow(vertices)-0:1,"SNPname1"] <- c("origin", "group1")
  vertices[nrow(vertices)-0:1,"pscale"] <- c(0, 0)
  
  
  
  
  
  #calculate the ANGLE of the labels
  vertices$id <- NA
  myleaves <- which(is.na( match(vertices$SNPname1, myedges$from) ))
  nleaves <- length(myleaves)
  vertices$id[ myleaves ] <- seq(1:nleaves)
  vertices$angle <- 90 - 360 * vertices$id / nleaves
  
  # calculate the alignment of labels: right or left
  # If I am on the left part of the plot, my labels have currently an angle < -90
  vertices$hjust <- ifelse( vertices$angle < -90, 1, 0)
  
  # flip angle BY to make them readable
  vertices$angle <- ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)
  
  
  vertices <- vertices %>% mutate(name = SNPname1) %>% select(-c("SNPname1","SNPname2"))
  
  mygraph <- igraph::graph_from_data_frame( myedges, vertices=vertices[c("name","pscale","id","angle","hjust")] )
  
  # The connection object must refer to the ids of the leaves:
  from  <-  match( connect$from, vertices$name)
  to  <-  match( connect$to, vertices$name)
  
  
  return(list(connect=connect, 
              mygraph=mygraph, 
              vertices=vertices, 
              from = from, to = to))
}



#### 1. First Stage - SNP screening ####

## Filter the SNPs with MAF < 0.05
KeepSet <- colMeans(genotype[,3:dim(genotype)[2]])/2>0.05 


## After filtering, keep 77848 SNPs
genotype_QCed <- genotype[,c(F,F,KeepSet)]


## prepare phenotype data ####

## Predict the TAstar
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



#### 1.1 GWAS study ####
# GWAS_cont <- apply(genotype_QCed, MARGIN = 2, FUN = function(SNP){
#   logitreg <- lm(Y1star~SNP)
#   Zscore <- coef(logitreg)[2]/sqrt(vcov(logitreg)[2, 2])
#   pvalue <- 2*(pnorm(-abs(Zscore)))
#   return(pvalue)
# })
# 
# GWAS_cont_Sort <- sort(GWAS_cont)
# 
# GWAS_bin <- apply(genotype_QCed, MARGIN = 2, FUN = function(SNP){
#   logitreg <- glm(Y2star~SNP, family = binomial())
#   Zscore <- coef(logitreg)[2]/sqrt(vcov(logitreg)[2, 2])
#   pvalue <- 2*(pnorm(-abs(Zscore)))
#   return(pvalue)
# })
# 
# GWAS_bin_Sort <- sort(GWAS_bin)
# 
# ## Select the top 25th genes appeared in each study
# Genelist <- unique(c(names(GWAS_cont_Sort)[1:50]))

# SNPlist <- c("id", "discard", Genelist)

SNPlist <- c("id","discard","rs45690064","rs27338905","rs32962338","rs33583459","rs224051056",
             "rs33217671","rs38916331","rs47869247","rs217439518","rs29477109",
             "rs252503010","rs265727287","rs246035173","rs231489766","rs46826545",
             "rs51809856","rs6279141","rs30535702","rs30201629","rs30549753",
             "rs30793692", "rs243608221", "rs27323259","rs31194230","rs51851360",
             "rs28316536","rs50703583","rs223979909","rs37116508","rs264716939","rs237368278",
             "rs33100460","rs48556900","rs37861542","rs46497021","rs49711091",
             "rs230308064","rs245357151","rs214901846","rs254351625","rs244502760","rs47083137","rs36724404",
             "cfw.9.114048825",
             "rs29426250", "rs29377917", "rs215894093", "cfw.2.10783680", "rs27100804", "rs48290901","rs215878200", "cfw.12.83303849", "cfw.4.7251197", "rs29473466", 
             "rs13462773", "rs47152068","rs29406933","rs33030063", "cfw.14.58744678", "rs32015836", "rs236499396", "rs245594080"
)

SNPlist <- unique(SNPlist)

genotypes <- genotype[,SNPlist]

save(genotypes,file="file/genotypes2.RData")




#### 2. Second Stage ####

## 2.1 Deciding covariates dependence structure ####

### 



### Standardize the candidate genotypes 

load("file/genotypes2.RData")

scaled.gene <- scale(genotypes[,3:dim(genotypes)[2]])

# keepindex <- colSums(cor(scaled.gene)>0.85)<=6
# 
# scaled.gene.QC <- scaled.gene[keepindex,keepindex]



varsel <- huge(scaled.gene, lambda = NULL, nlambda = 30, lambda.min.ratio = NULL, method = "mb",
               scr = F, scr.num = NULL, sym = "or", verbose = TRUE, cov.output =T)
out.select = huge.select(varsel, criterion = "stars", stars.thresh = 0.05,rep.num=30)


THETAgraph <- varsel$path[[out.select$opt.index]]

EdgeTrue <- NULL  
EdgeHash <- rep(T,dim(THETAgraph)[1])

for (i in 1:(dim(THETAgraph)[1]-1)){
  for (j in (i+1):dim(THETAgraph)[1]){
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
  colnames(DesMatrix)[dim(DesMatrix)[2]] <- paste(colnames(DesMatrix)[EdgeTrue[i,1]],'x',colnames(DesMatrix)[EdgeTrue[i,2]])
}


# DesMatrix <- scale(DesMatrix)


#### plot the selected graph
plot.selectgraph = function(x, ...){
  if(x$cov.input){
    cat("Model selection is not available when using the covariance matrix as input.")
  }
  if(!x$cov.input)
  {
    # par(mfrow=c(1,2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.5,0.5,0.5,0.5))
    par(mfrow=c(1,2))
    g = graph.adjacency(as.matrix(x$refit), mode="undirected", diag=FALSE)
    # layout.grid = layout.fruchterman.reingold(g)
    # layout.grid = layout.reingold.tilford(g)
    layout.grid = layout.kamada.kawai(g)
    # layout.grid = layout.lgl(g)
    
    plot(g, layout=layout.grid, edge.color='gray50',vertex.color="white", vertex.size=0.1, vertex.label=SNPlist[-1*1:2], vertex.label.dist=0.5,vertex.label.degree =pi)
    plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameter", ylab = "Sparsity Level", type = "l",xlim = rev(range(x$lambda)), main = "Solution path sparsity levels")
    lines(x$opt.lambda,x$opt.sparsity,type = "p")
  }
}

plot.selectgraph.ggplot2 = function(x, ...){
  if(x$cov.input){
    cat("Model selection is not available when using the covariance matrix as input.")
  }
  if(!x$cov.input)
  {
    # par(mfrow=c(1,2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.5,0.5,0.5,0.5))
    # par(mfrow=c(1,2))
    g = graph.adjacency(as.matrix(x$refit), mode="undirected", diag=FALSE)
    # layout.grid = layout.fruchterman.reingold(g)
    # layout.grid = layout.reingold.tilford(g)
    # layout.grid = layout.kamada.kawai(g)
    # layout.grid = layout.lgl(g)
    gmatrix <- as_edgelist(g)
    gdegree <- degree(g, mode ="all")
    
    SNPlistpure <- SNPlist[-c(1,2)]
    
    gdegreeMatrix <- data.frame(SNPlistpure, gdegree)
    
    
    edges <- as.data.frame(cbind(from_id = SNPlistpure[gmatrix[,1]], to_id = SNPlistpure[gmatrix[,2]]))
    vertices <- as.data.frame(SNPlistpure)
    MapObj <- fortify(as.edgedf(edges), vertices)
    
    MapObjd <- merge(MapObj, gdegreeMatrix, by.x = "from_id", by.y = "SNPlistpure" )
    
    MapObjd <- MapObjd %>% mutate(degree = sqrt(4 * gdegree + 4.5)) %>%
      mutate(highlight1 = gdegree > 0) %>%
      mutate(highlight2 = case_when(
        gdegree == 0 ~ 0,
        gdegree < 4 ~ 1,
        TRUE ~ 2))
    
    MapObjd$highlight2 <- factor(MapObjd$highlight2, levels = 0:2, labels = c("0","<=4",">4"))
    
    p1 <- ggplot(data = MapObjd, aes(from_id = from_id, to_id = to_id)) +
      geom_net(layout.alg = "kamadakawai", 
               aes(fontsize = degree, color = highlight2),
               size = 2, labelon = TRUE, vjust = -0.6, ecolour = "#C1AE8D",
               directed =FALSE, ealpha = 0.4) +
      scale_colour_manual(values = c("#7A989A", "#C1AE8D", "#C67052")) +
      xlim(c(-0.1, 1.05)) +
      labs(color = "degree") +
      theme_net() +
      theme(legend.position = "bottom")
    
    datalineplot <-  data.frame(lambda = x$lambda,sparsity = x$sparsity)
    
    p2 <- ggplot(data=datalineplot) +
      geom_line(aes(x=lambda, y=sparsity), color = "#849271",size = 1.2)+
      geom_point(x = x$opt.lambda, y = x$opt.sparsity, color = "#CF9546",size = 3) +
      labs(x = "Regularization Parameter", y = "Sparsity Level") + 
      theme_bw() +
      theme(text=element_text(size=12, family="mono"), axis.text = element_text(size = 14, family="mono"), panel.spacing = unit(1, "lines")) 
    
    
    p1 + p2 + plot_layout(widths = c(2, 1.2))
  }
}



pdf(file = "output/figdatacov2.pdf",height = 6, width = 14)
set.seed(2020)
plot.selectgraph.ggplot2(out.select)
dev.off()


## Mark the complete sample

# Response
Y1star <- data_GWAS$Y1star
Y2star <- data_GWAS$Y2star

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



Covariates1 <- Covariates2 <- as.matrix(cbind(Covariates_all[,! names(Covariates_all) %in% c("Y1star","Y2star","BW")])) 

CovMis1 <- as.matrix(data.frame(intercept=rep(1,nsample),cov=rep(1,nsample)))

CovMis2 <- as.matrix(data.frame(intercept=rep(1,nsample)))


completedata <- as.matrix(cbind(Covariates_all[,! names(Covariates_all) %in% c("BW")])) 

#### 3. Analysis  of 109 ####

## 3.1 Naive Analysis of 109 ####


## 3.2 Estimating Procedure ####

## 3.2.1 Naive model ####
naive.model1 <- lm(Y1star ~ -1+ ., data = data.frame(completedata[,-2]))
naive.model2 <- glm(Y2star ~ -1+ ., family = binomial(link = logit), data = data.frame(completedata[,-1]))

summary(naive.model1)
summary(naive.model2)


initial2 <- c(naive.model1$coefficients,naive.model2$coefficients,0.03,0)
cat("Initial value:", initial2, "\n")



## Display the most basic cases (10% misclassification, 0.56 measurment error)

naivebeta <- c(coef(naive.model1),coef(naive.model2))
naivesd <- c(sqrt(diag(vcov(naive.model1))),
             sqrt(diag(vcov(naive.model2))))

naiveZvalue <- naivebeta/naivesd
naivepvalue <- 2*(pnorm(-abs(naiveZvalue)))
naivepscale <- -log(2*(pnorm(-abs(naiveZvalue))), base = 10)
# pscale <- -log(2*(1-pnorm(abs(Zvalue))))

SNPnames <- names(naivebeta)
# SNPnames[(length(SNPnames)-1):length(SNPnames)] <- c("sigma","phi")

Table_numeric <- data.frame(SNPnames = SNPnames,
                            propobeta = naivebeta,
                            proposd = naivesd,
                            propoZ = naiveZvalue,
                            propop = naivepvalue,
                            pscale = naivepscale,
                            method = "naive",
                            type = rep(c("continuous", "discrete"),each = length(SNPnames)/2))

Tablemain <- SensiAnalysis(-2.197,0.60)  
Tablemains <- Tablemain %>%
  filter(!SNPnames %in% c("sigma", "phi")) %>%
  mutate(method = "Proposed") 
Tablemains <- Tablemains %>% 
  mutate(type = rep(c("continuous", "discrete"),each = dim(Tablemains)[1]/2))


## Create the table for ploting the results
Tablefirst <- rbind(Table_numeric,Tablemains)
Tablefirst$type <- factor(Tablefirst$type, levels=c("continuous", "discrete"))
Tablefirst$method <- factor(Tablefirst$method, levels=c("naive", "Proposed"),
                            labels = c("Naive", "Proposed"))



# print(TablefirstPara[c("SNPnames","propobeta","HLparaorder")], n=100)


pdf(file = "output/figresultsmain1.pdf",height = 7, width = 14)
set.seed(2019)

require("ggrepel")

TablefirstPara <- Tablefirst %>% 
  mutate(HLPara = ifelse(abs(propobeta)>1.5,"Highlight","Not Highlight")) %>%
  mutate(HLParalabel = ifelse(HLPara=="Highlight",SNPnames,"")) %>%
  mutate(HLParalabel = gsub("[.]x[.]", " x ", HLParalabel)) %>%
  mutate(HLParalabel = gsub("[.]", "-", HLParalabel)) 



p1 <- ggplot(TablefirstPara, aes(x=type,y=propobeta)) +
  geom_point(aes(color = HLPara), size =2) +
  facet_wrap(.~method) + 
  geom_text_repel(aes(label = HLParalabel), color="#C1AE8D", fontface=2, alpha = 0.9) +
  scale_color_manual(values = c("#C67052", "#7A989A")) +
  theme_bw() + 
  theme(text=element_text(size=12, family="mono"), axis.text = element_text(size = 14, family="mono"),
        legend.position = "none") +
  labs(x = "Data Type", y = "Parameter Estimate") +
  labs(title="(a)") +
  theme(strip.background =element_rect(fill="#4F534A",color="#4F534A"))+ # #535b44
  theme(strip.text = element_text(colour = 'white')) +
  theme(panel.border = element_rect(colour = "#4F534A")) 

p1
dev.off()


pdf(file = "output/figresultsmain2.pdf",height = 6, width = 14)
set.seed(2019)

TablefirstPvalue <- Tablefirst %>% 
  mutate(HLPvalue = ifelse(abs(pscale)>3.6,"Highlight","Not Highlight")) %>%
  mutate(HLPvaluelabel = ifelse(HLPvalue=="Highlight",SNPnames,""))  %>%
  mutate(HLPvaluelabel = gsub("[.]x[.]", " x ", HLPvaluelabel)) %>%
  mutate(HLPvaluelabel = gsub("[.]", "-", HLPvaluelabel)) 

p2 <- ggplot(TablefirstPvalue, aes(x=type,y=pscale)) +
  geom_point(aes(color = HLPvalue), size =2) +
  facet_wrap(.~method) + 
  geom_text_repel(aes(label = HLPvaluelabel), color="#C1AE8D", fontface=2, alpha = 0.9) +
  scale_color_manual(values = c("#C67052", "#7A989A")) +
  theme_bw() + 
  theme(text=element_text(size=12, family="mono"), axis.text = element_text(size = 14, family="mono"),
        legend.position = "none") +
  labs(x = "Data Type", y = bquote("-"~log[10]~"P-value")) +
  labs(title="(b)") +
  theme(strip.background =element_rect(fill="#4F534A",color="#4F534A"))+ # #535b44
  theme(strip.text = element_text(colour = 'white')) +
  theme(panel.border = element_rect(colour = "#4F534A"))

p2

dev.off()


# ### Create the foundation of the graph data.frame
# g = graph.adjacency(as.matrix(out.select$refit), mode="undirected", diag=FALSE)
# gmatrix <- as_edgelist(g)
# gdegree <- degree(g, mode ="all")
# 
# SNPlistpure <- SNPlist[-c(1,2)]
# gdegreeMatrix <- data.frame(SNPlistpure, gdegree)
# 
# edges <- as.data.frame(cbind(from_id = SNPlistpure[gmatrix[,1]], to_id = SNPlistpure[gmatrix[,2]]))
# vertices <- as.data.frame(SNPlistpure)
# MapObj <- fortify(as.edgedf(edges), vertices)
# 
# MapObjd <- merge(MapObj, gdegreeMatrix, by.x = "from_id", by.y = "SNPlistpure" )
# 
# MapObjd <- MapObjd %>% mutate(degree = sqrt(4 * gdegree + 4.5)) %>%
#   mutate(highlight1 = gdegree > 0) %>%
#   mutate(highlight2 = case_when(
#     gdegree == 0 ~ 0,
#     gdegree < 4 ~ 1,
#     TRUE ~ 2))
# 
# MapObjd$highlight2 <- factor(MapObjd$highlight2, levels = 0:2, labels = c("0","<=4",">4"))


SNPNameSplit <- str_split_fixed(Tablefirst$SNPnames, "[.]x[.]", 2)
colnames(SNPNameSplit) <- c("SNPname1", "SNPname2")

Tablefirstat <- cbind(Tablefirst, SNPNameSplit) 

# MapObjd_atv <- merge(MapObjd %>%
#                        mutate(to_id2=ifelse(is.na(to_id),"",to_id)) %>%
#                        select(c(-"to_id")),
#                      Tablefirstat %>%
#                        mutate(pscalev = pscale) %>%
#                        select(c("SNPname1","SNPname2","pscale","method","type")),
#                      by.x = c("from_id","to_id2"), by.y = c("SNPname1","SNPname2") )

Tablefirstat$SNPname1 <- gsub("[.]", "-", Tablefirstat$SNPname1)
Tablefirstat$SNPname2 <- gsub("[.]", "-", Tablefirstat$SNPname2)




pdf(file = "output/figresultsmain3.pdf",height = 6, width = 14)
set.seed(2018)

datasetg1 <- EdgeBundling(Tablefirstat, outcomearg = "continuous", methodarg = "Proposed")

edge.pscale <- datasetg1$connect$pscale

g1 = ggraph(datasetg1$mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, size= pscale), alpha = 0.8, color = "#7A989A") +
  geom_conn_bundle(data = get_con(from = datasetg1$from, to = datasetg1$to, edge.pscale=edge.pscale), aes(width = edge.pscale), colour= "#C67052", alpha=0.3) +
  geom_node_text(aes(x = x*1.25, y=y*1.25, filter = leaf, label=name, angle = angle, hjust=hjust), size=2.3, alpha=1, color = "#849271") +
  # scale_edge_colour_gradient(low="#CF9546", high="#C67052") +
  theme_void() +
  theme(
    legend.position="right",
    plot.margin=unit(c(0,0,0,0),"cm"),
    text=element_text(size=12, family="mono"),
  ) + 
  labs(title="(c) Continuous") +
  labs(size = bquote(atop("-"~log[10]~"P-value","of vertix")), edge_width = bquote(atop("-"~log[10]~"P-value","of edge"))) +
  guides(size = guide_legend(title.theme = element_text( size = 9 )),
         edge_width = guide_legend(title.theme = element_text( size = 9 ))) +
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))

datasetg2 <- EdgeBundling(Tablefirstat, outcomearg = "discrete", methodarg = "Proposed")

edge.pscale2 <- datasetg2$connect$pscale

g2 = ggraph(datasetg2$mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, size= pscale), alpha = 0.8, color = "#7A989A") +
  geom_conn_bundle(data = get_con(from = datasetg2$from, to = datasetg2$to, edge.pscale2=edge.pscale2), aes(width = edge.pscale2), colour= "#CF9546", alpha=0.3) +
  geom_node_text(aes(x = x*1.25, y=y*1.25, filter = leaf, label=name, angle = angle, hjust=hjust), size=2.3, alpha=1, color = "#849271") +
  # scale_edge_colour_gradient(low="#CF9546", high="#C67052") +
  theme_void() +
  theme(
    legend.position="right",
    plot.margin=unit(c(0,0,0,0),"cm"),
    text=element_text(size=12, family="mono"),
  ) + 
  labs(title="   Discrete") +
  labs(size = bquote(atop("-"~log[10]~"P-value","of vertix")), edge_width = bquote(atop("-"~log[10]~"P-value","of edge"))) +
  guides(size = guide_legend(title.theme = element_text( size = 9 )),
         edge_width = guide_legend(title.theme = element_text( size = 9 ))) +
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))

g1 + g2 +  plot_layout(widths = c(1, 1))
dev.off()


## Low Measurement Error
Tablelow0.52 <- SensiAnalysis(-2.944,0.52) 
Tablemedian0.52 <- SensiAnalysis(-2.197,0.52) 
Tablehigh0.52 <- SensiAnalysis( -1.386,0.52) 

Table1 <- integrate2(Tablelow0.52, Tablemedian0.52, Tablehigh0.52, title = c("5%","10%","20%"))

pdf(file = "output/figresults520.pdf",height = 12, width = 12)
plotintegrate2(Table1)
dev.off()


# Moderate Measurement Error

Tablelow0.60 <- SensiAnalysis(-2.944,0.60) 
Tablemedian0.60 <- SensiAnalysis(-2.197,0.60) 
Tablehigh0.60 <- SensiAnalysis( -1.386,0.60) 

Table2 <- integrate2(Tablelow0.60, Tablemedian0.60, Tablehigh0.60, title = c("5%","10%","20%"))
pdf(file = "output/figresults060.pdf",height = 12, width = 12)
plotintegrate2(Table2)
dev.off()


# Higher Measurement Error

Tablelow0.67 <- SensiAnalysis(-2.944,0.67) 
Tablemedian0.67 <- SensiAnalysis(-2.197,0.67) 
Tablehigh0.67 <- SensiAnalysis( -1.386,0.67) 

Table3 <- integrate2(Tablelow0.67, Tablemedian0.67, Tablehigh0.67, title = c("5%","10%","20%"))
pdf(file = "output/figresults067.pdf",height = 12, width = 12)
plotintegrate2(Table3)
dev.off()


