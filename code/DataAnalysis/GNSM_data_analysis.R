### Project GNMM
## As advised by referee, we consider a higher dimension of the SNP candidates.
## A sensitivity study is conducted to evaluate the performance of the method.

#### 0 Project Information and Configuration ####
ProjectName<- paste0("DataAnalysis")
set.seed(2019)

## 0.1 Set file path ####
# setwd("N:/Work/waterloo/2019/GNMM")
# WD=paste0(getwd(),"/results/EMsensi/")

## 0.2 Prepare package of GWAS toolbox ####
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

source("GNSM_data_analysis_fun.R")


#### 1. First Step - Data Preparation ####

## prepare phenotype data ####

phenotype<-read.csv(file="file/pheno.csv")

## Predict the TAstar

BMD90 <- quantile(phenotype$BMD,na.rm=T, probs = 0.9)
phenotype$BCstar <- ifelse(phenotype$BMD<BMD90,0,1)

## Association Data Frame
data_GWA_Pre <- data.frame(id = phenotype$id,Y1 = phenotype$tibia,Y2 = phenotype$BCstar,
                           discard = ifelse(phenotype$discard=="no",1,0),
                           mixup = ifelse(phenotype$mixup=="no",0,1), 
                           BMD=phenotype$BMD, Y1star=phenotype$tibia, Y2star=phenotype$abnormalbone)


## load proprocessed genotype data ####

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

# genotypes <- genotype[,SNPlist]
# 
# save(genotypes,file="file/genotypes2.RData")

load("file/genotypes2.RData")

data_GWAS_unsorted <- data_GWA_Pre[data_GWA_Pre$id %in% genotypes$id,]

data_GWAS <- data_GWAS_unsorted[order(data_GWAS_unsorted$id),]

# Response
Y1star <- data_GWAS$Y1star

Y2star <- data_GWAS$Y2star


#### 2. Second Stage ####

## 2.1 Deciding covariates dependence structure ####

### Standardize the candidate genotypes 

scaled.gene <- scale(genotypes[,3:dim(genotypes)[2]])

### Fit Gaussian Graphical Model

varsel <- huge(scaled.gene, lambda = NULL, nlambda = 30, lambda.min.ratio = NULL, method = "mb",
               scr = F, scr.num = NULL, sym = "or", verbose = TRUE, cov.output =T)
out.select = huge.select(varsel, criterion = "stars", stars.thresh = 0.05,rep.num=30)


## Save the optimally selected graph structure

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

### Adding the interaction terms into covariate data according to the optimal selected graph

for (i in 1:dim(EdgeTrue)[1]){
  DesMatrix <- cbind(DesMatrix,DesMatrix[,EdgeTrue[i,1]]*DesMatrix[,EdgeTrue[i,2]])
  colnames(DesMatrix)[dim(DesMatrix)[2]] <- paste(colnames(DesMatrix)[EdgeTrue[i,1]],'x',colnames(DesMatrix)[EdgeTrue[i,2]])
}

set.seed(2020)
figure3 <- plot.selectgraph.ggplot2(out.select)

pdf(file = "output/figdatacov21.pdf",height = 6, width = 6)
print(figure3[[1]])
dev.off()

pdf(file = "output/figdatacov22.pdf",height = 6, width = 7)
print(figure3[[2]])
dev.off()

### Data clearning: Only keep the sample with complete data

comp<-data.frame(Y1star,Y2star,DesMatrix)
index.notna<-complete.cases(comp)
index.notna <- ifelse(data_GWAS$discard == 1, index.notna , FALSE)

Y1star<-Y1star[index.notna]
Y2star<-Y2star[index.notna]

Covariates_all<-comp[index.notna,]
nsample<-dim(Covariates_all)[1]

Covariates1 <- Covariates2 <- as.matrix(cbind(Covariates_all[,! names(Covariates_all) %in% c("Y1star","Y2star","BW")])) 

CovMis1 <- as.matrix(data.frame(intercept=rep(1,nsample),cov=rep(1,nsample)))

CovMis2 <- as.matrix(data.frame(intercept=rep(1,nsample))) 

completedata <- as.matrix(cbind(Covariates_all[,! names(Covariates_all) %in% c("BW")])) 


#### 3. Analysis ####

## 3.1 Naive model ####
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

## 3.2 Proposed Methods ####


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




#### 4. Produce Visualization Results ####

## Produce figure 4(a)

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
  theme(strip.text = element_text(colour = 'white', size =16)) +
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




SNPNameSplit <- str_split_fixed(Tablefirst$SNPnames, "[.]x[.]", 2)
colnames(SNPNameSplit) <- c("SNPname1", "SNPname2")

Tablefirstat <- cbind(Tablefirst, SNPNameSplit) 


Tablefirstat$SNPname1 <- gsub("[.]", "-", Tablefirstat$SNPname1)
Tablefirstat$SNPname2 <- gsub("[.]", "-", Tablefirstat$SNPname2)


## Produce figure 4(d)


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


## Produce figure 4(c)


pdf(file = "output/figresultsmain4.pdf",height = 6, width = 14)
set.seed(2018)

datasetg1 <- EdgeBundling(Tablefirstat, outcomearg = "continuous", methodarg = "Naive")

edge.pscale <- datasetg1$connect$pscale

g1 = ggraph(datasetg1$mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, size= pscale), alpha = 0.8, color = "#C67052") +
  geom_conn_bundle(data = get_con(from = datasetg1$from, to = datasetg1$to, edge.pscale=edge.pscale), aes(width = edge.pscale), colour= "#7A989A", alpha=0.3) +
  geom_node_text(aes(x = x*1.25, y=y*1.25, filter = leaf, label=name, angle = angle, hjust=hjust), size=2.3, alpha=1, color = "#C1AE8D") +
  # scale_edge_colour_gradient(low="#CF9546", high="#C67052") +
  theme_void() +
  theme(
    legend.position="right",
    plot.margin=unit(c(0,0,0,0),"cm"),
    text=element_text(size=12, family="mono"),
  ) + 
  labs(title="(b) Continuous") +
  labs(size = bquote(atop("-"~log[10]~"P-value","of vertix")), edge_width = bquote(atop("-"~log[10]~"P-value","of edge"))) +
  guides(size = guide_legend(title.theme = element_text( size = 9 )),
         edge_width = guide_legend(title.theme = element_text( size = 9 ))) +
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))

datasetg2 <- EdgeBundling(Tablefirstat, outcomearg = "discrete", methodarg = "Naive")

edge.pscale2 <- datasetg2$connect$pscale

g2 = ggraph(datasetg2$mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, size= pscale), alpha = 0.8, color = "#C67052") +
  geom_conn_bundle(data = get_con(from = datasetg2$from, to = datasetg2$to, edge.pscale2=edge.pscale2), aes(width = edge.pscale2), colour= "#849271", alpha=0.3) +
  geom_node_text(aes(x = x*1.25, y=y*1.25, filter = leaf, label=name, angle = angle, hjust=hjust), size=2.3, alpha=1, color = "#C1AE8D") +
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



### Produce the results in the Supplementary Materials ###

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


