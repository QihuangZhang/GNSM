####
# This file is set up to summarize the result data of the simulation into Table2

#### 0 Simulation set up####
library(parallel)
library(scales)
library(xtable)

#### 1 Load Data ####

load(file="results/Table2/Table2_KNO.RData")
dataknow <- results_agg2

load(file="results/Table2/Table2_INV.RData")
dataintv <- results_agg2

load(file="results/Table2/Table2_EXV.RData")
dataextv <- results_agg2



#### 2 Combine the data ####

compdata <- cbind(dataknow,dataintv[,9:12],dataextv[,9:12])
colnames(compdata)[1:4] <- c("Type","graph","alpha","Sigma_e")

#### 3: Table 2  ####

Table2 <- compdata[compdata$Sigma_e==0.2,]
xtable(Table2, digits=3)

#### 4: Table 3  ####

Table3 <- compdata[compdata$Sigma_e==0.7,]
xtable(Table3, digits=3)



