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
  
  Table_numeric <- data.frame(SNPnames = c(SNPnames, "sigma", "phi"),
                              propobeta=betahat,
                              proposd = betaIsd_proposed,
                              propoZ = Zvalue,
                              propop = pvalue,
                              pscale = pscale)
}

#### plot the selected graph
plot.selectgraph = function(x, ...){
  if(x$cov.input){
    cat("Model selection is not available when using the covariance matrix as input.")
  }
  if(!x$cov.input)
  {
    par(mfrow=c(1,2))
    g = graph.adjacency(as.matrix(x$refit), mode="undirected", diag=FALSE)
    layout.grid = layout.kamada.kawai(g)
    
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
    g = graph.adjacency(as.matrix(x$refit), mode="undirected", diag=FALSE)
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
               size = 2, labelon = TRUE, vjust = -0.6, ecolour = "#7A989A",
               directed =FALSE, ealpha = 0.4) +
      scale_colour_manual(values = c("#626567","#7A989A", "#C67052")) +
      xlim(c(-0.1, 1.05)) +
      labs(color = "degree") +
      theme_net() +
      theme(legend.position = "bottom")
    
    datalineplot <-  data.frame(lambda = x$lambda,sparsity = x$sparsity)
    
    p2 <- ggplot(data=datalineplot) +
      geom_line(aes(x=lambda, y=sparsity), color = "#7A989A",size = 1.2)+
      geom_point(x = x$opt.lambda, y = x$opt.sparsity, color = "#C67052",size = 5) +
      labs(x = "Regularization Parameter", y = "Sparsity Level") + 
      theme_bw() +
      theme(text=element_text(size=12, family="mono"), axis.text = element_text(size = 14, family="mono"), panel.spacing = unit(1, "lines")) 
    
    return(list(p1=p1,p2=p2))
  }
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
    theme(strip.text = element_text(colour = 'white', size = 16)) +
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
    theme(strip.text = element_text(colour = 'white', size = 16)) +
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
