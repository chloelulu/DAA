## Script for extracting data 
pkg = c('dplyr','tidyr','reshape','RColorBrewer','ggplot2','ggpubr','stringr','scales','tibble','kableExtra','formattable','reactable','htmltools',"htmlwidgets","webshot2")
suppressPackageStartupMessages(sapply(pkg, require, character = T))
sub.methods = c('GMPR+Wilcox', 'GMPR+DESeq2','GMPR+edgeR',"corncob","eBay(Wilcox)",
                'Wrench+MSeq',"RAIDA","ANCOM-BC","DACOMP","LDM",
                "Omnibus",'MaAsLin2',"GMPR+glm","Aldex2(Wilcox)")
source('~/Documents/Mayo_project/Code/DailyCode.R')
source('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/code/SimulationEvaluation/func.R')

## D4
setwd('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationD4/')
files <- list.files(pattern = 'D4_summarymatrix1-')
length(files)
files[1:5]

load(paste0('D4_summarymatrix',1, '-Maaslin2.Rdata'))
nmss <- names(res_seqs)
nms <- gsub('Maaslin2','',nmss)
nms <- nms[grep('binary',nms)]
# nms <- c('L2low','L2high','L4low','L4high')

methods <- c('glmquassi2','Maaslin2','Rarefy', 'Wilcox','Wilcox.Wrench','Wilcox.gmpr','mbzinb','eBayW','eBayt','Aldex2we','Aldex2','RAIDA', 'DACOMP', 'LDM', 'glmquassi','BBinomial','ANCOMBC','DESeq2', 'DESeq2.Wrench', 'DESeq2.gmpr','edgeR', 'edgeR.Wrench', 'edgeR.gmpr','MSeq2','MSeq2.Wrench')
RES <- list()
for(nm in nms){
  for(i in 1:100){
    for(method in methods){
      try({
        load(paste0('D4_summarymatrix',i, '-',method,'.Rdata'))
        res <- res_seqs[[paste0(nm,method)]][,c('otu.id','fdr')]
        colnames(res)[2] <- method 
        head(res)
        RES[[paste0(method)]] <- res
      })
    }
    ## each rdata is each dataset with different method
    x = gsub('D3','',nm)
    x = gsub('loglinearSub_L1','',x)
    x = gsub('none','',x)
    x = gsub('nSam_L2L3','',x)
    save(RES, file = paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationD4/iter',i,'-',x,'.rdata'))
  }
}


for(nm in nms){
  for(i in 1:100){
    try({
      load(paste0('D4_summarymatrix',i, '-Wilcox.Rdata'))
      diff <- res_seqs[[paste0(nm,'Wilcox')]][,c('otu.id','diff.otu.ind')]
      x = gsub('D3','',nm)
      x = gsub('loglinearSub_L1','',x)
      x = gsub('none','',x)
      x = gsub('nSam_L2L3','',x)
      
      colnames(diff)[2] <- paste0(i,'_',x) 
    })
    save(diff, file = paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationD4/diff_iter',i,'-',x,'.rdata'))
  }
}



setwd('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationD4/')


nms <- list.files(pattern = 'rdata$')
nms <- gsub('.*-','',nms)
nms <- gsub('.rdata','',nms) %>% unique()

TT <- NULL
for(class in nms){
  files <- list.files(pattern = 'rdata$')
  files <- files[grep('^iter',files)]
  
  files <- files[grep(class, files)]
  
  
  TPRRs = FDRRs = NULL
  
  for(i in 1:length(files)){
    load(paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationD4/',files[i]))
    res <- RES[[1]];for(j in 2:length(RES)){res <- full_join(res, RES[[j]])}
    res <- res %>% column_to_rownames('otu.id')
    res[is.na(res)] <- 1
    dim(res)
    
    ## common findings
    colnames(res)[colnames(res) =='Wilcox'] = 'TSS+Wilcox'
    colnames(res)[colnames(res) =='Rarefy'] = 'Rarefy+Wilcox'
    colnames(res) = gsub('ANCOMBC','ANCOM-BC',colnames(res))
    colnames(res) = gsub('glmquassi2','GMPR+glm',colnames(res))
    colnames(res) = gsub('glmquassi','Wrench+glm',colnames(res))
    colnames(res) = gsub('DESeq2.gmpr','GMPR+DESeq2',colnames(res))
    colnames(res) = gsub('DESeq2.Wrench','Wrench+DESeq2',colnames(res))
    colnames(res) = gsub('edgeR.gmpr','GMPR+edgeR',colnames(res))
    colnames(res) = gsub('edgeR.Wrench','Wrench+edgeR',colnames(res))
    colnames(res) = gsub('MSeq2.Wrench','Wrench+MSeq',colnames(res))
    colnames(res) = gsub('MSeq2','MSeq',colnames(res))
    colnames(res) = gsub('eBayW','eBay(Wilcox)',colnames(res))
    colnames(res) = gsub('eBayt','eBay(t-test)',colnames(res))
    colnames(res) = gsub('BBinomial','corncob',colnames(res))
    colnames(res) = gsub('Aldex2we','Aldex2(t-test)',colnames(res))
    colnames(res) = gsub('^Aldex2$','Aldex2(Wilcox)',colnames(res))
    colnames(res) = gsub('^mbzinb','Omnibus',colnames(res))
    colnames(res) = gsub('^Wilcox.Wrench','Wrench+Wilcox',colnames(res))
    colnames(res) = gsub('^Wilcox.gmpr','GMPR+Wilcox',colnames(res))
    colnames(res) = gsub('Maaslin2','MaAsLin2',colnames(res))
    
    
    len <- sub.methods[!(sub.methods %in% colnames(res))]
    if(length(len) > 0){
      for(k in 1:length(len)){
        res[,len[k]] <- rep(1, nrow(res))
      }
    }
    cat(files[i],'_',len,'\n')
    
    # res[is.na(res)] <- 1
    # res[grep('NA',res)] <- 1
    
    res <- res[,sub.methods]
    
    load(paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationD4/','diff_',files[i]))
    colnames(diff)[2] = 'diff.otu.id'
    
    res <- full_join(res %>% rownames_to_column('otu.id'),diff) %>% column_to_rownames('otu.id')
    
    head(res)
    
    
    
    truth <- apply(res[,-ncol(res), drop = F],1, function(x) sum(x[!is.na(x)] <= 0.05)/(ncol(res)-1)) 
    truth <- as.data.frame(cbind(sapply(truth, function(x) x > 0.2),sapply(truth, function(x) x > 0.4),sapply(truth, function(x) x > 0.6),sapply(truth, function(x) x > 0.8)))
    colnames(truth) <- c('pct20','pct40','pct60','pct80')
    truth <- apply(truth, 2, function(x) ifelse(x==TRUE, 0, 1))
    
    res <- full_join(as.data.frame(truth) %>% rownames_to_column('otu.id'), res %>% rownames_to_column('otu.id')) %>% column_to_rownames('otu.id')
    
    
    pct.cal <- function(res, i){
      tp <- sum(res[,i] <= 0.05 & res$diff.otu.id ==TRUE, na.rm = T)
      tn <- sum(res[,i] > 0.05 &  res$diff.otu.id ==FALSE, na.rm = T)
      fn <- sum(res[,i] > 0.05 &  res$diff.otu.id ==TRUE, na.rm = T)
      fp <- sum(res[,i] <= 0.05 &  res$diff.otu.id ==FALSE, na.rm = T)
      tpr <- ifelse(!tp, 0, tp/(tp+fn))
      fdr <- ifelse(!fp, 0, fp/(tp+fp))
      return(c(tpr, fdr))
    }
    
    tt <- name <- NULL
    for(m in 1:(ncol(res)-1)){
      tt <- rbind(tt, pct.cal(res,m))
      name <- c(name, colnames(res)[m])
    }
    
    rownames(tt) <- name
    colnames(tt) <- c('tpr','fdr')
    tt <- as.data.frame(tt) %>% rownames_to_column('method')
    tt$class <- class
    tt$iter <- i
    TT <- rbind(TT, tt)
  }
}

save(TT, file = '~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationD4/RES_TT.RData')

load('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationD4/RES_TT.RData')
head(TT)
TT$diff.otu.pcts <- gsub('.*balanced','',TT$class)
TT$diff.otu.modes <- gsub('nOTU.*','',TT$class)

TT$nOTUs <- gsub('abundant','',TT$class)
TT$nOTUs <- gsub('rare','',TT$nOTUs)
TT$nOTUs <- gsub('binary.*','',TT$nOTUs)
TT$covariate.types <- 'binary'
TT$diff.otu.directs <- 'balanced'

TT <- TT %>% dplyr::select(-class)
head(TT)

TT1 <- TT %>% dplyr::select(-tpr) %>% mutate(measures ='FDR')
TT2 <- TT %>% dplyr::select(-fdr) %>% mutate(measures ='TPR')
colnames(TT1)[2] <- colnames(TT2)[2] <- 'value'
TT <- rbind(TT1, TT2)
colnames(TT)[1] = 'methods'
head(TT)

TT[,'nOTUs'] <- as.factor(TT[,'nOTUs'])
TT$diff.otu.pcts <- as.factor(TT$diff.otu.pcts)
TT$diff.otu.modes <- as.factor(TT$diff.otu.modes)
TT$covariate.types <- as.factor(TT$covariate.types)
TT$diff.otu.directs <- as.factor(TT$diff.otu.directs)
TT$measures <- as.factor(TT$measures)
str(TT)

levels(TT[,'nOTUs']) <- c('OTU=50','OTU=500')
levels(TT$diff.otu.modes) = str_to_title(levels(TT$diff.otu.modes))
levels(TT$diff.otu.pcts) = paste0(str_to_title(levels(TT$diff.otu.pcts)),' denisty')

head(TT)
res.df2 = data_summary(TT, as.formula(paste0('value ~ diff.otu.modes + nOTUs + diff.otu.pcts+ diff.otu.directs + methods + measures')))
unique(res.df2$methods)

save(res.df2, file = '~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/consensusD4.Rdata')









## C4
setwd('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationC4/')
files <- list.files(pattern = 'C4_summarymatrix1-')
length(files)
files[1:5]

load(paste0('C4_summarymatrix',1, '-Maaslin2.Rdata'))
nmss <- names(res_seqs)
nms <- gsub('Maaslin2','',nmss)
nms <- nms[grep('binary',nms)]
# nms <- c('L2low','L2high','L4low','L4high')

methods <- c('glmquassi2','Maaslin2','Rarefy', 'Wilcox','Wilcox.Wrench','Wilcox.gmpr','mbzinb','eBayW','eBayt','Aldex2we','Aldex2','RAIDA', 'DACOMP', 'LDM', 'glmquassi','BBinomial','ANCOMBC','DESeq2', 'DESeq2.Wrench', 'DESeq2.gmpr','edgeR', 'edgeR.Wrench', 'edgeR.gmpr','MSeq2','MSeq2.Wrench')
RES <- list()
for(nm in nms){
  for(i in 1:100){
    for(method in methods){
      try({
        load(paste0('C4_summarymatrix',i, '-',method,'.Rdata'))
        res <- res_seqs[[paste0(nm,method)]][,c('otu.id','fdr')]
        colnames(res)[2] <- method 
        head(res)
        RES[[paste0(method)]] <- res
      })
    }
    ## each rdata is each dataset with different method
    x = gsub('D3','',nm)
    x = gsub('loglinearSub_L1','',x)
    x = gsub('none','',x)
    x = gsub('nSam_L2L3','',x)
    save(RES, file = paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationC4/iter',i,'-',x,'.rdata'))
  }
}



for(nm in nms){
  for(i in 1:100){
    try({
      load(paste0('C4_summarymatrix',i, '-Wilcox.Rdata'))
      diff <- res_seqs[[paste0(nm,'Wilcox')]][,c('otu.id','diff.otu.ind')]
      x = gsub('D3','',nm)
      x = gsub('loglinearSub_L1','',x)
      x = gsub('none','',x)
      x = gsub('nSam_L2L3','',x)
      
      colnames(diff)[2] <- paste0(i,'_',x) 
    })
    save(diff, file = paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationC4/diff_iter',i,'-',x,'.rdata'))
  }
}


setwd('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationC4/')

nms <- list.files(pattern = 'rdata$')
nms <- gsub('.*-','',nms)
nms <- gsub('.rdata','',nms) %>% unique()

TT <- NULL
for(class in nms){
  files <- list.files(pattern = 'rdata$')
  files <- files[grep('^iter',files)]
  
  files <- files[grep(class, files)]
  
  
  TPRRs = FDRRs = NULL
  
  for(i in 1:length(files)){
    load(paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationC4/',files[i]))
    res <- RES[[1]];for(j in 2:length(RES)){res <- full_join(res, RES[[j]])}
    res <- res %>% column_to_rownames('otu.id')
    res[is.na(res)] <- 1
    dim(res)
    
    ## common findings
    colnames(res)[colnames(res) =='Wilcox'] = 'TSS+Wilcox'
    colnames(res)[colnames(res) =='Rarefy'] = 'Rarefy+Wilcox'
    colnames(res) = gsub('ANCOMBC','ANCOM-BC',colnames(res))
    colnames(res) = gsub('glmquassi2','GMPR+glm',colnames(res))
    colnames(res) = gsub('glmquassi','Wrench+glm',colnames(res))
    colnames(res) = gsub('DESeq2.gmpr','GMPR+DESeq2',colnames(res))
    colnames(res) = gsub('DESeq2.Wrench','Wrench+DESeq2',colnames(res))
    colnames(res) = gsub('edgeR.gmpr','GMPR+edgeR',colnames(res))
    colnames(res) = gsub('edgeR.Wrench','Wrench+edgeR',colnames(res))
    colnames(res) = gsub('MSeq2.Wrench','Wrench+MSeq',colnames(res))
    colnames(res) = gsub('MSeq2','MSeq',colnames(res))
    colnames(res) = gsub('eBayW','eBay(Wilcox)',colnames(res))
    colnames(res) = gsub('eBayt','eBay(t-test)',colnames(res))
    colnames(res) = gsub('BBinomial','corncob',colnames(res))
    colnames(res) = gsub('Aldex2we','Aldex2(t-test)',colnames(res))
    colnames(res) = gsub('^Aldex2$','Aldex2(Wilcox)',colnames(res))
    colnames(res) = gsub('^mbzinb','Omnibus',colnames(res))
    colnames(res) = gsub('^Wilcox.Wrench','Wrench+Wilcox',colnames(res))
    colnames(res) = gsub('^Wilcox.gmpr','GMPR+Wilcox',colnames(res))
    colnames(res) = gsub('Maaslin2','MaAsLin2',colnames(res))
    
    
    len <- sub.methods[!(sub.methods %in% colnames(res))]
    if(length(len) > 0){
      for(k in 1:length(len)){
        res[,len[k]] <- rep(1, nrow(res))
      }
    }
    cat(files[i],'_',len,'\n')
    
    # res[is.na(res)] <- 1
    # res[grep('NA',res)] <- 1
    
    res <- res[,sub.methods]
    
    load(paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationC4/','diff_',files[i]))
    colnames(diff)[2] = 'diff.otu.id'
    
    res <- full_join(res %>% rownames_to_column('otu.id'),diff) %>% column_to_rownames('otu.id')
    
    head(res)
    
    
    
    truth <- apply(res[,-ncol(res), drop = F],1, function(x) sum(x[!is.na(x)] <= 0.05)/(ncol(res)-1)) 
    truth <- as.data.frame(cbind(sapply(truth, function(x) x > 0.2),sapply(truth, function(x) x > 0.4),sapply(truth, function(x) x > 0.6),sapply(truth, function(x) x > 0.8)))
    colnames(truth) <- c('pct20','pct40','pct60','pct80')
    truth <- apply(truth, 2, function(x) ifelse(x==TRUE, 0, 1))
    
    res <- full_join(as.data.frame(truth) %>% rownames_to_column('otu.id'), res %>% rownames_to_column('otu.id')) %>% column_to_rownames('otu.id')
    
    
    pct.cal <- function(res, i){
      tp <- sum(res[,i] <= 0.05 & res$diff.otu.id ==TRUE, na.rm = T)
      tn <- sum(res[,i] > 0.05 &  res$diff.otu.id ==FALSE, na.rm = T)
      fn <- sum(res[,i] > 0.05 &  res$diff.otu.id ==TRUE, na.rm = T)
      fp <- sum(res[,i] <= 0.05 &  res$diff.otu.id ==FALSE, na.rm = T)
      tpr <- ifelse(!tp, 0, tp/(tp+fn))
      fdr <- ifelse(!fp, 0, fp/(tp+fp))
      return(c(tpr, fdr))
    }
    
    tt <- name <- NULL
    for(m in 1:(ncol(res)-1)){
      tt <- rbind(tt, pct.cal(res,m))
      name <- c(name, colnames(res)[m])
    }
    
    rownames(tt) <- name
    colnames(tt) <- c('tpr','fdr')
    tt <- as.data.frame(tt) %>% rownames_to_column('method')
    tt$class <- class
    tt$iter <- i
    TT <- rbind(TT, tt)
  }
}

save(TT, file = '~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationC4/RES_TT.RData')

load('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationC4/RES_TT.RData')
head(TT)
TT$diff.otu.pcts <- gsub('.*balanced','',TT$class)
TT$diff.otu.modes <- gsub('nOTU.*','',TT$class)

TT$nOTUs <- gsub('abundant','',TT$class)
TT$nOTUs <- gsub('rare','',TT$nOTUs)
TT$nOTUs <- gsub('binary.*','',TT$nOTUs)
TT$covariate.types <- 'binary'
TT$diff.otu.directs <- 'balanced'

TT <- TT %>% dplyr::select(-class)
head(TT)

TT1 <- TT %>% dplyr::select(-tpr) %>% mutate(measures ='FDR')
TT2 <- TT %>% dplyr::select(-fdr) %>% mutate(measures ='TPR')
colnames(TT1)[2] <- colnames(TT2)[2] <- 'value'
TT <- rbind(TT1, TT2)
colnames(TT)[1] = 'methods'
head(TT)

TT[,'nOTUs'] <- as.factor(TT[,'nOTUs'])
TT$diff.otu.pcts <- as.factor(TT$diff.otu.pcts)
TT$diff.otu.modes <- as.factor(TT$diff.otu.modes)
TT$covariate.types <- as.factor(TT$covariate.types)
TT$diff.otu.directs <- as.factor(TT$diff.otu.directs)
TT$measures <- as.factor(TT$measures)
str(TT)

levels(TT[,'nOTUs']) <- c('OTU=50','OTU=500')
levels(TT$diff.otu.modes) = str_to_title(levels(TT$diff.otu.modes))
levels(TT$diff.otu.pcts) = paste0(str_to_title(levels(TT$diff.otu.pcts)),' denisty')

head(TT)
res.df2 = data_summary(TT, as.formula(paste0('value ~ diff.otu.modes + nOTUs + diff.otu.pcts+ diff.otu.directs + methods + measures')))


save(res.df2, file = '~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/consensusC4.Rdata')




## C54
setwd('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationC54/')
files <- list.files(pattern = 'C54_summarymatrix1-')
length(files)
files[1:5]

load(paste0('C54_summarymatrix',1, '-Maaslin2.Rdata'))
nmss <- names(res_seqs)
nms <- gsub('Maaslin2','',nmss)
nms <- nms[grep('binary',nms)]
# nms <- c('L2low','L2high','L4low','L4high')

methods <- c('glmquassi2','Maaslin2','Rarefy', 'Wilcox','Wilcox.Wrench','Wilcox.gmpr','mbzinb','eBayW','eBayt','Aldex2we','Aldex2','RAIDA', 'DACOMP', 'LDM', 'glmquassi','BBinomial','ANCOMBC','DESeq2', 'DESeq2.Wrench', 'DESeq2.gmpr','edgeR', 'edgeR.Wrench', 'edgeR.gmpr','MSeq2','MSeq2.Wrench')
RES <- list()
for(nm in nms){
  for(i in 1:100){
    for(method in methods){
      try({
        load(paste0('C54_summarymatrix',i, '-',method,'.Rdata'))
        res <- res_seqs[[paste0(nm,method)]][,c('otu.id','fdr')]
        colnames(res)[2] <- method 
        head(res)
        RES[[paste0(method)]] <- res
      })
    }
    ## each rdata is each dataset with different method
    x = gsub('D3','',nm)
    x = gsub('loglinearSub_L1','',x)
    x = gsub('none','',x)
    x = gsub('nSam_L2L3','',x)
    save(RES, file = paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationC54/iter',i,'-',x,'.rdata'))
  }
}



for(nm in nms){
  for(i in 1:100){
    try({
      load(paste0('C54_summarymatrix',i, '-Wilcox.Rdata'))
      diff <- res_seqs[[paste0(nm,'Wilcox')]][,c('otu.id','diff.otu.ind')]
      x = gsub('D3','',nm)
      x = gsub('loglinearSub_L1','',x)
      x = gsub('none','',x)
      x = gsub('nSam_L2L3','',x)
      
      colnames(diff)[2] <- paste0(i,'_',x) 
    })
    save(diff, file = paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationC54/diff_iter',i,'-',x,'.rdata'))
  }
}

setwd('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationC54/')
nms <- list.files(pattern = 'rdata$')
nms <- gsub('.*-','',nms)
nms <- gsub('.rdata','',nms) %>% unique()

TT <- NULL
for(class in nms){
  files <- list.files(pattern = 'rdata$')
  files <- files[grep('^iter',files)]
  
  files <- files[grep(class, files)]
  
  
  TPRRs = FDRRs = NULL
  
  for(i in 1:length(files)){
    load(paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationC54/',files[i]))
    res <- RES[[1]];for(j in 2:length(RES)){res <- full_join(res, RES[[j]])}
    res <- res %>% column_to_rownames('otu.id')
    res[is.na(res)] <- 1
    dim(res)
    
    ## common findings
    colnames(res)[colnames(res) =='Wilcox'] = 'TSS+Wilcox'
    colnames(res)[colnames(res) =='Rarefy'] = 'Rarefy+Wilcox'
    colnames(res) = gsub('ANCOMBC','ANCOM-BC',colnames(res))
    colnames(res) = gsub('^glmquassi2$','GMPR+glm',colnames(res))
    colnames(res) = gsub('glmquassi','Wrench+glm',colnames(res))
    colnames(res) = gsub('DESeq2.gmpr','GMPR+DESeq2',colnames(res))
    colnames(res) = gsub('DESeq2.Wrench','Wrench+DESeq2',colnames(res))
    colnames(res) = gsub('edgeR.gmpr','GMPR+edgeR',colnames(res))
    colnames(res) = gsub('edgeR.Wrench','Wrench+edgeR',colnames(res))
    colnames(res) = gsub('MSeq2.Wrench','Wrench+MSeq',colnames(res))
    colnames(res) = gsub('MSeq2','MSeq',colnames(res))
    colnames(res) = gsub('eBayW','eBay(Wilcox)',colnames(res))
    colnames(res) = gsub('eBayt','eBay(t-test)',colnames(res))
    colnames(res) = gsub('BBinomial','corncob',colnames(res))
    colnames(res) = gsub('Aldex2we','Aldex2(t-test)',colnames(res))
    colnames(res) = gsub('^Aldex2$','Aldex2(Wilcox)',colnames(res))
    colnames(res) = gsub('^mbzinb','Omnibus',colnames(res))
    colnames(res) = gsub('^Wilcox.Wrench','Wrench+Wilcox',colnames(res))
    colnames(res) = gsub('^Wilcox.gmpr','GMPR+Wilcox',colnames(res))
    colnames(res) = gsub('Maaslin2','MaAsLin2',colnames(res))
    
    
    len <- sub.methods[!(sub.methods %in% colnames(res))]
    if(length(len) > 0){
      for(k in 1:length(len)){
        res[,len[k]] <- rep(1, nrow(res))
      }
    }
    cat(files[i],'_',len,'\n')
    
    # res[is.na(res)] <- 1
    # res[grep('NA',res)] <- 1
    
    res <- res[,sub.methods]
    
    load(paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationC54/','diff_',files[i]))
    colnames(diff)[2] = 'diff.otu.id'
    
    res <- full_join(res %>% rownames_to_column('otu.id'),diff) %>% column_to_rownames('otu.id')
    
    head(res)
    
    
    
    truth <- apply(res[,-ncol(res), drop = F],1, function(x) sum(x[!is.na(x)] <= 0.05)/(ncol(res)-1)) 
    truth <- as.data.frame(cbind(sapply(truth, function(x) x > 0.2),sapply(truth, function(x) x > 0.4),sapply(truth, function(x) x > 0.6),sapply(truth, function(x) x > 0.8)))
    colnames(truth) <- c('pct20','pct40','pct60','pct80')
    truth <- apply(truth, 2, function(x) ifelse(x==TRUE, 0, 1))
    
    res <- full_join(as.data.frame(truth) %>% rownames_to_column('otu.id'), res %>% rownames_to_column('otu.id')) %>% column_to_rownames('otu.id')
    
    
    pct.cal <- function(res, i){
      tp <- sum(res[,i] <= 0.05 & res$diff.otu.id ==TRUE, na.rm = T)
      tn <- sum(res[,i] > 0.05 &  res$diff.otu.id ==FALSE, na.rm = T)
      fn <- sum(res[,i] > 0.05 &  res$diff.otu.id ==TRUE, na.rm = T)
      fp <- sum(res[,i] <= 0.05 &  res$diff.otu.id ==FALSE, na.rm = T)
      tpr <- ifelse(!tp, 0, tp/(tp+fn))
      fdr <- ifelse(!fp, 0, fp/(tp+fp))
      return(c(tpr, fdr))
    }
    
    tt <- name <- NULL
    for(m in 1:(ncol(res)-1)){
      tt <- rbind(tt, pct.cal(res,m))
      name <- c(name, colnames(res)[m])
    }
    
    rownames(tt) <- name
    colnames(tt) <- c('tpr','fdr')
    tt <- as.data.frame(tt) %>% rownames_to_column('method')
    tt$class <- class
    tt$iter <- i
    TT <- rbind(TT, tt)
  }
}

save(TT, file = '~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationC54/RES_TT.RData')

load('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationC54/RES_TT.RData')
head(TT)
TT$diff.otu.pcts <- gsub('.*unbalanced','',TT$class)
TT$diff.otu.modes <- gsub('nOTU.*','',TT$class)

TT$nOTUs <- gsub('abundant','',TT$class)
TT$nOTUs <- gsub('rare','',TT$nOTUs)
TT$nOTUs <- gsub('binary.*','',TT$nOTUs)
TT$covariate.types <- 'binary'
TT$diff.otu.directs <- 'unbalanced'

TT <- TT %>% dplyr::select(-class)
head(TT)

TT1 <- TT %>% dplyr::select(-tpr) %>% mutate(measures ='FDR')
TT2 <- TT %>% dplyr::select(-fdr) %>% mutate(measures ='TPR')
colnames(TT1)[2] <- colnames(TT2)[2] <- 'value'
TT <- rbind(TT1, TT2)
colnames(TT)[1] = 'methods'
head(TT)

TT[,'nOTUs'] <- as.factor(TT[,'nOTUs'])
TT$diff.otu.pcts <- as.factor(TT$diff.otu.pcts)
TT$diff.otu.modes <- as.factor(TT$diff.otu.modes)
TT$covariate.types <- as.factor(TT$covariate.types)
TT$diff.otu.directs <- as.factor(TT$diff.otu.directs)
TT$measures <- as.factor(TT$measures)
str(TT)

levels(TT[,'nOTUs']) <- c('OTU=50','OTU=500')
levels(TT$diff.otu.modes) = str_to_title(levels(TT$diff.otu.modes))
levels(TT$diff.otu.pcts) = paste0(str_to_title(levels(TT$diff.otu.pcts)),' denisty')

head(TT)
res.df2 = data_summary(TT, as.formula(paste0('value ~ diff.otu.modes + nOTUs + diff.otu.pcts+ diff.otu.directs + methods + measures')))


save(res.df2, file = '~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/consensusC54.Rdata')







## D54
setwd('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationD54/')
files <- list.files(pattern = 'D54_summarymatrix1-')
length(files)
files[1:5]

load(paste0('D54_summarymatrix',1, '-Maaslin2.Rdata'))
nmss <- names(res_seqs)
nms <- gsub('Maaslin2','',nmss)
nms <- nms[grep('binary',nms)]

methods <- c('glmquassi2','Maaslin2','Rarefy', 'Wilcox','Wilcox.Wrench','Wilcox.gmpr','mbzinb','eBayW','eBayt','Aldex2we','Aldex2','RAIDA', 'DACOMP', 'LDM', 'glmquassi','BBinomial','ANCOMBC','DESeq2', 'DESeq2.Wrench', 'DESeq2.gmpr','edgeR', 'edgeR.Wrench', 'edgeR.gmpr','MSeq2','MSeq2.Wrench')
RES <- list()
for(nm in nms){
  for(i in 1:100){
    for(method in methods){
      try({
        load(paste0('D54_summarymatrix',i, '-',method,'.Rdata'))
        res <- res_seqs[[paste0(nm,method)]][,c('otu.id','fdr')]
        colnames(res)[2] <- method 
        head(res)
        RES[[paste0(method)]] <- res
      })
    }
    ## each rdata is each dataset with different method
    x = gsub('D3','',nm)
    x = gsub('loglinearSub_L1','',x)
    x = gsub('none','',x)
    x = gsub('nSam_L2L3','',x)
    save(RES, file = paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationD54/iter',i,'-',x,'.rdata'))
  }
}



for(nm in nms){
  for(i in 1:100){
    try({
      load(paste0('D54_summarymatrix',i, '-Wilcox.Rdata'))
      diff <- res_seqs[[paste0(nm,'Wilcox')]][,c('otu.id','diff.otu.ind')]
      x = gsub('D3','',nm)
      x = gsub('loglinearSub_L1','',x)
      x = gsub('none','',x)
      x = gsub('nSam_L2L3','',x)
      
      colnames(diff)[2] <- paste0(i,'_',x) 
    })
    save(diff, file = paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationD54/diff_iter',i,'-',x,'.rdata'))
  }
}

setwd('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationD54/')
nms <- list.files(pattern = 'rdata$')
nms <- gsub('.*-','',nms)
nms <- gsub('.rdata','',nms) %>% unique()

TT <- NULL
for(class in nms){
  files <- list.files(pattern = 'rdata$')
  files <- files[grep('^iter',files)]
  
  files <- files[grep(class, files)]
  
  
  TPRRs = FDRRs = NULL
  
  for(i in 1:length(files)){
    load(paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationD54/',files[i]))
    res <- RES[[1]];for(j in 2:length(RES)){res <- full_join(res, RES[[j]])}
    res <- res %>% column_to_rownames('otu.id')
    res[is.na(res)] <- 1
    dim(res)
    
    ## common findings
    colnames(res)[colnames(res) =='Wilcox'] = 'TSS+Wilcox'
    colnames(res)[colnames(res) =='Rarefy'] = 'Rarefy+Wilcox'
    colnames(res) = gsub('ANCOMBC','ANCOM-BC',colnames(res))
    colnames(res) = gsub('^glmquassi2$','GMPR+glm',colnames(res))
    colnames(res) = gsub('glmquassi','Wrench+glm',colnames(res))
    colnames(res) = gsub('DESeq2.gmpr','GMPR+DESeq2',colnames(res))
    colnames(res) = gsub('DESeq2.Wrench','Wrench+DESeq2',colnames(res))
    colnames(res) = gsub('edgeR.gmpr','GMPR+edgeR',colnames(res))
    colnames(res) = gsub('edgeR.Wrench','Wrench+edgeR',colnames(res))
    colnames(res) = gsub('MSeq2.Wrench','Wrench+MSeq',colnames(res))
    colnames(res) = gsub('MSeq2','MSeq',colnames(res))
    colnames(res) = gsub('eBayW','eBay(Wilcox)',colnames(res))
    colnames(res) = gsub('eBayt','eBay(t-test)',colnames(res))
    colnames(res) = gsub('BBinomial','corncob',colnames(res))
    colnames(res) = gsub('Aldex2we','Aldex2(t-test)',colnames(res))
    colnames(res) = gsub('^Aldex2$','Aldex2(Wilcox)',colnames(res))
    colnames(res) = gsub('^mbzinb','Omnibus',colnames(res))
    colnames(res) = gsub('^Wilcox.Wrench','Wrench+Wilcox',colnames(res))
    colnames(res) = gsub('^Wilcox.gmpr','GMPR+Wilcox',colnames(res))
    colnames(res) = gsub('Maaslin2','MaAsLin2',colnames(res))
    
    
    len <- sub.methods[!(sub.methods %in% colnames(res))]
    if(length(len) > 0){
      for(k in 1:length(len)){
        res[,len[k]] <- rep(1, nrow(res))
      }
    }
    cat(files[i],'_',len,'\n')
    
    # res[is.na(res)] <- 1
    # res[grep('NA',res)] <- 1
    
    res <- res[,sub.methods]
    
    load(paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationD54/','diff_',files[i]))
    colnames(diff)[2] = 'diff.otu.id'
    
    res <- full_join(res %>% rownames_to_column('otu.id'),diff) %>% column_to_rownames('otu.id')
    
    head(res)
    
    
    
    truth <- apply(res[,-ncol(res), drop = F],1, function(x) sum(x[!is.na(x)] <= 0.05)/(ncol(res)-1)) 
    truth <- as.data.frame(cbind(sapply(truth, function(x) x > 0.2),sapply(truth, function(x) x > 0.4),sapply(truth, function(x) x > 0.6),sapply(truth, function(x) x > 0.8)))
    colnames(truth) <- c('pct20','pct40','pct60','pct80')
    truth <- apply(truth, 2, function(x) ifelse(x==TRUE, 0, 1))
    
    res <- full_join(as.data.frame(truth) %>% rownames_to_column('otu.id'), res %>% rownames_to_column('otu.id')) %>% column_to_rownames('otu.id')
    
    
    pct.cal <- function(res, i){
      tp <- sum(res[,i] <= 0.05 & res$diff.otu.id ==TRUE, na.rm = T)
      tn <- sum(res[,i] > 0.05 &  res$diff.otu.id ==FALSE, na.rm = T)
      fn <- sum(res[,i] > 0.05 &  res$diff.otu.id ==TRUE, na.rm = T)
      fp <- sum(res[,i] <= 0.05 &  res$diff.otu.id ==FALSE, na.rm = T)
      tpr <- ifelse(!tp, 0, tp/(tp+fn))
      fdr <- ifelse(!fp, 0, fp/(tp+fp))
      return(c(tpr, fdr))
    }
    
    tt <- name <- NULL
    for(m in 1:(ncol(res)-1)){
      tt <- rbind(tt, pct.cal(res,m))
      name <- c(name, colnames(res)[m])
    }
    
    rownames(tt) <- name
    colnames(tt) <- c('tpr','fdr')
    tt <- as.data.frame(tt) %>% rownames_to_column('method')
    tt$class <- class
    tt$iter <- i
    TT <- rbind(TT, tt)
  }
}

save(TT, file = '~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationD54/RES_TT.RData')

load('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/simulationD54/RES_TT.RData')
head(TT)
TT$diff.otu.pcts <- gsub('.*unbalanced','',TT$class)
TT$diff.otu.modes <- gsub('nOTU.*','',TT$class)

TT$nOTUs <- gsub('abundant','',TT$class)
TT$nOTUs <- gsub('rare','',TT$nOTUs)
TT$nOTUs <- gsub('binary.*','',TT$nOTUs)
TT$covariate.types <- 'binary'
TT$diff.otu.directs <- 'unbalanced'

TT <- TT %>% dplyr::select(-class)
head(TT)

TT1 <- TT %>% dplyr::select(-tpr) %>% mutate(measures ='FDR')
TT2 <- TT %>% dplyr::select(-fdr) %>% mutate(measures ='TPR')
colnames(TT1)[2] <- colnames(TT2)[2] <- 'value'
TT <- rbind(TT1, TT2)
colnames(TT)[1] = 'methods'
head(TT)

TT[,'nOTUs'] <- as.factor(TT[,'nOTUs'])
TT$diff.otu.pcts <- as.factor(TT$diff.otu.pcts)
TT$diff.otu.modes <- as.factor(TT$diff.otu.modes)
TT$covariate.types <- as.factor(TT$covariate.types)
TT$diff.otu.directs <- as.factor(TT$diff.otu.directs)
TT$measures <- as.factor(TT$measures)
str(TT)

levels(TT[,'nOTUs']) <- c('OTU=50','OTU=500')
levels(TT$diff.otu.modes) = str_to_title(levels(TT$diff.otu.modes))
levels(TT$diff.otu.pcts) = paste0(str_to_title(levels(TT$diff.otu.pcts)),' denisty')

head(TT)
res.df2 = data_summary(TT, as.formula(paste0('value ~ diff.otu.modes + nOTUs + diff.otu.pcts+ diff.otu.directs + methods + measures')))


save(res.df2, file = '~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/consensusD54.Rdata')











