
edgeR.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  if(length(dat$covariates %>% unique())==2){
    cat('Covariate is 2-factor levels! \n')
    dat$covariates <- as.factor(dat$covariates)
  }else{
    dat$covariates <- dat$covariates
  }
  
  d <- DGEList(counts=dat$counts, group=dat$covariates)
  d <- edgeR::calcNormFactors(d) # TMM normalization
  design <- model.matrix(~dat$covariates)
  rownames(design) <- colnames(d)
  d <- estimateDisp(d, design)# also has robust = T mode
  
  if(length(table(dat$covariates))>2){ # see edgeR guide 
    # GLM likelihood ratio test
    fit <- glmFit(d,design) 
    res <- glmLRT(fit,coef=2)$table %>% mutate(fdr = p.adjust(PValue), pval = PValue) %>% dplyr::select(c('pval','fdr'))
    # rownames(res) = otu.name
    
    ## GLM quassi-likelihood ratio test
    # d=DGEList(counts=dat$counts, group=dat$covariates)
    # X <- ns(dat$covariates, df=3)
    # design <- model.matrix(~ X)
    # y <- estimateDisp(d, design)
    # fit <- glmQLFit(y, design, robust=TRUE) # negative binomial GLM for each tag and produces an object of class DGEGLM with some new components
    # res <- glmQLFTest(fit, coef=2:4)$table %>% #quasi-likelihood (QL) F-test
    #   rownames_to_column('otu.id') %>% 
    #   mutate(fdr = p.adjust(PValue), pval = PValue) %>% 
    #   dplyr::select(c('otu.id','pval','fdr')) %>% column_to_rownames('otu.id')
  }else{ # see edgeR guide Casestudy 4.1
    ## classic procedure comparing 2/+ grps
    fit=exactTest(d)
    res=topTags(fit, n=nrow(dat$counts))[[1]] %>% as.data.frame()
    colnames(res)[3:4] = c('pval','fdr')
  }
  
  if(sum(is.na(res$fdr))==nrow(dat$counts)){
    break
  }else{
    res <- res[,c('pval','fdr'), drop = F] %>% rownames_to_column('otu.id') 
    
    na <- sum(is.na(res$fdr))/nrow(res)
    sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
    time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
    
    if(FDR){
      res <- res %>% inner_join(dat$truth, by = 'otu.id')
      tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
      tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
      fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
      fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
      tpr <- ifelse(!tp, 0, tp/(tp+fn))
      fpr <- ifelse(!fp, 0, fp/(tn+fp))
      fdr <- ifelse(!fp, 0, fp/(tp+fp))
      F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
      MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
      return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
    }else{
      return(list(res = res, na = na, time = time, sig = sig))
    }
    
  }
  
  
  
}

edgeR.gmpr.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  if(length(dat$gmpr.covariates %>% unique())==2){
    cat('Covariate is 2-factor levels! \n')
    dat$covariates <- as.factor(dat$covariates)
  }else{
    dat$covariates <- dat$covariates
  }
  d=DGEList(counts=dat$counts, group=as.matrix(dat$meta.dat), norm.factors = dat$gmpr.size.factor)
  ## ignore normalization, none: the normalization factors are set to 1
  design <- model.matrix(~dat$covariates)
  d=estimateDisp(d, design)
  if(length(table(dat$covariates))>2){
    # GLM approach
    fit <- glmFit(d,design)
    res <- glmLRT(fit,coef=2)$table %>% mutate(fdr = p.adjust(PValue), pval = PValue) %>% dplyr::select(c('pval','fdr'))
    # rownames(res) = otu.name
  }else{
    ## the exact test in edgeR
    fit=exactTest(d)
    res=topTags(fit, n=nrow(dat$gmpr.counts))[[1]] %>% as.data.frame()
    colnames(res)[3:4] = c('pval','fdr')
  }
  
  if(sum(is.na(res$fdr))==nrow(dat$counts)){
    break
  }else{
    
  res <- res[,c('pval','fdr')] %>% rownames_to_column('otu.id')
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
}

edgeR.Wrench.wrapper <- function(dat, cutoff = 0.05, FDR = T) { 
  t1 = Sys.time()
  # see link: https://bioconductor.org/packages/devel/bioc/vignettes/Wrench/inst/doc/vignette.html
  d=edgeR::DGEList(counts=dat$counts,
                    group = as.matrix(dat$covariates),
                    norm.factors=dat$compositionalFactors)
  ## the classic approach in edgeR
    d=estimateDisp(d) # similar as estimateCommonDisp() + estimateTagWiseDisp()

  fit=exactTest(d)
  res=topTags(fit, n=nrow(dat$counts))[[1]] %>% as.data.frame()
  colnames(res)[3:4] = c('pval','fdr') 
  
  if(sum(is.na(res$fdr))==nrow(dat$counts)){
    break
  }else{
  res <- res[,c('pval','fdr')] %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
}

metagenomeSeq2.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  mgs = newMRexperiment(counts=dat$counts)
  mgs = metagenomeSeq::cumNorm(mgs) # calculate the scaling factors, run cumNorm
  pd = pData(mgs)
  mod = model.matrix(~1 + dat$covariates, data = pd)
  if(length(table(dat$covariates))>2){
    #fitZig 
    settings = zigControl(maxit=100,verbose=TRUE)
    fit = fitZig(obj = mgs,mod=mod,control=settings) #by default, the normalizing factors for obj=MGS are included in the model
    res <- MRcoefs(fit, coef=colnames(mod)[2],by=colnames(mod)[2], number = dim(fData(mgs))[1], group=2)
  }else{
    fit = fitFeatureModel(mgs,mod) # MSeq recommned this over the MSeq1, which uses fitZig()
    # res = MRcoefs(fit)
    res = MRfulltable(fit, number = nrow(assayData(mgs)$counts)) 
  }
  colnames(res)[colnames(res) == 'adjPvalues'] = 'fdr'
  colnames(res)[colnames(res) == 'pvalues'] = 'pval'
  
  if(sum(is.na(res$fdr))==nrow(dat$counts)){
    break
  }else{
  res <- res[,c('pval','fdr')] %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
}

metagenomeSeq2.Wrench.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  mgs = newMRexperiment(counts=dat$counts)
  metagenomeSeq::normFactors(mgs) <- dat$normalizationFactors
  mod = model.matrix(~dat$covariates)
  fit = fitFeatureModel(mgs,mod) # MSeq recommned this over the MSeq1, which uses fitZig()
  res = MRfulltable(fit, number = nrow(assayData(mgs)$counts))
  colnames(res)[colnames(res) == 'adjPvalues'] = 'fdr'
  colnames(res)[colnames(res) == 'pvalues'] = 'pval'
  
  if(sum(is.na(res$fdr))==nrow(dat$counts)){
    break
  }else{
  res <- res[,c('pval','fdr')] %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
}

metagenomeSeq2.gmpr.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  mgs = newMRexperiment(counts=dat$counts)
  metagenomeSeq::normFactors(mgs) <- dat$gmpr.size.factor
  mod = model.matrix(~dat$covariates)
  fit = fitFeatureModel(mgs,mod) # MSeq recommned this over the MSeq1, which uses fitZig()
  res = MRfulltable(fit, number = nrow(assayData(mgs)$counts))
  colnames(res)[colnames(res) == 'adjPvalues'] = 'fdr'
  colnames(res)[colnames(res) == 'pvalues'] = 'pval'
  
  if(sum(is.na(res$fdr))==nrow(dat$counts)){
    break
  }else{
    res <- res[,c('pval','fdr')] %>% rownames_to_column('otu.id') 
    
    na <- sum(is.na(res$fdr))/nrow(res)
    sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
    time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
    
    if(FDR){
      res <- res %>% inner_join(dat$truth, by = 'otu.id')
      tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
      tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
      fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
      fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
      tpr <- ifelse(!tp, 0, tp/(tp+fn))
      fpr <- ifelse(!fp, 0, fp/(tn+fp))
      fdr <- ifelse(!fp, 0, fp/(tp+fp))
      F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
      MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
      return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
    }else{
      return(list(res = res, na = na, time = time, sig = sig))
    }
  }
}

wilcox.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  if(length(table(dat$covariates))>2){
    pval <- apply(dat$prop, 1, function (y) cor.test(as.numeric(y), dat$covariates, method = 'spearman')$p.value)
  }else{
    cat('Covariate is 2-factor levels! \n')
    pval <- apply(dat$prop, 1,  function (y) wilcox.test(y ~ as.factor(dat$covariates))$p.value)
  }
  
  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  
  if(sum(is.na(fdr))==nrow(dat$counts)){
    break
  }else{
  res <- cbind(pval = pval, fdr = fdr) %>% as.data.frame() %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
}

wilcox.gmpr.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  if(length(table(dat$covariates))>2){
    pval <- apply(dat$counts, 1, function (y) cor.test(as.numeric(y), dat$covariates, method = 'spearman')$p.value)
  }else{
    cat('Covariate is 2-factor levels! \n')
    pval <- apply(dat$counts, 1,  function (y) wilcox.test(y ~ as.factor(dat$covariates))$p.value)
  }
  if(sum(is.na(pval))==nrow(dat$counts)){
    break
  }else{
    
  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  res <- cbind(pval = pval, fdr = fdr) %>% as.data.frame() %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
}

wilcox.Wrench.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  comm <- t(t(dat$counts) / dat$normalizationFactors)
  if(length(table(dat$covariates))>2){
    pval <- apply(comm, 1, function (y) cor.test(as.numeric(y), dat$covariates, method = 'spearman')$p.value)
  }else{
    cat('Covariate is 2-factor levels! \n')
    pval <- apply(comm, 1,  function (y) wilcox.test(y ~ as.factor(dat$covariates))$p.value)
  }
  if(sum(is.na(pval))==nrow(dat$counts)){
    break
  }else{
    
  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  res <- cbind(pval = pval, fdr = fdr) %>% as.data.frame() %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
}

Rarefy.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  prop <- dat$rff.prop
  if(length(table(dat$covariates))>2){
    pval <- apply(prop, 1, function (y) cor.test(as.numeric(y), dat$covariates, method = 'spearman')$p.value)
  }else{
    cat('Covariate is 2-factor levels! \n')
    pval <- apply(prop, 1,  function (y) wilcox.test(y ~ as.factor(dat$covariates))$p.value)
  }
  
  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  if(sum(is.na(fdr))==nrow(dat$counts)){
    break
  }else{
    
    
  res <- as.data.frame(cbind(pval = pval, fdr = fdr)) %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
}

DESeq.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  comm = dat$counts + 1L # ADD 11/05/2020, for avoiding ERROR: every gene contains at least one zero, cannot compute log geometric means
  design <- as.formula(paste('~', names(dat$meta.dat[1])))
  dds <- DESeqDataSetFromMatrix(countData = comm, colData = dat$meta.dat, design= design)
  dds <- DESeq2::DESeq(dds)
  res <- as.data.frame(as.matrix(results(dds)))
  pval <- res$pvalue
  names(pval) <- rownames(res)
  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  if(sum(is.na(fdr))==nrow(dat$counts)){
    break
  }else{
    
  res <- as.data.frame(cbind(pval = pval, fdr = fdr)) %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
}

DESeq.gmpr.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()

  design <- as.formula(paste('~', names(dat$meta.dat[1])))
  dds <- DESeqDataSetFromMatrix(countData = dat$counts, colData = dat$meta.dat, design = design)
  sizeFactors(dds) <- dat$gmpr.size.factor
  dds <- DESeq(dds)
  res <- as.data.frame(as.matrix(results(dds)))
  pval <- res$pvalue
  names(pval) <- rownames(res)
  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  if(sum(is.na(fdr))==nrow(dat$counts)){
    break
  }else{
    
  res <- as.data.frame(cbind(pval = pval, fdr = fdr)) %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
}

DESeq.Wrench.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  # see link: https://bioconductor.org/packages/devel/bioc/vignettes/Wrench/inst/doc/vignette.html
  design <- as.formula(paste('~', names(dat$meta.dat[1])))
  dds <- DESeqDataSetFromMatrix(countData =dat$counts, colData = dat$meta.dat, design= design)
  DESeq2::sizeFactors(dds) <- dat$normalizationFactors
  dds <- DESeq(dds)
  res <- as.data.frame(as.matrix(results(dds)))
  pval <- res$pvalue
  names(pval) <- rownames(res)
  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  if(sum(is.na(fdr))==nrow(dat$counts)){
    break
  }else{
    
  res <- as.data.frame(cbind(pval = pval, fdr = fdr)) %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
}

Aldex2.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  comm <- as.data.frame(dat$counts)
  if(length(table(dat$covariates))>2){
    design <- model.matrix(as.formula('~ grp'), dat$meta.dat)
    x <- aldex.clr(comm, design)
    res.obj <- aldex.glm(x, design)
    res <- as.data.frame(res.obj[,c(8,10)])%>% rownames_to_column('otu.id')
    colnames(res)[2:3] = c('pval','fdr')
    # covariate <- dat$covariates
    # aldex.obj <- aldex.clr(comm, covariate)
    # corr.test <- aldex.corr(aldex.obj, covariate)
    ## for package version 1.22.0
    # pval <- corr.test$spearman.ep
    # fdr <- corr.test[,'spearman.eBH'] 
    # names(fdr) <- names(pval)<- rownames(corr.test)
    ## for package version 1.18.0
    # res <- as.data.frame(corr.test[,c('p','BH')])%>% rownames_to_column('otu.id') 
  }else{
    covariate <- as.factor(dat$covariates)
    aldex.obj <- aldex(comm, covariate)
    pval <- aldex.obj$wi.ep
    fdr <- aldex.obj$wi.eBH # wilcox method with BH adjustment
    names(fdr) <- names(pval) <- rownames(aldex.obj)
    res <- as.data.frame(cbind(pval, fdr)) %>% rownames_to_column('otu.id') 
  }
  
  if(sum(is.na(res$fdr))==nrow(dat$counts)){
    break
  }else{
    
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
}

mbzinb.wrapper <- function(dat, cutoff = 0.05, FDR = T) { # grp.name can be categorical and numerical
  t1 = Sys.time()
  
  mbzinb.data <- mbzinb.dataset(dat$counts, dat$meta.dat)
  mbzinb.data <- mbzinb::norm.factors(mbzinb.data,  norm.factors = dat$gmpr.size.factor)
  mbzinb.obj <- mbzinb.test(mbzinb.data, group = 'grp',filter.if.unfiltered =F)
  res <- mbzinb.results(mbzinb.obj,nreturn = nrow(dat$counts)) %>% rownames_to_column('otu.id') %>% mutate(pval = PValue, fdr = Padj) %>% dplyr::select(otu.id, pval, fdr)
  if(sum(is.na(res$fdr))==nrow(dat$counts)){
    break
  }else{
    
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
}

ANCOMBC.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  dat$meta.dat$new <- dat$meta.dat$grp
  sam <- phyloseq::sample_data(dat$meta.dat)
  otu.tab <- otu_table(dat$counts, taxa_are_rows = T)
  phyloseq <- merge_phyloseq(otu.tab, sam)
  
  if(length(dat$covariates %>% unique())==2){
    out = ancombc(phyloseq = phyloseq, formula = names(dat$meta.dat[1]), p_adj_method = "fdr", lib_cut = 0, zero_cut = 1,group = names(dat$meta.dat[1]), struc_zero =T, conserve = T)
    #out = ancombc(phyloseq = phyloseq, formula = names(dat$meta.dat)[1], p_adj_method = "fdr", lib_cut = 800, group = names(dat$meta.dat[1]))
  }else{
    out = ancombc(phyloseq = phyloseq, formula = names(dat$meta.dat)[1], p_adj_method = "fdr",lib_cut = 0, zero_cut = 1, group = NULL, struc_zero = FALSE, neg_lb = FALSE,conserve = T)
  }
  res <- as.data.frame(as.matrix(cbind(pval = out$res$p_val$grp,fdr = out$res$q_val$grp))) # Caution: the code adjusted since the cluster makes errors
  colnames(res) = c('pval','fdr')
  res$otu.id <- rownames(out$res$p_val)
  if(sum(is.na(res$fdr))==nrow(dat$counts)){
    break
  }else{
    
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
}

ZIBB.wrapper <- function(dat, FDR = T, cutoff = 0.05){
  # https://cran.r-project.org/web/packages/ZIBBSeqDiscovery/vignettes/vignettes.html
  t1 = Sys.time()
  x <- dat$counts
  y <- as.matrix(cbind(1, dat$meta.dat))
  z <- as.matrix(cbind(1, log(colSums(x))))
  out.free <- fitZIBB(x, # otu * sample
                      y, # metadata, 1st column is intercept = 1, 2nd column is covariate
                      z, # 1st column is intercept =1, 2nd column is log(depth) 
                      mode="free")
  out.constrained <- fitZIBB(x, y, z, mode="constrained", 
                             gn=3, betastart=out.free$betahat, 
                             psi.start=out.free$psi, eta.start=out.free$zeroCoef)
  # out.constrained.mcc <- mcc.adj(out.constrained, kostic.x, kostic.y, kostic.z, K=4) # default adjusted pvalue 
  if(sum(is.na(out.constrained$p))==nrow(dat$counts)){
    break
  }else{
    
  res <- as.data.frame(as.matrix(out.constrained$p))
  colnames(res) = 'pval'
  res <- res %>% mutate(fdr = p.adjust(pval, method = 'fdr'))  %>% rownames_to_column('otu.id')# Caution: the code adjusted since the cluster makes errors
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
}

ZicoSeq.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  grp.name <- names(dat$meta.dat)[1]
  
  zicoseq.obj <- ZicoSeq.plus(meta.dat = dat$meta.dat, comm = dat$counts, grp.name = grp.name, adj.name = NULL,
                              prev.filter = 0, abund.filter = 0, nonzero.filter = 0, min.prop = 0, 
                              is.winsor = TRUE, winsor.qt = 0.97,
                              is.prior = TRUE,  post.sample.no = 25, 
                              link.func = list(function (x) x^0.5), stats.combine.func = max, 
                              perm.no = 99,  strata = NULL, 
                              stage.no = 6, excl.pct = 0.20,
                              is.fwer = FALSE, verbose = TRUE, return.comm = FALSE, return.perm.F = FALSE)
  
  if(sum(is.na(zicoseq.obj$p.adj.fdr))==nrow(dat$counts)){
    break
  }else{
    fdr <- zicoseq.obj$p.adj.fdr
    res <- as.data.frame(cbind(pval = zicoseq.obj$p.raw,fdr = fdr)) %>% rownames_to_column('otu.id')
      na <- sum(is.na(res$fdr))/nrow(res)
      sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
      time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
    
    
    if(FDR){
      res <- res %>% inner_join(dat$truth, by = 'otu.id')
      tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
      tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
      fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
      fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
      tpr <- ifelse(!tp, 0, tp/(tp+fn))
      fpr <- ifelse(!fp, 0, fp/(tn+fp))
      fdr <- ifelse(!fp, 0, fp/(tp+fp))
      F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
      MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
      return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
    }else{
      return(list(res = res, na = na, time = time, sig = sig))
    }
    
  }
}

RAIDA.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  df.raida <- as.data.frame(dat$meta.dat) %>%
    rownames_to_column('otu.id') %>%
    inner_join(as.data.frame(t(dat$counts))  %>% rownames_to_column('otu.id'))
  df.raida <- df.raida[order(-df.raida$grp),]
  rownames(df.raida) <- NULL
  n.lib <- c(as.numeric(table(df.raida$grp)[1]),as.numeric(table(df.raida$grp)[2]))
  c.data <- df.raida %>% dplyr::select(-'grp') %>% column_to_rownames('otu.id') %>% t(.) %>% as.data.frame(.)
  
  # c.data <- as.data.frame(dat$counts)
  # n.lib <- c(as.numeric(table(dat$meta.dat)[1]),as.numeric(table(dat$meta.dat)[2]))
  # raida.obj <- raida(c.data, n.lib, show.ref.features = F, show.all.features = T)
  
  raida.obj <- raida(c.data, n.lib, show.ref.features = T, show.all.features = T)
  res <- raida.obj$result %>%
    rownames_to_column('otu.id') %>%
    mutate(pval = p, fdr = p.adj) %>%
    dplyr::select(otu.id, pval, fdr)
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  # add for references features are regarded as non-DAs, thus it is not wrong theoratically
  res <- rbind(res,as.data.frame(cbind(otu.id = raida.obj$reference.features, pval = rep(1,length(raida.obj$reference.features)),
                            fdr = rep(1,length(raida.obj$reference.features)))))
  res$pval <- as.numeric(as.character(res$pval))
  res$fdr <- as.numeric(as.character(res$fdr))
  if(sum(is.na(res$fdr))==nrow(dat$counts)){
    break
  }else{
    
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
}

dacomp.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  dacomp.otu <- t(dat$counts) 
  dacomp.meta <- as.vector(dat$covariates)
  otu.name <- colnames(dacomp.otu)
    
  references = dacomp.select_references(X = dacomp.otu, median_SD_threshold = 0.6, verbose = T)
  
  if(length(table(dat$covariates))>2){
    dacomp.obj <- dacomp.test(X = dacomp.otu, y = dacomp.meta, ind_reference_taxa = references,
                              nr_perm = 1000, test = DACOMP.TEST.NAME.SPEARMAN, verbose = T,q = cutoff,
                              disable_DSFDR=T)
  }else{
    dacomp.obj <- dacomp.test(X = dacomp.otu, y = dacomp.meta, ind_reference_taxa = references,
                              nr_perm = 1000, test = DACOMP.TEST.NAME.WILCOXON, verbose = T,q = cutoff,
                              disable_DSFDR=T)
  }
  pval <- dacomp.obj$p.values.test
  names(pval) <- otu.name
  qval <- dacomp.obj$p.values.test.adjusted
  names(qval) <- otu.name
  ref.otu = colnames(dacomp.otu)[references$selected_references]
  # replace reference otu with 1, while NA is still NA
  pval[names(pval) %in% ref.otu] = 1
  qval[names(qval) %in% ref.otu] = 1
  res <- as.data.frame(cbind(pval, fdr=qval)) %>% rownames_to_column('otu.id')
  if(sum(is.na(res$fdr))==nrow(dat$counts)){
    break
  }else{
    
    
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
}

ldm.wrapper <- function(dat, cutoff = 0.05, FDR = T) {#https://github.com/yijuanhu/LDM/
  t1 = Sys.time()
  ldm.otu <<- t(dat$counts) %>% as.data.frame()
  ldm.obj <- ldm(formula=ldm.otu ~ grp, data = dat$meta.dat, fdr.nominal = cutoff) # otu: row sample, col taxa
  res <- as.data.frame(t(ldm.obj$p.otu.omni)) %>% # only select p-values for the OTU-specific tests based on the omnibus statistics
    rownames_to_column('otu.id') %>%
    dplyr::rename(pval = V1) %>% mutate(fdr = p.adjust(pval, 'fdr'))
  
  if(sum(is.na(res$fdr))==nrow(dat$counts)){
    break
  }else{
    
    
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
}

glmquassi.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  comm <- dat$counts[,names(dat$gmpr.size.factor)]
  covariate <- dat$meta.dat[names(dat$gmpr.size.factor),]
  pvs <- NULL
  tryCatch({
    for(i in 1:nrow(comm)){
    m1.op <- glm(comm[i,] ~ covariate, offset=log(dat$gmpr.size.factor),family=quasipoisson)
    pv.op <- wald.test(b = coef(m1.op), Sigma = vcov(m1.op), Terms = 2)$result$chi2['P'] #coef(summary(m1.op))[2,4] 
    pvs <- c(pvs,pv.op)
  }}, error =function(e){cat(paste0(' ERROR : '),conditionMessage(e), "\n")})
  
  res <- as.data.frame(pvs)
  res$otu.id <- rownames(comm)
  colnames(res) <- c('pval','otu.id')
  res$fdr = p.adjust(res$pval, 'fdr') 
  if(sum(is.na(res$fdr))==nrow(dat$counts)){
    break
  }else{
    
    
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
}

Maaslin2.wrapper <- function(dat, cutoff = 0.05, FDR = T,output){
  t1 = Sys.time()
  otu <- t(dat$counts)
  meta <- dat$meta.dat 
  if(length(unique(meta$grp))==2){meta$grp <- as.factor(meta$grp)}
  
fit = Maaslin2(input_data = otu, input_metadata = meta, output = output, 
    fixed_effects = "grp",plot_heatmap = F, plot_scatter = F)
  head(fit$results)
  res <- fit$results[,c('feature','pval')] %>% dplyr::rename(otu.id = feature)
  res$fdr <- p.adjust(res$pval, method = 'fdr')
  if(sum(is.na(res$fdr))==nrow(dat$counts)){
    break
  }else{
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
  
}

BBinomial.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  sam <- phyloseq::sample_data(dat$meta.dat)
  otu.tab <- otu_table(dat$counts, taxa_are_rows = T)
  phyloseq <- merge_phyloseq(otu.tab, sam)
  # Scenario: Testing for differential abundance across Day, without controlling for anything else:
  corncob <- differentialTest(formula = ~ grp,
                   phi.formula = ~ 1,
                   formula_null = ~ 1,
                   phi.formula_null = ~ 1,
                   test = "Wald", boot = FALSE,
                   fdr_cutoff = cutoff,
                   data = phyloseq)
  res <- as.data.frame(cbind(corncob$p,corncob$p_fdr)) %>% rownames_to_column('otu.id')
  colnames(res)[2:3] <- c('pval','fdr')
  # res <- res[res$otu.id %in% otu.name,]
  if(sum(is.na(res$fdr))==nrow(dat$counts)){
    break
  }else{
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
}

eBayW.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  comm <- as.data.frame(t(dat$counts))
  covariate <- as.factor(dat$covariates)
  levels(covariate) <- c(0,1) # package require the case = 1, control = 0. If use other value, will cause fail
  
  eBay.res <- eBay(otu.data=comm, group=covariate, cutf=cutoff, test.methods='wilcoxon', adj.m='BH')
  res <- as.data.frame(cbind(pval = rep(1,length(eBay.res$final.p)),fdr = eBay.res$final.p)) %>% rownames_to_column('otu.id')
  if(sum(is.na(res$fdr))==nrow(dat$counts)){
    break
  }else{
    
    
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  }
}





