## for binary + covariate 

ZicoSeq.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  grp.name <- names(dat$meta.dat)[1]
  
  zicoseq.obj <- ZicoSeq(meta.dat = dat$meta.dat, comm = dat$counts, grp.name = grp.name, adj.name = 'covariate',
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

ZicoSeq3.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  grp.name <- names(dat$meta.dat)[1]
  
  zicoseq.obj <- ZicoSeq(meta.dat = dat$meta.dat, comm = dat$counts, grp.name = grp.name, adj.name = 'covariate',
                         prev.filter = 0, abund.filter = 0, nonzero.filter = 0, min.prop = 0,
                         is.winsor = TRUE, winsor.qt = 0.97,
                         is.prior = TRUE,  post.sample.no = 25,
                         link.func = list(function (x) x^0.5), stats.combine.func = max,
                         perm.no = 99,  strata = NULL,
                         stage.no = 6, excl.pct = 0.10, ref.pct = 0.40,
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


glmquassi2.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  comm <- dat$counts[,names(dat$gmpr.size.factor)]
  grp <- dat$meta.dat[names(dat$gmpr.size.factor),]$grp
  covariate <- dat$meta.dat[names(dat$gmpr.size.factor),]$covariate
  pvs <- NULL
  tryCatch({
    for(i in 1:nrow(comm)){
      m1.op <- glm(comm[i,] ~ grp + covariate, offset=log(dat$normalizationFactors),family=quasipoisson) # with Wrench normalized factors
      pv.op <- wald.test(b = coef(m1.op), Sigma = vcov(m1.op), Terms = 2)$result$chi2['P']
      pvs <- c(pvs, pv.op)
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

glmquassi.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  comm <- dat$counts[,names(dat$gmpr.size.factor)]
  grp <- dat$meta.dat[names(dat$gmpr.size.factor),]$grp
  covariate <- dat$meta.dat[names(dat$gmpr.size.factor),]$covariate
  pvs <- NULL
  tryCatch({
    for(i in 1:nrow(comm)){
      m1.op <- glm(comm[i,] ~ grp + covariate, offset=log(dat$gmpr.size.factor),family=quasipoisson) # with GMPR factors
      pv.op <- wald.test(b = coef(m1.op), Sigma = vcov(m1.op), Terms = 2)$result$chi2['P']
      pvs <- c(pvs, pv.op)
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


Maaslin2.wrapper <- function(dat, cutoff = 0.05, FDR = T, output){
  t1 = Sys.time()
  otu <- t(dat$counts)
  meta <- dat$meta.dat 
  if(length(unique(meta$grp))==2){meta$grp <- as.factor(meta$grp)}
  if(length(colnames(meta))==1){
    fit = Maaslin2(input_data = otu, input_metadata = meta, output = output, 
                   min_abundance = 0,min_prevalence = 0, min_variance = 0,
                   fixed_effects = c('grp'), plot_heatmap = F, plot_scatter = F)
  }else{
    fit = Maaslin2(input_data = otu, input_metadata = meta, output = output, 
                   min_abundance = 0,min_prevalence = 0, min_variance = 0,
                   fixed_effects = c('grp','covariate'), plot_heatmap = F, plot_scatter = F)
  }
  
  res <- fit$results[,c('feature','pval','metadata')]%>%filter(metadata=='grp') %>% dplyr::select(-metadata) %>%dplyr::rename(otu.id = feature)
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

DESeq.Wrench.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  design <- as.formula(paste('~', names(dat$meta.dat[2]),'+',names(dat$meta.dat[1])))
  dds <- DESeqDataSetFromMatrix(countData = dat$counts, colData = dat$meta.dat, design = design)
  sizeFactors(dds) <- dat$normalizationFactors
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

ANCOMBC.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  dat$meta.dat$new <- dat$meta.dat$grp
  sam <- phyloseq::sample_data(dat$meta.dat)
  otu.tab <- otu_table(dat$counts, taxa_are_rows = T)
  phyloseq <- merge_phyloseq(otu.tab, sam)
  
  if(length(dat$covariates %>% unique())==2){
    out = ancombc(phyloseq = phyloseq, formula = paste0(names(dat$meta.dat[1]), '+', names(dat$meta.dat[2])), p_adj_method = "fdr", lib_cut = 0, zero_cut = 1,group = names(dat$meta.dat[1]), struc_zero =T, conserve = T)
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

Aldex2.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  comm <- as.data.frame(dat$counts)
  design <- model.matrix(as.formula('~ grp + covariate'), dat$meta.dat)
  x <- aldex.clr(comm, design)
  res.obj <- aldex.glm(x, design)
  res <- as.data.frame(res.obj[,c(8,14)])%>% rownames_to_column('otu.id')
  colnames(res)[2:3] = c('pval','fdr')

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
  ldm.obj <- ldm(formula=ldm.otu | covariate ~ grp, data = dat$meta.dat, fdr.nominal = cutoff) # otu: row sample, col taxa
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

BBinomial.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  sam <- phyloseq::sample_data(dat$meta.dat)
  otu.tab <- otu_table(dat$counts, taxa_are_rows = T)
  phyloseq <- merge_phyloseq(otu.tab, sam)
  # Scenario: Testing for differential abundance across Day, without controlling for anything else:
  corncob <- differentialTest(formula = ~ grp + covariate,
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

edgeR.Wrench.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  if(length(dat$covariates %>% unique())==2){
    cat('Covariate is 2-factor levels! \n')
    dat$covariates <- as.factor(dat$covariates)
  }else{
    dat$covariates <- dat$covariates
  }
  
  d <- DGEList(counts=dat$counts, group=dat$covariates,norm.factors=dat$compositionalFactors)
  # d <- edgeR::calcNormFactors(d) # TMM normalization
  design <- model.matrix(~dat$meta.dat$grp + dat$meta.dat$covariate)
  rownames(design) <- colnames(d)
  d <- estimateDisp(d, design)# also has robust = T mode
  
  fit=exactTest(d)
  res=topTags(fit, n=nrow(dat$gmpr.counts))[[1]] %>% as.data.frame()
  colnames(res)[3:4] = c('pval','fdr')
  res <- res[,c('pval','fdr')] %>% rownames_to_column('otu.id')
  
  ## glmFit
  # fit <- glmFit(d,design) 
  # res <- glmLRT(fit,coef=2:3)$table %>% rownames_to_column('otu.id') %>% mutate(fdr = p.adjust(PValue), pval = PValue)

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

edgeR.gmpr.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  if(length(dat$covariates %>% unique())==2){
    cat('Covariate is 2-factor levels! \n')
    dat$covariates <- as.factor(dat$covariates)
  }else{
    dat$covariates <- dat$covariates
  }
  
  d <- DGEList(counts=dat$counts, group=dat$covariates,norm.factors=dat$gmpr.size.factor)
  # d <- edgeR::calcNormFactors(d) # TMM normalization
  design <- model.matrix(~dat$meta.dat$grp + dat$meta.dat$covariate)
  rownames(design) <- colnames(d)
  d <- estimateDisp(d, design)# also has robust = T mode
  
  fit=exactTest(d)
  res=topTags(fit, n=nrow(dat$gmpr.counts))[[1]] %>% as.data.frame()
  colnames(res)[3:4] = c('pval','fdr')
  res <- res[,c('pval','fdr')] %>% rownames_to_column('otu.id')
  
  ## glmFit
  # fit <- glmFit(d,design) 
  # res <- glmLRT(fit,coef=2:3)$table %>% rownames_to_column('otu.id') %>% mutate(fdr = p.adjust(PValue), pval = PValue)
  
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

DESeq.gmpr.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  design <- as.formula(paste('~', names(dat$meta.dat[2]),'+',names(dat$meta.dat[1])))
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






## add 11/18/2021 for most recently developed methods
fastANCOM.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  library(fastANCOM)
  t1 = Sys.time()
  data <- t(dat$counts)
  group <- dat$meta.dat$grp
  z <- dat$meta.dat$covariate
  fit <- fastANCOM(Y=data, x=group, z = z)
  final_fit <- fit$results$final
  pval <- final_fit$log2FC.pval
  names(pval) <- rownames(final_fit)
  fdr <- p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  
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

linda.wrapper <- function(dat, cutoff = 0.05, FDR = T){
  t1 = Sys.time()
  feature.dat <- dat$counts
  meta.dat <- dat$meta.dat
  linda.obj <- MicrobiomeStat::linda(feature.dat, meta.dat, formula = '~ grp+covariate', alpha = cutoff, feature.dat.type ='count')
  res <- as.data.frame(linda.obj$output$grp) %>% rownames_to_column('otu.id') %>% mutate(pval = pvalue,fdr = padj) 
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

ZINQ.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  library(ZINQ)
  t1 = Sys.time()
  otu.tab <- dat$counts
  meta <- dat$meta.dat[colnames(otu.tab),,drop =F]
  tab <- merge(t(otu.tab),meta, by = 0)
  pval <- otu <- NULL
  for(i in rownames(otu.tab)){
    if(unique(meta$grp)==2){
      result <- ZINQ_tests(formula.logistic=as.formula(paste0(i,'~grp + covariate')), formula.quantile=as.formula(paste0(i,'~grp+ covariate')), C="grp", y_CorD = "D", data=tab)
    }else{
      result <- ZINQ_tests(formula.logistic=as.formula(paste0(i,'~grp+ covariate')), formula.quantile=as.formula(paste0(i,'~grp+ covariate')), C="grp", y_CorD = "C", data=tab)
    }
    P <- ZINQ_combination(result)
    otu <- c(otu, i)
    pval <- c(pval,P)
  }
  
  fdr <- p.adjust(pval, 'fdr')
  res <- as.data.frame(cbind(otu,pval, fdr)) %>% column_to_rownames('otu')
  res$pval <- as.numeric(as.character(res$pval));res$fdr <- as.numeric(as.character(res$fdr))
  
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

