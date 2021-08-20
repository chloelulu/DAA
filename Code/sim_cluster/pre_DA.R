func <- function(part, method, paras) {
  
  set.seed(part)
  methods <- paras$methods
  model.paras <- paras$model.paras
  
  nSams <- paras$nSams
  nOTUs <- paras$nOTUs
  diff.otu.pcts <- paras$diff.otu.pcts
  diff.otu.directs <- paras$diff.otu.directs
  diff.otu.modes <- paras$diff.otu.modes
  covariate.types <- paras$covariate.types
  covariate.eff.means <- paras$covariate.eff.means
  confounder.types <- paras$confounder.types
  depth.mus <- paras$depth.mus
  depth.conf.factors <- paras$depth.conf.factors
  models <- paras$models
  nSubs <- paras$nSubs
  
  include.top.otu <- paras$include.top.otu
  
  measures <- c('TPR', 'FDR', 'FP','FPR', 'na','time','MCC','F1')
  
  resdir <- paras$resdir
  prefix <- paras$prefix
  resdir <- gsub('/$', '', resdir)
  
  # 23 methods totally
  methods_funs <- list('ZicoSeq'= "ZicoSeq.wrapper", 'Rarefy'= "Rarefy.wrapper", 'Aldex2'= "Aldex2.wrapper", 'Omnibus'="Omnibus.wrapper", 'ANCOM2'="ANCOM2.wrapper", 'ZicoSeq0'= "ZicoSeq0.wrapper",'ZicoSeq1'= "ZicoSeq1.wrapper",'ZicoSeq2'= "ZicoSeq2.wrapper",'ZicoSeq3'= "ZicoSeq3.wrapper",
                       'RAIDA'="RAIDA.wrapper", 'DACOMP' = "dacomp.wrapper", 'LDM'="ldm.wrapper", 'RioNorm2'="RioNorm.wrapper",'glmquassi'='glmquassi.wrapper','Maaslin2'="Maaslin2.wrapper",'glmquassi2'="glmquassi2.wrapper",
                       'BBinomial' = 'BBinomial.wrapper','ANCOMBC'="ANCOMBC.wrapper", 'MSeq'="metagenomeSeq.wrapper",'ZIBB' = "ZIBB.wrapper",'RDB'="RDB.wrapper",
                       'DESeq2'="DESeq.wrapper", 'DESeq2.Wrench'="DESeq.Wrench.wrapper", 'DESeq2.gmpr'="DESeq.gmpr.wrapper",
                       'Wilcox'='wilcox.wrapper' , 'Wilcox.Wrench'='wilcox.Wrench.wrapper' , 'Wilcox.gmpr'='wilcox.gmpr.wrapper',
                       'edgeR'="edgeR.wrapper", 'edgeR.gmpr'="edgeR.gmpr.wrapper", 'edgeR.Wrench'="edgeR.Wrench.wrapper",
                       'MSeq2.Wrench'="metagenomeSeq2.Wrench.wrapper", 'MSeq2'="metagenomeSeq2.wrapper",'MSeq2.gmpr'="metagenomeSeq2.gmpr.wrapper",
                       'linda'="linda.wrapper",'eBayt'="eBayt.wrapper",'eBayW'="eBayW.wrapper",'Aldex2we'="Aldex2we.wrapper",
                       'CLRBC' ='CLRBC.wrapper','ttest.gmpr' = 'ttest.gmpr.wrapper','ttest.Wrench' = 'ttest.Wrench.wrapper','mbzinb' = "mbzinb.wrapper","Rarefyttest"="Rarefyttest.wrapper")
  
  res <- array(NA, c(length(depth.mus), length(diff.otu.modes), length(models), length(nSubs), length(nOTUs), length(covariate.types), length(depth.conf.factors),
                     length(diff.otu.directs), length(confounder.types), length(nSams), length(covariate.eff.means), length(diff.otu.pcts),
                     length(measures), length(method)),
               dimnames=list(depth.mus, diff.otu.modes, models, nSubs, nOTUs, covariate.types, depth.conf.factors,
                             diff.otu.directs, confounder.types, nSams, covariate.eff.means, diff.otu.pcts,
                             measures, method))
  
  sink(file.path(resdir, paste(prefix, "_",  part, '_',method,".log", sep="")))
  cat(date(), '\n')
 path0 = paste0(resdir,'/demo/',part)
if(!dir.exists(path0)){dir.create(path0)}
 
  res_seqs <- list()
  for (diff.otu.pct in diff.otu.pcts) {
    for (diff.otu.mode in diff.otu.modes) {
      for (confounder.type in confounder.types) {
        for (depth.conf.factor in depth.conf.factors){
          for (covariate.type in covariate.types) {
            for (nSam in nSams) {
              for (nOTU in nOTUs) {
                for (depth.mu in depth.mus) {
                  for (model in models){
                    for (nSub in nSubs){
                      for (covariate.eff.mean in covariate.eff.means){
                        for (diff.otu.direct in diff.otu.directs){
                          cat(getwd())
                          load(paste0('../iter',part,'/',depth.mu, diff.otu.mode, model, nSub, nOTU, covariate.type, depth.conf.factor, diff.otu.direct, confounder.type, nSam, covariate.eff.mean, diff.otu.pct,'.Rdata'))
                          output = paste0(resdir,'/demo/',part,'/', depth.mu, diff.otu.mode, model, nSub, nOTU, covariate.type, depth.conf.factor, diff.otu.direct, confounder.type, nSam, covariate.eff.mean, diff.otu.pct)
			  if(!dir.exists(output)){dir.create(output)}
                          # apply each DA method
                          MCC <- F1 <- fpr <- fdr <- tpr <- fp <- na <- time <- NA;res_seq <- NULL
                          wrapper <- match.fun(methods_funs[[method]])
                          tryCatch({
                            cat(paste(method, '\n'))
                            if(method=='Maaslin2'){
                              wrapper.obj <- wrapper(dat,cutoff=0.05,FDR=T, output = output)
                            }else{
                              wrapper.obj <- wrapper(dat,cutoff=0.05,FDR=T)
                            }
                            fdr <- wrapper.obj$fdr
                            tpr <- wrapper.obj$tpr
                            fp <- wrapper.obj$fp
                            fpr <- wrapper.obj$fpr
                            MCC <- wrapper.obj$MCC
                            F1 <- wrapper.obj$F1
                            res_seq <- as.data.frame(wrapper.obj$res)
                            na <- wrapper.obj$na
                            time <- wrapper.obj$time
                            res[depth.mu, diff.otu.mode, model, nSub, nOTU, covariate.type, depth.conf.factor, diff.otu.direct, confounder.type, nSam, covariate.eff.mean, diff.otu.pct,'TPR', method] <- tpr
                            res[depth.mu, diff.otu.mode, model, nSub, nOTU, covariate.type, depth.conf.factor, diff.otu.direct, confounder.type, nSam, covariate.eff.mean, diff.otu.pct, 'FDR', method] <-  fdr
                            res[depth.mu, diff.otu.mode, model, nSub, nOTU, covariate.type, depth.conf.factor, diff.otu.direct, confounder.type, nSam, covariate.eff.mean, diff.otu.pct, 'FP', method] <-  fp
                            res[depth.mu, diff.otu.mode, model, nSub, nOTU, covariate.type, depth.conf.factor, diff.otu.direct, confounder.type, nSam, covariate.eff.mean, diff.otu.pct, 'FPR', method] <-  fpr
                            res[depth.mu, diff.otu.mode, model, nSub, nOTU, covariate.type, depth.conf.factor, diff.otu.direct, confounder.type, nSam, covariate.eff.mean, diff.otu.pct, 'na', method] <-  na
                            res[depth.mu, diff.otu.mode, model, nSub, nOTU, covariate.type, depth.conf.factor, diff.otu.direct, confounder.type, nSam, covariate.eff.mean, diff.otu.pct, 'time', method] <-  time
                            res[depth.mu, diff.otu.mode, model, nSub, nOTU, covariate.type, depth.conf.factor, diff.otu.direct, confounder.type, nSam, covariate.eff.mean, diff.otu.pct, 'MCC', method] <-  MCC
                            res[depth.mu, diff.otu.mode, model, nSub, nOTU, covariate.type, depth.conf.factor, diff.otu.direct, confounder.type, nSam, covariate.eff.mean, diff.otu.pct, 'F1', method] <-  F1
                            res_seqs[[paste0(depth.mu, diff.otu.mode, model, nSub, nOTU, covariate.type, depth.conf.factor, diff.otu.direct, confounder.type, nSam, covariate.eff.mean, diff.otu.pct, method)]] <- res_seq
                          }, error =function(e){cat(paste0(paste0(depth.mu, diff.otu.mode, model, nSub, nOTU, covariate.type, depth.conf.factor, diff.otu.direct, confounder.type, nSam, covariate.eff.mean, diff.otu.pct), method, ' ERROR : '),conditionMessage(e), "\n")})
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  warnings()
  cat('\n', date(), '\n', method,'Finished!')
  sink()
  save(res_seqs,file=file.path(resdir, paste(prefix, "_summarymatrix",  part, '-', method, ".Rdata", sep="")))
  #save(res, file=file.path(resdir, paste(prefix, "_res",  part, '-', method, ".Rdata", sep="")))
  return(res)
}

