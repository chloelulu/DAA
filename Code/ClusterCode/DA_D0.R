prefix <- 'D0'
resdir <-  file.path(paste0("/research/bsi/projects/staff_analysis/m216453/SemiSim/", prefix))
source("/research/bsi/projects/staff_analysis/m216453/SemiSim/code/Cluster_mayo.R")
temp <- load('VaginalIntroitus_V35.RData',envir=.GlobalEnv)
temp0 <- load('VaginalIntroitus_V35_dirmult.RData',envir=.GlobalEnv)
#otu.tab <- VaginalIntroitus_V35
source('/research/bsi/projects/staff_analysis/m216453/SemiSim/code/pre_DA.R')

paras <- list()
paras$model.paras <- NULL
paras$resdir <- resdir
paras$prefix <- prefix
paras$methods <- c('ZicoSeq','glmquassi2','Rarefy', 'Wilcox','Wilcox.Wrench','Wilcox.gmpr','mbzinb','eBayW','eBayt','Aldex2we','Aldex2','RAIDA', 'DACOMP', 'LDM', 'glmquassi','BBinomial','ANCOMBC','DESeq2', 'DESeq2.Wrench', 'DESeq2.gmpr','edgeR', 'edgeR.Wrench', 'edgeR.gmpr','MSeq2','MSeq2.Wrench','fastANCOM','linda','ZicoSeq3','ZINQ')
paras$nOTUs = c('nOTU_L1','nOTU_L5')
paras$nSams = c('nSam_L1','nSam_L4')
paras$diff.otu.pcts= c('none')
paras$diff.otu.modes = c('abundant')
paras$diff.otu.directs = c('balanced')
paras$covariate.types = c('binary')
paras$covariate.eff.means = c('none')
paras$confounder.types = c('none')
paras$depth.mus = c('D3')
paras$depth.conf.factors = c('none','DL3')
paras$include.top.otu = FALSE
paras$models = 'loglinear'
paras$nSubs = 'Sub_L1'

setwd(resdir)
res <- clsapply(1:1000, func, paras, queque='1-hour', tempdir=file.path(resdir, 'tmpC'))
setwd(resdir)
save(res, file=file.path(resdir, paste(prefix, "_res.Rdata", sep="")))
