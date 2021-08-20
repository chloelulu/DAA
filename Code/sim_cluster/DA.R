prefix <- 'D4'
resdir <-  file.path(paste0("/research/bsi/projects/staff_analysis/m216453/SemiSim/", prefix))
source(file.path(paste0(resdir, "/Cluster_mayo.R")))
temp <- load('UrogenitalTract_V35.RData',envir=.GlobalEnv)
temp0 <- load('UrogenitalTract_V35_dirmult.RData',envir=.GlobalEnv)
source('pre_DA.R')


rm(.Random.seed)
paras <- list()
paras$model.paras <- NULL
paras$resdir <- resdir
paras$prefix <- prefix
paras$methods <- c('Maaslin2','ZicoSeq','glmquassi2','Rarefy', 'Wilcox','Wilcox.Wrench','Wilcox.gmpr','mbzinb','eBayW','eBayt','Aldex2we','Aldex2','RAIDA', 'DACOMP', 'LDM', 'glmquassi','BBinomial','ANCOMBC','DESeq2', 'DESeq2.Wrench', 'DESeq2.gmpr','edgeR', 'edgeR.Wrench', 'edgeR.gmpr','MSeq2','MSeq2.Wrench')
paras$nOTUs = c('nOTU_L1','nOTU_L5')
paras$nSams = c('nSam_L1','nSam_L2','nSam_L4')
paras$diff.otu.pcts= c('low', 'medium', 'high')
paras$diff.otu.modes = c('rare','abundant')
paras$diff.otu.directs = c('balanced','unbalanced')
paras$covariate.types = "binary"
paras$covariate.eff.means = 'L3'
paras$confounder.types = c('none','continuous')
paras$depth.mus = 'D3'
paras$depth.conf.factors = 'none'
paras$include.top.otu = FALSE
paras$models = 'loglinear'
paras$nSubs = 'Sub_L1'


setwd(resdir)
res <- clsapply(1:100, func, paras, queque='1-day', tempdir=file.path(resdir, 'tmpC'))
setwd(resdir)
save(res, file=file.path(resdir, paste(prefix, "_res.Rdata", sep="")))
