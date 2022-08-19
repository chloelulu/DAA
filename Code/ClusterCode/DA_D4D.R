prefix <- 'D4D'
resdir <-  file.path(paste0("/research/bsi/projects/staff_analysis/m216453/SemiSim/", prefix))
source("/research/bsi/projects/staff_analysis/m216453/SemiSim/code/Cluster_mayo1.R")
temp <- load('VaginalIntroitus_V35.RData',envir=.GlobalEnv)
temp0 <- load('VaginalIntroitus_V35_dirmult.RData',envir=.GlobalEnv)
#otu.tab <- VaginalIntroitus_V35
source('/research/bsi/projects/staff_analysis/m216453/SemiSim/code/pre_DA1.R')

rm(.Random.seed)
paras <- list()
paras$model.paras <- NULL
paras$resdir <- resdir
paras$prefix <- prefix
paras$methods <- c('ZicoSeq','Aldex2', 'LDM', 'glmquassi2','BBinomial','ANCOMBC','DESeq2.Wrench', 'edgeR.Wrench','Maaslin2','linda',"fastANCOM",'ZINQ','ZicoSeq3')
paras$nOTUs = c('nOTU_L5')
paras$nSams = c('nSam_L2')
paras$diff.otu.pcts= c('low', 'medium', 'high')
paras$diff.otu.modes = c('abundant', 'rare')
paras$diff.otu.directs = c('balanced')
paras$covariate.types = c("binary")
paras$covariate.eff.means = 'L3'
paras$confounder.types = c('continuous')
paras$depth.mus = c('D3')
paras$depth.conf.factors = c('none')
paras$include.top.otu = FALSE
paras$models = 'loglinear'
paras$nSubs = 'Sub_L1'


setwd(resdir)
res <- clsapply(1:100, func, paras, queque='1-day', tempdir=file.path(resdir, 'tmpC'))
setwd(resdir)
save(res, file=file.path(resdir, paste(prefix, "_res.Rdata", sep="")))
