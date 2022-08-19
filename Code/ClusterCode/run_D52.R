prefix <- 'D52'
resdir <-  file.path(paste0("/research/bsi/projects/staff_analysis/m216453/SemiSim/", prefix))
temp <- load('VaginalIntroitus_V35.RData',envir=.GlobalEnv)
temp0 <- load('VaginalIntroitus_V35_dirmult.RData',envir=.GlobalEnv)
#otu.tab <- Stool_V35
source('pre_run.R')

paras <- list()
paras$model.paras <- NULL
paras$resdir <- resdir
paras$prefix <- prefix
paras$nOTUs = c('nOTU_L5')#,'nOTU_L2','nOTU_L3','nOTU_L4','nOTU_L5')
paras$nSams = c('nSam_L1','nSam_L4')
paras$diff.otu.pcts= c('low', 'medium', 'high')
paras$diff.otu.modes = c('abundant', 'rare')
paras$diff.otu.directs = c('unbalanced')
paras$covariate.types = c("continuous","binary")
paras$covariate.eff.means = c('L3')
paras$confounder.types = c('none')
paras$depth.mus = c('D3')
paras$depth.conf.factors = c('none')
paras$include.top.otu = FALSE
paras$models = 'loglinear'
paras$nSubs = 'Sub_L1'


setwd(resdir)
clsapply(1:100, func, paras, queque='1-hour', tempdir=file.path(resdir, 'tmpC'))

