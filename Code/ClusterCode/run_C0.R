prefix <- 'C0'
resdir <-  file.path(paste0("/research/bsi/projects/staff_analysis/m216453/SemiSim/", prefix))
temp <- load('Stool_V35.RData',envir=.GlobalEnv)
temp0 <- load('Stool_V35_dirmult.RData',envir=.GlobalEnv)
#otu.tab <- Stool_V35
source('pre_run.R')

paras <- list()
paras$model.paras <- NULL
paras$resdir <- resdir
paras$prefix <- prefix
paras$nOTUs = c('nOTU_L1','nOTU_L5')
paras$nSams = c('nSam_L1','nSam_L4')
paras$diff.otu.pcts= c('none')
paras$diff.otu.modes = 'abundant'
paras$diff.otu.directs = c('balanced')
paras$covariate.types = c("binary")
paras$covariate.eff.means = c('none')
paras$confounder.types = c('none')
paras$depth.mus = c('D3')
paras$depth.conf.factors = c('none')
paras$include.top.otu = FALSE
paras$models = 'loglinear'
paras$nSubs = 'Sub_L1'


setwd(resdir)
clsapply(1:1000, func, paras, queque='1-day', tempdir=file.path(resdir, 'tmpC'))

