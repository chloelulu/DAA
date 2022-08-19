### Figure S26: add the most recently developed methods
pkg = c('dplyr','tidyr','reshape','RColorBrewer','ggplot2','ggpubr','stringr','scales','tibble','kableExtra','formattable','reactable','htmltools',"htmlwidgets","webshot2")
suppressPackageStartupMessages(sapply(pkg, require, character = T))

root <- '/Users/m216453'
setwd(paste0(root,'/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/SemiSimulation/Submit/DAA/'))
wd = getwd()
source(paste0(root,'/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Code/Function.R'))
output = paste0(wd,'/Result/SimulationEvaluation07022022/')

##############################
## Figure 3
sub.methods = c('ZicoSeq','ZINQ','fastANCOM','LinDA')
data.C1 <- clean_null(name = 'C0', type.name = 'Stool', sub.methods = sub.methods, root = root)
data.D1 <- clean_null(name = 'D0', type.name = 'Vaginal', sub.methods = sub.methods, root = root)

data <- rbind(data.C1 %>% dplyr::filter(depth.conf.factors=='None') %>% dplyr::select(-depth.conf.factors) %>% mutate(depth.conf.factors ='Stool'),
              data.D1 %>% dplyr::filter(depth.conf.factors=='None') %>% dplyr::select(-depth.conf.factors) %>%mutate(depth.conf.factors ='Vaginal'))

head(data);dim(data)
plot_FDR_kable(data, covariate.type='binary',kable.methods=sub.methods, type.name = 'Stool',output = output)

##############################
covariate.type = 'binary'
output = paste0(wd,'/Result/SimulationEvaluation07022022/Figure4/');if(!(dir.exists(output))){dir.create(output)}
##  balanced change: moderate: stool/vaginal
# Stool
res.obj1 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output, filename = 'C4', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L5',
                       sub.level.rename = "OTU=500",sub.methods = sub.methods)
plot.reactable8a(res.df2= res.obj1$res.df2, factor ='nOTUs',diff.otu.mode =c('Abundant','Rare'), output =output, name = 'A.Moderate.Stool')
# Vaginal
res.obj2 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'D4', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L5',
                       sub.level.rename =  "OTU=500",sub.methods = sub.methods)
plot.reactable8a(res.df2= res.obj2$res.df2, factor ='nOTUs',diff.otu.mode =c('Abundant','Rare'),output =output, name = 'B.Moderate.Vaginal')



## Figure 5. Unbalanced change: moderate: stool/vaginal
output = paste0(wd,'/Result/SimulationEvaluation07022022/Figure5/');if(!(dir.exists(output))){dir.create(output)}
# Stool
res.obj1 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'C54', factor ='nOTUs',
                       covariate.type =covariate.type, sub.level ='nOTU_L5',
                       sub.level.rename = 'OTU=500',sub.methods = sub.methods)
plot.reactable7b(res.df2= res.obj1$res.df2, factor ='nOTUs',diff.otu.mode =c('Abundant'),output =output, name = 'A.Moderate.Stool')

# Vaginal
res.obj2 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'D54', factor ='nOTUs',
                       covariate.type =covariate.type, sub.level ='nOTU_L5',
                       sub.level.rename = 'OTU=500',sub.methods = sub.methods)
plot.reactable7b(res.df2= res.obj2$res.df2, factor ='nOTUs',diff.otu.mode =c('Abundant'),output =output, name = 'B.Moderate.Vaginal')


## Supplementary Figure 5
output = paste0(wd,'/Result/SimulationEvaluation07022022/SupplementaryFigure5/');if(!(dir.exists(output))){dir.create(output)}
## ------combine stool and vaginal: small sample; OTU = 500; sample =50; effectsize = L3
res.obj1 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'C2', factor ='nSams',
                       covariate.type =covariate.type, sub.level ='nSam_L1',
                       sub.level.rename = 'sample=50',sub.methods = sub.methods)
plot.reactable8a(res.df2= res.obj1$res.df2, factor ='nSams',  diff.otu.mode =c('Abundant','Rare'),output =output, name = 'A.SmallSample.Stool')

res.obj2 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'D2', factor ='nSams',
                       covariate.type =covariate.type, sub.level ='nSam_L1',
                       sub.level.rename = 'sample=50',sub.methods = sub.methods)
plot.reactable8a(res.df2= res.obj2$res.df2, factor ='nSams', diff.otu.mode =c('Abundant','Rare'),output =output, name = 'B.SmallSample.Vaginal')

## ------combine stool and vaginal: small sample; OTU = 500; sample =50; effectsize = L3
res.obj1 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'C52', factor ='nSams',
                       covariate.type =covariate.type, sub.level ='nSam_L1',
                       sub.level.rename = 'sample=50',sub.methods = sub.methods)
plot.reactable7b(res.df2= res.obj1$res.df2, factor ='nSams',  diff.otu.mode =c('Abundant'),output =output, name = 'C.SmallSample.Stool')

res.obj2 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'D52', factor ='nSams',
                       covariate.type =covariate.type, sub.level ='nSam_L1',
                       sub.level.rename = 'sample=50',sub.methods = sub.methods)
plot.reactable7b(res.df2= res.obj2$res.df2, factor ='nSams',  diff.otu.mode =c('Abundant'),output =output, name = 'D.SmallSample.Vaginal')



## Supplementary Figure 6
output = paste0(wd,'/Result/SimulationEvaluation07022022/SupplementaryFigure6/');if(!(dir.exists(output))){dir.create(output)}
## ------combine stool and vaginal: small OTU; OTU = 50; sample =100; effectsize = L3
res.obj1 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'C4', factor ='nOTUs',
                       covariate.type =covariate.type, sub.level ='nOTU_L1',
                       sub.level.rename = 'OTU=50',sub.methods = sub.methods)
plot.reactable8a(res.df2= res.obj1$res.df2, factor ='nOTUs', diff.otu.mode =c('Abundant','Rare'),output =output, name = 'A.SmallOTU.Stool')

res.obj2 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'D4', factor ='nOTUs',
                       covariate.type =covariate.type, sub.level ='nOTU_L1',
                       sub.level.rename = 'OTU=50',sub.methods = sub.methods)
plot.reactable8a(res.df2= res.obj2$res.df2, factor ='nOTUs', diff.otu.mode =c('Abundant','Rare'),output =output, name = 'B.SmallOTU.Vaginal')

## ------combine stool and vaginal: small OTU; OTU = 50; sample =100; effectsize = L3
res.obj1 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'C54', factor ='nOTUs',
                       covariate.type =covariate.type, sub.level ='nOTU_L1',
                       sub.level.rename = 'OTU=50',sub.methods = sub.methods)
plot.reactable7b(res.df2= res.obj1$res.df2, factor ='nOTUs', diff.otu.mode =c('Abundant'),output =output, name = 'C.SmallOTU.Stool')

res.obj2 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'D54', factor ='nOTUs',
                       covariate.type =covariate.type, sub.level ='nOTU_L1',
                       sub.level.rename = 'OTU=50',sub.methods = sub.methods)
plot.reactable7b(res.df2= res.obj2$res.df2, factor ='nOTUs', diff.otu.mode =c('Abundant'),output =output, name = 'D.SmallOTU.Vaginal')








##############################
output = paste0(wd,'/Result/SimulationEvaluation07022022/Figure8/');if(!(dir.exists(output))){dir.create(output)}
sub.methods1 = sub.methods
## binary + covariate: non-compositional
# Stool
res.obj1 <- clean_data(dir = paste0(wd,'/Data/SimulationEvaluation/'),
                       output = output, filename = 'C4C', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L5',
                       sub.level.rename = "OTU=500",sub.methods = sub.methods1)
plot.reactable8a(res.df2= res.obj1$res.df2, factor ='nOTUs',diff.otu.mode =c('Abundant','Rare'), output =output, name = 'A.Moderate.Binary+covariate.Stool')
# Vaginal
res.obj2 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'D4D', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L5',
                       sub.level.rename =  "OTU=500",sub.methods = sub.methods1)
plot.reactable8a(res.df2= res.obj2$res.df2, factor ='nOTUs',diff.otu.mode =c('Abundant','Rare'),output =output, name = 'B.Moderate.Binary+covariate.Vaginal')


## binary + covariate: compositional
# Stool
res.obj1 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'C54C', factor ='nOTUs',
                       covariate.type =covariate.type, sub.level ='nOTU_L5',
                       sub.level.rename = 'OTU=500',sub.methods = sub.methods1)
plot.reactable7b(res.df2= res.obj1$res.df2, factor ='nOTUs',diff.otu.mode =c('Abundant'),output =output, name = 'C.Moderate.Binary+covariate.Stool')

# Vaginal
res.obj2 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'D54D', factor ='nOTUs',
                       covariate.type =covariate.type, sub.level ='nOTU_L5',
                       sub.level.rename = 'OTU=500',sub.methods = sub.methods1)
plot.reactable7b(res.df2= res.obj2$res.df2, factor ='nOTUs',diff.otu.mode =c('Abundant'),output =output, name = 'D.Moderate.Binary+covariate.Vaginal')




##############################
##Figure 9: Depth Confounding
output = paste0(wd,'/Result/SimulationEvaluation07022022/Figure9/');if(!(dir.exists(output))){dir.create(output)}
res.obj1 <- clean_data(dir = paste0(wd,'/Data/SimulationEvaluation/'),
                       output = output ,filename = 'C7', factor ='depth.conf.factors',
                       covariate.type ='binary', sub.level ='DL1',
                       sub.level.rename = 'Depth confounding +',sub.methods = sub.methods)
plot.reactable8a(res.df2= res.obj1$res.df2, factor ='depth.conf.factors', diff.otu.mode =c('Abundant','Rare'),output =output, name = 'Figure9.DepthCounfounding.Stool+')

output = paste0(wd,'/Result/SimulationEvaluation07022022/SupplementaryFigure7/');if(!(dir.exists(output))){dir.create(output)}
res.obj1 <- clean_data(dir = paste0(wd,'/Data/SimulationEvaluation/'),
                       output = output ,filename = 'C7', factor ='depth.conf.factors',
                       covariate.type ='binary', sub.level ='DL3',
                       sub.level.rename = 'Depth confounding ++',sub.methods = sub.methods)

plot.reactable8a(res.df2= res.obj1$res.df2, factor ='depth.conf.factors', diff.otu.mode =c('Abundant','Rare'),output =output, name = 'SupplementaryFigure7.DepthCounfounding.Stool++')








