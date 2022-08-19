pkg = c('dplyr','tidyr','reshape','RColorBrewer','ggplot2','ggpubr','stringr','scales','tibble','kableExtra','formattable','reactable','htmltools',"htmlwidgets","webshot2")
suppressPackageStartupMessages(sapply(pkg, require, character = T))

root <- '/Users/m216453'
setwd(paste0(root,'/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/SemiSimulation/Submit/DAA/'))
wd = getwd()
source(paste0(root,'/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Code/Function.R'))
output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/SemiSimulation/Submit/Revision1/SupplementaryTablesFigs/'

them <- theme(axis.text.x = element_blank(),
              axis.text.y = element_text(color="black", size = 20),
              axis.title = element_text(color="black", size = 20),
              axis.ticks.x = element_blank(),
              title = element_text(size = 20),
              legend.position = 'top',
              plot.margin = margin(1, 1, 1, 1, "cm"),
              legend.text = element_text(size = 20),
              strip.text = element_text(size = 20)) 

## ------- Fig. S10 (Figure 1 supplementary )-------
sub.methods = c("Rarefy+Wilcox","TSS+Wilcox",'GMPR+Wilcox', 'GMPR+DESeq2','GMPR+edgeR','Wrench+MSeq',
                "RAIDA","ANCOM-BC","DACOMP","LDM","Omnibus",'MaAsLin2',#'ZINQ','fastANCOM',
                "Aldex2(Wilcox)","GMPR+glm","corncob","eBay(Wilcox)")
data.C1 <- clean_null(name = 'C0', type.name = 'Stool', sub.methods = sub.methods, root = root)
data.D1 <- clean_null(name = 'D0', type.name = 'Vaginal', sub.methods = sub.methods, root = root)

data <- rbind(data.C1 %>% dplyr::filter(depth.conf.factors=='None') %>% dplyr::select(-depth.conf.factors) %>% mutate(depth.conf.factors ='Stool'),
              data.D1 %>% dplyr::filter(depth.conf.factors=='None') %>% dplyr::select(-depth.conf.factors) %>%mutate(depth.conf.factors ='Vaginal'))
df <- data %>% dplyr::filter(covariate.types==covariate.type) %>% dplyr::select(depth.conf.factors,nOTUs, nSams, methods, value) 
df[df$value>0,'ct'] = 1
df[df$value==0,'ct'] = 0
df = df %>% dplyr::select(-value) 
formula = paste0('ct ~ depth.conf.factors + nOTUs + nSams + methods')
kb <- data_summary(df, as.formula(formula))
ord <- kb[kb$depth.conf.factors == 'Stool' & kb$nOTUs == 'OTU=50' & kb$nSams =='sample=50',]
kb <- within(kb, methods <- factor(methods, levels = ord[order(ord$ct, decreasing = T),'methods']))
kb$nOTUs <- gsub('OTU=','taxa number = ', kb$nOTUs)
kb$nSams <- gsub('sample=','sample size = ', kb$nSams)
unique(kb$nSams)
kb <- within(kb, nOTUs <- factor(nOTUs, levels = c(unique(kb$nOTUs))))
kb <- within(kb, nSams <- factor(nSams, levels = c("sample size = 50", "sample size = 200")))


p1 <- ggplot(kb[kb$depth.conf.factors == 'Stool',],aes(x = methods, y = ct, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(nOTUs ~ nSams, scales = 'free_y') +
  scale_fill_manual(values = cols[sub.methods]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Stool')

p2 <- ggplot(kb[kb$depth.conf.factors != 'Stool',],aes(x = methods, y = ct, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(nOTUs ~ nSams, scales = 'free_y') +
  scale_fill_manual(values = cols[sub.methods]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them+ 
  ggtitle('Vaginal')
p12 <- ggarrange(p1, p2, nrow = 1, common.legend = T)
annotate_figure(p12, fig.lab = 'Fig.S10',fig.lab.pos = "top.left",
                fig.lab.size = 24, fig.lab.face = 'bold')

ggsave(file = paste0(output, 'Fig.S10.pdf'), width = 15, height =8)




## ------- Fig S12ab.(Figure 2 supplementary) -------
sub.methods = c("Rarefy+Wilcox","TSS+Wilcox",'GMPR+Wilcox', 
                'GMPR+DESeq2','GMPR+edgeR','Wrench+MSeq',
                "RAIDA","ANCOM-BC","DACOMP",#'ZINQ','fastANCOM',
                "LDM","Omnibus",'MaAsLin2',
                "GMPR+glm","Aldex2(Wilcox)","GMPR+glm",
                "corncob","eBay(Wilcox)")
covariate.type = 'binary'
##  balanced change: moderate: stool/vaginal
# Stool
res.obj1 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output, filename = 'C4', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L5', 
                       sub.level.rename = "OTU=500",sub.methods = sub.methods)
factor ='nOTUs';diff.otu.mode =c('Abundant','Rare');name = 'A.Moderate.Stool'
df <- res.obj1$res2 %>% dplyr::filter(diff.otu.modes %in% diff.otu.mode & measures %in% c('FDR','TPR')) %>% 
  dplyr::select(as.symbol(factor),diff.otu.pcts, diff.otu.modes, methods,measures,value)
kb <- data_summary(df, as.formula(paste0('value ~ measures + diff.otu.pcts + diff.otu.modes + methods+ ', factor)))
ord <- kb[kb$measures == 'FDR' & kb$diff.otu.pcts == 'Low density' & kb$diff.otu.modes =='Abundant',]
kb <- within(kb, methods <- factor(methods, levels = ord[order(ord$value, decreasing = T),'methods']))
kb$diff.otu.pcts <- gsub('\\ density','',kb$diff.otu.pcts)
kb <- within(kb, diff.otu.pcts <- factor(diff.otu.pcts, levels = c('Low','Medium','High')))
p1 <- ggplot(kb[kb$measures == 'FDR',],aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[sub.methods]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Stool')

p2 <- ggplot(kb[kb$measures != 'FDR',],aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[sub.methods]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Stool')



# Vaginal
res.obj1 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output, filename = 'D4', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L5', 
                       sub.level.rename = "OTU=500",sub.methods = sub.methods)
factor ='nOTUs';diff.otu.mode =c('Abundant','Rare')
df <- res.obj1$res2 %>% dplyr::filter(diff.otu.modes %in% diff.otu.mode & measures %in% c('FDR','TPR')) %>% 
  dplyr::select(as.symbol(factor),diff.otu.pcts, diff.otu.modes, methods,measures,value)
kb <- data_summary(df, as.formula(paste0('value ~ measures + diff.otu.pcts + diff.otu.modes + methods+ ', factor)))
ord <- kb[kb$measures == 'FDR' & kb$diff.otu.pcts == 'Low density' & kb$diff.otu.modes =='Abundant',]
kb <- within(kb, methods <- factor(methods, levels = ord[order(ord$value, decreasing = T),'methods']))
kb$diff.otu.pcts <- gsub('\\ density','',kb$diff.otu.pcts)
kb <- within(kb, diff.otu.pcts <- factor(diff.otu.pcts, levels = c('Low','Medium','High')))
kb <- within(kb, diff.otu.modes <- factor(diff.otu.modes, levels = c('Abundant','Rare')))

p3 <- ggplot(kb[kb$measures == 'FDR',],aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[sub.methods]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Vaginal')

p4 <- ggplot(kb[kb$measures != 'FDR',],aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts, scales = 'free_y') +
  scale_fill_manual(values = cols[sub.methods]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Vaginal')




## ------- Fig S12cd. (Figure 2 supplementary) -------
##  balanced change: moderate: stool/vaginal
res.obj1 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'C54', factor ='nOTUs', 
                       covariate.type =covariate.type, sub.level ='nOTU_L5', 
                       sub.level.rename = 'OTU=500',sub.methods = sub.methods)
factor ='nOTUs';diff.otu.mode =c('Abundant')
df1 <- res.obj1$res2 %>% dplyr::filter(diff.otu.modes %in% diff.otu.mode & measures %in% c('FDR','TPR')) %>% 
  dplyr::select(as.symbol(factor),diff.otu.pcts, diff.otu.modes, methods,measures,value)

res.obj1 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'D54', factor ='nOTUs', 
                       covariate.type =covariate.type, sub.level ='nOTU_L5', 
                       sub.level.rename = 'OTU=500',sub.methods = sub.methods)
factor ='nOTUs';diff.otu.mode =c('Abundant')
df2 <- res.obj1$res2 %>% dplyr::filter(diff.otu.modes %in% diff.otu.mode & measures %in% c('FDR','TPR')) %>% 
  dplyr::select(as.symbol(factor),diff.otu.pcts, diff.otu.modes, methods,measures,value)

df <- rbind(df1 %>% mutate(type = 'Stool'), df2 %>% mutate(type = 'Vaginal'))

kb <- data_summary(df, as.formula(paste0('value ~ measures + diff.otu.pcts + diff.otu.modes + methods+ type + ', factor)))
head(kb)
ord <- kb[kb$measures == 'FDR' & kb$diff.otu.pcts == 'Low density' & kb$type=='Stool',]
kb <- within(kb, methods <- factor(methods, levels = ord[order(ord$value, decreasing = T),'methods']))
kb$diff.otu.pcts <- gsub('\\ density','',kb$diff.otu.pcts)
kb <- within(kb, diff.otu.pcts <- factor(diff.otu.pcts, levels = c('Low','Medium','High')))
kb <- within(kb, measures <- factor(measures, levels = c('FDR','TPR')))

p51 <- ggplot(kb[kb$type == 'Stool' & kb$measures == 'FDR',],aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[sub.methods]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Stool')
p52 <- ggplot(kb[kb$type == 'Stool' & kb$measures == 'TPR',],aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[sub.methods]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Stool')

p61 <- ggplot(kb[kb$type == 'Vaginal' & kb$measures == 'FDR',],aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[sub.methods]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Vaginal')
p62 <- ggplot(kb[kb$type == 'Vaginal' & kb$measures == 'TPR',],aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[sub.methods]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Vaginal')
p123456 <- ggarrange(p1, p2, p3, p4, p51, p52, p61, p62, 
                     labels = c('a','','b','','c','','d',''),
                     font.label = list(size = 24),
                     ncol = 2, nrow = 4, common.legend = T,
                     heights = c(1,1,0.7, 0.7))
annotate_figure(p123456, fig.lab = 'Fig.S12',fig.lab.pos = "top.left",
                fig.lab.size = 24, fig.lab.face = 'bold')
ggsave(file = paste0(output, 'Fig.S12.pdf'), width = 14, height = 18)






## ------------ Fig. S13 (Figure 3 barplot)  ------------ 
sub.methods = c("Rarefy+Wilcox","TSS+Wilcox",'GMPR+Wilcox', 'GMPR+DESeq2','GMPR+edgeR','Wrench+MSeq','ZicoSeq',#'fastANCOM','ZINQ',
                "RAIDA","ANCOM-BC","DACOMP","LDM","Omnibus",'MaAsLin2',"GMPR+glm","Aldex2(Wilcox)","GMPR+glm","corncob","eBay(Wilcox)")
# Stool - binary --Balanced
res.obj1 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output, filename = 'C4', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L5', 
                       sub.level.rename = "OTU=500",sub.methods = sub.methods)
best1 <- cal.best(res.df2= res.obj1$res.df2 %>% filter(methods != 'ZicoSeq'), factor ='nOTUs',diff.otu.mode =c('Abundant','Rare'))
res1 <- res.obj1$res2
# Vaginal - binary --Balanced
res.obj2 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'D4', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L5', 
                       sub.level.rename = "OTU=500",sub.methods = sub.methods)
best2 <- cal.best(res.df2= res.obj2$res.df2%>% filter(methods != 'ZicoSeq'), factor ='nOTUs',diff.otu.mode =c('Abundant','Rare'))
res2 <- res.obj2$res2
### Unbalanced
### moderate
# Stool - binary --unBalanced
res.obj3 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'C54', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L5', 
                       sub.level.rename = "OTU=500",sub.methods = sub.methods)
best3 <- cal.best1(res.df2= res.obj3$res.df2 %>% filter(methods != 'ZicoSeq'), factor ='nOTUs',diff.otu.mode =c('Abundant'))
res3 <- res.obj3$res2
# Vaginal - binary --unBalanced
res.obj4 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'D54', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L5', 
                       sub.level.rename = "OTU=500",sub.methods = sub.methods)
best4 <- cal.best1(res.df2= res.obj4$res.df2%>% filter(methods != 'ZicoSeq'), factor ='nOTUs',diff.otu.mode =c('Abundant'))
res4 <- res.obj4$res2

res <- rbind(res1 %>% mutate(type = 'Stool', grp = 'Balanced') %>% dplyr::filter(methods %in% c(best1,'ZicoSeq')), 
             res2 %>% mutate(type = 'Vaginal', grp = 'Balanced') %>% dplyr::filter(methods %in% c(best2,'ZicoSeq')), 
             res3 %>% mutate(type = 'Stool', grp = 'Unbalanced') %>% dplyr::filter(methods %in% c(best3,'ZicoSeq')), 
             res4 %>% mutate(type = 'Vaginal', grp = 'Unbalanced') %>% dplyr::filter(methods %in% c(best4,'ZicoSeq'))) %>%
  dplyr::filter(measures %in% c('TPR','FDR')) %>% droplevels()
head(res)
res$diff.otu.pcts <- gsub('\\ density','',res$diff.otu.pcts)
res <- within(res, diff.otu.pcts <- factor(diff.otu.pcts, levels = c('Low','Medium','High')))
summary(res)
res0 <- data_summary(res, value ~ diff.otu.modes + diff.otu.pcts + methods + measures + type + grp)


p1 <- ggplot(res0 %>% dplyr::filter(type == 'Stool' & grp == 'Balanced' & measures =='FDR'),
       aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best1,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Stool:Balanced')

p2 <- ggplot(res0 %>% dplyr::filter(type == 'Stool' & grp == 'Balanced' & measures =='TPR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best1,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  theme(legend.position = 'top') + 
  ggtitle('Stool:Balanced')

p3 <- ggplot(res0 %>% dplyr::filter(type == 'Vaginal' & grp == 'Balanced' & measures =='FDR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best2,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Vaginal:Balanced')

p4 <- ggplot(res0 %>% dplyr::filter(type == 'Vaginal' & grp == 'Balanced' & measures =='TPR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best2,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Vaginal:Balanced')


p51 <- ggplot(res0 %>% dplyr::filter(type == 'Stool' & grp == 'Unbalanced' & diff.otu.modes == 'Abundant' & measures == 'FDR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best3,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Stool:Unbalanced')

p52 <- ggplot(res0 %>% dplyr::filter(type == 'Stool' & grp == 'Unbalanced' & diff.otu.modes == 'Abundant' & measures == 'TPR'),
              aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best3,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Stool:Unbalanced')

p61 <- ggplot(res0 %>% dplyr::filter(type == 'Vaginal' & grp == 'Unbalanced' & diff.otu.modes == 'Abundant' & measures == 'FDR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best4,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them+ 
  ggtitle('Vaginal:Unbalanced')
p62 <- ggplot(res0 %>% dplyr::filter(type == 'Vaginal' & grp == 'Unbalanced' & diff.otu.modes == 'Abundant' & measures == 'TPR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best4,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  them+ 
  ggtitle('Vaginal:Unbalanced')


p12 <- ggarrange(p1, p2, common.legend = T, nrow = 1, labels = 'a',font.label = list(size = 24, color = "black", face = "bold", family = NULL))
p34 <- ggarrange(p3, p4, common.legend = T, nrow = 1, labels = '',font.label = list(size = 24, color = "black", face = "bold", family = NULL))
p5 <- ggarrange(p51, p52,common.legend =T, nrow = 1, ncol = 2, labels = c('b'),font.label = list(size = 24, color = "black", face = "bold", family = NULL))
p6 <- ggarrange(p61,p62, common.legend = T, nrow = 1, ncol = 2, labels = c(''),font.label = list(size = 24, color = "black", face = "bold", family = NULL))

pdf(paste0(output,'Fig.S16ab.pdf'), width = 12, height = 18)
p.ab <- ggarrange(p12, p34, p5, p6, ncol = 1, nrow =4, heights = c(1,1,0.7,0.7))
annotate_figure(p.ab, fig.lab = 'Fig.S16',fig.lab.pos = "top.left",fig.lab.size = 24, fig.lab.face = 'bold')
dev.off()




## small sample  
# Stool - binary --Balanced
res.obj1 <- clean_data(dir =paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'C2', factor ='nSams', covariate.type =covariate.type, sub.level ='nSam_L1', 
                       sub.level.rename = "sample=50",sub.methods = sub.methods)
best1 <- cal.best(res.df2= res.obj1$res.df2 %>% filter(methods != 'ZicoSeq'), factor ='nSams',diff.otu.mode =c('Abundant','Rare'))
res1 <- res.obj1$res2
# Vaginal - binary --Balanced
res.obj2 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'D2', factor ='nSams', covariate.type =covariate.type, sub.level ='nSam_L1', 
                       sub.level.rename = "sample=50",sub.methods = sub.methods)
best2 <- cal.best(res.df2= res.obj2$res.df2%>% filter(methods != 'ZicoSeq'), factor ='nSams',diff.otu.mode =c('Abundant','Rare'))
res2 <- res.obj2$res2
### Unbalanced
### moderate
# Stool - binary --unBalanced
res.obj3 <- clean_data(dir =paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'C52', factor ='nSams', covariate.type =covariate.type, sub.level ='nSam_L1', 
                       sub.level.rename = "sample=50",sub.methods = sub.methods)
best3 <- cal.best1(res.df2= res.obj3$res.df2 %>% filter(methods != 'ZicoSeq'), factor ='nSams',diff.otu.mode =c('Abundant'))
res3 <- res.obj3$res2
# Vaginal - binary --unBalanced
res.obj4 <- clean_data(dir =paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'D52', factor ='nSams', covariate.type =covariate.type, sub.level ='nSam_L1', 
                       sub.level.rename = "sample=50",sub.methods = sub.methods)
best4 <- cal.best1(res.df2= res.obj4$res.df2%>% filter(methods != 'ZicoSeq'), factor ='nSams',diff.otu.mode =c('Abundant'))
res4 <- res.obj4$res2
# arrange
res.cd <- rbind(res1 %>% mutate(type = 'Stool', grp = 'Balanced') %>% dplyr::filter(methods %in% c(best1,'ZicoSeq')), 
                       res2 %>% mutate(type = 'Vaginal', grp = 'Balanced') %>% dplyr::filter(methods %in% c(best2,'ZicoSeq')), 
                       res3 %>% mutate(type = 'Stool', grp = 'Unbalanced') %>% dplyr::filter(methods %in% c(best3,'ZicoSeq')), 
                       res4 %>% mutate(type = 'Vaginal', grp = 'Unbalanced') %>% dplyr::filter(methods %in% c(best4,'ZicoSeq'))) %>%
  dplyr::filter(measures %in% c('TPR','FDR')) %>% droplevels()

res.cd$diff.otu.pcts <- gsub('\\ density','',res.cd$diff.otu.pcts)
res.cd <- within(res.cd, diff.otu.pcts <- factor(diff.otu.pcts, levels = c('Low','Medium','High')))
summary(res)
res.cd0 <- data_summary(res.cd, value ~ diff.otu.modes + diff.otu.pcts + methods + measures + type + grp)
p1 <- ggplot(res.cd0 %>% dplyr::filter(type == 'Stool' & grp == 'Balanced' & measures =='FDR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts ,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best1,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them+ 
  ggtitle('Stool: Balanced')

p2 <- ggplot(res.cd0 %>% dplyr::filter(type == 'Stool' & grp == 'Balanced' & measures =='TPR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts ,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best1,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Stool: Balanced')

p3 <- ggplot(res.cd0 %>% dplyr::filter(type == 'Vaginal' & grp == 'Balanced' & measures =='FDR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts ,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best2,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them+ 
  ggtitle('Vaginal: Balanced')

p4 <- ggplot(res.cd0 %>% dplyr::filter(type == 'Vaginal' & grp == 'Balanced' & measures =='TPR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts ,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best2,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Vaginal: Balanced')


p51 <- ggplot(res.cd0 %>% dplyr::filter(type == 'Stool' & grp == 'Unbalanced' & diff.otu.modes == 'Abundant' & measures == 'FDR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts ,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best3,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Stool: Unbalanced')

p52 <- ggplot(res.cd0 %>% dplyr::filter(type == 'Stool' & grp == 'Unbalanced' & diff.otu.modes == 'Abundant' & measures == 'TPR'),
              aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts ,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best3,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Stool: Unbalanced')

p61 <- ggplot(res.cd0 %>% dplyr::filter(type == 'Vaginal' & grp == 'Unbalanced' & diff.otu.modes == 'Abundant' & measures == 'FDR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts ,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best4,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Vaginal: Unbalanced')

p62 <- ggplot(res.cd0 %>% dplyr::filter(type == 'Vaginal' & grp == 'Unbalanced' & diff.otu.modes == 'Abundant' & measures == 'TPR'),
              aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts ,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best4,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Vaginal: Unbalanced')


p12 <- ggarrange(p1, p2, common.legend = T, nrow = 1, labels = 'c',font.label = list(size = 24, color = "black", face = "bold", family = NULL))
p34 <- ggarrange(p3, p4, common.legend = T, nrow = 1, labels = '',font.label = list(size = 24, color = "black", face = "bold", family = NULL))
p5 <- ggarrange(p51, p52,common.legend =T, nrow = 1, ncol = 2, labels = c('d'),font.label = list(size = 24, color = "black", face = "bold", family = NULL))
p6 <- ggarrange(p61, p62, common.legend = T, nrow = 1, ncol = 2, labels = c(''),font.label = list(size = 24, color = "black", face = "bold", family = NULL))

pdf(paste0(output,'Fig.S16cd.pdf'), width = 12, height = 18)
p.ab <- ggarrange(p12, p34, p5, p6, ncol = 1, nrow =4, heights = c(1,1,0.7,0.7))
annotate_figure(p.ab, fig.lab = 'Fig.S16',fig.lab.pos = "top.left",fig.lab.size = 24, fig.lab.face = 'bold')
dev.off()

### small OTU
# Stool - binary --Balanced
res.obj1 <- clean_data(dir =paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'C4', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L1', 
                       sub.level.rename = "OTU=50",sub.methods = sub.methods)
best1 <- cal.best(res.df2= res.obj1$res.df2 %>% filter(methods != 'ZicoSeq'), factor ='nOTUs',diff.otu.mode =c('Abundant','Rare'))
res1 <- res.obj1$res2
# Vaginal - binary --Balanced
res.obj2 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'D4', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L1', 
                       sub.level.rename = "OTU=50",sub.methods = sub.methods)
best2 <- cal.best(res.df2= res.obj2$res.df2%>% filter(methods != 'ZicoSeq'), factor ='nOTUs',diff.otu.mode =c('Abundant','Rare'))
res2 <- res.obj2$res2
### Unbalanced
### moderate
# Stool - binary --unBalanced
res.obj3 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'C54', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L1', 
                       sub.level.rename = "OTU=50",sub.methods = sub.methods)
best3 <- cal.best1(res.df2= res.obj3$res.df2 %>% filter(methods != 'ZicoSeq'), factor ='nOTUs',diff.otu.mode =c('Abundant'))
res3 <- res.obj3$res2
# Vaginal - binary --unBalanced
res.obj4 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'D54', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L1', 
                       sub.level.rename = "OTU=50",sub.methods = sub.methods)
best4 <- cal.best1(res.df2= res.obj4$res.df2%>% filter(methods != 'ZicoSeq'), factor ='nOTUs',diff.otu.mode =c('Abundant'))
res4 <- res.obj4$res2

# arrange
res.ef <- rbind(res1 %>% mutate(type = 'Stool', grp = 'Balanced') %>% dplyr::filter(methods %in% c(best1,'ZicoSeq')), 
             res2 %>% mutate(type = 'Vaginal', grp = 'Balanced') %>% dplyr::filter(methods %in% c(best2,'ZicoSeq')), 
             res3 %>% mutate(type = 'Stool', grp = 'Unbalanced') %>% dplyr::filter(methods %in% c(best3,'ZicoSeq')), 
             res4 %>% mutate(type = 'Vaginal', grp = 'Unbalanced') %>% dplyr::filter(methods %in% c(best4,'ZicoSeq'))) %>%
  dplyr::filter(measures %in% c('TPR','FDR')) %>% droplevels()
head(res.ef)
res.ef$diff.otu.pcts <- gsub('\\ density','',res.ef$diff.otu.pcts)
res.ef <- within(res.ef, diff.otu.pcts <- factor(diff.otu.pcts, levels = c('Low','Medium','High')))
summary(res)
res.ef0 <- data_summary(res.ef, value ~ diff.otu.modes + diff.otu.pcts + methods + measures + type + grp)


p1 <- ggplot(res.ef0 %>% dplyr::filter(type == 'Stool' & grp == 'Balanced' & measures =='FDR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts ,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best1,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Stool: Balanced')

p2 <- ggplot(res.ef0 %>% dplyr::filter(type == 'Stool' & grp == 'Balanced' & measures =='TPR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts ,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best1,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Stool: Balanced')

p3 <- ggplot(res.ef0 %>% dplyr::filter(type == 'Vaginal' & grp == 'Balanced' & measures =='FDR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts ,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best2,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Vaginal: Balanced')

p4 <- ggplot(res.ef0 %>% dplyr::filter(type == 'Vaginal' & grp == 'Balanced' & measures =='TPR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts ,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best2,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Vaginal: Balanced')


p51 <- ggplot(res.ef0 %>% dplyr::filter(type == 'Stool' & grp == 'Unbalanced' & diff.otu.modes == 'Abundant' & measures == 'FDR'),
              aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts ,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best3,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Stool: Unbalanced')

p52 <- ggplot(res.ef0 %>% dplyr::filter(type == 'Stool' & grp == 'Unbalanced' & diff.otu.modes == 'Abundant' & measures == 'TPR'),
              aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts ,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best3,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Stool: Unbalanced')

p61 <- ggplot(res.ef0 %>% dplyr::filter(type == 'Vaginal' & grp == 'Unbalanced' & diff.otu.modes == 'Abundant' & measures == 'FDR'),
              aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts ,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best4,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Vaginal: Unbalanced')

p62 <- ggplot(res.ef0 %>% dplyr::filter(type == 'Vaginal' & grp == 'Unbalanced' & diff.otu.modes == 'Abundant' & measures == 'TPR'),
              aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts ,  scales = 'free_y') +
  scale_fill_manual(values = cols[c(best4,'ZicoSeq')]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Vaginal: Unbalanced')


p12 <- ggarrange(p1, p2, common.legend = T, nrow = 1, labels = 'e',font.label = list(size = 24, color = "black", face = "bold", family = NULL))
p34 <- ggarrange(p3, p4, common.legend = T, nrow = 1, labels = '',font.label = list(size = 24, color = "black", face = "bold", family = NULL))
p5 <- ggarrange(p51, p52,common.legend =T, nrow = 1, ncol = 2, labels = c('f'),font.label = list(size = 24, color = "black", face = "bold", family = NULL))
p6 <- ggarrange(p61,p62, common.legend = T, nrow = 1, ncol = 2, labels = c(''),font.label = list(size = 24, color = "black", face = "bold", family = NULL))

pdf(paste0(output,'Fig.S16ef.pdf'), width = 12, height = 18)
p.ab <- ggarrange(p12, p34, p5, p6, ncol = 1, nrow =4, heights = c(1,1,0.7,0.7))
annotate_figure(p.ab, fig.lab = 'Fig.S16',fig.lab.pos = "top.left",fig.lab.size = 24, fig.lab.face = 'bold')
dev.off()






## ------ Figure S18: (coresponding to: Figure 4)------
sub.methods1 = c('GMPR+DESeq2','GMPR+edgeR',"ANCOM-BC","LDM","Aldex2(glm)","GMPR+glm",'MaAsLin2',"corncob",'ZicoSeq')

## binary + covariate: non-compositional
# Stool
res.obj1 <- clean_data(dir = paste0(wd,'/Data/SimulationEvaluation/old/'),
                       output = output, filename = 'C4C', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L5', 
                       sub.level.rename = "OTU=500",sub.methods = sub.methods1)
res1 <- res.obj1$res2
# Vaginal
res.obj2 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/old/'),
                       output = output ,filename = 'D4D', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L5',
                       sub.level.rename =  "OTU=500",sub.methods = sub.methods1)
res2 <- res.obj2$res2
## binary + covariate: compositional
# Stool
res.obj1 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/old/'),
                       output = output ,filename = 'C54C', factor ='nOTUs', 
                       covariate.type =covariate.type, sub.level ='nOTU_L5', 
                       sub.level.rename = 'OTU=500',sub.methods = sub.methods1)
res3 <- res.obj1$res2
# Vaginal
res.obj2 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/old/'),
                       output = output ,filename = 'D54D', factor ='nOTUs', 
                       covariate.type =covariate.type, sub.level ='nOTU_L5', 
                       sub.level.rename = 'OTU=500',sub.methods = sub.methods1)
res4 <- res.obj2$res2

# arrange
res.confounder <- rbind(res1 %>% mutate(type = 'Stool', grp = 'Balanced') %>% dplyr::filter(methods %in% sub.methods1), 
                res2 %>% mutate(type = 'Vaginal', grp = 'Balanced') %>% dplyr::filter(methods %in% sub.methods1), 
                res3 %>% mutate(type = 'Stool', grp = 'Unbalanced') %>% dplyr::filter(methods %in% sub.methods1), 
                res4 %>% mutate(type = 'Vaginal', grp = 'Unbalanced') %>% dplyr::filter(methods %in% sub.methods1)) %>%
  dplyr::filter(measures %in% c('TPR','FDR')) %>% droplevels()
head(res.confounder)
res.confounder$diff.otu.pcts <- gsub('\\ density','',res.confounder$diff.otu.pcts)
res.confounder <- within(res.confounder, diff.otu.pcts <- factor(diff.otu.pcts, levels = c('Low','Medium','High')))
res.confounder0 <- data_summary(res.confounder, value ~ diff.otu.modes + diff.otu.pcts + methods + measures + type + grp)

p1 <- ggplot(res.confounder0 %>% dplyr::filter(type == 'Stool' & grp == 'Balanced' & measures =='FDR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[sub.methods1]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Stool: Balanced')

p2 <- ggplot(res.confounder0 %>% dplyr::filter(type == 'Stool' & grp == 'Balanced' & measures =='TPR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[sub.methods1]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Stool: Balanced')

p3 <- ggplot(res.confounder0 %>% dplyr::filter(type == 'Vaginal' & grp == 'Balanced' & measures =='FDR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[sub.methods1]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Vaginal: Balanced')

p4 <- ggplot(res.confounder0 %>% dplyr::filter(type == 'Vaginal' & grp == 'Balanced' & measures =='TPR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[sub.methods1]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Vaginal: Balanced')


p51 <- ggplot(res.confounder0 %>% dplyr::filter(type == 'Stool' & grp == 'Unbalanced' & diff.otu.modes == 'Abundant' & measures == 'FDR'),
              aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts ,  scales = 'free_y') +
  scale_fill_manual(values = cols[sub.methods1]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Stool: Unbalanced')

p52 <- ggplot(res.confounder0 %>% dplyr::filter(type == 'Stool' & grp == 'Unbalanced' & diff.otu.modes == 'Abundant' & measures == 'TPR'),
              aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts ,  scales = 'free_y') +
  scale_fill_manual(values = cols[sub.methods1]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Stool: Unbalanced')

p61 <- ggplot(res.confounder0 %>% dplyr::filter(type == 'Vaginal' & grp == 'Unbalanced' & diff.otu.modes == 'Abundant' & measures == 'FDR'),
              aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts ,  scales = 'free_y') +
  scale_fill_manual(values = cols[sub.methods1]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Vaginal: Unbalanced')

p62 <- ggplot(res.confounder0 %>% dplyr::filter(type == 'Vaginal' & grp == 'Unbalanced' & diff.otu.modes == 'Abundant' & measures == 'TPR'),
              aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts ,  scales = 'free_y') +
  scale_fill_manual(values = cols[sub.methods1]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('Vaginal: Unbalanced')


graphics.off()
p.ab <- ggarrange(p1, p2, p3, p4, p51, p52,p61,p62, common.legend = T, nrow = 4, ncol = 2, 
          labels = c('a','','b','','c','','d',''), heights = c(1,1,0.7,0.7),
          font.label = list(size = 24, color = "black", face = "bold", 
                            family = NULL)) + 
  theme(plot.margin = margin(1,0.2,0.1,0.1, "cm")) 
annotate_figure(p.ab, fig.lab = 'Fig.S18',fig.lab.pos = "top.left",fig.lab.size = 24, fig.lab.face = 'bold')
ggsave(paste0(output,'Fig.S18.pdf'), width = 12, height = 18)









##---------Fig.S19: barplot (corresponding to Fig 5 Depth Confounding)------
res.obj1 <- clean_data(dir = paste0(wd,'/Data/SimulationEvaluation/'),
                       output = output ,filename = 'C7', factor ='depth.conf.factors',
                       covariate.type ='binary', sub.level ='DL1',
                       sub.level.rename = 'Depth confounding +',sub.methods = sub.methods)
res1 <- res.obj1$res2 %>% dplyr::filter(diff.otu.mode %in% c('Abundant','Rare'))
res1$diff.otu.pcts <- gsub('\\ density','',res1$diff.otu.pcts)
res1 <- within(res1, diff.otu.pcts <- factor(diff.otu.pcts, levels = c('Low','Medium','High')))
head(res1)
res10 <- data_summary(res1, value ~ diff.otu.modes + diff.otu.pcts + methods + measures)

p1 <- ggplot(res10 %>% dplyr::filter(measures =='FDR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[sub.methods]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('')
p2 <- ggplot(res10 %>% dplyr::filter(measures =='TPR'),
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[sub.methods]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  them + 
  ggtitle('')
p.ab <- ggarrange(p1, p2, nrow = 1, common.legend = T) + 
  theme(plot.margin = margin(1,0.2,0.1,0.1, "cm")) 
annotate_figure(p.ab, fig.lab = 'Fig.S19',fig.lab.pos = "top.left",fig.lab.size = 24, fig.lab.face = 'bold')
ggsave(paste0(output,'Fig.S19.pdf'), width = 12, height =7)






## add in revision1 for sample size = 5000
library(dplyr)
root <- '/Users/m216453'
source(paste0(root,'/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Code/Function.R'))

load('/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/SimulationEvaluation/largesample_C_res.Rdata')
res0 <- reshape2::melt(res)
colnames(res0) <- c('depth.mus', 'diff.otu.modes', 'model', 'nSub', 'nOTUs', 'covariate.types', 'depth.conf.factors', 'diff.otu.directs', 'confounder.types', 'nSams', 'covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters')
res.C <- data_summary(res0, value ~ diff.otu.directs + diff.otu.modes + diff.otu.pcts + nSams + covariate.eff.means + measures + methods)
head(res.C)


load('/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/SimulationEvaluation/largesample_D_res1.Rdata')
res0 <- reshape2::melt(res)
colnames(res0) <- c('depth.mus', 'diff.otu.modes', 'model', 'nSub', 'nOTUs', 'covariate.types', 'depth.conf.factors', 'diff.otu.directs', 'confounder.types', 'nSams', 'covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters')
res.DB <- data_summary(res0, value ~ diff.otu.directs + diff.otu.modes + diff.otu.pcts + nSams + covariate.eff.means + measures + methods)


load('/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/SimulationEvaluation/largesample_CB_res1.Rdata')
res0 <- reshape2::melt(res)
colnames(res0) <- c('depth.mus', 'diff.otu.modes', 'model', 'nSub', 'nOTUs', 'covariate.types', 'depth.conf.factors', 'diff.otu.directs', 'confounder.types', 'nSams', 'covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters')
res.CB <- data_summary(res0, value ~ diff.otu.directs + diff.otu.modes + diff.otu.pcts + nSams + covariate.eff.means + measures + methods)

res.D <- rbind(res.CB, res.DB)

res.CD <- rbind(res.C %>% mutate(type = 'Stool') %>% dplyr::filter(!(diff.otu.directs=='unbalanced' & diff.otu.modes == 'rare')), 
                res.D %>% mutate(type = 'Vaginal')) %>% dplyr::filter(nSams == 'nSam_L7') %>% droplevels()

Stool <- res.CD %>% dplyr::filter(type =='Stool' & covariate.eff.means =='L2') %>% droplevels()
Vaginal <- res.CD %>% dplyr::filter(type !='Stool' & covariate.eff.means =='L1') %>% droplevels()

p1 <- ggplot(Stool[Stool$measures == 'TPR',], aes(x = diff.otu.pcts, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_wrap(diff.otu.directs ~ diff.otu.modes, nrow = 1) +
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  scale_fill_manual(values = brewer.pal(9,'Paired')[5]) +
  labs(y = 'TPR', x = '',  fill = '') + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 1, size = 20, colour = 'black'),
        axis.text.y = element_text(color="black", size = 20),
        axis.title = element_text(color="black", size = 20),
        axis.ticks.x = element_blank(),
        title = element_text(size = 20),
        legend.position = 'top',
        plot.margin = margin(1, 1, 1, 1, "cm"),
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 20)) +
  ggtitle('a. Stool(n=5000)')

p2 <- ggplot(Stool[Stool$measures == 'FDR',], aes(x = diff.otu.pcts, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_wrap(diff.otu.directs ~ diff.otu.modes, nrow = 1) +
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  scale_fill_manual(values = brewer.pal(9,'Paired')[5]) +
  labs(y = 'FDR', x = '', fill = '') + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 1, size = 20, colour = 'black'),
        axis.text.y = element_text(color="black", size = 20),
        axis.title = element_text(color="black", size = 20),
        axis.ticks.x = element_blank(),
        title = element_text(size = 20),
        legend.position = 'top',
        plot.margin = margin(1, 1, 1, 1, "cm"),
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 20)) +
  ggtitle('')

p3 <- ggplot(Vaginal[Vaginal$measures == 'TPR',], aes(x = diff.otu.pcts, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_wrap(diff.otu.directs ~ diff.otu.modes, nrow = 1) +
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  scale_fill_manual(values = brewer.pal(9,'Paired')[5]) +
  labs(y = 'TPR', x = '',  fill = '') + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 1, size = 20, colour = 'black'),
        axis.text.y = element_text(color="black", size = 20),
        axis.title = element_text(color="black", size = 20),
        axis.ticks.x = element_blank(),
        title = element_text(size = 20),
        legend.position = 'top',
        plot.margin = margin(1, 1, 1, 1, "cm"),
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 20)) +
  ggtitle('b. Vaginal(n=5000)')

p4 <- ggplot(Vaginal[Vaginal$measures == 'FDR',], aes(x = diff.otu.pcts, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_wrap(diff.otu.directs ~ diff.otu.modes, nrow = 1) +
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  scale_fill_manual(values = brewer.pal(9,'Paired')[5]) +
  labs(y = 'FDR', x = '', fill = '') + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 1, size = 20, colour = 'black'),
        axis.text.y = element_text(color="black", size = 20),
        axis.title = element_text(color="black", size = 20),
        axis.ticks.x = element_blank(),
        title = element_text(size = 20),
        legend.position = 'top',
        plot.margin = margin(1, 1, 1, 1, "cm"),
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 20)) +
  ggtitle('')

ggarrange(p1, p2, p3, p4, common.legend = T, nrow = 2, ncol =2)
getwd()
ggsave(file = paste0('Result/SimulationEvaluation03122022/ZicoSeq5000.pdf'), width = 13, height = 10)





## Fig.S6 (Reviewer's suggestion in changing ref.cutoff of ZicoSeq)
root <- '/Users/m216453'
setwd(paste0(root,'/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/SemiSimulation/Submit/DAA/'))
wd <- paste0(root,'/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/SemiSimulation/Submit/DAA/',"Data/SimulationEvaluation/zicoseq4010ref/")
sub.methods = c("ZicoSeq","ZicoSeq3")
covariate.type = 'binary'



##  balanced change: moderate: stool/vaginal
# Stool
res.obj1 <- clean_data(dir = wd,
                       output = output, filename = 'C4', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L5', 
                       sub.level.rename = "OTU=500",sub.methods = sub.methods)
# Vaginal
res.obj2 <- clean_data(dir = wd,
                       output = output ,filename = 'D4', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L5',
                       sub.level.rename = "OTU=500",sub.methods = sub.methods)

plot.ref.balanced(stool.df = res.obj1$res.df2, vaginal.df = res.obj1$res.df2, 
                  fig.lab = 'Balanced(moderate)', output = output)


## Figure 5. Unbalanced change: moderate: stool/vaginal
# Stool
res.obj1 <- clean_data(dir = wd,
                       output = output ,filename = 'C54', factor ='nOTUs', 
                       covariate.type =covariate.type, sub.level ='nOTU_L5', 
                       sub.level.rename = 'OTU=500',sub.methods = sub.methods)

# Vaginal
res.obj2 <- clean_data(dir = wd,
                       output = output ,filename = 'D54', factor ='nOTUs', 
                       covariate.type =covariate.type, sub.level ='nOTU_L5', 
                       sub.level.rename = 'OTU=500',sub.methods = sub.methods)

plot.ref.balanced(stool.df = res.obj1$res.df2, vaginal.df = res.obj1$res.df2, 
                  balanced = F,
                  fig.lab = 'Unbalanced(moderate)', output = output)

## ------combine stool and vaginal: small sample; OTU = 500; sample =50; effectsize = L3
res.obj1 <- clean_data(dir = wd,
                       output = output ,filename = 'C2', factor ='nSams', 
                       covariate.type =covariate.type, sub.level ='nSam_L1', 
                       sub.level.rename = 'sample=50',sub.methods = sub.methods)

res.obj2 <- clean_data(dir = wd,
                       output = output ,filename = 'D2', factor ='nSams', 
                       covariate.type =covariate.type, sub.level ='nSam_L1', 
                       sub.level.rename = 'sample=50',sub.methods = sub.methods)
plot.ref.balanced(stool.df = res.obj1$res.df2, vaginal.df = res.obj1$res.df2, 
                  fig.lab = 'Balanced(sample=50)', output = output)



## ------combine stool and vaginal: small sample; OTU = 500; sample =50; effectsize = L3
res.obj1 <- clean_data(dir = wd,
                       output = output ,filename = 'C52', factor ='nSams', 
                       covariate.type =covariate.type, sub.level ='nSam_L1', 
                       sub.level.rename = 'sample=50',sub.methods = sub.methods)

res.obj2 <- clean_data(dir = wd,
                       output = output ,filename = 'D52', factor ='nSams', 
                       covariate.type =covariate.type, sub.level ='nSam_L1', 
                       sub.level.rename = 'sample=50',sub.methods = sub.methods)
plot.ref.balanced(stool.df = res.obj1$res.df2, vaginal.df = res.obj1$res.df2, 
                  balanced = F, fig.lab = 'Unbalanced(sample=50)', output = output)






## ------combine stool and vaginal: small OTU; OTU = 50; sample =100; effectsize = L3
res.obj1 <- clean_data(dir = wd,
                       output = output ,filename = 'C4', factor ='nOTUs', 
                       covariate.type =covariate.type, sub.level ='nOTU_L1', 
                       sub.level.rename = 'OTU=50',sub.methods = sub.methods)

res.obj2 <- clean_data(dir = wd,
                       output = output ,filename = 'D4', factor ='nOTUs', 
                       covariate.type =covariate.type, sub.level ='nOTU_L1', 
                       sub.level.rename = 'OTU=50',sub.methods = sub.methods)

plot.ref.balanced(stool.df = res.obj1$res.df2, vaginal.df = res.obj1$res.df2, 
                  fig.lab = 'Balanced(taxa number=50)',output = output)



## ------combine stool and vaginal: small OTU; OTU = 50; sample =100; effectsize = L3
res.obj1 <- clean_data(dir = wd,
                       output = output ,filename = 'C54', factor ='nOTUs', 
                       covariate.type =covariate.type, sub.level ='nOTU_L1', 
                       sub.level.rename = 'OTU=50',sub.methods = sub.methods)

res.obj2 <- clean_data(dir = wd,
                       output = output ,filename = 'D54', factor ='nOTUs', 
                       covariate.type =covariate.type, sub.level ='nOTU_L1', 
                       sub.level.rename = 'OTU=50',sub.methods = sub.methods)
plot.ref.balanced(stool.df = res.obj1$res.df2, vaginal.df = res.obj1$res.df2, 
                  balanced = F, fig.lab = 'Unbalanced(taxa number=50)', output = output)




##############################
## binary + covariate: non-compositional
# Stool
res.obj1 <- clean_data(dir = wd,
                       output = output, filename = 'C4C', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L5', 
                       sub.level.rename = "OTU=500",sub.methods = sub.methods)
# Vaginal
res.obj2 <- clean_data(dir = wd,
                       output = output ,filename = 'D4D', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L5',
                       sub.level.rename =  "OTU=500",sub.methods = sub.methods)
plot.ref.balanced(stool.df = res.obj1$res.df2, vaginal.df = res.obj1$res.df2, 
                  fig.lab = 'Balanced(moderate) binary+confouder',output = output)

## binary + covariate: compositional
# Stool
res.obj1 <- clean_data(dir = wd,
                       output = output, filename = 'C54C', factor ='nOTUs', 
                       covariate.type =covariate.type, sub.level ='nOTU_L5', 
                       sub.level.rename = 'OTU=500',sub.methods = sub.methods)

# Vaginal
res.obj2 <- clean_data(dir =wd,
                       output = output ,filename = 'D54D', factor ='nOTUs', 
                       covariate.type =covariate.type, sub.level ='nOTU_L5', 
                       sub.level.rename = 'OTU=500',sub.methods = sub.methods)
plot.ref.balanced(stool.df = res.obj1$res.df2, vaginal.df = res.obj1$res.df2, 
                  balanced = F,
                  fig.lab = 'Unbalanced(moderate) binary+confouder',output = output)



##############################
## Fig. S19: (corresponding to Fig 5: Depth Confounding)
res.obj1 <- clean_data(dir = wd,
                       output = output ,filename = 'C7', factor ='depth.conf.factors',
                       covariate.type ='binary', sub.level ='DL1',
                       sub.level.rename = 'Depth confounding +',sub.methods = sub.methods)
res.obj2 <- clean_data(dir = wd,
                       output = output ,filename = 'C7', factor ='depth.conf.factors',
                       covariate.type ='binary', sub.level ='DL3',
                       sub.level.rename = 'Depth confounding ++',sub.methods = sub.methods)










## Fig S17: large sample size 100 for ZicoSeq to reviewer's suggestion
source(paste0(root,'/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Code/Function.R'))
setwd('/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/SimulationEvaluation/largesample1000')
load('largesample_C1000B_res.Rdata')
res0 <- reshape2::melt(res)
colnames(res0) <- c('depth.mus', 'diff.otu.modes', 'model', 'nSub', 'nOTUs', 'covariate.types', 'depth.conf.factors', 'diff.otu.directs', 'confounder.types', 'nSams', 'covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters')
res.CB <- data_summary(res0, value ~ diff.otu.directs + diff.otu.modes + diff.otu.pcts + nSams + covariate.eff.means + measures + methods)
res.CB$diff.otu.modes <- str_to_title(res.CB$diff.otu.modes)
res.CB$diff.otu.pcts <- str_to_title(res.CB$diff.otu.pcts)
res.CB <- within(res.CB, diff.otu.pcts <- factor(diff.otu.pcts, levels = c('Low','Medium','High')))

load('largesample_C1000U_res.Rdata')
res0 <- reshape2::melt(res)
colnames(res0) <- c('depth.mus', 'diff.otu.modes', 'model', 'nSub', 'nOTUs', 'covariate.types', 'depth.conf.factors', 'diff.otu.directs', 'confounder.types', 'nSams', 'covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters')
res.CU <- data_summary(res0, value ~ diff.otu.directs + diff.otu.modes + diff.otu.pcts + nSams + covariate.eff.means + measures + methods)
res.CU$methods <- gsub('MSeq2.Wrench','Wrench+MSeq',res.CU$methods)
res.CU$diff.otu.modes <- str_to_title(res.CU$diff.otu.modes)
res.CU$diff.otu.pcts <- str_to_title(res.CU$diff.otu.pcts)
res.CU <- within(res.CU, diff.otu.pcts <- factor(diff.otu.pcts, levels = c('Low','Medium','High')))


load('largesample_D1000B_res.Rdata')
res0 <- reshape2::melt(res)
colnames(res0) <- c('depth.mus', 'diff.otu.modes', 'model', 'nSub', 'nOTUs', 'covariate.types', 'depth.conf.factors', 'diff.otu.directs', 'confounder.types', 'nSams', 'covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters')
res.DB <- data_summary(res0, value ~ diff.otu.directs + diff.otu.modes + diff.otu.pcts + nSams + covariate.eff.means + measures + methods)
res.DB$diff.otu.modes <- str_to_title(res.DB$diff.otu.modes)
res.DB$diff.otu.pcts <- str_to_title(res.DB$diff.otu.pcts)
res.DB <- within(res.DB, diff.otu.pcts <- factor(diff.otu.pcts, levels = c('Low','Medium','High')))
res.DB$methods <- gsub('MSeq2.Wrench','Wrench+MSeq',res.DB$methods)


load('largesample_D1000U_res.Rdata')
res0 <- reshape2::melt(res)
colnames(res0) <- c('depth.mus', 'diff.otu.modes', 'model', 'nSub', 'nOTUs', 'covariate.types', 'depth.conf.factors', 'diff.otu.directs', 'confounder.types', 'nSams', 'covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters')
res.DU <- data_summary(res0, value ~ diff.otu.directs + diff.otu.modes + diff.otu.pcts + nSams + covariate.eff.means + measures + methods)
res.DU$methods <- gsub('MSeq2.Wrench','Wrench+MSeq',res.DU$methods)
res.DU$diff.otu.modes <- str_to_title(res.DU$diff.otu.modes)
res.DU$diff.otu.pcts <- str_to_title(res.DU$diff.otu.pcts)
res.DU <- within(res.DU, diff.otu.pcts <- factor(diff.otu.pcts, levels = c('Low','Medium','High')))


them1 <- theme(axis.text.x = element_blank(),
      axis.text.y = element_text(color="black", size = 20),
      axis.title = element_text(color="black", size = 20),
      axis.ticks.x = element_blank(),
      title = element_text(size = 20),
      legend.position = 'top',
      plot.margin = margin(1, 1, 1, 1, "cm"),
      legend.text = element_text(size = 20),
      strip.text = element_text(size = 20)) 

p1 <- ggplot(res.CB %>% dplyr::filter(measures == 'TPR'), 
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  scale_y_continuous(limits = c(0,1)) + 
  scale_fill_manual(values = cols[as.character(unique(res.CB$methods))]) + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  ggtitle('Stool:Balanced (n = 1000)') + them1
p2 <- ggplot(res.DB %>% dplyr::filter(measures == 'TPR'), 
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  scale_y_continuous(limits = c(0,1)) + 
  scale_fill_manual(values = cols[as.character(unique(res.DB$methods))]) + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  ggtitle('Vaginal:Balanced (n = 1000)') + them1
p3 <- ggplot(res.CU %>% dplyr::filter(diff.otu.modes == 'Abundant') %>% dplyr::filter(measures == 'TPR'), 
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  scale_y_continuous(limits = c(0,1)) + 
  scale_fill_manual(values = cols[as.character(unique(res.CU$methods))]) + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  ggtitle('Stool:Unbalanced (n = 1000)') + them1
p4 <- ggplot(res.DU %>% dplyr::filter(diff.otu.modes == 'Abundant') %>% dplyr::filter(measures == 'TPR'), 
             aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  scale_y_continuous(limits = c(0,1)) + 
  scale_fill_manual(values = cols[as.character(unique(res.DU$methods))]) + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'TPR', x = '', color = "", fill = '') + theme_bw() + 
  ggtitle('Vaginal:Unbalanced (n = 1000)') + them1
p11 <- ggplot(res.CB %>% dplyr::filter(measures == 'FDR'), 
              aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[as.character(unique(res.CB$methods))]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  ggtitle('Stool:Balanced (n = 1000)') + them1
p22 <- ggplot(res.DB %>% dplyr::filter(measures == 'FDR'), 
              aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[as.character(unique(res.DB$methods))]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  ggtitle('Vaginal:Balanced (n = 1000)') + them1
p33 <- ggplot(res.CU %>% dplyr::filter(diff.otu.modes == 'Abundant') %>% dplyr::filter(measures == 'FDR'), 
              aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[as.character(unique(res.CU$methods))]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  ggtitle('Stool:Unbalanced (n = 1000)') + them1
p44 <- ggplot(res.DU %>% dplyr::filter(diff.otu.modes == 'Abundant') %>% dplyr::filter(measures == 'FDR'), 
              aes(x = methods, y = value, fill = methods)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(diff.otu.modes ~ diff.otu.pcts,  scales = 'free_y') +
  scale_fill_manual(values = cols[as.character(unique(res.DU$methods))]) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,
                position=position_dodge2(.7, preserve = "single")) +
  labs(y = 'FDR', x = '', color = "", fill = '') + theme_bw() + 
  ggtitle('Vaginal:Unbalanced (n = 1000)') + them1

# p.ab <- ggarrange(p11, p1,p22, p2,p33, p3,p44, p4,nrow = 4, ncol = 2, common.legend = F, 
#                   labels = c('a','','','','b','','',''), heights = c(1,1,0.7,0.7),
#                   font.label = list(size = 24, color = "black", face = "bold", 
#                                     family = NULL)) + 
#   theme(plot.margin = margin(1,0.2,0.1,0.1, "cm"))
# annotate_figure(p.ab, fig.lab = 'Fig.S17',fig.lab.pos = "top.left",fig.lab.size = 24, fig.lab.face = 'bold')
# ggsave(file = paste0(output,'FigS17.pdf'),width = 15, height = 18)

p12 <- ggarrange(p11, p1, common.legend = T, nrow = 1, labels = 'a',font.label = list(size = 24, color = "black", face = "bold", family = NULL))
p34 <- ggarrange(p22, p2, common.legend = T, nrow = 1, labels = '',font.label = list(size = 24, color = "black", face = "bold", family = NULL))
p56 <- ggarrange(p33, p3,p44, p4, common.legend = T, nrow = 2, ncol = 2, labels = c('b'),font.label = list(size = 24, color = "black", face = "bold", family = NULL))

pdf(paste0(output,'Fig.S17.pdf'), width = 12, height = 18)
p.ab <- ggarrange(p12, p34, p56, ncol = 1, nrow =3, heights = c(1,1,1.5))
annotate_figure(p.ab, fig.lab = 'Fig.S17',fig.lab.pos = "top.left",fig.lab.size = 24, fig.lab.face = 'bold')
dev.off()





# Fig. S22 run time for ZicoSeq n=1000
setwd('/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/SimulationEvaluation/')
load('largesample_C_res.Rdata')
res0 <- reshape2::melt(res)
colnames(res0) <- c('depth.mus', 'diff.otu.modes', 'model', 'nSub', 'nOTUs', 'covariate.types', 'depth.conf.factors', 'diff.otu.directs', 'confounder.types', 'nSams', 'covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters')
res0 <- res0 %>% dplyr::filter(measures == 'time' & diff.otu.modes =='abundant' & 
                                 diff.otu.directs =='balanced'& covariate.eff.means =='L1' & diff.otu.pcts=='low') %>% 
  droplevels()
dim(res0)
res.time <- data_summary(res0, value ~ nSams + methods)
(res.time$value)/60
res0$nSams <- gsub('nSam_L6','n=1000', res0$nSams)
res0$nSams <- gsub('nSam_L7','n=5000', res0$nSams)

ggplot(res0, aes(x = nSams, y = value/60, fill = nSams)) + 
  geom_boxplot() + 
  scale_fill_brewer(palette = 'Set2') + 
  theme_bw() + 
  labs(x = 'sample size', y = 'time(mins)', fill = '') +
  theme(axis.text.x = element_text(color="black", size = 20),
        axis.text.y = element_text(color="black", size = 20),
        axis.title = element_text(color="black", size = 20),
        axis.ticks.x = element_blank(),
        title = element_text(size = 20),
        legend.position = 'none',
        plot.margin = margin(1, 1, 1, 1, "cm"),
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 20)) +
  ggtitle('Fig.S22')
ggsave(file = paste0('/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/SemiSimulation/Submit/Revision1/SupplementaryTablesFigs/','Fig.S22.pdf'),
       width = 5, height = 5)





## Fig.S23: stability evaluated by speraman correlation of pvalues between filter and no filter
pkg =  c('ggpubr','microbiome',"eBay","modeest","ANCOMBC","aod","phyloseq","reshape","corncob","MASS","readr","DESeq2", "ALDEx2", "metagenomeSeq", "edgeR", "GUniFrac", "grDevices", "dirmult", "exactRankTests","nlme", "dplyr", "magrittr", "tidyr", "protoclust", "ggplot2", "compositions","rmutil","tibble","reticulate","dacomp","LDM","Wrench")
sapply(pkg, require, character = TRUE)
# setwd('~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/stability/')

# on cluster folder /research/bsi/projects/staff_analysis/m216453/SemiSim/filter; and /research/bsi/projects/staff_analysis/m216453/SemiSim/nofilter
setwd('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/stability/') # 03/22/2021 new add
methods = c('ZicoSeq','glmquassi2','Rarefy', 'Wilcox','Wilcox.Wrench','Wilcox.gmpr','mbzinb','eBayW','eBayt','Aldex2we','Aldex2','RAIDA', 'DACOMP', 'LDM', 'glmquassi','BBinomial','ANCOMBC','DESeq2', 'DESeq2.Wrench', 'DESeq2.gmpr','edgeR', 'edgeR.Wrench', 'edgeR.gmpr','MSeq2','MSeq2.Wrench','Maaslin2')

P = EST = DF = list()
for(method in methods){
  try({
    files = NULL
    for(i in 1:100){
      file = paste0("FILTER_summarymatrix",i,"-",method,".Rdata")
      files = c(files, file)
    }
    p = est = NULL
    df2 = NULL
    
    for(file in files){
      load(paste0('NOFILTER/NO',file))
      res_seqs0 = res_seqs
      names0 = names(res_seqs0)
      names0 = names0[grep('binary',names0)] %>% na.omit()
      
      load(paste0('FILTER/',file))
      res_seqs40 = res_seqs
      names40 = names(res_seqs40)
      names40 = names40[grep('binary',names40)] %>% na.omit()
      names = intersect(names0, names40)
      names = names[grep('nOTU_L5',names)]
      for(name in names){
        df0 = res_seqs0[[name]]
        df40 = res_seqs40[[name]]
        cat(nrow(df0),':', nrow(df40),'\n')
        df = inner_join(df0 %>% dplyr::select(otu.id, fdr) %>% dplyr::rename(nofilter = fdr),
                        df40 %>% dplyr::select(otu.id, fdr) %>% dplyr::rename(filter = fdr))
        df1 = df %>% mutate(data = name)
        df2 = rbind(df2, df1)
        cor = cor.test(df$nofilter,df$filter, method = 'spearman')
        p = c(p,cor$p.value)
        est = c(est,cor$estimate)
      }
    }
  })
  
  DF[[method]] = df2
  P[[method]] = p
  EST[[method]] = est
}

save(P,file= 'correlationP3.Rdata') # ended with 3 means the 3rd version.
save(EST,file= 'correlationR3.Rdata')
save(DF,file= 'correlationDataframe3.Rdata')
save(P, EST, DF, file ='stability3.Rdata')

load('~/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/SimulationEvaluation/stability3.Rdata')
corr= NULL
for(i in 1:length(EST)){
  cor = as.data.frame(EST[[i]])
  colnames(cor) = 'R'
  cor$method = names(EST)[i]
  corr = rbind(corr, cor)
}
head(corr)

unique(corr$method)


corr.sum <- corr %>% filter(method != 'Omnibus')
colnames(corr.sum) = c('Stability','methods')
corr.sum$methods = as.character(corr.sum$methods)
corr.sum$methods = gsub('^Wilcox$','TSS+Wilcox',corr.sum$methods)
corr.sum$methods = gsub('^Rarefy$','Rarefy+Wilcox',corr.sum$methods)
corr.sum$methods = gsub('ANCOMBC','ANCOM-BC',corr.sum$methods)
corr.sum$methods = gsub('glmquassi2','Wrench+glm',corr.sum$methods)
corr.sum$methods = gsub('glmquassi','GMPR+glm',corr.sum$methods)
corr.sum$methods = gsub('DESeq2.gmpr','GMPR+DESeq2',corr.sum$methods)
corr.sum$methods = gsub('DESeq2.Wrench','Wrench+DESeq2',corr.sum$methods)
corr.sum$methods = gsub('edgeR.gmpr','GMPR+edgeR',corr.sum$methods)
corr.sum$methods = gsub('edgeR.Wrench','Wrench+edgeR',corr.sum$methods)
corr.sum$methods = gsub('MSeq2.Wrench','Wrench+MSeq',corr.sum$methods)
corr.sum$methods = gsub('^MSeq2$','metagenomeSeq',corr.sum$methods)
corr.sum$methods = gsub('eBayW','eBay(Wilcox)',corr.sum$methods)
corr.sum$methods = gsub('eBayt','eBay(t-test)',corr.sum$methods)
corr.sum$methods = gsub('BBinomial','corncob',corr.sum$methods)
corr.sum$methods = gsub('Aldex2we','Aldex2(t-test)',corr.sum$methods)
corr.sum$methods = gsub('^Aldex2$','Aldex2(Wilcox)',corr.sum$methods)
corr.sum$methods = gsub('^ttest$',"TSS+t-test",corr.sum$methods)
corr.sum$methods = gsub('^Rarefyttest$','Rarefy+t-test',corr.sum$methods)
corr.sum$methods = gsub('linda',"LinDA",corr.sum$methods)
corr.sum$methods = gsub('mbzinb',"Omnibus",corr.sum$methods)
corr.sum$methods = gsub('Wilcox.gmpr',"GMPR+Wilcox",corr.sum$methods)
corr.sum$methods = gsub('Wilcox.Wrench',"Wrench+Wilcox",corr.sum$methods)
corr.sum$methods = gsub('Maaslin2',"MaAsLin2",corr.sum$methods)

head(corr.sum)
corr.sum.mean = data_summary(corr.sum, Stability ~ methods)
corr.sum.mean = corr.sum.mean[order(-corr.sum.mean$median),]

corr.sum <- within(corr.sum, methods <- factor(methods, levels=corr.sum.mean$methods))

ggplot(corr.sum, aes(x = methods, y = Stability, color = methods)) + 
  theme_bw() +
  stat_boxplot(geom='errorbar', linetype=1, width=0.5)+
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(aes(colour = methods), width = 0.2, cex = 0.6, alpha = 0.4) +
  scale_color_manual(values = cols) +
  # geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.4, size = 0.5,position=position_dodge2(.7, preserve = "single")) +
  # scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
  labs(y = 'time(seconds)', x = '', color = "", fill = '') +
  theme(axis.text.x = element_text(color="black", size = 20, angle = 90, vjust = 1, hjust = 1),
        axis.text.y = element_text(color="black", size = 26),
        axis.title = element_text(color="black", size = 24),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 30),
        plot.margin = margin(1,1,1,1,'cm'),
        strip.background = element_rect(fill="white",color = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size= 1),
        legend.position = 'none')
ggsave(paste0(wd,'/Result/SimulationEvaluation/SupplementaryFigure23.pdf'), width = 9,height = 6)







##############################
## Figure S26: add the most recent developed methods after the submission of paper
sub.methods = c('ZicoSeq','ZINQ','fastANCOM','LinDA')
data.C1 <- clean_null(name = 'C0', type.name = 'Stool', sub.methods = sub.methods, root = root)
data.D1 <- clean_null(name = 'D0', type.name = 'Vaginal', sub.methods = sub.methods, root = root)

data <- rbind(data.C1 %>% dplyr::filter(depth.conf.factors=='None') %>% dplyr::select(-depth.conf.factors) %>% mutate(depth.conf.factors ='Stool'),
              data.D1 %>% dplyr::filter(depth.conf.factors=='None') %>% dplyr::select(-depth.conf.factors) %>%mutate(depth.conf.factors ='Vaginal'))

head(data);dim(data)
plot_FDR_kable(data, covariate.type='binary',kable.methods=sub.methods, type.name = 'Stool',output = output)

##############################
covariate.type = 'binary'
output = paste0(wd,'/Result/SimulationEvaluation07022022/FigureS26ab/');if(!(dir.exists(output))){dir.create(output)}
## FigureS26ab balanced change: moderate: stool/vaginal
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



## FigureS26cd. Unbalanced change: moderate: stool/vaginal
output = paste0(wd,'/Result/SimulationEvaluation07022022/FigureS26cd/');if(!(dir.exists(output))){dir.create(output)}
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


