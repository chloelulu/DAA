### code for produce Figure 3-10, Supplementary Figure 4-9 in the manuscript
pkg = c('dplyr','tidyr','reshape','RColorBrewer','ggplot2','ggpubr','stringr','scales','tibble','kableExtra','formattable','reactable','htmltools',"htmlwidgets","webshot2")
suppressPackageStartupMessages(sapply(pkg, require, character = T))

root <- '/Users/m216453'
setwd(paste0(root,'/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/SemiSimulation/Submit/DAA/'))
wd = getwd()
source(paste0(root,'/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Code/Function.R'))
output = paste0(wd,'/Result/SimulationEvaluation07192022/')

##############################
## Fig. 1
sub.methods = c("Rarefy+Wilcox","TSS+Wilcox",'GMPR+Wilcox', 'GMPR+DESeq2','GMPR+edgeR','Wrench+MSeq',
                "RAIDA","ANCOM-BC","DACOMP","LDM","Omnibus",'MaAsLin2','ZINQ','fastANCOM','IFAA','RDB',
                "Aldex2(Wilcox)","GMPR+glm","corncob","eBay(Wilcox)")
data.C1 <- clean_null(name = 'C0', type.name = 'Stool', sub.methods = sub.methods, root = root)
data.D1 <- clean_null(name = 'D0', type.name = 'Vaginal', sub.methods = sub.methods, root = root)

data <- rbind(data.C1 %>% dplyr::filter(depth.conf.factors=='None') %>% dplyr::select(-depth.conf.factors) %>% mutate(depth.conf.factors ='Stool'),
              data.D1 %>% dplyr::filter(depth.conf.factors=='None') %>% dplyr::select(-depth.conf.factors) %>%mutate(depth.conf.factors ='Vaginal'))

head(data);dim(data)
plot_FDR_kable(data, covariate.type='binary',kable.methods=sub.methods, type.name = 'Stool',output = output)

##############################
#### Figure S11 : wrench vs default for DESeq2 and edgeR
sub.methods = c('GMPR+DESeq2','DESeq2','edgeR','GMPR+edgeR')
data.C1 <- clean_null(name = 'C0', type.name = 'Stool', sub.methods = sub.methods, root = root)
data.D1 <- clean_null(name = 'D0', type.name = 'Vaginal', sub.methods = sub.methods, root = root)

data <- rbind(data.C1 %>% dplyr::filter(depth.conf.factors=='None') %>% dplyr::select(-depth.conf.factors) %>% mutate(depth.conf.factors ='Stool'),
              data.D1 %>% dplyr::filter(depth.conf.factors=='None') %>% dplyr::select(-depth.conf.factors) %>%mutate(depth.conf.factors ='Vaginal')) %>% 
  filter(covariate.types=='binary', methods %in% sub.methods) %>% dplyr::select(-covariate.types)
head(data)
data[data$value>0,'ct'] = 1
data[data$value==0,'ct'] = 0
data = data %>% dplyr::select(-value)

df <- data_summary(data, ct~ depth.conf.factors + nOTUs + nSams + methods)
df$methods <- as.factor(df$methods)
df <- within(df, methods <- factor(methods, levels=c("DESeq2","GMPR+DESeq2","edgeR",'GMPR+edgeR')))
head(df)

ggplot(df, aes(x = nOTUs, y = ct,  fill = methods)) +
  theme_bw() +
  geom_bar(position = position_dodge2(width = 0.7, preserve = "single"), stat="identity", width = 0.7) +
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,position=position_dodge2(.7, preserve = "single")) +
  scale_fill_manual(values = cols[sub.methods]) +
  # scale_y_continuous(limits = c(0,1),expand = c(0.1, 0, 0, 0)) +
  facet_grid(nSams~depth.conf.factors, scales = 'free_y')+
  labs(y = 'FDR', x = '', color = "", fill = '') +
  theme(axis.text.x = element_text(color="black", size =26),
        axis.text.y = element_text(color="black", size = 26),
        axis.title = element_text(color="black", size = 30),
        strip.text = element_text(size = 30),
        strip.background = element_rect(fill="white",color = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size= 1),
        legend.title=element_text(size=30),
        legend.text = element_text(size=20),
        plot.title = element_text(size=22))
ggsave(paste0(output,'NULL/','SupplementaryFigure11.pdf'), width =12, height = 7)




##############################
sub.methods = c("Rarefy+Wilcox","TSS+Wilcox",'GMPR+Wilcox', 
                'GMPR+DESeq2','GMPR+edgeR','Wrench+MSeq',
                "RAIDA","ANCOM-BC","DACOMP",
                "LDM","Omnibus",'MaAsLin2',
                "GMPR+glm","Aldex2(Wilcox)","GMPR+glm",
                "corncob","eBay(Wilcox)")
covariate.type = 'binary'
output = paste0(wd,'/Result/SimulationEvaluation07192022/Figure4/');if(!(dir.exists(output))){dir.create(output)}
##  Fig 2ab. balanced change: moderate: stool/vaginal
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



## Fig 2cd. Unbalanced change: moderate: stool/vaginal
output = paste0(wd,'/Result/SimulationEvaluation07192022/Figure5/');if(!(dir.exists(output))){dir.create(output)}
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


## Fig S13
output = paste0(wd,'/Result/SimulationEvaluation07192022/SupplementaryFigure5/');if(!(dir.exists(output))){dir.create(output)}
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



## Fig S14
output = paste0(wd,'/Result/SimulationEvaluation07192022/SupplementaryFigure6/');if(!(dir.exists(output))){dir.create(output)}
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
## FigS15: ZicoSeq NULL setting
sub.methods = c('ZicoSeq')
data.C1 <- clean_null(name = 'C0', type.name = 'Stool', sub.methods = sub.methods, root = root)
data.D1 <- clean_null(name = 'D0', type.name = 'Vaginal', sub.methods = sub.methods, root = root)

df <- rbind(data.C1 %>% dplyr::filter(depth.conf.factors=='None') %>% dplyr::select(-depth.conf.factors) %>% mutate(depth.conf.factors ='Stool'),
              data.D1 %>% dplyr::filter(depth.conf.factors=='None') %>% dplyr::select(-depth.conf.factors) %>%mutate(depth.conf.factors ='Vaginal')) %>% 
  filter(covariate.types=='binary', methods %in% sub.methods) %>% dplyr::select(-covariate.types)
df[df$value>0,'ct'] = 1
df[df$value==0,'ct'] = 0
df = df %>% dplyr::select(-value)
df <- data_summary(df, ct ~ depth.conf.factors + nOTUs + nSams + methods)
head(df)
ggplot(df, aes(x = nOTUs, y = ct,  fill = methods)) +
  theme_bw() +
  geom_bar(position = position_dodge2(width = 0.7, preserve = "single"), stat="identity", width = 0.7) +
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.5, size = 0.2,position=position_dodge2(.7, preserve = "single")) +
  scale_fill_manual(values = cols) +
  # scale_y_continuous(limits = c(0,1),expand = c(0.1, 0, 0, 0)) +
  facet_grid(nSams~depth.conf.factors, scales = 'free_y')+
  labs(y = 'FDR', x = '', color = "", fill = '') +
  theme(axis.text.x = element_text(color="black", size =26),
        axis.text.y = element_text(color="black", size = 26),
        axis.title = element_text(color="black", size = 30),
        strip.text = element_text(size = 30),
        strip.background = element_rect(fill="white",color = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size= 1),
        legend.title=element_text(size=30),
        legend.text = element_text(size=20),legend.position = 'none',
        plot.title = element_text(size=22))
ggsave(paste0(root,"/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Result/SimulationEvaluation07192022/NULL/",'FigS15.pdf'), width =9, height = 7)


##############################
output = paste0(wd,'/Result/SimulationEvaluation07192022/Figure7/');if(!(dir.exists(output))){dir.create(output)}
## Fig3: ZicoSeq compare to others
sub.methods = c("Rarefy+Wilcox","TSS+Wilcox",'GMPR+Wilcox', 'GMPR+DESeq2','GMPR+edgeR','Wrench+MSeq','ZicoSeq',
                "RAIDA","ANCOM-BC","DACOMP","LDM","Omnibus",'MaAsLin2',"GMPR+glm","Aldex2(Wilcox)","GMPR+glm","corncob","eBay(Wilcox)")
# Stool - binary --Balanced
res.obj1 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output, filename = 'C4', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L5', 
                       sub.level.rename = "OTU=500",sub.methods = sub.methods)
best1 <- cal.best(res.df2= res.obj1$res.df2 %>% filter(methods != 'ZicoSeq'), factor ='nOTUs',diff.otu.mode =c('Abundant','Rare'))
# Vaginal - binary --Balanced
res.obj2 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'D4', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L5', 
                       sub.level.rename = "OTU=500",sub.methods = sub.methods)
best2 <- cal.best(res.df2= res.obj2$res.df2%>% filter(methods != 'ZicoSeq'), factor ='nOTUs',diff.otu.mode =c('Abundant','Rare'))

### Unbalanced
### moderate
# Stool - binary --unBalanced
res.obj3 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'C54', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L5', 
                       sub.level.rename = "OTU=500",sub.methods = sub.methods)
best3 <- cal.best1(res.df2= res.obj3$res.df2 %>% filter(methods != 'ZicoSeq'), factor ='nOTUs',diff.otu.mode =c('Abundant'))
# Vaginal - binary --unBalanced
res.obj4 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'D54', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L5', 
                       sub.level.rename = "OTU=500",sub.methods = sub.methods)
best4 <- cal.best1(res.df2= res.obj4$res.df2%>% filter(methods != 'ZicoSeq'), factor ='nOTUs',diff.otu.mode =c('Abundant'))

# arrange
res.obj <- list(df1 = res.obj1$res.df2, df2 = res.obj2$res.df2,df3 = res.obj3$res.df2, df4 = res.obj4$res.df2)

best <- list(best1, best2, best3, best4)

plot.reactable0222(res.obj, best, output =output, factor = 'nOTUs',name = paste0('Moderate.',covariate.type))





### small OTU
# Stool - binary --Balanced
res.obj1 <- clean_data(dir =paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'C4', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L1', 
                       sub.level.rename = "OTU=50",sub.methods = sub.methods)
best1 <- cal.best(res.df2= res.obj1$res.df2 %>% filter(methods != 'ZicoSeq'), factor ='nOTUs',diff.otu.mode =c('Abundant','Rare'))
# Vaginal - binary --Balanced
res.obj2 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'D4', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L1', 
                       sub.level.rename = "OTU=50",sub.methods = sub.methods)
best2 <- cal.best(res.df2= res.obj2$res.df2%>% filter(methods != 'ZicoSeq'), factor ='nOTUs',diff.otu.mode =c('Abundant','Rare'))

### Unbalanced
### moderate
# Stool - binary --unBalanced
res.obj3 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'C54', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L1', 
                       sub.level.rename = "OTU=50",sub.methods = sub.methods)
best3 <- cal.best1(res.df2= res.obj3$res.df2 %>% filter(methods != 'ZicoSeq'), factor ='nOTUs',diff.otu.mode =c('Abundant'))
# Vaginal - binary --unBalanced
res.obj4 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'D54', factor ='nOTUs', covariate.type =covariate.type, sub.level ='nOTU_L1', 
                       sub.level.rename = "OTU=50",sub.methods = sub.methods)
best4 <- cal.best1(res.df2= res.obj4$res.df2%>% filter(methods != 'ZicoSeq'), factor ='nOTUs',diff.otu.mode =c('Abundant'))

# arrange
res.obj <- list(df1 = res.obj1$res.df2, df2 = res.obj2$res.df2,df3 = res.obj3$res.df2, df4 = res.obj4$res.df2)

best <- list(best1, best2, best3, best4)

plot.reactable0222(res.obj, best, output =output, factor = 'nOTUs',name = paste0('SmallOTU.',covariate.type))


## small sample  
# Stool - binary --Balanced
res.obj1 <- clean_data(dir =paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'C2', factor ='nSams', covariate.type =covariate.type, sub.level ='nSam_L1', 
                       sub.level.rename = "sample=50",sub.methods = sub.methods)
best1 <- cal.best(res.df2= res.obj1$res.df2 %>% filter(methods != 'ZicoSeq'), factor ='nSams',diff.otu.mode =c('Abundant','Rare'))
# Vaginal - binary --Balanced
res.obj2 <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'D2', factor ='nSams', covariate.type =covariate.type, sub.level ='nSam_L1', 
                       sub.level.rename = "sample=50",sub.methods = sub.methods)
best2 <- cal.best(res.df2= res.obj2$res.df2%>% filter(methods != 'ZicoSeq'), factor ='nSams',diff.otu.mode =c('Abundant','Rare'))

### Unbalanced
### moderate
# Stool - binary --unBalanced
res.obj3 <- clean_data(dir =paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'C52', factor ='nSams', covariate.type =covariate.type, sub.level ='nSam_L1', 
                       sub.level.rename = "sample=50",sub.methods = sub.methods)
best3 <- cal.best1(res.df2= res.obj3$res.df2 %>% filter(methods != 'ZicoSeq'), factor ='nSams',diff.otu.mode =c('Abundant'))
# Vaginal - binary --unBalanced
res.obj4 <- clean_data(dir =paste0(wd,'/data/SimulationEvaluation/'),
                       output = output ,filename = 'D52', factor ='nSams', covariate.type =covariate.type, sub.level ='nSam_L1', 
                       sub.level.rename = "sample=50",sub.methods = sub.methods)
best4 <- cal.best1(res.df2= res.obj4$res.df2%>% filter(methods != 'ZicoSeq'), factor ='nSams',diff.otu.mode =c('Abundant'))

# arrange
res.obj <- list(df1 = res.obj1$res.df2, df2 = res.obj2$res.df2,df3 = res.obj3$res.df2, df4 = res.obj4$res.df2)

best <- list(best1, best2, best3, best4)

plot.reactable0222(res.obj, best, output =output, factor = 'nSams',name = paste0('SmallSample.',covariate.type))







##############################
output = paste0(wd,'/Result/SimulationEvaluation07192022/Figure8/');if(!(dir.exists(output))){dir.create(output)}
sub.methods1 = c('GMPR+DESeq2','GMPR+edgeR',"ANCOM-BC","LDM","Aldex2(glm)","GMPR+glm",'MaAsLin2',"corncob",'ZicoSeq')
## Fig 4: binary + covariate: non-compositional
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
##Fig 5: Depth Confounding
output = paste0(wd,'/Result/SimulationEvaluation07192022/Figure9/');if(!(dir.exists(output))){dir.create(output)}
res.obj1 <- clean_data(dir = paste0(wd,'/Data/SimulationEvaluation/'),
                       output = output ,filename = 'C7', factor ='depth.conf.factors',
                       covariate.type ='binary', sub.level ='DL1',
                       sub.level.rename = 'Depth confounding +',sub.methods = sub.methods)
plot.reactable8a(res.df2= res.obj1$res.df2, factor ='depth.conf.factors', diff.otu.mode =c('Abundant','Rare'),output =output, name = 'Fig5.DepthCounfounding.Stool+')

output = paste0(wd,'/Result/SimulationEvaluation07192022/FigureS20/');if(!(dir.exists(output))){dir.create(output)}
res.obj1 <- clean_data(dir = paste0(wd,'/Data/SimulationEvaluation/'),
                       output = output ,filename = 'C7', factor ='depth.conf.factors',
                       covariate.type ='binary', sub.level ='DL3',
                       sub.level.rename = 'Depth confounding ++',sub.methods = sub.methods)

plot.reactable8a(res.df2= res.obj1$res.df2, factor ='depth.conf.factors', diff.otu.mode =c('Abundant','Rare'),output =output, name = 'SupplementaryFigureS20.DepthCounfounding.Stool++')

##############################
##Supplementary Figure S23 
output = paste0(wd,'/Result/SimulationEvaluation07192022/SupplementaryFigure9/');if(!(dir.exists(output))){dir.create(output)}
methods <- c('mbzinb','Wilcox.gmpr','Rarefy','Wilcox',
             'eBayW', 'Aldex2', 'Maaslin2',
             'RAIDA', 'DACOMP', 'LDM', 
             'glmquassi2','BBinomial','ZicoSeq',
             'ANCOMBC', 'DESeq2.gmpr', 'edgeR.gmpr','MSeq2.Wrench')
Methods <- c('Omnibus','GMPR+Wilcox','Rarefy+Wilcox','Wilcox',
             'eBay(Wilcox)', 'Aldex2(Wilcox)', 'MaAsLin2',
             'RAIDA', 'DACOMP', 'LDM', 
             'GMPR+glm','corncob','ZicoSeq',
             'ANCOM-BC','GMPR+DESeq2','GMPR+edgeR', 'Wrench+MSeq')
load(paste0(wd,'/Data/SimulationEvaluation/stability3.Rdata'))
sink(paste0(output,'stability.txt'))
cat('correlation of p-values between filter vs nofilter \n')
for(i in 1:length(methods)){
  method = methods[i]
  data = DF[[method]]
  data$data = gsub('loglinearSub_L1|binarynone|nSam_L2|L3|^D3|none|unbalanced','',data$data)
  # each has 18 facets
  unq = unique(data$data)
  data$diff.otu.modes = 'rare'
  data[grep('abundant',data$data),'diff.otu.modes']= 'abundant'
  
  data$signaldensity = 'low'
  data[grep('medium',data$data),'signaldensity']= 'medium'
  data[grep('high',data$data),'signaldensity']= 'high'
  data$data = gsub(method,'',data$data)
  data$data = gsub('low|medium|high|abundant|rare','',data$data)
  data = data %>% dplyr::select(otu.id, nofilter, filter, diff.otu.modes, signaldensity, data) 
  
  # subset one level to plot, for saving figure size
  #data = data %>% filter(diff.otu.modes=='abundant' & signaldensity =='high' & data =='nOTU_L5')
  cor <- cor.test(data$nofilter, data$filter, method = 'spearman')$estimate
  cat(Methods[i],': ',cor,'\n' )
  # graphics.off()
  try({
    p1 = ggplot(data, aes(x = nofilter, y = filter)) +
      geom_point(size = 0.1) +theme_bw() +
      # facet_wrap(data ~ signaldensity, ncol = 5)+
      labs(x = 'No filter', y = 'Filter by 40% prevelance') +
      theme(axis.text.x = element_text(color="black", size =26),
            axis.text.y = element_text(color="black", size = 26),
            axis.title = element_text(color="black", size = 30),
            strip.text = element_text(size = 30),
            strip.background = element_rect(fill="white",color = "black", size = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black",size= 1),
            legend.position = 'none',
            legend.title=element_text(size=30),
            legend.text = element_text(size=30),
            plot.title = element_text(size=22))+ ggtitle(Methods[i])
    ggsave(file = paste0(output,'/',Methods[i],'_stability.pdf'), width = 10, height = 10, dpi = 100)
  })
}
sink()

##############################
##Fig. 6: heatmap
output = paste0(wd,'/Result/SimulationEvaluation07192022/Figure10/');if(!(dir.exists(output))){dir.create(output)}
sub.methods = c("Rarefy+Wilcox","TSS+Wilcox",'GMPR+Wilcox', 'GMPR+DESeq2','GMPR+edgeR','Wrench+MSeq',
                "RAIDA","ANCOM-BC","DACOMP","LDM","Omnibus",'MaAsLin2','ZicoSeq',
                "Aldex2(Wilcox)","GMPR+glm","corncob","eBay(Wilcox)")

# Stool
## no compositional
c2b <- clean_data(dir =paste0(wd,'/data/SimulationEvaluation/'),filename = 'C2', factor = 'nSams', covariate.type =covariate.type, sub.level = c('nSam_L1','nSam_L4'), 
                  sub.level.rename = c('sample=50', 'sample=200'),sub.methods =  sub.methods)

c4b <- clean_data(dir =paste0(wd,'/data/SimulationEvaluation/'),filename = 'C4', factor ='nOTUs', covariate.type =covariate.type, sub.level = c('nOTU_L1','nOTU_L5'), 
                  sub.level.rename = c('OTU=50','OTU=500'),sub.methods = sub.methods)

c5b <- clean_data(dir =paste0(wd,'/data/SimulationEvaluation/'),filename = 'C7', factor ='depth.conf.factors', covariate.type =covariate.type, sub.level = c('DL1','DL3'), 
                  sub.level.rename = c('Depth confounding +','Depth confounding ++'),sub.methods = sub.methods)
## compositional
c7b <- clean_data(dir =paste0(wd,'/data/SimulationEvaluation/'),filename = 'C52', factor ='nSams', covariate.type =covariate.type, sub.level = c('nSam_L1','nSam_L4'), diff.otu.mode = 'abundant',
                  sub.level.rename = c('sample=50', 'sample=200'),sub.methods = sub.methods)
c9b <- clean_data(dir =paste0(wd,'/data/SimulationEvaluation/'),filename = 'C54', factor ='nOTUs', covariate.type =covariate.type, sub.level = c('nOTU_L1','nOTU_L5'), diff.otu.mode = 'abundant',
                  sub.level.rename = c('OTU=50','OTU=500'),sub.methods = sub.methods)

# save(c2b, c4b, c5b, c7b, c9b, file = paste0(wd,'/data/SimulationEvaluation/Stool_cleaned.Rdata'))
## no compositional
d2b <- clean_data(dir =paste0(wd,'/data/SimulationEvaluation/'),filename = 'D2', factor = 'nSams', covariate.type =covariate.type, sub.level = c('nSam_L1','nSam_L4'), 
                  sub.level.rename = c('sample=50', 'sample=200'),sub.methods =  sub.methods)
d4b <- clean_data(dir =paste0(wd,'/data/SimulationEvaluation/'),filename = 'D4', factor ='nOTUs', covariate.type =covariate.type, sub.level = c('nOTU_L1','nOTU_L5'), 
                  sub.level.rename = c('OTU=50','OTU=500'),sub.methods = sub.methods)
## compositional
d7b <- clean_data(dir =paste0(wd,'/data/SimulationEvaluation/'),filename = 'D52', factor ='nSams', covariate.type =covariate.type, sub.level = c('nSam_L1','nSam_L4'), diff.otu.mode = 'abundant',
                  sub.level.rename = c('sample=50', 'sample=200'),sub.methods = sub.methods)
d9b <- clean_data(dir =paste0(wd,'/data/SimulationEvaluation/'),filename = 'D54', factor ='nOTUs', covariate.type =covariate.type, sub.level = c('nOTU_L1','nOTU_L5'), diff.otu.mode = 'abundant',
                  sub.level.rename = c('OTU=50','OTU=500'),sub.methods = sub.methods)

# save(d2b, d4b, d7b,  d9b, file = paste0(wd,'/data/SimulationEvaluation/Vaginal_cleaned.Rdata'))



load(paste0(wd,'/data/SimulationEvaluation/Vaginal_cleaned.Rdata'))
load(paste0(wd,'/data/SimulationEvaluation/Stool_cleaned.Rdata'))
## baseline : stool/vaginal
stool.baseline.fdr <- fdr.class(data = c4b$res.df2, grp = 'nOTUs', level='OTU=500', name = '', setting ='Basic setting:')
stool.baseline.tpr <- tpr.class(data = c4b$tpr.rank, patn = 'OTU\\=500',setting = 'Basic setting:', type = 'stool', name = '')
stool.baseline <- full_join(stool.baseline.fdr, stool.baseline.tpr)
colnames(stool.baseline) <- gsub('\\_','', colnames(stool.baseline))

vaginal.baseline.fdr <- fdr.class(data = d4b$res.df2, grp = 'nOTUs', level='OTU=500', type = 'vaginal',name = '', setting ='Basic setting:')
vaginal.baseline.tpr <- tpr.class(data = d4b$tpr.rank, patn = 'OTU\\=500',setting = 'Basic setting:', type = 'vaginal', name = '')
vaginal.baseline <- full_join(vaginal.baseline.fdr, vaginal.baseline.tpr)
colnames(vaginal.baseline) <- gsub('\\_','', colnames(vaginal.baseline))


## compostional: baseline : stool/vaginal (OTU=500, diff.otu.pcts = 3 levels, diff.otu.modes =)
compositional.stool.baseline.fdr <- fdr.class(data = c9b$res.df2, grp = 'nOTUs', level='OTU=500', name = '', setting ='Basic setting:')
compositional.stool.baseline.tpr <- tpr.class(data = c9b$tpr.rank, patn = 'OTU\\=500',setting = 'Basic setting:', type = 'stool', name = '')
compositional.stool.baseline <- full_join(compositional.stool.baseline.fdr, compositional.stool.baseline.tpr)
colnames(compositional.stool.baseline) <- gsub('\\_','', colnames(compositional.stool.baseline))
colnames(compositional.stool.baseline) <- gsub('Basic setting\\:','Basic setting:compositional:', colnames(compositional.stool.baseline))


compositional.vaginal.baseline.fdr <- fdr.class(data = d9b$res.df2, grp = 'nOTUs', level='OTU=500', type = 'vaginal',name = '', setting ='Basic setting:')
compositional.vaginal.baseline.tpr <- tpr.class(data = d9b$tpr.rank, patn = 'OTU\\=500',setting = 'Basic setting:', type = 'vaginal', name = '')
compositional.vaginal.baseline <- full_join(compositional.vaginal.baseline.fdr, compositional.vaginal.baseline.tpr)
colnames(compositional.vaginal.baseline) <- gsub('\\_','', colnames(compositional.vaginal.baseline))
colnames(compositional.vaginal.baseline) <- gsub('Basic setting\\:','Basic setting:compositional:', colnames(compositional.vaginal.baseline))

# save(stool.baseline, vaginal.baseline,compositional.stool.baseline,compositional.vaginal.baseline, file = paste0(wd,'/data/SimulationEvaluation/baseline.Rdata'))

## Challenge setting
##  small sample size : sample=50
SampleSizeStool = fdr.class(data = c2b$res.df2, grp = 'nSams', level='sample=50', name = 'SmallSample')
SampleSizeVaginal = fdr.class(data = c2b$res.df2, grp = 'nSams', level='sample=50', name = 'SmallSample', type = 'Vaginal')

SampleSizeStool.tpr = tpr.class(data = c2b$tpr.rank, grp = 'nSams', patn='sample\\=50', name = 'SmallSample')
SampleSizeVaginal.tpr = tpr.class(data = c2b$tpr.rank, grp = 'nSams', patn='sample\\=50', name = 'SmallSample', type = 'Vaginal')


##  small OTU: nOTU = 50
OTUStool = fdr.class(data = c4b$res.df2, grp = 'nOTUs', level='OTU=50', name = 'SmallOTU')
OTUVaginal = fdr.class(data = d4b$res.df2, grp = 'nOTUs', level='OTU=50', name = 'SmallOTU', type = 'Vaginal')

OTUStool.tpr = tpr.class(data = c4b$tpr.rank, grp = 'nOTUs', patn='OTU\\=50$', name = 'SmallOTU')
OTUVaginal.tpr =tpr.class(data = d4b$tpr.rank, grp = 'nOTUs', patn='OTU\\=50$', name = 'SmallOTU', type = 'Vaginal')

## Depth confounding
DepthConfoundingStool = fdr.class(data = c5b$res.df2, grp = 'depth.conf.factors', level=c('Depth confounding +','Depth confounding ++'), name = 'Depth confounding')
DepthConfoundingStool.tpr = tpr.class(data = c5b$tpr.rank, grp = 'depth.conf.factors', patn='Depth confounding', name = 'Depth confounding')

## Caution: compostional setting: diff.otu.modes=='abundant'
## compositional - small sample size : sample=50
CompositionalSampleSizeStool = fdr.class(data = c7b$res.df2, grp = 'nSams', level='sample=50', name = 'CompositionalSmallSample',diff.otu.mode= c('Abundant'))
CompositionalSampleSizeVaginal = fdr.class(data = c7b$res.df2, grp = 'nSams', level='sample=50', name = 'CompositionalSmallSample',diff.otu.mode= c('Abundant'),type = 'Vaginal')

CompositionalSampleSizeStool.tpr = tpr.class(data = c7b$tpr.rank, grp = 'nSams', patn='sample\\=50', name = 'CompositionalSmallSample',patn2 = 'Abundant')
CompositionalSampleSizeVaginal.tpr = tpr.class(data = c7b$tpr.rank, grp = 'nSams', patn='sample\\=50', name = 'CompositionalSmallSample',patn2 = 'Abundant', type = 'Vaginal')


## compositional - small OTU: nOTU = 50
CompositionalOTUStool = fdr.class(data = c9b$res.df2, grp = 'nOTUs', level='OTU=50', name = 'CompositionalSmallOTU',diff.otu.mode= c('Abundant'))
CompositionalOTUVaginal = fdr.class(data = d9b$res.df2, grp = 'nOTUs', level='OTU=50', name = 'CompositionalSmallOTU',diff.otu.mode= c('Abundant'), type = 'Vaginal')

CompositionalOTUStool.tpr = tpr.class(data = c9b$tpr.rank, grp = 'nOTUs', patn='OTU\\=50$', name = 'CompositionalSmallOTU',patn2 = 'Abundant')
CompositionalOTUVaginal.tpr =tpr.class(data = d9b$tpr.rank, grp = 'nOTUs', patn='OTU\\=50$', name = 'CompositionalSmallOTU',patn2 = 'Abundant', type = 'Vaginal')


stool.fdr <- full_join(OTUStool, SampleSizeStool)  %>% full_join(DepthConfoundingStool) %>% 
  full_join(CompositionalSampleSizeStool) %>% full_join(CompositionalOTUStool)

vaginal.fdr  <- full_join(OTUVaginal, SampleSizeVaginal)  %>% 
  full_join(CompositionalSampleSizeVaginal) %>% full_join(CompositionalOTUVaginal)

stool.tpr <- full_join(OTUStool.tpr, SampleSizeStool.tpr)  %>% full_join(DepthConfoundingStool.tpr) %>% 
  full_join(CompositionalSampleSizeStool.tpr) %>% full_join(CompositionalOTUStool.tpr)

vaginal.tpr  <- full_join(OTUVaginal.tpr, SampleSizeVaginal.tpr)  %>% 
  full_join(CompositionalSampleSizeVaginal.tpr) %>% full_join(CompositionalOTUVaginal.tpr)

baseline <- full_join(full_join(stool.baseline, vaginal.baseline), full_join(compositional.stool.baseline, compositional.vaginal.baseline))

stool <- full_join(stool.baseline,stool.fdr) %>% full_join(stool.tpr) %>% full_join(compositional.stool.baseline)
vaginal <- full_join(vaginal.baseline,vaginal.fdr) %>% full_join(vaginal.tpr)%>% full_join(compositional.vaginal.baseline)



ord1 <- colnames(stool)[grep('ompositional',colnames(stool))]
ord2 <- colnames(stool)[-grep('ompositional',colnames(stool))]
stool <- stool %>% dplyr::select(c('methods',ord1[grep('medianFDR',ord1)],ord1[grep('MaxMedianFDR',ord1)],ord1[grep('medianTPR',ord1)],
                                   ord2[grep('medianFDR',ord2)],ord2[grep('MaxMedianFDR',ord2)],ord2[grep('medianTPR',ord2)]))


ord1 <- colnames(vaginal)[grep('ompositional',colnames(vaginal))]
ord2 <- colnames(vaginal)[-grep('ompositional',colnames(vaginal))]
vaginal <- vaginal%>% dplyr::select(c('methods',ord1[grep('medianFDR',ord1)],ord1[grep('MaxMedianFDR',ord1)],ord1[grep('medianTPR',ord1)],
                                      ord2[grep('medianFDR',ord2)],ord2[grep('MaxMedianFDR',ord2)],ord2[grep('medianTPR',ord2)]))


# save(vaginal, stool, file = paste0(wd,'/data/SimulationEvaluation/baselineChanllenge.Rdata'))

## time
# load(paste0(wd,'/data/SimulationEvaluation/Stool_cleaned.Rdata'))
# load(paste0(wd,'/data/SimulationEvaluation/Vaginal_cleaned.Rdata'))
# 
# time = c4b$res.df2 %>% filter(measures =='time' & nOTUs=='OTU=500') %>% group_by(methods) %>% summarize(Speed = mean(value, na.rm = TRUE)) %>% column_to_rownames('methods')
# time2 = d4b$res.df2 %>% filter(measures =='time' & nOTUs=='OTU=500') %>% group_by(methods) %>% summarize(Speed = mean(value, na.rm = TRUE))
# time3 = c9b$res.df2 %>% filter(measures =='time' & nOTUs=='OTU=500') %>% group_by(methods) %>% summarize(Speed = mean(value, na.rm = TRUE)) #%>% column_to_rownames('methods')
# time4 = d9b$res.df2 %>% filter(measures =='time' & nOTUs=='OTU=500') %>% group_by(methods) %>% summarize(Speed = mean(value, na.rm = TRUE))
# colnames(time1)[2] <- 'a';colnames(time2)[2] <- 'b';colnames(time3)[2] <- 'c';colnames(time4)[2] <- 'd'
# time <- full_join(time1, time2) %>% full_join(time3)  %>% full_join(time4)  %>% column_to_rownames('methods')
# time <- rowMeans(time) %>% as.data.frame()
# colnames(time) <- 'Speed'
time <- read.csv(paste0(wd,'/Data/SimulationEvaluation/time.csv')) %>% dplyr::select(methods, value) %>% column_to_rownames('methods')
time.level <- apply(time, 2, function(x) ifelse(x < 60, 'Good', ifelse(x>600,'Poor','Intermediate'))) %>% as.data.frame()
colnames(time.level) = 'Speed'
time.level 
dim(time.level)



## Stability
# load(paste0(wd,'/data/SimulationEvaluation/correlationR3.Rdata'))
corr.sum = corr.sum #%>% rownames_to_column('methods')
corr.sum = aggregate(Stability~methods, corr.sum, function(x) mean(x))
corr.sum  = apply(corr.sum %>% column_to_rownames('methods'), 2, function(x) ifelse(x >0.9, 'Good', ifelse(x<0.7,'Poor','Intermediate'))) %>% as.data.frame() 

## bias DAs 
load(paste0(wd,'/Data/SimulationEvaluation/old/StoolbinaryBiasDAs.Rdata'))
Stool.bias = aggregate(ct ~ methods+depth.conf.factors, data = kb, function(x) mean(x))
Stool.bias.binary.none.median = Stool.bias[Stool.bias$depth.conf.factors=='None',] 
Stool.bias.binary.none.median$score = 'Intermediate'
Stool.bias.binary.none.median[Stool.bias.binary.none.median$ct<=0.05,'score'] = 'Good'
Stool.bias.binary.none.median[Stool.bias.binary.none.median$ct>0.1,'score'] = 'Poor'
Stool.bias.binary.none.median =Stool.bias.binary.none.median[,c('methods','score')]
colnames(Stool.bias.binary.none.median)[2] = 'NULLsetting:medianFDR.stool'
dim(Stool.bias.binary.none.median)

Stool.bias.max = aggregate(ct ~ methods+depth.conf.factors, data = kb, function(x) max(x))
Stool.bias.binary.none.max = Stool.bias[Stool.bias.max$depth.conf.factors=='None',] 
Stool.bias.binary.none.max$score = 'Intermediate'
Stool.bias.binary.none.max[Stool.bias.binary.none.max$ct<=0.05,'score'] = 'Good'
Stool.bias.binary.none.max[Stool.bias.binary.none.max$ct>0.1,'score'] = 'Poor'
Stool.bias.binary.none.max =Stool.bias.binary.none.max[,c('methods','score')]
colnames(Stool.bias.binary.none.max)[2] = 'NULLsetting:maxFDR.stool'
dim(Stool.bias.binary.none.max)

Stool.bias.binary.none = Stool.bias.binary.none.median %>% full_join(Stool.bias.binary.none.max) %>% filter(methods %in% sub.methods)


## *binaryBiasDAs.Rdata are generated by plot_NULL0.R script.
load(paste0(wd,'/data/SimulationEvaluation/VaginalbinaryBiasDAs.Rdata'))
Vaginal.bias = aggregate(ct ~ methods+depth.conf.factors, data = kb, function(x) mean(x))
Vaginal.bias.binary.none.median = Vaginal.bias[Vaginal.bias$depth.conf.factors=='None',] 
Vaginal.bias.binary.none.median$score = 'Intermediate'
Vaginal.bias.binary.none.median[Vaginal.bias.binary.none.median$ct<=0.05,'score'] = 'Good'
Vaginal.bias.binary.none.median[Vaginal.bias.binary.none.median$ct>0.1,'score'] = 'Poor'
Vaginal.bias.binary.none.median =Vaginal.bias.binary.none.median[,c('methods','score')]
colnames(Vaginal.bias.binary.none.median)[2] = 'NULLsetting:medianFDR.vaginal'
dim(Vaginal.bias.binary.none.median)

Vaginal.bias.max = aggregate(ct ~ methods+depth.conf.factors, data = kb, function(x) max(x))
Vaginal.bias.binary.none.max = Vaginal.bias[Vaginal.bias.max$depth.conf.factors=='None',] 
Vaginal.bias.binary.none.max$score = 'Intermediate'
Vaginal.bias.binary.none.max[Vaginal.bias.binary.none.max$ct<=0.05,'score'] = 'Good'
Vaginal.bias.binary.none.max[Vaginal.bias.binary.none.max$ct>0.1,'score'] = 'Poor'
Vaginal.bias.binary.none.max =Vaginal.bias.binary.none.max[,c('methods','score')]
colnames(Vaginal.bias.binary.none.max)[2] = 'NULLsetting:maxFDR.vaginal'
dim(Vaginal.bias.binary.none.max)

Vaginal.bias.binary.none = Vaginal.bias.binary.none.median %>% full_join(Vaginal.bias.binary.none.max)%>% filter(methods %in% sub.methods)


sub.methods1 = c("Rarefy+Wilcox","TSS+Wilcox",'GMPR+Wilcox', 'GMPR+DESeq2','GMPR+edgeR',
                 'Wrench+MSeq',"RAIDA","ANCOM-BC","DACOMP","LDM",
                 "Omnibus",'MaAsLin2',"GMPR+glm","Aldex2(Wilcox)",
                 "corncob","eBay(Wilcox)")
## 03/21/2021 add 
load(paste0(wd,'/Data/SimulationEvaluation/old/baselineChanllenge.Rdata'))
head(vaginal);head(stool)

heat.binary = full_join(vaginal, stool) %>%
  full_join(time.level%>% rownames_to_column('methods')) %>%# full_join(as.data.frame(na.binary.level) %>% rownames_to_column('methods')) %>%
  full_join(Vaginal.bias.binary.none) %>% 
  full_join(Stool.bias.binary.none) %>% 
  full_join(corr.sum %>% rownames_to_column('methods'))%>%
  filter(methods %in% sub.methods) %>%
  column_to_rownames('methods') 
dim(heat.binary )


colnames(heat.binary)


heat.binary$`Complex design` = NA
heat.binary[rownames(heat.binary) %in% c('GMPR+DESeq2','GMPR+edgeR','LDM','ANCOM-BC','ZicoSeq','Aldex2(Wilcox)','Aldex2(t-test)','MaAsLin2',
                                         'GMPR+glm','Omnibus','corncob','DESeq2','edgeR','GMPR+DESeq2','GMPR+edgeR'),'Complex design'] = 'Good'
heat.binary[rownames(heat.binary) %in% c('DACOMP'),'Complex design'] = 'Intermediate'
heat.binary[rownames(heat.binary) %in% c('metagenomeSeq','eBay(t-test)', 'eBay(Wilcox)','RAIDA', 
                                         'Rarefy+Wilcox','GMPR+Wilcox','TSS+Wilcox','Wrench+MSeq'),'Complex design'] = 'Poor'
# Note: Aldex2 can perform covariate adjustment with glm and correlation with cor.test

heat = apply(heat.binary,2, function(x) as.factor(x))
rownames(heat) = rownames(heat.binary)

keep.n1 <- colnames(heat)[-grep('MaxMedian',colnames(heat))]
heat <- heat[,keep.n1, drop =F]
keep.n2 <- colnames(heat)[-grep('maxFDR',colnames(heat))]
heat <- heat[,keep.n2, drop =F]
keep.n3 <- colnames(heat)[-grep('Depth\\ confounding_Challenge:medianTPR.stool',colnames(heat))]
heat <- heat[,keep.n3, drop =F]
dim(heat )
colnames(heat) <- gsub('median','',colnames(heat))



## add 06/02/2021, use ComplexHeatmap package to draw heatmap
library(ComplexHeatmap)
complxheat <- as.data.frame(t(heat))
class(complxheat)
rownames(complxheat)


d <- complxheat %>% rownames_to_column('level')
d$level <- gsub('Vaginal','vaginal', d$level)


idx.fdr <- d$level[grep('FDR',d$level)]
idx1 <- idx.fdr[grep('NULL',idx.fdr)]
idx2 <- idx.fdr[-grep('NULL',idx.fdr)]
idx3 <- d$level[grep('TPR',d$level)]
idx4 <- c("Speed","Complex design","Stability")
length(idx1) + length(idx2) + length(idx3) + length(idx4) ==nrow(complxheat)

module1.idx <- c(idx1, idx2) # FDR + NULL + depthconfounding
module2.idx <- idx3 # TPR
module3.idx <- idx4 # Other


d$class <- NA
d[which(d$level%in% module1.idx[grep('NULL',module1.idx)]),'class'] <- 'A:NULL'

d[which(d$level%in% module1.idx[-grep('NULL',module1.idx)]),'class'] <- 'B:FDR'

d[which(d$level%in% module2.idx),'class'] <- 'C:TPR'
d[which(d$level%in% module3.idx),'class'] <- 'D:Other'


d$type <- 'Vaginal'
d[grep('stool',d$level),'type'] <- 'Stool'
d[grep('Speed|Complex design|Stability',d$level),'type'] <- 'Other'

d$comp <- 'Balanced'
d[grep('ompositional',d$level),'comp'] <- 'Unbalanced'
d[grep('Speed|Complex design|Stability',d$level),'comp'] <- 'Other'


d1 = d
dim(d1)
ord <- c(idx1, 
         ## balanced FDR-vaginal
         "Basic setting:compositional:FDR.vaginal","CompositionalSmallOTU_Challenge:FDR.vaginal","CompositionalSmallSample_Challenge:FDR.vaginal",
         ## balanced FDR-stool
         "Basic setting:compositional:FDR.stool","CompositionalSmallOTU_Challenge:FDR.stool","CompositionalSmallSample_Challenge:FDR.stool",
         ## unbalanced FDR - vaginal
         "Basic setting:FDR.vaginal","SmallOTU_Challenge:FDR.vaginal","SmallSample_Challenge:FDR.vaginal",
         ## unbalanced FDR - stool
         "Basic setting:FDR.stool","SmallOTU_Challenge:FDR.stool","SmallSample_Challenge:FDR.stool",
         
         "Depth confounding_Challenge:FDR.stool",
         
         ##nbalanced TPR-vaginal
         "Basic setting:compositional:TPR.vaginal","CompositionalSmallOTU_Challenge:TPR.vaginal","CompositionalSmallSample_Challenge:TPR.vaginal",
         ## balanced TPR-stool
         "Basic setting:compositional:TPR.stool","CompositionalSmallOTU_Challenge:TPR.stool","CompositionalSmallSample_Challenge:TPR.stool",
         ## unbalanced TPR - vaginal
         "Basic setting:TPR.vaginal","SmallOTU_Challenge:TPR.vaginal","SmallSample_Challenge:TPR.vaginal",
         ## unbalanced TPR - stool
         "Basic setting:TPR.stool","SmallOTU_Challenge:TPR.stool","SmallSample_Challenge:TPR.stool",
         
         module3.idx)
d1 <- d1 %>% column_to_rownames('level')

d1 <- d1[ord,,drop =F];head(d1)
head(d1)
d2 <- as.data.frame(d1[,-c(18:20)]) 
head(d2)

dd = apply(d2, 2, function(x)(as.numeric(as.factor(x)))) %>% as.data.frame()
rownames(dd) = rownames(d2)
ddd = dist(t(dd), method = "euclidean")
dd1 = hclust(ddd, method="ward.D2")
ord.cl = colnames(dd)[dd1$order]

d2 <- d2[,ord.cl]


cols1 = c("Poor"=brewer.pal(8,'Set1')[1],"Good"=brewer.pal(8,'Set1')[3],"Intermediate"=brewer.pal(8,'Set2')[6])
rowsplit = c(rep(as.character(unique(d1$class)[1]),table(d1$class)[1]),rep(as.character(unique(d1$class)[2]),table(d1$class)[2]),
             rep(as.character(unique(d1$class)[3]),table(d1$class)[3]),rep(as.character(unique(d1$class)[4]),table(d1$class)[4]))
rowann = d1[,c('type','comp'),drop =F]
ann_name <- rownames(d2)
ann_name[grep('Depth confounding',ann_name)] <- 'Depth confounding'
ann_name[grep('SmallSample',ann_name)] <- 'Small sample'
ann_name[grep('SmallOTU',ann_name)] <- 'Small number of taxa'
ann_name[grep('Basic setting',ann_name)] <- 'Basic setting'
ann_name[grep('NULLsetting',ann_name)] <- 'Global NULL'

row_ha = rowAnnotation(df = rowann, foo = anno_text(ann_name, location =0, just = "left"),simple_anno_size = unit(0.3, "cm"), 
                       col = list(type = c("Other" = brewer.pal(8,'Dark2')[3], "Stool" = brewer.pal(8,'Dark2')[2],'Vaginal' = brewer.pal(8,'Dark2')[1]),
                                  comp = c("Unbalanced" = brewer.pal(8,'Pastel2')[1], "Balanced" = brewer.pal(8,'Pastel2')[2],"Other" = brewer.pal(8,'Pastel2')[3])))
pdf(paste0(output,'Figure10.pdf'), width = 8, height = 9)
ht = Heatmap(as.matrix(d2), col = cols1,rect_gp = gpar(col = "white", lwd = 0.5), 
             show_row_names =F,
             border = TRUE,row_split = rowsplit,
             right_annotation = row_ha)
draw(ht, padding = unit(c(2, 2, 2, 2), "cm"))
dev.off()



###################
### Supplementary  Figure  S25 : consensus ensemble 
sub.methods = c("pct20","pct40","pct60","pct80")
output = paste0(wd,'/Result/SimulationEvaluation07192022/SupplementaryFigure15/');if(!(dir.exists(output))){dir.create(output)}
load(paste0(wd,'/Data/SimulationEvaluation/old/consensusD4.Rdata'))
plot.reactable8a(res.df2 %>% filter(methods %in% sub.methods), factor ='nOTUs',diff.otu.mode =c('Abundant','Rare'),output =output, name = 'Moderate.Vaginal.balanced.Ensemble')
load(paste0(wd,'/Data/SimulationEvaluation/old/consensusC4.Rdata'))
plot.reactable8a(res.df2 %>% filter(methods %in% sub.methods), factor ='nOTUs',diff.otu.mode =c('Abundant','Rare'),output =output, name = 'Moderate.Stool.balanced.Ensemble')

load(paste0(wd,'/Data/SimulationEvaluation/old/consensusD54.Rdata'))
plot.reactable7b(res.df2=res.df2 %>% filter(diff.otu.modes=='Abundant' &nOTUs=='OTU=500'&methods %in% sub.methods), factor ='nOTUs',diff.otu.mode ='Abundant',output =output, name = 'Moderate.Vaginal.unbalanced.Ensemble')
load(paste0(wd,'/data/SimulationEvaluation/old/consensusC54.Rdata'))
plot.reactable7b(res.df2 = res.df2 %>% filter(diff.otu.modes=='Abundant' &nOTUs=='OTU=500'&methods %in% sub.methods), factor ='nOTUs',diff.otu.mode ='Abundant',output =output, name = 'Moderate.Stool.unbalanced.Ensemble')



################
## plot for cumputation time for n=100, m = 500, here choose NOFILTER setting, which includes n=100, m = 500
load(paste0(wd,'/Data/SimulationEvaluation/old/NOFILTER_res.Rdata'))
sub.methods = c("Rarefy+Wilcox","TSS+Wilcox",'GMPR+Wilcox', 'GMPR+DESeq2','GMPR+edgeR','Wrench+MSeq',
                "RAIDA","ANCOM-BC","DACOMP","LDM","Omnibus",'MaAsLin2','ZicoSeq',
                "Aldex2(Wilcox)","GMPR+glm","corncob","eBay(Wilcox)")
res.obj <- clean_data(dir = paste0(wd,'/data/SimulationEvaluation/old/'),
                       output = output ,filename = 'NOFILTER', factor ='nOTUs', covariate.type ='binary', sub.level ='nOTU_L5', 
                       sub.level.rename = "OTU=500",sub.methods = sub.methods)
time <- res.obj$res2 %>% filter(measures =='time' & diff.otu.pcts=='Low density')
head(time)
dim(time)
mn.time <- data_summary(value ~ methods, data = time)
mn.time <- mn.time[order(-mn.time$median),]
write.csv(mn.time, file = paste0(wd,'/Result/SimulationEvaluation07192022/time.csv'), row.names = F)
time <- within(time, methods <- factor(methods, levels=mn.time$methods))
ggplot(time, aes(x = methods, y = value, color = methods)) + 
  theme_bw() +
  stat_boxplot(geom='errorbar', linetype=1, width=0.5)+
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(aes(colour = methods), width = 0.2, cex = 0.6, alpha = 0.4) +
  scale_color_manual(values = cols) +
  labs(y = 'time(seconds)', x = '', color = "", fill = '') +
  theme(axis.text.x = element_text(color="black", size = 24, angle = 90, vjust = 0.4, hjust = 1),
        axis.text.y = element_text(color="black", size = 24),
        axis.title = element_text(color="black", size = 24),
        strip.text = element_text(size = 20),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        strip.background = element_rect(fill="white",color = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size= 1),
        legend.position = 'none')
ggsave(paste0(wd,'/Result/SimulationEvaluation07192022/Fig.S21.pdf'), width = 7,height = 7)










#########################
## From plotStability.R script
load(paste0(wd,'/Data/SimulationEvaluation/old/stability3.Rdata'))
corr= NULL
for(i in 1:length(EST)){
  cor = as.data.frame(EST[[i]])
  colnames(cor) = 'R'
  cor$method = names(EST)[i]
  corr = rbind(corr, cor)
}
head(corr)
unique(corr$method)
dim(corr)
corr.sum = corr
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
corr.sum$methods = gsub('^MSeq2$','MSeq',corr.sum$methods)
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

sink(paste0(wd,'/Data/SimulationEvaluation/stablity.txt'))
cat('All methods:\n')
summary(corr.sum$Stability)
cat('\n')
corr.sum.mean = data_summary(corr.sum, Stability ~ methods)
corr.sum.mean = corr.sum.mean[order(-corr.sum.mean$median),]
corr.sum.mean
cat('\n')
dacomp <- corr.sum[corr.sum$methods%in%c('DACOMP'),]
cat('DACOMP: \n')
summary(dacomp$Stability)
raida <- corr.sum[corr.sum$methods%in%c('RAIDA'),]
cat('RAIDA: \n')
summary(raida$Stability)
sink()

# write.csv(corr.sum.mean, file = paste0(wd,'/Data/SimulationEvaluation/stablity_summary.csv'),row.names = F)

corr.sum <- corr.sum %>% filter(methods %in% sub.methods)
corr.sum.mean = data_summary(corr.sum, Stability ~ methods) 
corr.sum.mean = corr.sum.mean[order(-corr.sum.mean$median),]

corr.sum <- within(corr.sum, methods <- factor(methods, levels=corr.sum.mean$methods))
range <- inner_join(aggregate(Stability ~ methods, corr.sum, function(x) round(min(x),2)) %>% dplyr::rename(min = Stability),
                    aggregate(Stability ~ methods, corr.sum, function(x) round(max(x),2)) %>% dplyr::rename(max = Stability)) %>% 
  inner_join(aggregate(Stability ~ methods, corr.sum, function(x) round(mean(x),2)) %>% dplyr::rename(mean = Stability))
range[order(range$mean),]


ggplot(corr.sum, aes(x = methods, y = Stability, color = methods)) + 
  theme_bw() +
  stat_boxplot(geom='errorbar', linetype=1, width=0.5)+
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(aes(colour = methods), width = 0.2, cex = 0.6, alpha = 0.4) +
  scale_color_manual(values = cols) +
  labs(y = 'Stability', x = '', color = "", fill = '') +
  theme(axis.text.x = element_text(color="black", size = 24, angle = 90, vjust = 0.4, hjust = 1),
        axis.text.y = element_text(color="black", size = 24),
        axis.title = element_text(color="black", size = 24),
        # axis.ticks = element_blank(),
        strip.text = element_text(size = 24),
        plot.margin = margin(1,1,1,1,'cm'),
        strip.background = element_rect(fill="white",color = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size= 1),
        legend.position = 'none')
ggsave(paste0(wd,'/Result/SimulationEvaluation07192022/SupplementaryFigure9.pdf'), width = 7,height = 7)






