cols <- c('ZicoSeq'=brewer.pal(9,'Paired')[5], 'Omnibus' = brewer.pal(9,'PuBuGn')[7], 'mbzinb' = brewer.pal(9,'PuBuGn')[8], 
          'GMPR+glm' = brewer.pal(9,'PuRd')[3],'LinDA'=brewer.pal(9,'Set1')[8],'MaAsLin2'=brewer.pal(9,'Set1')[7],
          'RAIDA'=brewer.pal(11,'BrBG')[7],'RioNorm2' = brewer.pal(11,'BrBG')[8],
          'eBay(Wilcox)' = brewer.pal(8,'RdPu')[4],'eBay(t-test)' = brewer.pal(8,'RdPu')[6], 'ANCOM-BC' = brewer.pal(8,'Set2')[7],
          'LDM'=brewer.pal(9,'Set1')[1], 'DACOMP'=brewer.pal(8,'Dark2')[6],'Aldex2(Wilcox)'=brewer.pal(8,'YlGn')[3],'Aldex2(t-test)'=brewer.pal(8,'YlGn')[2],
          'corncob'=brewer.pal(11,'BrBG')[2],'DESeq2'=brewer.pal(9,'Greys')[7], 'Wrench+DESeq2'= brewer.pal(9,'Greys')[3], 'GMPR+DESeq2'=brewer.pal(9,'Greys')[5],
          'TSS+Wilcox'=brewer.pal(9,'Purples')[7],'Rarefy+Wilcox'= brewer.pal(9,'Set1')[5], 'GMPR+Wilcox'= brewer.pal(9,'Purples')[5], 'Wrench+Wilcox'= brewer.pal(9,'Purples')[3], 
          'TSS+Spearman'=brewer.pal(9,'Purples')[7],'Rarefy+Spearman'=brewer.pal(9,'Set1')[5], 'GMPR+Spearman'= brewer.pal(9,'Purples')[5], 
          'TSS+t-test'=brewer.pal(9,'GnBu')[7],'Rarefy+t-test'= brewer.pal(9,'GnBu')[5], 'GMPR+t-test'= brewer.pal(9,'GnBu')[5], 'Wrench+t-test'= brewer.pal(9,'GnBu')[3], 
          'edgeR'=brewer.pal(9,'Greens')[7], 'Wrench+edgeR'=brewer.pal(9,'Greens')[3], 'GMPR+edgeR'=brewer.pal(9,'Greens')[5], 
          'MSeq'=brewer.pal(9,'Blues')[7], 'Wrench+MSeq'=brewer.pal(9,'Blues')[5])

data_summary <- function(data, formula){
  ## formula: value ~ variables to be aggregated
  error.info <- aggregate(as.formula(formula), data, function(x) mean(is.na(x)))
  m <- aggregate(as.formula(formula), data, function(x) mean(x[!is.na(x)]))
  med <- aggregate(as.formula(formula), data, function(x) median(x[!is.na(x)]))
  se <- aggregate(as.formula(formula), data, function(x) {
    ind <- !is.na(x)
    sd(x[ind]) / sqrt(length(x[ind]))})
  sd <- aggregate(as.formula(formula), data, function(x) {
    ind <- !is.na(x)
    sd(x[ind])})
  ymin <- m[, ncol(m)] - 1.96 * se[, ncol(m)]
  ymin <- ifelse(ymin > 0, ymin, 0)
  ymax <- m[, ncol(m)] + 1.96 * se[, ncol(m)]
  sum <- cbind(m, median = med[, ncol(med)],SD = sd[, ncol(sd)], ymax=ymax, ymin=ymin,SE=se[, ncol(se)], ErrRate=error.info[, ncol(m)]) 
  return(sum)
}


clean_null <- function(name, type.name, sub.methods){
  # load(paste0('/Users/m216453/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/SimulationEvaluation/',name,'_res.Rdata'))
  # res0 <- melt(res)
  # colnames(res0) <- c('depth.mus', 'diff.otu.modes', 'model', 'nSub', 'nOTUs',
  #                     'covariate.types', 'depth.conf.factors', 'diff.otu.directs', 'confounder.types', 'nSams',
  #                     'covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters')
  # save(res0,file = paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/',name,'_melt.Rdata'))
  load(paste0('~/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/SimulationEvaluation/',name,'_melt.Rdata'))
  res1 = res0 %>% filter(nOTUs %in% c('nOTU_L1','nOTU_L5') & nSams %in% c('nSam_L1','nSam_L4') & measures =='FP') %>% 
    dplyr::select(nOTUs,nSams, covariate.types, depth.conf.factors, methods,value) %>% na.omit() %>%droplevels() 
  
  res1$depth.conf.factors = gsub('DL3','Depth confounding',res1$depth.conf.factors)
  res1$depth.conf.factors = gsub('none','None',res1$depth.conf.factors)
  
  levels(res1$nSams) <- c('sample=50','sample=200')
  levels(res1$nOTUs) <- c('OTU=50','OTU=500')
  levels(res1$depth.conf.factors) <- c('None','Depth confounding')
  
  # na = res1[is.na(res1$value),]
  res1$methods = as.character(res1$methods)
  res1$methods[res1$methods =='Wilcox.gmpr' & res1$covariate.types =='binary'] = 'GMPR+Wilcox'
  res1$methods[res1$methods =='Wilcox.Wrench' & res1$covariate.types =='binary'] = 'Wrench+Wilcox'
  res1$methods[res1$methods =='Wilcox' & res1$covariate.types =='binary'] = 'TSS+Wilcox'
  res1$methods[res1$methods =='Rarefy' & res1$covariate.types =='binary'] = 'Rarefy+Wilcox'
  res1$methods[res1$methods =='Wilcox.gmpr' & res1$covariate.types =='continuous'] = 'GMPR+Spearman'
  res1$methods[res1$methods =='Wilcox' & res1$covariate.types =='continuous'] = 'TSS+Spearman'
  res1$methods[res1$methods =='Rarefy' & res1$covariate.types =='continuous'] = 'Rarefy+Spearman'
  res1$methods[res1$methods =='Aldex2' & res1$covariate.types =='continuous'] = 'Aldex2(glm)'
  res1$methods[res1$methods =='Aldex2' & res1$covariate.types =='binary'] = 'Aldex2(Wilcox)'
  res1$methods = gsub('ANCOMBC','ANCOM-BC',res1$methods)
  res1$methods = gsub('^glmquassi2$','Wrench+glm',res1$methods)
  res1$methods = gsub('^glmquassi$','GMPR+glm',res1$methods)
  res1$methods = gsub('DESeq2.gmpr','GMPR+DESeq2',res1$methods)
  res1$methods = gsub('DESeq2.Wrench','Wrench+DESeq2',res1$methods)
  res1$methods = gsub('edgeR.gmpr','GMPR+edgeR',res1$methods)
  res1$methods = gsub('edgeR.Wrench','Wrench+edgeR',res1$methods)
  res1$methods = gsub('MSeq2.Wrench','Wrench+MSeq',res1$methods)
  res1$methods = gsub('^MSeq2$','MSeq',res1$methods)
  res1$methods = gsub('eBayW','eBay(Wilcox)',res1$methods)
  res1$methods = gsub('eBayt','eBay(t-test)',res1$methods)
  res1$methods = gsub('BBinomial','corncob',res1$methods)
  res1$methods = gsub('Aldex2we','Aldex2(t-test)',res1$methods)
  res1$methods = gsub('^linda$','LinDA',res1$methods)
  res1$methods <- gsub('^mbzinb$','Omnibus', res1$methods)
  res1$methods <- gsub('Maaslin2','MaAsLin2', res1$methods)
  
  res1 = res1 %>% filter(methods %in% sub.methods) %>% dplyr::select(covariate.types,depth.conf.factors,nOTUs, nSams, methods, value) 
  
  return(res= res1)
}

plot_FDR_kable <- function(data, covariate.type,kable.methods, type.name, 
                           output = '/Users/m216453/Documents/Mayo_Research/SemiSimulation/Submit/DAA/result/SimulationEvaluation/'){
  ## For FDIs
  df = data %>% filter(covariate.types==covariate.type) %>% dplyr::select(depth.conf.factors,nOTUs, nSams, methods, value) 
  df[df$value>0,'ct'] = 1
  df[df$value==0,'ct'] = 0
  df = df %>% dplyr::select(-value) 
  formula = paste0('ct~ depth.conf.factors + nOTUs + nSams + methods')
  kb <- data_summary(df, as.formula(formula))
  kb <- kb[,colnames(kb)[c(1:5,9)]]
  colnames(kb) <- gsub('ymin','Low.CI',colnames(kb))
  head(kb)
  
  kb[kb$ct > 0.2,'shape'] = 'x'
  kb[kb$ct <= 0.2 & kb$ct>0.1,'shape'] = '*'
  kb[kb$ct <= 0.1 & kb$ct>0.05,'shape'] = '**'
  # kb[kb$Low.CI <= 0.05,'shape'] = '***'
  kb[kb$ct <= 0.05,'shape'] = '***'
  
  x1 = kb %>% dplyr::select(depth.conf.factors,nOTUs, nSams, methods, shape) %>%filter(methods %in% kable.methods)
  x1$nOTUs = factor(x1$nOTUs, levels = c('OTU=50','OTU=500'))
  x1$nSams = factor(x1$nSams, levels = c('sample=50','sample=200'))
  x1$depth.conf.factors = factor(x1$depth.conf.factors, levels = c('Stool','Vaginal'))
  x2 = x1 %>% unite('grp',c('depth.conf.factors','nOTUs','nSams'))
  x3 = x2 %>% spread(c('grp'), shape)
  head(x3)
  x4 = x3[,c('methods',
             "Stool_OTU=50_sample=50","Stool_OTU=50_sample=200",
             "Stool_OTU=500_sample=50","Stool_OTU=500_sample=200",
             "Vaginal_OTU=50_sample=50","Vaginal_OTU=50_sample=200",
             "Vaginal_OTU=500_sample=50","Vaginal_OTU=500_sample=200")]
  colnames(x4)[1] = 'Sample size'
  colnames(x4) = gsub('.*sample=','',colnames(x4))
  head(x4)
  x4 = x4 %>% column_to_rownames('Sample size')
  
  x41 = x4[,c(1:4)] 
  x41$score1 = apply(x41, 1, function(x) str_count(x, "\\*") %>% sum())
  x42 = x4[,c(5:8)]
  x42$score2 = apply(x42, 1, function(x) str_count(x, "\\*") %>% sum())
  x4 = merge(x41, x42, by = 0)
  colnames(x4) = gsub('\\..*','',colnames(x4))
  colnames(x4)[1] = 'Sample size'
  # x4$Overall = x41$score + x42$score 
  x4 = x4[order(x4$score1, decreasing = T),]
  colnames(x4)[grep('score',colnames(x4))] = c(" "," ")
  rownames(x4) = NULL
  
  count.ct = function(data = x4,j = 2){
    ct = NULL
    for(i in 1:nrow(data)){
      ct=c(ct,str_count(x4[i,j], "\\*") %>% sum())
    }
    ct = gsub('^3$',brewer.pal(11,'PiYG')[9],ct)
    ct = gsub('^2$',brewer.pal(8,'Set2')[6],ct)
    ct = gsub('^1$',brewer.pal(8,'RdGy')[2],ct)
    ct = gsub('^0$',brewer.pal(12,'Set3')[9],ct)
    return(ct)
  } 
  
  count.score = function(data = x4,j = 11, color1 = brewer.pal(8,'Paired')[1], color2 = 'black', value = 12){
    ct = x4[,j]
    ct[ct==value] = color1
    ct[-grep('^\\#',ct)] =color2
    return(ct)
  }
  
  ## order the methods, while makes ZicoSeq top 1
  ord <- x4[,c(1,6,11), drop =F]
  ord$ord <- ord[,2] + ord[,3]
  ord <- ord[order(-ord$ord),]
  # if(ord[1,1] !='ZicoSeq' & ord[ord$`Sample size`=='ZicoSeq',4] == max(ord$ord)){
  #   ord1 <- ord$`Sample size`[!(ord$`Sample size` %in% 'ZicoSeq')]
  #   ord <- c('ZicoSeq',ord1)
  # }else{
  ord <- ord$`Sample size`
  # }
  
  x4 <- ((x4 %>% column_to_rownames('Sample size'))[ord,,drop = F]) %>% rownames_to_column('Sample size')
  colnames(x4) <- gsub('\\..*','',colnames(x4))
  
  
  
  kbl(x4, align = 'c') %>%
    kable_classic_2(full_width = F) %>% 
    add_header_above(c("Taxa number" = 1, "50" =2,"500" =2,"Score" = 1,"50" =2,"500" =2,"Score" = 1), bold = T) %>%
    add_header_above(c(" " = 1, "Stool" =4," " = 1, "Vaginal" = 4," " = 1), bold = T) %>%
    row_spec(0,bold=TRUE)  %>%
    column_spec(2, background= count.ct(data = x4,2),bold = T) %>%
    column_spec(3, background= count.ct(data = x4,3),bold = T) %>%
    column_spec(4, background= count.ct(data = x4,4),bold = T) %>%
    column_spec(5, background= count.ct(data = x4,5),bold = T) %>% 
    column_spec(7, background= count.ct(data = x4,7),bold = T) %>% 
    column_spec(8, background= count.ct(data = x4,8),bold = T) %>% 
    column_spec(9, background= count.ct(data = x4,9),bold = T) %>% 
    column_spec(10, background= count.ct(data = x4,10),bold = T) %>% 
    column_spec(11, color = count.score(data = x4, j = 11,color1 = 'white', color2 = 'black'),
                background = count.score(data = x4, j = 11,color1 = brewer.pal(8,'Dark2')[1], color2 = 'white')) %>%
    column_spec(6, color = count.score(data = x4, j = 6,color1 = 'white', color2 = 'black'),
                background = count.score(data = x4, j = 6,color1 = brewer.pal(8,'Dark2')[1], color2 = 'white')) %>%
    save_kable(file = paste0(output,'NULL/',covariate.type,'_FDI_Null.pdf'))
  
}



clean_data <- function(dir = '~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/', 
                       output ='~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/',
                       filename = 'C1', factor ='covariate.eff.means', covariate.type ='binary',sub.level = c('L2','L4'),
                       sub.level.rename = c('Weak effect','Strong effect'), diff.otu.mode = c('abundant','rare'),
                       sub.methods = c("Rarefy+Wilcox","TSS+Wilcox",'Wrench+DESeq2','Wrench+edgeR','Wrench+MSeq',
                                       "RAIDA","ANCOM-BC","DACOMP","LDM","Omnibus","Aldex2(Wilcox)","Wrench+glm","corncob","eBay(Wilcox)")
){
  load(paste0(dir,filename,'_res.Rdata'))
  res0 <- reshape2::melt(res)
  colnames(res0) <- c('depth.mus', 'diff.otu.modes', 'model', 'nSub', 'nOTUs', 'covariate.types', 'depth.conf.factors', 'diff.otu.directs', 'confounder.types', 'nSams', 'covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters')
  res0$methods = as.character(res0$methods)
  res0$methods[res0$methods =='Aldex2' & res0$confounder.types !='none'] = 'Aldex2glm'
  
  if(length(grep('smallsample',filename))>0){
    res1 <- res0 %>% dplyr::select(c('diff.otu.modes','covariate.types','covariate.eff.means', factor, 'diff.otu.pcts','measures', 'methods', 'value','iters')) %>% 
      filter(diff.otu.modes %in% diff.otu.mode & !!as.name(factor) %in% sub.level & covariate.eff.means %in% 'L3') %>% 
      filter(methods != 'Omnibus')  %>% droplevels() 
    
  }else{
    res1 <- res0 %>% dplyr::select(c('diff.otu.modes','covariate.types',factor, 'diff.otu.pcts','measures', 'methods', 'value','iters')) %>% 
      filter(diff.otu.modes %in% diff.otu.mode & !!as.name(factor) %in% sub.level) %>% 
      filter(methods != 'Omnibus')  %>% droplevels() 
    
  }
  levels(res1[,factor]) <- sub.level.rename 
  levels(res1$diff.otu.modes) = str_to_title(levels(res1$diff.otu.modes))
  levels(res1$diff.otu.pcts) = paste0(str_to_title(unique(res1$diff.otu.pcts)),' density')
  formula = paste0('value ~ diff.otu.modes+',factor,'+ diff.otu.pcts+ methods + measures')
  
  ## Rename methods
  res1$methods = as.character(res1$methods)
  res1$methods[res1$methods =='Wilcox.gmpr' & res1$covariate.types =='binary'] = 'GMPR+Wilcox'
  res1$methods[res1$methods =='Wilcox.Wrench' & res1$covariate.types =='binary'] = 'Wrench+Wilcox'
  res1$methods[res1$methods =='Wilcox' & res1$covariate.types =='binary'] = 'TSS+Wilcox'
  res1$methods[res1$methods =='Rarefy' & res1$covariate.types =='binary'] = 'Rarefy+Wilcox'
  res1$methods[res1$methods =='Wilcox.gmpr' & res1$covariate.types =='continuous'] = 'GMPR+Spearman'
  res1$methods[res1$methods =='Wilcox' & res1$covariate.types =='continuous'] = 'TSS+Spearman'
  res1$methods[res1$methods =='Rarefy' & res1$covariate.types =='continuous'] = 'Rarefy+Spearman'
  res1$methods[res1$methods =='Rarefyttest' & res1$covariate.types =='binary'] = 'Rarefy+t-test'
  res1$methods[res1$methods =='ttest' & res1$covariate.types =='binary'] = 'TSS+t-test'
  res1$methods[res1$methods =='ttest.gmpr' & res1$covariate.types =='binary'] = 'GMPR+t-test'
  res1$methods[res1$methods =='ttest.Wrench' & res1$covariate.types =='binary'] = 'Wrench+t-test'
  res1$methods[res1$methods =='Aldex2' & res1$covariate.types =='continuous'] = 'Aldex2(glm)'
  res1$methods[res1$methods =='Aldex2' & res1$covariate.types =='binary'] = 'Aldex2(Wilcox)'
  res1$methods[res1$methods =='Aldex2glm'] = 'Aldex2(glm)'
  res1$methods = gsub('ANCOMBC','ANCOM-BC',res1$methods)
  res1$methods = gsub('glmquassi2','Wrench+glm',res1$methods)
  res1$methods = gsub('glmquassi','GMPR+glm',res1$methods)
  res1$methods = gsub('DESeq2.gmpr','GMPR+DESeq2',res1$methods)
  res1$methods = gsub('DESeq2.Wrench','Wrench+DESeq2',res1$methods)
  res1$methods = gsub('edgeR.gmpr','GMPR+edgeR',res1$methods)
  res1$methods = gsub('edgeR.Wrench','Wrench+edgeR',res1$methods)
  res1$methods = gsub('MSeq2.Wrench','Wrench+MSeq',res1$methods)
  res1$methods = gsub('MSeq2','MSeq',res1$methods)
  res1$methods = gsub('eBayW','eBay(Wilcox)',res1$methods)
  res1$methods = gsub('eBayt','eBay(t-test)',res1$methods)
  res1$methods = gsub('BBinomial','corncob',res1$methods)
  res1$methods = gsub('Aldex2we','Aldex2(t-test)',res1$methods)
  res1$methods = gsub('^linda$','LinDA',res1$methods)
  res1$methods = gsub('^mbzinb$','Omnibus', res1$methods)
  res1$methods = gsub('^Maaslin2$','MaAsLin2', res1$methods)
  
  res2 <- dplyr::filter(res1, methods %in% sub.methods & covariate.types == covariate.type) %>% dplyr::select(-covariate.types)%>% 
    separate(iters, c('iters','m'),'-') %>% dplyr::select(-m) %>%
    tidyr::complete(diff.otu.modes,!!as.name(factor),diff.otu.pcts, methods, measures, iters) %>% na.omit() %>% droplevels()
  # sub1 = res2[res2$measures=='FDR' & res2$diff.otu.modes=='Abundant' &res2[,factor]  =='Weak effect' &res2$diff.otu.pcts =='Low density',]
  res20 = res2
  
  res.df20 <- data_summary(data= res2, formula = formula)
  #res2[res2$measures =='FDR','value'] = (res2[res2$measures =='FDR','value'])#/0.05 # change FDR to FDR inflation
  # summary <- data_summary(data= res2, formula ='value ~ methods + measures')
  res.df2 <- data_summary(data= res2, formula = formula)
  ord = res.df2 %>% filter(!!as.name(factor) ==names(table(res.df2[factor]))[1] & diff.otu.pcts =='Low density' & measures =='TPR'& diff.otu.modes =='Abundant')
  ord = ord[order(ord$value),]$methods
  res.df2$methods = factor(res.df2$methods, levels = ord)
  res2$methods = factor(res2$methods, levels = ord)
  
  letters1 = c(letters,rev(toupper(letters)))
  for(i in 1:length(unique(res.df2$methods))){
    res.df2[(res.df2$methods ==ord[i]),'label'] = letters1[i]
    res2[(res2$methods ==ord[i]),'label'] = letters1[i]
  }
  
  res.df2$legend = paste0(res.df2$label,':',res.df2$methods)
  res2$legend = paste0(res2$label,':',res2$methods)
  
  for(i in 1:length(unique(res.df2$methods))){
    res.df2[is.na(res.df2$value),'label'] = NA
  }
  
  tpr = aggregate(value~., data = res.df2 %>% filter(measures =='TPR')%>% dplyr::select(c(factor,'diff.otu.pcts', 'diff.otu.modes', 'methods', 'value')), function(x)  mean(x[!is.na(x)])) %>% 
    unite('grp',c('diff.otu.modes','diff.otu.pcts',factor)) %>% 
    spread(c('grp'), value) %>% column_to_rownames('methods')
  tpr.rank = apply(tpr, 2, rank)
  
  return(list(res.df2=res.df2, res2=res2, res20 = res20, res.df20=res.df20, tpr.rank=tpr.rank))
}



plot.reactable8a <- function(res.df2, factor = 'depth.conf.factors', diff.otu.mode =c('Abundant','Rare'),name='',
                             output ='~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/Depth/'){
  fdr = aggregate(value~., data= res.df2 %>% filter(measures =='FDR' & diff.otu.modes %in% diff.otu.mode) %>% dplyr::select(as.symbol(factor),diff.otu.pcts, diff.otu.modes, methods,value) , function(x)  mean(x[!is.na(x)])) %>%
    unite('grp',c('diff.otu.modes','diff.otu.pcts',factor)) %>% spread(c('grp'), value)
  fdr[is.na(fdr)] = 100
  colnames(fdr)
  
  names <- c('methods',
             paste0(sort(levels(res.df2$diff.otu.modes))[1],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[factor]])[1]),
             paste0(sort(levels(res.df2$diff.otu.modes))[1],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[factor]])[1]),
             paste0(sort(levels(res.df2$diff.otu.modes))[1],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[factor]])[1]),
             
             paste0(sort(levels(res.df2$diff.otu.modes))[2],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[factor]])[1]),
             paste0(sort(levels(res.df2$diff.otu.modes))[2],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[factor]])[1]),
             paste0(sort(levels(res.df2$diff.otu.modes))[2],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[factor]])[1]))
  
  fdr = fdr[,names] %>% column_to_rownames('methods')
  
  # Lower CI: fdr.est - 1.96 * sqrt(fdr.est * (1 - fdr.est)/ num.iter))
  fdr.CI = res.df2 %>% filter(measures =='FDR' & diff.otu.modes %in% diff.otu.mode) %>% dplyr::select(as.symbol(factor),diff.otu.pcts, diff.otu.modes, methods, ymin) %>% unite('grp',c('diff.otu.modes','diff.otu.pcts',factor)) %>% spread(c('grp'), ymin)
  fdr.CI[is.na(fdr.CI)] = 100
  fdr.CI = fdr.CI[,names] %>% column_to_rownames('methods')
  dim(fdr);dim(fdr.CI);fdr[1:3,1:3];fdr.CI[1:3,1:3]
  test <- function(v1, v2) ifelse(v1<=0.05,"***", ifelse(v2<=0.1, "**", ifelse(v2<=0.2, "*", 'x')))
  label.fdr <- mapply(test,v1 = fdr.CI, v2 = fdr)
  colnames(label.fdr) = colnames(fdr);rownames(label.fdr) = rownames(fdr)
  # Weak effect szie 
  fdr = label.fdr
  fdr = as.data.frame(fdr[,c(1:6)])
  fdr$score = apply(fdr, 1, function(x) str_count(x, "\\*") %>% sum())
  colnames(fdr) = gsub('\\..*','',colnames(fdr))
  fdr = (fdr) %>% rownames_to_column('Signal density')
  fdr$score <- rank(fdr$score)
  
  ## TPR
  tpr = aggregate(value~., data = res.df2 %>% filter(measures =='TPR' & diff.otu.modes %in% diff.otu.mode) %>% dplyr::select(factor, diff.otu.pcts, diff.otu.modes, methods, value), function(x)  mean(x[!is.na(x)])) %>% 
    unite('grp',c('diff.otu.modes','diff.otu.pcts',factor)) %>% 
    spread(c('grp'), value) %>% dplyr::select(names) %>% column_to_rownames('methods')
  tpr.rank = apply(tpr, 2, rank)
  # add score
  tpr[,'FDR score'] = fdr[,8]
  tpr[,'TPR score'] = rank(apply(tpr.rank[,c(1:6)], 1, function(x) round(mean(x), digits = 1)))
  
  # reorder the columns
  tpr = tpr  %>% rownames_to_column('Signal density')
  
  ord <- tpr[,c(1,8:9),drop=F]
  ord$ord <- ord[,2] + ord[,3]
  ord <- ord[order(-ord$ord),]
  
  # if(ord[1,1] !='ZicoSeq' & ord[ord$`Signal density`=='ZicoSeq',4] == max(ord$ord)){
  #   ord1 <- ord$`Signal density`[!(ord$`Signal density` %in% 'ZicoSeq')]
  #   ord <- c('ZicoSeq',ord1)
  # }else{
  #   ord <- ord$`Signal density`
  # }
    # order by sum(TPR + FDR score )
  ord <- ord$`Signal density`

  tpr = (tpr %>% column_to_rownames('Signal density'))[ord,,drop = F] %>% rownames_to_column('Signal density')
  fdr = fdr %>% column_to_rownames('Signal density')
  fdr = fdr[tpr$`Signal density`,] %>% rownames_to_column('Signal density')
  
  for(i in c(2:7)){
    tpr[,i] = format(round(tpr[,i], digits = 2), nsmall = 2)
  }
  
  bar_chart <- function(label, width = "100%", height = "10px", fill = "forestgreen", background = NULL) {
    bar <- div(style = list(background = fill, width = width, height = height))
    chart <- div(style = list(flexGrow = 1, marginLeft = "1px", background = background), bar)
    div(style = list(display = "flex", alignItems = "center"), label, chart)
  }
  
  get_color <- function(fdr){
    if(fdr =='***'){
      col = brewer.pal(8,'Dark2')[5]
    }else if(fdr =='**'){
      col = brewer.pal(8,'Set2')[6]
    }else if(fdr =='*'){
      col = brewer.pal(8,'RdGy')[2]
    }else{
      col = 'grey'
    }
    return(col)
  }
  
  # making table 
  fdr.dt <- fdr[,-c(8),drop = F] # for colddef.list
  tpr[,c(2:7)] = sapply(tpr[,c(2:7)], function(x) round(as.numeric(x), digits = 2), USE.NAMES = F)
  tpr[,c(8)] = sapply(tpr[,c(8)], function(x) round(as.numeric(x), digits = 1), USE.NAMES = F)
  
  
  n1 <- gsub('\\ density.*','', colnames(tpr))
  n2 <- gsub('.*\\_','', n1)
  NM <- n2[grep('Low|Medium|High',n2)]
  
  coldefs.list <- function(numcols){
    coldefs_list = NULL
    color_list = list()
    for (idx in 1:length(numcols)){
      fdr_col = (fdr.dt)[, (idx+1)] # column 8 is fdr score column, which does not match tpr table, need to be deleted; +1: fdr first column is methods
      colors <- sapply(fdr_col, function(x) get_color(x), USE.NAMES = F)
      name = numcols[idx]
      color_list[[name]] = colors
      
      cell.func <- function(value, index, name) {
        width <- value
        bar_chart(format(round(value, digits = 2), nsmall = 2), width = width*80, background = brewer.pal(12,'Set3')[9],
                  fill = color_list[[name]][index])
      }
      
      style.func <- function(value, index, name) {
        color <- color_list[[name]][index]
        list(fontWeight = 800, fontSize = 16,color = color)
      }
      
      coldefs <- list(
        reactable::colDef(style = style.func, cell=cell.func, name = NM[idx], align = 'center')
      )
      
      coldefs_list = c(coldefs_list,coldefs)
    }
    
    # change names
    names(coldefs_list) <- numcols
    
    return (list(coldefs_list, color_list))
  }
  
  numcols <- colnames(tpr)[c(2:7)]
  
  coldefs_list <- coldefs.list(numcols)[[1]]
  color_list <- coldefs.list(numcols)[[2]]
  
  minW = 39
  coldefs_list[[colnames(tpr)[1]]] = colDef(minWidth = 120, align = 'center', cell = function(value, index) {
    name = colnames(tpr)[1]
    name <- tpr[,1][index]
    tagList(
      div(style = list(fontWeight = 600, fontSize = 16), name)
    )
  })
  coldefs_list[[colnames(tpr)[9]]] = colDef(name = "TPR",minWidth = minW, align = 'center', cell = function(value, index) {
    name = colnames(tpr)[9]
    name <- tpr[,9][index]
    tagList(
      div(style = list(fontWeight = 600, fontSize = 16), name)
    )
  })
  coldefs_list[[colnames(tpr)[8]]] = colDef(name = "FDR",minWidth = minW, align = 'center', cell = function(value, index) {
    name = colnames(tpr)[8]
    name <- tpr[,8][index]
    tagList(
      div(style = list(fontWeight = 600, fontSize = 16), name)
    )
  })
  
  tpr[,1] <- sub("Wrench\\+MSeq",'Wrench+MSeq',tpr[,1])
  
  # tpr[,'TPR score'] <- rank(tpr[,'TPR score'])
  html_file <- paste0(output,name,'_NonCompostional.html')
  tb <- reactable(tpr, pagination=FALSE,resizable = FALSE, wrap = FALSE, bordered = F,
                  style = list(fontFamily = "Work Sans, sans-serif", fontSize = "16px"),
                  theme = reactableTheme(
                    headerStyle = list(
                      "&:hover[aria-sort]" = list(background = "#f7f7f8"),
                      "&[aria-sort='ascending'], &[aria-sort='descending']" = list(background = "hsl(0, 0%, 98%)"),
                      borderColor = "#555"
                    )
                  ),
                  columnGroups = list(
                    colGroup(name = '', columns = colnames(tpr)[1]),
                    colGroup(name = gsub('_.*','',colnames(tpr)[2]), columns = colnames(tpr)[2:4]),
                    colGroup(name = gsub('_.*','',colnames(tpr)[5]), columns = colnames(tpr)[5:7]),
                    colGroup(name = "Rank", columns = colnames(tpr)[8:9])
                  ),
                  columns = coldefs_list
  )
  saveWidget(widget = tb, file = html_file, selfcontained = TRUE)
  
  
  img_file <- paste0(output,name,'_NonCompostional.png')
  webshot(url = html_file, file = img_file, vwidth = 1300, vheight = 1000, zoom = 3)
}




plot.reactable7b <- function(res.df2, factor = 'depth.conf.factors', diff.otu.mode =c('Abundant'), name = '', 
                             output ='~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/Depth/'){
  fdr = aggregate(value~., data= res.df2 %>% filter(measures =='FDR' & diff.otu.modes %in% diff.otu.mode) %>% droplevels()%>% dplyr::select(diff.otu.pcts,  methods,value) , function(x)  mean(x[!is.na(x)])) %>%
    unite('grp',c('diff.otu.pcts')) %>% spread(c('grp'), value)
  fdr[is.na(fdr)] = 100
  colnames(fdr)
  
  names <- c('methods',
             paste0(levels(res.df2$diff.otu.pcts)[1]),
             paste0(levels(res.df2$diff.otu.pcts)[2]),
             paste0(levels(res.df2$diff.otu.pcts)[3]))
  
  fdr = fdr[,names] %>% column_to_rownames('methods')
  
  # Lower CI: fdr.est - 1.96 * sqrt(fdr.est * (1 - fdr.est)/ num.iter))
  fdr.CI = res.df2 %>% filter(measures =='FDR' & diff.otu.modes %in% diff.otu.mode) %>% dplyr::select(diff.otu.pcts, methods, ymin) %>% unite('grp',c('diff.otu.pcts')) %>% spread(c('grp'), ymin)
  fdr.CI[is.na(fdr.CI)] = 100
  fdr.CI = fdr.CI[,names] %>% column_to_rownames('methods')
  dim(fdr);dim(fdr.CI);fdr[1:3,1:3];fdr.CI[1:3,1:3]
  test <- function(v1, v2) ifelse(v1<=0.05,"***", ifelse(v2<=0.1, "**", ifelse(v2<=0.2, "*", 'x')))
  label.fdr <- mapply(test,v1 = fdr.CI, v2 = fdr)
  colnames(label.fdr) = colnames(fdr);rownames(label.fdr) = rownames(fdr)
  # Weak effect szie 
  fdr = label.fdr
  fdr = as.data.frame(fdr[,c(1:3)])
  fdr$score = apply(fdr, 1, function(x) str_count(x, "\\*") %>% sum())
  colnames(fdr) = gsub('\\..*','',colnames(fdr))
  fdr = fdr %>% rownames_to_column('Signal density')
  fdr$score <- rank(fdr$score)
  
  ## TPR
  tpr = aggregate(value~., data = res.df2 %>% filter(measures =='TPR' & diff.otu.modes %in% diff.otu.mode)%>% dplyr::select(diff.otu.pcts, methods, value), function(x)  mean(x[!is.na(x)])) %>% 
    unite('grp',c('diff.otu.pcts')) %>% 
    spread(c('grp'), value) %>% dplyr::select(names) %>% column_to_rownames('methods')
  tpr = tpr[fdr[,'Signal density'],]
  tpr.rank = apply(tpr, 2, rank)
  # add score
  tpr[,'FDR score'] = fdr[,5]
  tpr[,'TPR score'] = rank(apply(tpr.rank[,c(1:3)], 1, function(x) round(mean(x), digits = 1)))
  
  # reorder the columns
  tpr = tpr  %>% rownames_to_column('Signal density')
  
  ord <- tpr[,c(1,5:6),drop=F]
  ord$ord <- ord[,2] + ord[,3]
  ord <- ord[order(-ord$ord),]
  # if(ord[1,1] !='ZicoSeq' & ord[ord$`Signal density`=='ZicoSeq',4] == max(ord$ord)){
  #   ord1 <- ord$`Signal density`[!(ord$`Signal density` %in% 'ZicoSeq')]
  #   ord <- c('ZicoSeq',ord1)
  # }else{
  # ord <- ord$`Signal density`
  # }
  
    # order by sum(TPR + FDR score )
  ord <- ord$`Signal density`

  tpr = (tpr %>% column_to_rownames('Signal density'))[ord,,drop = F] %>% rownames_to_column('Signal density')
  fdr = fdr %>% column_to_rownames('Signal density')
  fdr = fdr[tpr$`Signal density`,] %>% rownames_to_column('Signal density')
  
  for(i in c(2:4)){
    tpr[,i] = format(round(tpr[,i], digits = 2), nsmall = 2)
  }
  
  bar_chart <- function(label, width = "100%", height = "10px", fill = "forestgreen", background = NULL) {
    bar <- div(style = list(background = fill, width = width, height = height))
    chart <- div(style = list(flexGrow = 1, marginLeft = "1px", background = background), bar)
    div(style = list(display = "flex", alignItems = "center"), label, chart)
  }
  
  get_color <- function(fdr){
    if(fdr =='***'){
      col = brewer.pal(8,'Dark2')[5]
    }else if(fdr =='**'){
      col = brewer.pal(8,'Set2')[6]
    }else if(fdr =='*'){
      col = brewer.pal(8,'RdGy')[2]
    }else{
      col = 'grey'
    }
    return(col)
  }
  
  # making table 
  fdr.dt <- fdr[,-c(5),drop = F]
  tpr[,c(2:4)] = sapply(tpr[,c(2:4)], function(x) round(as.numeric(x), digits = 2), USE.NAMES = F)
  tpr[,c(5:6)] = sapply(tpr[,c(5:6)], function(x) round(as.numeric(x), digits = 1), USE.NAMES = F)
  
  n1 <- gsub('\\ density.*','', colnames(tpr))
  n2 <- gsub('.*\\_','', n1)
  NM <- n2[grep('Low|Medium|High',n2)]
  
  coldefs.list <- function(numcols){
    coldefs_list = NULL
    color_list = list()
    for (idx in 1:length(numcols)){
      fdr_col = (fdr.dt)[, (idx+1)] # column 8 is fdr score column, which does not match tpr table, need to be deleted; +1: fdr first column is methods
      colors <- sapply(fdr_col, function(x) get_color(x), USE.NAMES = F)
      name = numcols[idx]
      color_list[[name]] = colors
      
      cell.func <- function(value, index, name) {
        width <- value
        bar_chart(format(round(value, digits = 2), nsmall = 2), width = width*80, background = brewer.pal(12,'Set3')[9],
                  fill = color_list[[name]][index])
      }
      
      style.func <- function(value, index, name) {
        color <- color_list[[name]][index]
        list(fontWeight = 800, fontSize = 16,color = color)
      }
      
      coldefs <- list(
        reactable::colDef(style = style.func, cell=cell.func, name = NM[idx], align = 'center')
      )
      
      coldefs_list = c(coldefs_list,coldefs)
    }
    
    # change names
    names(coldefs_list) <- numcols
    
    return (list(coldefs_list, color_list))
  }
  
  numcols <- colnames(tpr)[c(2:4)]
  
  coldefs_list <- coldefs.list(numcols)[[1]]
  color_list <- coldefs.list(numcols)[[2]]
  
  minW = 39
  coldefs_list[[colnames(tpr)[1]]] = colDef(minWidth = 120, align = 'center', cell = function(value, index) {
    name = colnames(tpr)[1]
    name <- tpr[,1][index]
    tagList(
      div(style = list(fontWeight = 600, fontSize = 16), name)
    )
  })
  coldefs_list[[colnames(tpr)[5]]] = colDef(name = "FDR",minWidth = minW, align = 'center', cell = function(value, index) {
    name = colnames(tpr)[5]
    name <- tpr[,5][index]
    tagList(
      div(style = list(fontWeight = 600, fontSize = 16), name)
    )
  })
  coldefs_list[[colnames(tpr)[6]]] = colDef(name = "TPR",minWidth = minW, align = 'center', cell = function(value, index) {
    name = colnames(tpr)[6]
    name <- tpr[,6][index]
    tagList(
      div(style = list(fontWeight = 600, fontSize = 16), name)
    )
  })
  
  tpr[,1] <- sub("Wrench\\+MSeq",'Wrench+MSeq',tpr[,1])
  # tpr[,'TPR score'] <- rank(tpr[,'TPR score'])
  html_file <- paste0(output,name,'_Compostional.html')
  tb <- reactable(tpr, pagination=FALSE,resizable = FALSE, wrap = FALSE, bordered = F,
                  style = list(fontFamily = "Work Sans, sans-serif", fontSize = "16px"),
                  theme = reactableTheme(
                    headerStyle = list(
                      "&:hover[aria-sort]" = list(background = "#f7f7f8"),
                      "&[aria-sort='ascending'], &[aria-sort='descending']" = list(background = "hsl(0, 0%, 98%)"),
                      borderColor = "#555"
                    )
                  ),
                  columnGroups = list(
                    colGroup(name = '', columns = colnames(tpr)[1]),
                    colGroup(name = gsub('.*\\.','', name), columns = colnames(tpr)[2:4]),
                    colGroup(name = "Rank", columns = colnames(tpr)[5:6])),
                  columns = coldefs_list
  )
  saveWidget(widget = tb, file = html_file, selfcontained = TRUE)
  
  
  img_file <- paste0(output,name,'_Compostional.png')
  webshot(url = html_file, file = img_file, vwidth = 800, vheight = 1000, zoom = 3)
  
}




cal.best <- function(res.df2, factor = 'depth.conf.factors', diff.otu.mode =c('Abundant','Rare'),name=''){
  fdr = aggregate(value~., data= res.df2 %>% filter(measures =='FDR' & diff.otu.modes %in% diff.otu.mode) %>% dplyr::select(as.symbol(factor),diff.otu.pcts, diff.otu.modes, methods,value) , function(x)  mean(x[!is.na(x)])) %>%
    unite('grp',c('diff.otu.modes','diff.otu.pcts',factor)) %>% spread(c('grp'), value) 
  fdr[is.na(fdr)] = 100
  colnames(fdr)
  
  names <- c('methods',
             paste0(sort(levels(res.df2$diff.otu.modes))[1],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[factor]])[1]),
             paste0(sort(levels(res.df2$diff.otu.modes))[1],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[factor]])[1]),
             paste0(sort(levels(res.df2$diff.otu.modes))[1],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[factor]])[1]),
             
             paste0(sort(levels(res.df2$diff.otu.modes))[2],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[factor]])[1]),
             paste0(sort(levels(res.df2$diff.otu.modes))[2],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[factor]])[1]),
             paste0(sort(levels(res.df2$diff.otu.modes))[2],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[factor]])[1]))
  
  fdr = fdr[,names] %>% column_to_rownames('methods')
  
  # Lower CI: fdr.est - 1.96 * sqrt(fdr.est * (1 - fdr.est)/ num.iter))
  fdr.CI = res.df2 %>% filter(measures =='FDR' & diff.otu.modes %in% diff.otu.mode) %>% dplyr::select(as.symbol(factor),diff.otu.pcts, diff.otu.modes, methods, ymin) %>% unite('grp',c('diff.otu.modes','diff.otu.pcts',factor)) %>% spread(c('grp'), ymin)
  fdr.CI[is.na(fdr.CI)] = 100
  fdr.CI = fdr.CI[,names] %>% column_to_rownames('methods')
  dim(fdr);dim(fdr.CI);fdr[1:3,1:3];fdr.CI[1:3,1:3]
  test <- function(v1, v2) ifelse(v1<=0.05,"***", ifelse(v2<=0.1, "**", ifelse(v2<=0.2, "*", 'x')))
  label.fdr <- mapply(test,v1 = fdr.CI, v2 = fdr)
  colnames(label.fdr) = colnames(fdr);rownames(label.fdr) = rownames(fdr)
  # Weak effect szie 
  fdr = label.fdr
  fdr = as.data.frame(fdr[,c(1:6)])
  fdr$score = apply(fdr, 1, function(x) str_count(x, "\\*") %>% sum())
  colnames(fdr) = gsub('\\..*','',colnames(fdr))
  fdr = (fdr) %>% rownames_to_column('Signal density')
  fdr$score <- rank(fdr$score)
  
  
  ## TPR
  tpr = aggregate(value~., data = res.df2 %>% filter(measures =='TPR' & diff.otu.modes %in% diff.otu.mode) %>% dplyr::select(factor, diff.otu.pcts, diff.otu.modes, methods, value), function(x)  mean(x[!is.na(x)])) %>% 
    unite('grp',c('diff.otu.modes','diff.otu.pcts',factor)) %>% 
    spread(c('grp'), value) %>% dplyr::select(names) %>% column_to_rownames('methods')
  tpr.rank = apply(tpr, 2, rank)
  # add score
  tpr[,'FDR score'] = fdr[,8]
  tpr[,'TPR score'] = rank(apply(tpr.rank[,c(1:6)], 1, function(x) round(mean(x), digits = 1)))
  
  # reorder the columns
  tpr = tpr  %>% rownames_to_column('Signal density')
  
  ord <- tpr[,c(1,8:9),drop=F]
  # ord$`TPR score` <- rank(ord$`TPR score`)
  
  ord$ord <- ord[,2] + ord[,3]
  ord <- ord[order(-ord$ord),]
  best <- ord[1,1]
  return(best)
}


cal.best1 <- function(res.df2, factor = 'depth.conf.factors', diff.otu.mode =c('Abundant'),name=''){
  fdr = aggregate(value~., data= res.df2 %>% filter(measures =='FDR' & diff.otu.modes %in% diff.otu.mode) %>% dplyr::select(diff.otu.pcts, methods,value) , function(x)  mean(x[!is.na(x)])) %>%
    unite('grp',c('diff.otu.pcts')) %>% spread(c('grp'), value) 
  fdr[is.na(fdr)] = 100
  
  names <- c('methods', "Low density", "Medium density", "High density")
  
  fdr = fdr[,names] %>% column_to_rownames('methods')
  
  # Lower CI: fdr.est - 1.96 * sqrt(fdr.est * (1 - fdr.est)/ num.iter))
  fdr.CI = res.df2 %>% filter(measures =='FDR' & diff.otu.modes %in% diff.otu.mode) %>% dplyr::select(diff.otu.pcts, methods, ymin) %>% unite('grp',c('diff.otu.pcts')) %>% spread(c('grp'), ymin)
  fdr.CI[is.na(fdr.CI)] = 100
  fdr.CI = fdr.CI[,names] %>% column_to_rownames('methods')
  dim(fdr);dim(fdr.CI);fdr[1:3,1:3];fdr.CI[1:3,1:3]
  test <- function(v1, v2) ifelse(v1<=0.05,"***", ifelse(v2<=0.1, "**", ifelse(v2<=0.2, "*", 'x')))
  label.fdr <- mapply(test,v1 = fdr.CI, v2 = fdr)
  colnames(label.fdr) = colnames(fdr);rownames(label.fdr) = rownames(fdr)
  # Weak effect szie 
  fdr = label.fdr
  fdr = as.data.frame(fdr[,c(1:3)])
  fdr$score = apply(fdr, 1, function(x) str_count(x, "\\*") %>% sum())
  colnames(fdr) = gsub('\\..*','',colnames(fdr))
  fdr = (fdr) %>% rownames_to_column('methods')
  fdr$score <- rank(fdr$score)
  
  
  ## TPR
  tpr = aggregate(value~., data = res.df2 %>% filter(measures =='TPR' & diff.otu.modes %in% diff.otu.mode) %>% dplyr::select( diff.otu.pcts,methods, value), function(x)  mean(x[!is.na(x)])) %>% 
    unite('grp',c('diff.otu.pcts')) %>% 
    spread(c('grp'), value) %>% dplyr::select(names) %>% column_to_rownames('methods')
  tpr.rank = apply(tpr, 2, rank)
  # add score
  tpr[,'FDR score'] = fdr[,5]
  tpr[,'TPR score'] = rank(apply(tpr.rank[,c(1:3)], 1, function(x) round(mean(x), digits = 2)))
  
  # reorder the columns
  tpr = tpr  %>% rownames_to_column('Signal density')
  
  ord <- tpr[,c(1,5:6),drop=F]
  # ord$`TPR score` <- rank(ord$`TPR score`)
  ord$ord <- (ord[,2] + ord[,3])
  
  ord <- ord[order(-ord$ord),]
  best <- ord[1,1]
  return(best)
}









plot.reactable0222 <- function(res.obj, best, output, name, target.method = 'ZicoSeq', factor = 'nOTUs'){
  ## ---------------non-compositional---------------
  names <- c('methods','scenario',"Abundant_Low density", "Abundant_Medium density","Abundant_High density" ,"Rare_Low density" ,"Rare_Medium density","Rare_High density")
  fdr = rbind(aggregate(value~., data= res.obj$df1 %>% filter(measures =='FDR'&methods %in% c(target.method,best[1])) %>% dplyr::select(diff.otu.pcts, diff.otu.modes, methods,value) , function(x)  mean(x[!is.na(x)])) %>%
                unite('grp',c('diff.otu.modes','diff.otu.pcts')) %>% spread(c('grp'), value) %>% mutate(scenario ='stool_balanced'),
              aggregate(value~., data= res.obj$df2 %>% filter(measures =='FDR'&methods %in% c(target.method,best[2])) %>% dplyr::select(diff.otu.pcts, diff.otu.modes, methods,value) , function(x)  mean(x[!is.na(x)])) %>%
                unite('grp',c('diff.otu.modes','diff.otu.pcts')) %>% spread(c('grp'), value) %>% mutate(scenario ='vaginal_balanced'))
  idx <- which(fdr$methods ==target.method)
  idx1 <- which(fdr$methods !=target.method)
  idxx <- c(idx[1],idx1[1],idx[2],idx1[2])
  fdr = fdr[idxx,]
  fdr = fdr[,names] %>% unite('methods',c('scenario','methods'), sep = ': ') 
  rownames(fdr) <- NULL
  fdr = fdr %>% column_to_rownames('methods')
  head(fdr)
  
  # Lower CI: fdr.est - 1.96 * sqrt(fdr.est * (1 - fdr.est)/ num.iter))
  fdr.CI = rbind(res.obj$df1 %>% filter(measures =='FDR' &methods %in% c(target.method,best[1])) %>% dplyr::select(diff.otu.pcts, diff.otu.modes, methods, ymin) %>% unite('grp',c('diff.otu.modes','diff.otu.pcts')) %>% spread(c('grp'), ymin) %>% mutate(scenario ='stool_balanced'),
                 res.obj$df2 %>% filter(measures =='FDR' &methods %in% c(target.method,best[2])) %>% dplyr::select(diff.otu.pcts, diff.otu.modes, methods, ymin) %>% unite('grp',c('diff.otu.modes','diff.otu.pcts')) %>% spread(c('grp'), ymin) %>% mutate(scenario ='vaginal_balanced'))
  idx <- which(fdr.CI$methods ==target.method)
  idx1 <- which(fdr.CI$methods !=target.method)
  idxx <- c(idx[1],idx1[1],idx[2],idx1[2])
  fdr.CI = fdr.CI[idxx,]
  
  fdr.CI = fdr.CI[,names] %>% unite('methods',c('scenario','methods'), sep = ': ') 
  rownames(fdr.CI) <- NULL
  fdr.CI = fdr.CI %>% column_to_rownames('methods')
  
  test <- function(v1, v2) ifelse(v1<=0.05,"***", ifelse(v2<=0.1, "**", ifelse(v2<=0.2, "*", 'x')))
  label.fdr <- mapply(test,v1 = fdr.CI, v2 = fdr)
  colnames(label.fdr) = colnames(fdr);rownames(label.fdr) = rownames(fdr)
  
  # Weak effect szie 
  fdr = label.fdr
  # fdr = as.data.frame(fdr[,c(1:6)])
  # fdr$score = apply(fdr, 1, function(x) str_count(x, "\\*") %>% sum())
  # colnames(fdr) = gsub('\\..*','',colnames(fdr))
  fdr = as.data.frame(fdr) %>% rownames_to_column('methods')
  
  
  data1 <- rbind(res.obj$df1 %>% mutate(scenario = 'stool_balanced') %>% dplyr::select(-factor) %>% filter(methods %in% c(target.method,best[1])),
                 res.obj$df2 %>% mutate(scenario = 'vaginal_balanced')%>% dplyr::select(-factor) %>% filter(methods %in% c(target.method,best[2])))
  
  tpr <- data1 %>% filter(measures %in% 'TPR') %>% dplyr::select(scenario, diff.otu.pcts, diff.otu.modes, methods,value)  %>%
    unite('grp',c('diff.otu.modes','diff.otu.pcts')) %>% spread(c('grp'), value)  %>% unite('methods',c('scenario','methods'), sep = ': ')
  tpr <- tpr %>% column_to_rownames('methods')
  tpr <- tpr[fdr$methods,,drop =F] %>% rownames_to_column('methods')
  fdr;tpr
  colnames(tpr)[1] <- colnames(fdr)[1] <- 'Signal density'
  
  
  for(i in c(2:7)){
    tpr[,i] = format(round(tpr[,i], digits = 2), nsmall = 2)
  }
  tpr1 <- tpr; fdr1 <- fdr
  
  ## ---------------compostional ---------------
  names <- c('methods','scenario',"Low density", "Medium density","High density")
  fdr = rbind(aggregate(value~., data= res.obj$df3 %>% filter(measures =='FDR'&methods %in% c(target.method,best[3]) & diff.otu.modes=='Abundant') %>% dplyr::select(diff.otu.pcts,methods,value) , function(x)  mean(x[!is.na(x)])) %>%
                unite('grp',c('diff.otu.pcts')) %>% spread(c('grp'), value) %>% mutate(scenario ='stool_unbalanced'),
              aggregate(value~., data= res.obj$df4 %>% filter(measures =='FDR'&methods %in% c(target.method,best[4])& diff.otu.modes=='Abundant') %>% dplyr::select(diff.otu.pcts,  methods,value) , function(x)  mean(x[!is.na(x)])) %>%
                unite('grp',c('diff.otu.pcts')) %>% spread(c('grp'), value) %>% mutate(scenario ='vaginal_unbalanced'))
  idx <- which(fdr$methods ==target.method)
  idx1 <- which(fdr$methods !=target.method)
  idxx <- c(idx[1],idx1[1],idx[2],idx1[2])
  fdr = fdr[idxx,]
  
  fdr = fdr[,names] %>% unite('methods',c('scenario','methods'), sep = ': ') 
  rownames(fdr) <- NULL
  fdr = fdr %>% column_to_rownames('methods')
  
  
  # Lower CI: fdr.est - 1.96 * sqrt(fdr.est * (1 - fdr.est)/ num.iter))
  fdr.CI = rbind(res.obj$df3 %>% filter(measures =='FDR' &methods %in% c(target.method,best[3])& diff.otu.modes=='Abundant') %>% dplyr::select(diff.otu.pcts,  methods, ymin) %>% unite('grp',c('diff.otu.pcts')) %>% spread(c('grp'), ymin) %>% mutate(scenario ='stool_unbalanced'),
                 res.obj$df4 %>% filter(measures =='FDR' &methods %in% c(target.method,best[4])& diff.otu.modes=='Abundant') %>% dplyr::select(diff.otu.pcts,  methods, ymin) %>% unite('grp',c('diff.otu.pcts')) %>% spread(c('grp'), ymin) %>% mutate(scenario ='vaginal_unbalanced'))
  idx <- which(fdr.CI$methods ==target.method)
  idx1 <- which(fdr.CI$methods !=target.method)
  idxx <- c(idx[1],idx1[1],idx[2],idx1[2])
  fdr.CI = fdr.CI[idxx,]
  
  fdr.CI = fdr.CI[,names] %>% unite('methods',c('scenario','methods'), sep = ': ') 
  rownames(fdr.CI) <- NULL
  fdr.CI = fdr.CI %>% column_to_rownames('methods')
  
  dim(fdr);dim(fdr.CI);fdr[1:3,1:3];fdr.CI[1:3,1:3]
  test <- function(v1, v2) ifelse(v1<=0.05,"***", ifelse(v2<=0.1, "**", ifelse(v2<=0.2, "*", 'x')))
  label.fdr <- mapply(test,v1 = fdr.CI, v2 = fdr)
  colnames(label.fdr) = colnames(fdr);rownames(label.fdr) = rownames(fdr)
  
  # Weak effect szie 
  fdr = label.fdr
  # fdr = as.data.frame(fdr[,c(1:3)])
  # fdr$score = apply(fdr, 1, function(x) str_count(x, "\\*") %>% sum())
  # colnames(fdr) = gsub('\\..*','',colnames(fdr))
  fdr = as.data.frame(fdr) %>% rownames_to_column('methods')
  
  
  data1 <- rbind(res.obj$df3 %>% mutate(scenario = 'stool_unbalanced') %>% dplyr::select(-factor) %>% filter(methods %in% c(target.method,best[3])& diff.otu.modes=='Abundant'),
                 res.obj$df4 %>% mutate(scenario = 'vaginal_unbalanced')%>% dplyr::select(-factor) %>% filter(methods %in% c(target.method,best[4])& diff.otu.modes=='Abundant'))
  
  tpr <- data1 %>% filter(measures %in% 'TPR') %>% dplyr::select(scenario, diff.otu.pcts, methods,value)  %>%
    unite('grp',c('diff.otu.pcts')) %>% spread(c('grp'), value)  %>% unite('methods',c('scenario','methods'), sep = ': ')
  tpr <- tpr %>% column_to_rownames('methods')
  tpr <- tpr[fdr$methods,,drop =F] %>% rownames_to_column('methods')
  fdr;tpr
  colnames(tpr)[1] <- colnames(fdr)[1] <- 'Signal density'
  
  for(i in c(2:4)){
    tpr[,i] = format(round(tpr[,i], digits = 2), nsmall = 2)
  }
  
  tpr2 <- tpr; fdr2 <- fdr
  
  ## --------combine fdr and tpr--------
  tpr <- cbind(tpr1, tpr2); fdr <- cbind(fdr1, fdr2)
  colnames(tpr)[8] <- colnames(fdr)[8] <- "Signal density "
  
  tpr$Data <- str_to_title(gsub('_.*','',tpr$`Signal density`))
  fdr$Data <- str_to_title(gsub('_.*','',fdr$`Signal density`))
  tpr$Methods <- gsub('.*: ','',tpr$`Signal density`)
  tpr$`Methods ` <- gsub('.*: ','',tpr$`Signal density `)
  fdr$Methods <- gsub('.*: ','',fdr$`Signal density`)
  fdr$`Methods ` <- gsub('.*: ','',fdr$`Signal density `)
  
  names <- c("Data",'Methods',"Abundant_Low density", "Abundant_Medium density","Abundant_High density" ,
             "Rare_Low density" ,"Rare_Medium density","Rare_High density",
             "Methods ","Low density", "Medium density","High density")
  tpr <- tpr[,names]
  fdr <- fdr[,names]
  colnames(tpr)[1] <-colnames(fdr)[1] <- ' '
  colnames(tpr) <- gsub('Methods', 'Signal density',colnames(tpr))
  colnames(tpr) <- gsub('Methods ', 'Signal density ',colnames(tpr))
  colnames(fdr) <- gsub('Methods', 'Signal density',colnames(fdr))
  colnames(fdr) <- gsub('Methods ', 'Signal density ',colnames(fdr))
  
  ##-------- design table --------
  bar_chart <- function(label, width = "100%", height = "10px", fill = "forestgreen", background = NULL) {
    bar <- div(style = list(background = fill, width = width, height = height))
    chart <- div(style = list(flexGrow = 1, marginLeft = "1px", background = background), bar)
    div(style = list(display = "flex", alignItems = "center"), label, chart)
  }
  
  get_color <- function(fdr){
    if(fdr =='***'){
      col = brewer.pal(8,'Dark2')[5]
    }else if(fdr =='**'){
      col = brewer.pal(8,'Set2')[6]
    }else if(fdr =='*'){
      col = brewer.pal(8,'RdGy')[2]
    }else{
      col = 'grey'
    }
    return(col)
  }
  
  
  # making table 
  # fdr.dt <- fdr[,-c(8),drop = F] # for colddef.list
  tpr[,c(3:8,10:12)] = sapply(tpr[,c(3:8,10:12)], function(x) round(as.numeric(x), digits = 2), USE.NAMES = F)
  # tpr[,c(8)] = sapply(tpr[,c(8)], function(x) round(as.numeric(x), digits = 1), USE.NAMES = F)
  
  # colnames(tpr)[1] <- colnames(fdr)[1] <- 'Signal density1'
  # names <- c("Signal density","Abundant_Low density", "Abundant_Medium density","Abundant_High density" ,
  #            "Rare_Low density" ,"Rare_Medium density","Rare_High density",
  #            "Signal density ","Low density", "Medium density","High density")
  # 
  # tpr <- tpr[,names];fdr <- fdr[,names]
  n1 <- gsub('\\ density.*','', colnames(tpr))
  n2 <- gsub('.*\\_','', n1)
  NM1 <- n2[grep('Low|Medium|High',n2)][1:6]
  NM2 <- n2[grep('Low|Medium|High',n2)][7:9]
  NM <- c(NM1, '',NM2)
  
  coldefs.list <- function(numcols){
    coldefs_list = NULL
    color_list = list()
    for (idx in c(3:8,10:12)){
      if(idx %in% c(3:8)){
        fdr_col = fdr[, (idx)] # column 8 is fdr score column, which does not match tpr table, need to be deleted; +1: fdr first column is methods
        colors <- sapply(fdr_col, function(x) get_color(x), USE.NAMES = F)
        name = numcols[idx-2]
        color_list[[name]] = colors
      }else{
        fdr_col = fdr[, (idx)] # column 8 is fdr score column, which does not match tpr table, need to be deleted; +1: fdr first column is methods
        colors <- sapply(fdr_col, function(x) get_color(x), USE.NAMES = F)
        name = numcols[idx-3]
        color_list[[name]] = colors
      }
      
      cell.func <- function(value, index, name) {
        width <- value
        bar_chart(format(round(value, digits = 2), nsmall = 2), width = width*80, background = brewer.pal(12,'Set3')[9],
                  fill = color_list[[name]][index])
      }
      
      style.func <- function(value, index, name) {
        color <- color_list[[name]][index]
        list(fontWeight = 800, fontSize = 16,color = color)
      }
      
      coldefs <- list(
        reactable::colDef(style = style.func, cell=cell.func, name = NM[idx-2], align = 'center')
      )
      
      coldefs_list = c(coldefs_list,coldefs)
    }
    
    # change names
    names(coldefs_list) <- numcols
    
    return (list(coldefs_list, color_list))
  }
  
  numcols <- colnames(tpr)[c(3:8,10:12)]
  
  coldefs_list <- coldefs.list(numcols)[[1]]
  color_list <- coldefs.list(numcols)[[2]]
  
  minW = 39
  # coldefs_list[[colnames(tpr)[1]]] = colDef(
  #   minWidth = 140, align = 'center', cell = function(value, index) {
  #   name = colnames(tpr)[1]
  #   name <- tpr[,1][index]
  #   tagList(
  #     div(style = list(fontWeight = 600, fontSize = 16), name)
  #   )
  # })
  
  coldefs_list[[colnames(tpr)[1]]] = colDef(style = JS("function(rowInfo, colInfo, state) {
                                var firstSorted = state.sorted[0]
                                // Merge cells if unsorted or sorting by school
                                if (!firstSorted || firstSorted.id === ' ') {
                                  var prevRow = state.pageRows[rowInfo.viewIndex - 1]
                                  if (prevRow && rowInfo.row[' '] === prevRow[' ']) {
                                    return { visibility: 'hidden' }
                                  }
                                }
                              }"),
                                            minWidth = 90, align = 'center', cell = function(value, index) {
                                              name = colnames(tpr)[1]
                                              name <- tpr[,1][index]
                                              tagList(
                                                div(style = list(fontWeight = 600, fontSize = 16), name)
                                              )
                                            })
  
  coldefs_list[[colnames(tpr)[2]]] = colDef(minWidth = 100, align = 'center', cell = function(value, index) {
    name = colnames(tpr)[2]
    name <- tpr[,2][index]
    tagList(
      div(style = list(fontWeight = 600, fontSize = 16), name)
    )
  })
  
  coldefs_list[[colnames(tpr)[9]]] = colDef(minWidth = 100, align = 'center', cell = function(value, index) {
    name = colnames(tpr)[9]
    name <- tpr[,9][index]
    tagList(
      div(style = list(fontWeight = 600, fontSize = 16), name)
    )
  })
  
  
  html_file <- paste0(output,target.method,'_',name,'.html')
  tb <- reactable(tpr, 
                  pagination=FALSE,
                  resizable = FALSE, wrap = FALSE, bordered = F,
                  style = list(fontFamily = "Work Sans, sans-serif", fontSize = "16px"),
                  rowStyle = JS("
                function(rowInfo, state) {

                var nextRow = state.pageRows[rowInfo.viewIndex + 1]

                if (nextRow && rowInfo.row[' '] !== nextRow[' ']) {
                      // Use box-shadow to add a 2px border without taking extra space
                      return { boxShadow: 'inset 0 -1px 0 rgba(0, 0, 0, 0.5)' }
                    }
                  }
              "),
                  
                  theme = reactableTheme(
                    headerStyle = list(
                      "&:hover[aria-sort]" = list(background = "#f7f7f8"),
                      "&[aria-sort='ascending'], &[aria-sort='descending']" = list(background = "hsl(0, 0%, 98%)"),
                      borderColor = "#555"
                    )
                  ),
                  
                  columnGroups = list(
                    colGroup(name = '', columns = colnames(tpr)[1]),
                    colGroup(name = 'Balanced', columns = colnames(tpr)[2]),
                    colGroup(name = gsub('_.*','',colnames(tpr)[3]), columns = colnames(tpr)[3:5]),
                    colGroup(name = gsub('_.*','',colnames(tpr)[6]), columns = colnames(tpr)[6:8]),
                    colGroup(name = 'Unbalanced', columns = colnames(tpr)[9]),
                    colGroup(name = gsub('_.*','',colnames(tpr)[3]), columns = colnames(tpr)[10:12])
                  ),
                  
                  columns = coldefs_list
  )  
  saveWidget(widget = tb, file = html_file, selfcontained = TRUE)
  
  
  
  img_file <- paste0(output,target.method,'_',name,'.png')
  webshot(url = html_file, file = img_file, vwidth = 1800, vheight = 1000, zoom = 3)
}






## mean(fdr) of all iters -- median(6 box) -- max(6 box)
fdr.class <- function(data, i = 2, grp, level, name, measure = 'FDR', diff.otu.mode= c('Rare','Abundant'), 
                      setting = 'Challenge:', type = 'stool'){
  colnames(data)[i] = 'grp'
  df <- filter(data, grp ==level & measures ==measure & diff.otu.modes %in% diff.otu.mode) %>% 
    dplyr::select(methods, value) %>% 
    group_by(methods) %>% 
    summarise(medianOFmedian=mean(value), maxOFmedian = max(value))
  colnames(df)[2:3] = paste0(name,'_',setting,c('medianFDR.','MaxMedianFDR.'),type)
  df[,2] <- apply(df[,2], 2, function(x) ifelse(x <= 0.05, 'Good', ifelse(x>0.1,'Poor','Intermediate')))%>% as.data.frame()
  df[,3] <- apply(df[,3], 2, function(x) ifelse(x <= 0.05, 'Good', ifelse(x>0.1,'Poor','Intermediate')))%>% as.data.frame()
  return(challeng.df = df)
}
## rank dataframe(kable shows) -- choose specific columns(settings) -- median(selected columns)
tpr.class <- function(data, grp, level, name, patn = 'Medium\\ effect', patn2 = NULL, setting = 'Challenge:', type = 'stool'){
  df <- as.data.frame(data)
  col.sel <- colnames(df)
  df <- df %>% dplyr::select(col.sel[grep(patn,col.sel)])
  if(!is.null(patn2)){
    col.sel <- colnames(df)
    df <- df %>% dplyr::select(col.sel[grep(patn2,col.sel)])
  }
  df <- apply(df, 1, function(x) mean(x)) %>% as.data.frame() %>% rownames_to_column('methods')
  colnames(df)[2] = paste0(name,'_',setting,'medianTPR.', type)
  df <- df %>% column_to_rownames('methods')
  df <- apply(df, 2, function(x) ifelse(x > 10, 'Good', ifelse(x<=5,'Poor','Intermediate')))%>% as.data.frame()
  df <- df %>% rownames_to_column('methods')
  return(df)
}




