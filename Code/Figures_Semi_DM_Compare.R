#### Below code for FigureS7-8-9 in the manuscript for comparison of DM and SemiSimulation framework

pkg =  c('plyr','scales','ggpubr','microbiome',"eBay","modeest","ANCOMBC","aod","phyloseq","reshape","corncob","MASS","readr","DESeq2", "ALDEx2", "metagenomeSeq", "edgeR", "GUniFrac", "grDevices", "dirmult", "exactRankTests","nlme", "dplyr", "magrittr", "tidyr", "protoclust", "ggplot2", "compositions","rmutil","tibble","reticulate","dacomp","LDM","Wrench")
sapply(pkg, require, character = TRUE)

setwd("/Users/m216453/Documents/Mayo_Research/SemiSimulation/Submit/DAA/")
wd = getwd()


source(paste0(wd,'/code/Func.R'))
QC <- function(otu.tab, prev = 0.1, minp = 0.002){
  prop = t(t(otu.tab)/colSums(otu.tab))
  prop <- prop[rowSums(prop!=0) > prev * ncol(prop), , drop=FALSE]
  otu.tab <- otu.tab[rownames(prop), , drop=FALSE]

  prop <- prop[rowMaxs(prop) > minp, , drop=FALSE]
  otu.tab <- otu.tab[rownames(prop), , drop=FALSE]

  idx.f <- which(colSums(otu.tab)>1000)
  # idx.f <- apply(otu.tab, 1, function(x) sum(x > 2)) > round(ncol(otu.tab) * 0.05)
  otu.tab <- otu.tab[,idx.f]
  return(otu.tab = otu.tab)
}


load(paste0(wd,'/data/CompareDMandSemiparametric/Stool_V35.RData')) # OTU dataset
load(paste0(wd,'/data/CompareDMandSemiparametric/Stool_V35_dirmult.Rdata')) # For saving computational time, already prepared dirmult.paras
otu.tab = Stool_V35
name = 'Stool'
otu.tab <- QC(otu.tab = otu.tab)

name = 'UrogenitalTract'
load(paste0(wd,'/data/CompareDMandSemiparametric/','Semi_',name,'.Rdata'))
load(paste0(wd,'/data/CompareDMandSemiparametric/','DM_',name,'.Rdata'))

## SemiSimulation
nOTU = 'nOTU_L5'
nSam = 'nSam_L1'
diff.otu.pct= 'low'
diff.otu.mode = 'abundant'
diff.otu.direct = 'balanced'
covariate.type = 'binary'
covariate.eff.mean = 'none'
confounder.type = 'none'
depth.mu = 'D3'
depth.conf.factor = 'none'
include.top.otu = FALSE
model = 'loglinear'
nSub = 'Sub_L1'
model.paras = NULL
par(mfrow = c(2,2))
Sim.obj <- SimulateSeq(otu.tab = otu.tab,
                       nOTU = nOTU, diff.otu.pct = diff.otu.pct, diff.otu.direct = diff.otu.direct, diff.otu.mode = diff.otu.mode,
                       covariate.type = covariate.type, covariate.eff.mean = covariate.eff.mean, confounder.type = confounder.type, depth.mu =depth.mu, depth.conf.factor = depth.conf.factor,
                       model = model, nSub = nSub,include.top.otu = include.top.otu)
otu.tab.dir <- Sim.obj$otu.tab.sim
save(Sim.obj, file = paste0('Semi_',name,'.Rdata'))

## Dirichlet multinomial Simulation
source('code/Func_DM.R')
data <- SimulateMSeq(otu.tab = otu.tab, diff.otu.mode = 'abundant', diff.otu.direct = 'balanced', covariate.type = 'binary',
                     covariate.eff.mean = 0, confounder.eff.mean = 0, depth.mu = 10000, zinfl.otu.pct = 0)
otu.tab.DM <- data$otu.tab.sim
save(data, file = paste0('DM_',name,'.Rdata'))






for(name in c('Stool','UrogenitalTract')){
  load(paste0(wd,'/data/CompareDMandSemiparametric/','Semi_overfitting_',name,'.Rdata'))
  load(paste0(wd,'/data/CompareDMandSemiparametric/','DM_overfitting_',name,'.Rdata'))
  load(paste0(wd,'/data/CompareDMandSemiparametric/','sample_overfitting_',name,'.Rdata'))
  
  otu.tab = otu.nonsample
  otu.tab.dir <- Sim.obj$otu.tab.sim
  otu.tab.DM <- data$otu.tab.sim
  dim(otu.tab);dim(otu.tab.dir);dim(otu.tab.DM)
  if((nrow(otu.tab)==nrow(otu.tab.dir)) &nrow(otu.tab.DM)==nrow(otu.tab.dir)){
    ################################  Data rarefraction  #########################################
    set.seed(123)
    rary2 = GUniFrac::Rarefy(t(otu.tab), depth = 1000)$otu.tab.rff %>% t()
    set.seed(123)
    rary3 = GUniFrac::Rarefy(t(otu.tab.DM), depth = 1000)$otu.tab.rff %>% t()
    set.seed(123)
    rary4 = GUniFrac::Rarefy(t(otu.tab.dir), depth = 1000)$otu.tab.rff %>% t()
    
    
    rary2.r = t(t(rary2)/colSums(rary2))
    rary3.r = t(t(rary3)/colSums(rary3))
    rary4.r = t(t(rary4)/colSums(rary4))
    
    ################################  OTU distribution  #########################################
    # OTU prevalence
    sum(rary2 ==0);sum(rary3 ==0);sum(rary4 ==0)
    zero2= apply(rary2, 2, function(x) sum(x ==0)/nrow(rary2))
    zero3= apply(rary3, 2, function(x) sum(x ==0)/nrow(rary3))
    zero4= apply(rary4, 2, function(x) sum(x ==0)/nrow(rary4))
    sink(paste0(wd,'/Result/CompareDMandSemiparametric/',name,'_Sparsity_Overfitting.txt'))
    cat('Raw median, max, min \n')
    round(median(zero2),2);cat(paste0('(',round(min(zero2),2),'-',round(max(zero2),2),')'),'\n');
    cat('DM \n')
    round(median(zero3),2);cat(paste0('(',round(min(zero3),2),'-',round(max(zero3),2),')'),'\n');
    cat('Semi \n')
    round(median(zero4),2);cat(paste0('(',round(min(zero4),2),'-',round(max(zero4),2),')'),'\n');
    cat('raw vs DM:',wilcox.test(zero2,zero3)$p.value,'\n')
    cat('raw vs Semi:',wilcox.test(zero2,zero4)$p.value,'\n')
    sink()
    cols = c(`DM` = '#3db7e4', Semiparametric = '#ff8849', `Real`= '#69be28')
    
    bind.df = cbind(`DM` = zero3, `Real` = zero2, `Semiparametric` = zero4) %>% melt() 
    colnames(bind.df) <- c('X1','X2','value')
    
    bind.df$X2 <- factor(bind.df$X2, levels=c('Real','DM','Semiparametric'))
    
    bind.df %>% ggplot(aes(x = X2, y = value, fill = X2)) + 
      geom_violin() +
      geom_boxplot(width=0.1,outlier.size = 1, outlier.colour = 'grey30', fill = 'white') +
      geom_hline(yintercept = median(zero2), colour = '#69be28', linetype = 'dashed', size = 0.5) +
      geom_hline(yintercept = max(zero2), colour ='#69be28', linetype = 'dashed', size = 0.5) +
      geom_hline(yintercept = min(zero2), colour = '#69be28', linetype = 'dashed', size = 0.5) +
      ylab('Sparsity') +
      xlab('') +
      labs(fill = '') +theme_bw() +
      scale_fill_manual(values = cols) + 
      scale_y_continuous(breaks=c(round(min(zero2), digits = 2), round(median(zero2), digits = 2), round(max(zero2), digits = 2))) +
      theme(axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_text(color="black", size = 26),
            axis.title = element_text(color="black", size = 30),
            strip.text = element_text(size = 30),
            strip.background = element_rect(fill="white",color = "black", size = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black",size= 1),
            legend.title=element_text(size=30),
            legend.text = element_text(size=30),
            plot.title = element_text(size=22)) +
      ggtitle(name)
    ggsave(paste0(wd,'/result/CompareDMandSemiparametric/',name,'_Sparsity_violin_Overfitting.pdf'),  width = 8, height =3.5, dpi = 60, units = 'in')
    
    pre2 = as.data.frame(prevalence(rary2,count=FALSE)) %>% rownames_to_column('OTU')
    colnames(pre2)[2] = 'raw'
    pre3 = as.data.frame(prevalence(rary3,count=FALSE)) %>% rownames_to_column('OTU')
    colnames(pre3)[2] = 'DM'
    pre4 = as.data.frame(prevalence(rary4,count=FALSE)) %>% rownames_to_column('OTU')
    colnames(pre4)[2] = 'dirchlet'
    sink(paste0(wd,'/result/CompareDMandSemiparametric/',name,'_Prevelance_Overfitting.txt'))
    cat('Raw median, max, min \n')
    round(median(pre2[,2]),2);cat(paste0('(',round(min(pre2[,2]),2),'-',round(max(pre2[,2]),2),')'),'\n');
    cat('DM \n')
    round(median(pre3[,2]),2);cat(paste0('(',round(min(pre3[,2]),2),'-',round(max(pre3[,2]),2),')'),'\n');
    cat('Semi \n')
    round(median(pre4[,2]),2);cat(paste0('(',round(min(pre4[,2]),2),'-',round(max(pre4[,2]),2),')'),'\n');
    cat('raw vs DM:',wilcox.test(pre2[,2],pre3[,2])$p.value,'\n')
    cat('raw vs Semi:',wilcox.test(pre2[,2],pre4[,2])$p.value,'\n')
    sink()
    
    qq = inner_join(pre3, pre2) %>% inner_join(pre4);head(qq)
    qq1 = cbind(`DM` = zero3, `Real` = zero2, `Semiparametric` = zero4) %>% as.data.frame() %>% rownames_to_column('OTU')
    
    head(qq1)
    
    colnames(qq)[2:4] = c('DM','Real','Semiparametric')
    qq.m <- melt(qq)
    head(qq.m)
    qq.m$variable <- factor(qq.m$variable, levels=c('Real','DM','Semiparametric'))
    
    ggplot(qq.m,aes(x = variable, y = value, fill = variable)) + 
      geom_violin() +
      geom_boxplot(width=0.1,outlier.size = 0, outlier.colour = 'grey30', fill = 'white') +
      geom_hline(yintercept = median(pre2$raw), colour = '#69be28', linetype = 'dashed', size = 0.5) +
      geom_hline(yintercept = max(pre2$raw), colour ='#69be28', linetype = 'dashed', size = 0.5) +
      geom_hline(yintercept = min(pre2$raw), colour = '#69be28', linetype = 'dashed', size = 0.5) +
      ylab('Prevelance') +
      xlab('') +
      labs(fill = '') +theme_bw() +
      scale_fill_manual(values = cols) + 
      scale_y_continuous(breaks=c(round(min(pre2$raw), digits = 2), round(median(pre2$raw), digits = 2), round(max(pre2$raw), digits = 2))) +
      theme(axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_text(color="black", size = 26),
            axis.title = element_text(color="black", size = 30),
            strip.text = element_text(size = 30),
            strip.background = element_rect(fill="white",color = "black", size = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black",size= 1),
            legend.title=element_text(size=30),
            legend.text = element_text(size=30),
            plot.title = element_text(size=22)) +
      ggtitle(name)
    ggsave(paste0(wd,'/result/CompareDMandSemiparametric/',name,'_PrevelanceBoxplot_Overfitting.pdf'), width = 8, height =5, dpi = 60, units = 'in')
    
    
    
    qqplot3 <- function(qq, title = title, digits = 2,name = name, QQ= T,
                        sqrt= F, # sqrt transform of the dataset, then plot density plot
                        sqrtqq = F, # sqrt transform of the dataset, then plot qqplot
                        xlab = 'Prevelance (Real)',ylab = 'Prevelance (Compared)',xdensity = ''){
      colnames(qq)[2:4] = c('DM','Real','Semiparametric')
      cols = c(`DM` = '#3db7e4',Semiparametric = '#ff8849', `Real`= '#69be28')
      
      qq_23 = as.data.frame(stats::qqplot(qq[,'Real'], qq[,'DM'], plot.it=FALSE)) %>% dplyr::rename(`DM` = y, `Real` = x)
      qq_23 <- unique(qq_23)
      qq_24 = as.data.frame(stats::qqplot(qq[,'Real'], qq[,'Semiparametric'], plot.it=FALSE)) %>% dplyr::rename(Semiparametric = y, `Real` = x)
      qq_24 <- unique(qq_24)
      df = rbind(qq_23 %>% mutate(method ='DM') %>% dplyr::rename(value = `DM`),
                 qq_24 %>%mutate(method ='Semiparametric') %>% dplyr::rename(value = Semiparametric))
      head(df)
      
      fmt_dcimals <- function(decimals=0){
        # return a function responpsible for formatting the 
        # axis labels with a given number of decimals 
        function(x) as.character(round(x,decimals))
      }
      df$id = paste0('id',1:nrow(df))
      head(df)
      df$method <- factor(df$method, levels=c('Real','DM','Semiparametric'))
      if(QQ){
        if(sqrtqq){
          qq234 = df %>% dplyr::mutate(`Real` = sqrt(`Real`), value = sqrt(value))  %>% ggplot() + 
            geom_abline(intercept = 0, color = 'grey10', cex = 1.5, alpha = 0.7) +
            geom_point(aes(x=`Real`, y=value, color = method), size = 0.4) +
            scale_x_continuous(labels =fmt_dcimals(digits)) +
            scale_y_continuous(labels =fmt_dcimals(digits)) +
            scale_color_manual(values = cols) + 
            theme_bw()+
            theme(axis.text = element_text(color="black", size = 16),
                  plot.margin = unit(c(1,1,1,1), "cm"),
                  legend.text= element_text(color="black", size = 16),
                  legend.title = element_text(color="black", size = 16),
                  axis.title = element_text(color="black", size = 16),
                  legend.position = 'none')+
            labs(x=xlab,y=ylab)
        } else{
          qq234 = ggplot(df) + 
            geom_abline(intercept = 0, color = 'grey10', cex = 1.5, alpha = 0.7) +
            geom_point(aes(x=`Real`, y=value, color = method), size = 0.4) +
            scale_x_continuous(labels =fmt_dcimals(digits)) +
            scale_y_continuous(labels =fmt_dcimals(digits)) +
            scale_color_manual(values = cols) + 
            theme_bw()+
            theme(axis.text = element_text(color="black", size = 16),
                  plot.margin = unit(c(1,1,1,1), "cm"),
                  legend.text= element_text(color="black", size = 16),
                  legend.title = element_text(color="black", size = 16),
                  axis.title = element_text(color="black", size = 16),
                  legend.position = 'none') +
            labs(x=xlab,y=ylab)
        } 
      }
      
      
      
      # qq23 = ggplot(qq_23) + 
      #   geom_abline(intercept = 0, color = 'grey10', cex = 1.5, alpha = 0.7) +
      #   geom_point(aes(x=`Real`, y=`DM`), color = brewer.pal(8,'Dark2')[3]) +
      #   theme_minimal() + theme_classic() +
      #   theme_bw()+
      #   theme(axis.text.x = element_text(color="black", size = 18,vjust = 0.5),
      #         axis.text.y = element_text(color="black", size = 18),
      #         axis.title = element_text(color="black", size = 18),
      #         legend.position = 'none')+
      #   labs(x='Real',y='DM')
      # qq24 = ggplot(qq_24) + 
      #   geom_abline(intercept = 0, color = 'grey10', cex = 1.5, alpha = 0.7) +
      #   geom_point(aes(x=`Real`, y= Semiparametric), color = 'forestgreen') +
      #   # geom_jitter(aes(x=raw, y= dirchlet), color = 'forestgreen') +
      #   theme_minimal() + theme_classic() +
      #   theme_bw()+
      #   theme(axis.text.x = element_text(color="black", size = 18,vjust = 0.5),
      #         axis.text.y = element_text(color="black", size = 18),
      #         axis.title = element_text(color="black", size = 18),
      #         legend.position = 'none')+
      #   labs(x='Real',y='Semiparametric')
      # qq1234 = ggarrange(qq23, qq24, nrow = 2)
      
      prev = melt(qq)
      prev$variable <- factor(prev$variable, levels=c('Real','DM','Semiparametric'))
      if(sqrt){
        prev = prev %>% dplyr::mutate(value = sqrt(value))
      }
      
      # mu <- aggregate(value~ variable, data = prev, function(x) mean(x))
      ## histogram plot
      # p1 = ggplot(prev, aes(x=value, color=variable)) +
      #   geom_histogram(fill="white", position="dodge")+
      #   geom_vline(data=mu, aes(xintercept=grp.mean, color=variable),linetype="dashed", cex = 1)+
      #   scale_color_brewer(palette="Dark2") +
      #   theme_minimal() +theme_classic() +
      #   theme(axis.text.x = element_text(color="black", size = 18,vjust = 0.5),
      #         axis.text.y = element_text(color="black", size = 18),
      #         axis.title = element_text(color="black", size = 18),
      #         legend.position = 'right')+
      #   labs(x='',y='count', title = title)
      p2 = ggplot(prev, aes(x=value, fill=variable)) +
        geom_density(alpha=0.5)+
        scale_fill_manual(values = cols) +
        theme_minimal() + theme_classic() +
        theme(axis.text = element_text(color="black", size = 16,vjust = 0.5),
              plot.margin = unit(c(1,1,1,1), "cm"),
              legend.text= element_text(color="black", size = 16),
              legend.title = element_text(color="black", size = 16),
              axis.title = element_text(color="black", size = 16),
              legend.position = 'right')+
        labs(x = xdensity, y = 'Density', fill = '') +
        scale_y_continuous(trans = sqrt_trans(),
                           breaks = trans_breaks("sqrt", function(x) x^2),
                           labels = trans_format("sqrt", math_format(.x^2)))
      # prev1230 = ggarrange(p1, p2, nrow = 2)
      plt = ggarrange(p2, qq234, nrow = 1, common.legend = TRUE)# %>% annotate_figure(top = name)
      return(plt)
    }
    
    qqplot3(qq, title = 'OTU prevalence',xdensity = 'OTU prevalence', sqrt = F, name = name)
    ggsave(paste0(wd,'/result/CompareDMandSemiparametric/',name,'_OTU prevalence_Overfitting.pdf'), width = 8, height =4, dpi = 60, units = 'in')
    qqplot3(qq1, title = 'Sparsity',xdensity = 'Sparsity', sqrt = F, xlab = 'Sparsity(Real)',ylab = 'Sparsity(Compared)', name = name)
    ggsave(paste0(wd,'/result/CompareDMandSemiparametric/',name,'_Sparsity_Overfitting.pdf'), width = 8, height =4, dpi = 60, units = 'in')
    
    # OTU mean abundance
    mn2 = as.data.frame(rowMeans(rary2.r)) %>% rownames_to_column('OTU')
    colnames(mn2)[2] ='raw'
    mn3 = as.data.frame(rowMeans(rary3.r)) %>% rownames_to_column('OTU')
    colnames(mn3)[2] = 'DM'
    mn4 = as.data.frame(rowMeans(rary4.r)) %>% rownames_to_column('OTU')
    colnames(mn4)[2] = 'dirchlet'
    mn_abund = inner_join(mn3, mn2) %>%inner_join(mn4)
    # mn_abund[,-1] = sqrt(mn_abund[,-1])
    head(mn_abund)
    
    tran2 = asin(sqrt(rary2.r))
    tran3 = asin(sqrt(rary3.r))
    tran4 = asin(sqrt(rary4.r))
    
    
    var2 = apply(asin(sqrt(rary2.r)), 1, function(x) var(x))
    var3 = apply(asin(sqrt(rary3.r)), 1, function(x) var(x))
    var4 = apply(asin(sqrt(rary4.r)), 1, function(x) var(x))
    var = 
      as.data.frame(var2) %>% rownames_to_column('OTU') %>%
      inner_join(as.data.frame(var3) %>% rownames_to_column('OTU')) %>% 
      inner_join(as.data.frame(var4) %>% rownames_to_column('OTU')) %>% 
      dplyr::rename(raw = var2, DM = var3, dirchlet = var4)
    head(var)
    
    ## DM var > raw var
    var$diff = var$DM - var$raw
    id = var[var$diff > 0,]$OTU # DM> raw
    id1 = var[var$diff < 0,]$OTU
    
    df <- melt(mn_abund %>% dplyr::rename(`DM`=DM,`Real`=raw,Semiparametric = dirchlet) %>% dplyr::filter(OTU %in% id1))
    df$variable <- factor(df$variable, levels=c('Real','DM','Semiparametric'))
    ggplot(df,
           aes(x = variable, y = value, fill = variable)) + 
      geom_violin() +
      geom_boxplot(width = 0.1, fill = 'white', outlier.colour = NA) +
      ylab('log(Mean abundance)') +
      xlab('') +
      labs(fill = '') +theme_bw() +
      scale_fill_manual(values = cols) + 
      theme(axis.text.x = element_blank(),
            plot.margin = unit(c(1,1,1,1), "cm"),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(color="black", size = 26),
            axis.title = element_text(color="black", size = 30),
            strip.text = element_text(size = 30),
            strip.background = element_rect(fill="white",color = "black", size = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black",size= 1),
            legend.title=element_text(size=30),
            legend.text = element_text(size=30),
            plot.title = element_text(size=22)) +
      scale_y_continuous(trans = sqrt_trans(),
                         breaks = trans_breaks("sqrt", function(x) x^2),
                         labels = trans_format("sqrt", math_format(.x^2))) +
      geom_hline(yintercept = median(mn_abund$raw), colour = 'red', linetype = 'dashed', size = 0.5) 
    
    ggsave(paste0(wd,'/result/CompareDMandSemiparametric/',name,'_OTU mean abundance violin_Overfitting.pdf'), width = 10, height = 5, dpi = 60, units = 'in')
    
    ggplot(mn_abund %>% filter(OTU %in% id))+
      geom_abline(intercept = 0, color = 'grey30', cex = 1.5, alpha = 0.7) +
      geom_point(aes(x=sqrt(raw), y=sqrt(DM)), color = '#3db7e4') +
      geom_point(aes(x=sqrt(raw), y=sqrt(dirchlet)), color ='#ff8849') +
      ylab('Mean abundance (compared)') +
      xlab('Mean abundance (Real)') +
      labs(fill = '') +theme_bw() +
      theme(axis.text = element_text(color="black", size = 26),
            plot.margin = unit(c(1,1,1,1), "cm"),
            axis.title = element_text(color="black", size = 30),
            strip.text = element_text(size = 30),
            strip.background = element_rect(fill="white",color = "black", size = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black",size= 1),
            legend.title=element_text(size=30),
            legend.text = element_text(size=30),
            plot.title = element_text(size=22))
    ggsave(paste0(wd,'/result/CompareDMandSemiparametric/',name,'_OTU mean abundance scatter_Overfitting.pdf'),width = 10, height = 7, dpi = 60, units = 'in')
    
    # id1: DM > raw variation; id: raw > DM variation
    qqplot3(qq= mn_abund %>% filter(OTU %in% id),sqrt = T,sqrtqq = T,
            title = 'OTU mean abundance',digits=3, xlab = 'sqrt(Mean abundance (Real))', 
            ylab =  'sqrt(Abundance (Compared))', xdensity = 'sqrt(Mean abundance)', name = name)
    ggsave(paste0(wd,'/result/CompareDMandSemiparametric/',name,'_OTU mean abundance_lowVarInRaw_Overfitting.pdf'),  width = 10, height = 5, dpi = 60, units = 'in')
    qqplot3(qq= mn_abund %>% filter(OTU %in% id1),sqrt = T,sqrtqq = T,
            title = 'OTU mean abundance',digits=3, xlab = 'sqrt(Mean abundance (Real))', 
            ylab =  'sqrt(Abundance (Compared))', xdensity = 'sqrt(Mean abundance)', name = name)
    ggsave(paste0(wd,'/result/CompareDMandSemiparametric/',name,'_OTU mean abundance_highVarInRaw_Overfitting.pdf'),  width = 10, height = 5, dpi = 60, units = 'in')
    
    qqplot3(qq= mn_abund,sqrt = T,sqrtqq = T,
            title = 'OTU mean abundance',digits=3, xlab = 'sqrt(Mean abundance (Real))', 
            ylab =  'sqrt(Abundance (Compared))', xdensity = 'sqrt(Mean abundance)', name = name)
    ggsave(paste0(wd,'/result/CompareDMandSemiparametric/',name,'_OTU mean abundance_Overfitting.pdf'),  width = 10, height = 5, dpi = 60, units = 'in')
    
    # OTU abundance variation (sd / mu or arcsin(sqrt))
    
    tran2 = asin(sqrt(rary2.r))
    tran3 = asin(sqrt(rary3.r))
    tran4 = asin(sqrt(rary4.r))
    
    
    var2 = apply(tran2, 1, function(x) var(x))
    var3 = apply(tran3, 1, function(x) var(x))
    var4 = apply(tran4, 1, function(x) var(x))
    sqrt(mean(var2));sqrt(mean(var3));sqrt(mean(var4))
    sqrt(median(var2));sqrt(median(var3));sqrt(median(var4))
    var = 
      as.data.frame(var3) %>% rownames_to_column('OTU') %>%
      inner_join(as.data.frame(var2) %>% rownames_to_column('OTU')) %>% 
      inner_join(as.data.frame(var4) %>% rownames_to_column('OTU')) %>% 
      dplyr::rename(DM = var3, raw = var2, dirchlet = var4)
    # DM var > raw var
    var$diff = var$DM - var$raw
    id = var[var$diff > 0,]$OTU # DM> raw
    id1 = var[var$diff < 0,]$OTU
    var = var %>% dplyr::select(-diff)
    df = melt(var %>% dplyr::rename(`DM`=DM,`Real`=raw,Semiparametric = dirchlet) %>% filter(OTU %in% id1))
    df$variable <- factor(df$variable, levels=c('Real','DM','Semiparametric'))
    ggplot(df,
           aes(x = variable, y = value, fill = variable)) + 
      geom_violin() +
      geom_boxplot(width = 0.1, fill = 'white', outlier.colour = NA) +
      ylab('Variation of abundance') +
      xlab('') +
      labs(fill = '') +theme_bw() +
      scale_fill_manual(values = cols) + 
      theme(axis.text.x = element_blank(),
            plot.margin = unit(c(1,1,1,1), "cm"),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(color="black", size = 26),
            axis.title = element_text(color="black", size = 30),
            strip.text = element_text(size = 30),
            strip.background = element_rect(fill="white",color = "black", size = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black",size= 1),
            legend.title=element_text(size=30),
            legend.text = element_text(size=30),
            plot.title = element_text(size=22)) +
      scale_y_continuous(trans = sqrt_trans(),
                         breaks = trans_breaks("sqrt", function(x) x^2),
                         labels = trans_format("sqrt", math_format(.x^2))) +
      geom_hline(yintercept = median((var %>% dplyr::rename(`DM`=DM,`Real`=raw,Semiparametric = dirchlet) %>% filter(OTU %in% id1))[,'Real']), colour = 'red', linetype = 'dashed', size = 0.5) 
    ggsave(paste0(wd,'/result/CompareDMandSemiparametric/',name,'_OTU abundance variation violin lowDM_Overfitting.pdf'),  width = 10, height = 6, dpi = 100, units = 'in')
    
    df = melt(var %>% dplyr::rename(`DM`=DM,`Real`=raw,Semiparametric = dirchlet) %>% filter(OTU %in% id))
    df$variable <- factor(df$variable, levels=c('Real','DM','Semiparametric'))
    
    ggplot(df,
           aes(x = variable, y = value, fill = variable)) + 
      geom_violin() +
      geom_boxplot(width = 0.1, fill = 'white', outlier.colour = NA) +
      ylab('Variation of abundance') +
      xlab('') +
      labs(fill = '') +theme_bw() +
      scale_fill_manual(values = cols) + 
      theme(axis.text.x = element_blank(),
            plot.margin = unit(c(1,1,1,1), "cm"),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(color="black", size = 26),
            axis.title = element_text(color="black", size = 30),
            strip.text = element_text(size = 30),
            strip.background = element_rect(fill="white",color = "black", size = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black",size= 1),
            legend.title=element_text(size=30),
            legend.text = element_text(size=30),
            plot.title = element_text(size=22)) +
      scale_y_continuous(trans = sqrt_trans(),
                         breaks = trans_breaks("sqrt", function(x) x^2),
                         labels = trans_format("sqrt", math_format(.x^2))) +
      geom_hline(yintercept = median((var %>% dplyr::rename(`DM`=DM,`Real`=raw,Semiparametric = dirchlet) %>% filter(OTU %in% id))[,'Real']), colour = 'red', linetype = 'dashed', size = 0.5) 
    ggsave(paste0(wd,'/result/CompareDMandSemiparametric/',name,'_OTU abundance variation violin highDM_Overfitting.pdf'),  width = 10, height = 6, dpi = 100, units = 'in')
    
    
    qqplot3(var%>% filter(OTU %in% id1), title = 'OTU abundance variation',digits=3, sqrt = F,sqrtqq = F,
            xlab = 'sqrt(Variance of abundance (Real))', ylab =  'sqrt(Variance of abundance (Compared))', xdensity = 'sqrt(Variance of abundance)', name = name)
    ggsave(paste0(wd,'/result/CompareDMandSemiparametric/',name,'_OTU abundance variation_Overfitting.pdf'),  width = 10, height = 5, dpi = 100, units = 'in')
    
    
    
    ########### correlation structure ########################
    # taxa-taxa correlation
    clean.data <- function(data){
      data = data[rowSums(data!=0) > 0 * ncol(data), , drop=FALSE]
      data = data[rowMaxs(data) > 0, , drop=FALSE]
      return(data)
    }
    
    tax2 = clean.data(data = rary2.r)
    # tax2 = rary2.r
    r.rary2.r = t(tax2)
    is.matrix(r.rary2.r)
    cor2 = cor(r.rary2.r, method = 'spearman')
    diag(cor2)=NA
    mcor2 = melt(cor2) %>% na.omit()
    mcor2 = mcor2 %>% tidyr::unite('relation',1:2,sep = '_', remove = T)%>% dplyr::rename(raw = value)
    
    tax3 = clean.data(data = rary3.r)
    # tax3 = rary3.r
    r.rary3.r = t(tax3)
    is.matrix(r.rary3.r)
    cor3 = cor(r.rary3.r, method = 'spearman')
    diag(cor3)=NA
    mcor3 = melt(cor3) %>% na.omit() %>% tidyr::unite('relation',1:2,sep = '_', remove = T)%>% dplyr::rename(DM = value)
    
    tax4 = clean.data(data = rary4.r)
    # tax4 = rary4.r
    r.rary4.r = t(tax4)
    is.matrix(r.rary4.r)
    cor4 = cor(r.rary4.r, method = 'spearman')
    diag(cor4)=NA
    mcor4 = melt(cor4) %>% na.omit() %>% tidyr::unite('relation',1:2,sep = '_', remove = T)%>% dplyr::rename(dirchlet = value)
    
    cor = inner_join(mcor3, mcor2) %>% inner_join(mcor4) 
    head(cor);dim(cor)
    
    sink(paste0(wd,'/result/CompareDMandSemiparametric/',name,'_correlation_Overfitting.txt'))
    summary(cor)
    cat(round(median(cor$DM)),'\n')
    cat(paste0('DM: ','(',round(min(cor$DM),2),'-',round(max(cor$DM),2),')','\n'))
    cat(paste0('raw: ','(',round(min(cor$raw),2),'-',round(max(cor$raw),2),')','\n'))
    cat(paste0('Semi: ','(',round(min(cor$dirchlet),2),'-',round(max(cor$dirchlet),2),')','\n'))
    sink()
    
    
    
    qqplot3(cor, title = 'Taxa-Taxa(proportion 0.1% > 5% samples) correlation structure', sqrtqq = F, sqrt = F, QQ = T,
            xlab = 'Real', 
            ylab =  'Compared', 
            xdensity = 'Taxa-taxa correlation')
    ggsave(paste0(wd,'/result/CompareDMandSemiparametric/',name,'_taxa-tax-correlation_filtered_Overfitting.pdf'),  width = 8, height =4, dpi = 30,units = 'in')
    
    
    ################################  Sample distribution  #########################################
    # !!! shall be based on observed proportion
    # (nonzero  distribution at the same depth <alpha-diversity distribution>, pairwise distance distribution)
    div2 = vegan::diversity(rary2, index = "shannon", MARGIN = 1, base = exp(1)) %>% as.matrix()%>% as.data.frame() %>% rownames_to_column('sampleid')%>% dplyr::rename(raw = V1)
    div3 = vegan::diversity(rary3, index = "shannon", MARGIN = 1, base = exp(1)) %>% as.matrix()%>% as.data.frame() %>% rownames_to_column('sampleid')%>% dplyr::rename(DM = V1)
    div4 = vegan::diversity(rary4, index = "shannon", MARGIN = 1, base = exp(1)) %>% as.matrix()%>% as.data.frame() %>% rownames_to_column('sampleid')%>% dplyr::rename(dirchlet = V1)
    div0 = inner_join(div3, div2) %>% inner_join(div4) ;head(div0)
    
    sink(paste0(wd,'/result/CompareDMandSemiparametric/',name,'_Shannon_Overfitting.txt'))
    cat('Raw median, max, min \n')
    round(median(div2[,2]),2);cat(paste0('(',round(min(div2[,2]),2),'-',round(max(div2[,2]),2),')'),'\n');
    cat('DM \n')
    round(median(div3[,2]),2);cat(paste0('(',round(min(div3[,2]),2),'-',round(max(div3[,2]),2),')'),'\n');
    cat('Semi \n')
    round(median(div4[,2]),2);cat(paste0('(',round(min(div4[,2]),2),'-',round(max(div4[,2]),2),')'),'\n');
    cat('raw vs DM:',wilcox.test(div2[,2],div3[,2])$p.value,'\n')
    cat('raw vs Semi:',wilcox.test(div2[,2],div4[,2])$p.value,'\n')
    sink()
    
    
    df = melt(div0%>% dplyr::rename(`DM`=DM,`Real`=raw,Semiparametric = dirchlet))
    df$variable <- factor(df$variable, levels=c('Real','DM','Semiparametric'))
    
    ggplot(df,
           aes(x = variable, y = value, fill = variable)) + 
      geom_violin() +
      geom_boxplot(width=0.1,outlier.size = 0, outlier.colour = 'grey30', fill = 'white') +
      geom_hline(yintercept = median(div2[,2]), colour = '#69be28', linetype = 'dashed', size = 0.5) +
      geom_hline(yintercept = max(div2[,2]), colour ='#69be28', linetype = 'dashed', size = 0.5) +
      geom_hline(yintercept = min(div2[,2]), colour = '#69be28', linetype = 'dashed', size = 0.5) +
      ylab('Shannon Diversity') +
      xlab('') +
      labs(fill = '') +theme_bw() +
      scale_fill_manual(values = cols) + 
      scale_y_continuous(breaks=c(round(min(div2[,2]), digits = 2), round(median(div2[,2]), digits = 2), round(max(div2[,2]), digits = 2))) +
      theme(axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_text(color="black", size = 26),
            axis.title = element_text(color="black", size = 30),
            strip.text = element_text(size = 30),
            strip.background = element_rect(fill="white",color = "black", size = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black",size= 1),
            legend.title=element_text(size=30),
            legend.text = element_text(size=30),
            plot.title = element_text(size=22)) +
      ggtitle(name)
    ggsave(paste0(wd,'/result/CompareDMandSemiparametric/',name,'_alphaBoxplot_Overfitting.pdf'),  width = 8, height =5, dpi = 60, units = 'in')
    
    
    qqplot3(div0, title = 'Shannon Diversity', xlab = 'Shannon diversity  (Real)', ylab =  'Shannon diversity (Compared)', xdensity = 'Shannon diversity', name = name)
    ggsave(paste0(wd,'/result/CompareDMandSemiparametric/',name,'_alphadiversity1_Overfitting.pdf'),  width = 10, height = 5, dpi = 60, units = 'in')
    
    
    
    ## PCoA compare sample dist
    a2 = otu.tab
    colnames(a2) = paste0('raw.', colnames(a2))
    a3 = otu.tab.DM
    colnames(a3) = paste0('DM.', colnames(a3))
    a4 = otu.tab.dir
    colnames(a4) = paste0('dirchlet.', colnames(a4))
    a2[1:5,1:5];a3[1:5,1:5];a4[1:5,1:5]
    intersect(rownames(a2),rownames(a3))
    intersect(rownames(a2),rownames(a4))
    dim(a2);dim(a3);dim(a4)
    a = inner_join(rownames_to_column(as.data.frame(a3),'otu'),rownames_to_column(as.data.frame(a2),'otu')) %>% 
      inner_join(rownames_to_column(as.data.frame(a4),'otu')) %>% column_to_rownames('otu')
    set.seed(123)
    rary = Rarefy(t(a), depth = 1000)$otu.tab.rff %>% t()
    head(a)
    # anosim(a, grouping, permutations = 999, distance = "bray", strata = NULL,
    #        parallel = getOption("mc.cores"))
    
    ta = t(a)
    dim(ta)
    ta.DM = ta[grep('^raw|^DM',rownames(ta)),]
    grp = gsub('\\..*','',rownames(ta.DM))
    ta.dist <- vegdist(ta.DM)
    set.seed(23)
    ano <- adonis2(ta.dist ~ grp, permutations = 999, method="bray")
    sink(paste0(wd,'/result/CompareDMandSemiparametric/',name,'_ANOSIM_DM_Overfitting.txt'))
    ano
    sink()
    
    
    ta.dirchlet = ta[grep('^raw|^dirchlet',rownames(ta)),]
    grp = gsub('\\..*','',rownames(ta.dirchlet))
    ta.dist <- vegdist(ta.dirchlet, method = 'bray')
    ano <- adonis2(ta.dist ~ grp, permutations = 999, method="bray")
    sink(paste0(wd,'/result/CompareDMandSemiparametric/',name,'_ANOSIM_Semiparametric_Overfitting.txt'))
    ano
    sink()
    
    
    
    m2 = colnames(a2) %>% as.data.frame() %>% mutate(grp ='Real')
    m3 = colnames(a3) %>% as.data.frame() %>% mutate(grp ='DM')
    m4 = colnames(a4) %>% as.data.frame() %>% mutate(grp ='Semiparametric')
    m = rbind(m2, m3, m4)
    colnames(m)[1] = 'SampleID'
    
    metadata <- sample_data(m %>% mutate(sampleid = SampleID) %>% column_to_rownames('sampleid'))
    otutable <- otu_table(rary, taxa_are_rows = T)
    rary <- merge_phyloseq(otutable,metadata) 
    
    pslog <- transform_sample_counts(rary, function(x) log(1 + x))
    out.pcoa.log <- ordinate(pslog,  method = "PCoA", distance = "bray")
    evals <- out.pcoa.log$values[,1]
    eng <- as.data.frame(out.pcoa.log$vectors)[,c(1:2)] %>% rownames_to_column('SampleID') %>% inner_join(as.data.frame(as.matrix(sample_data(rary))) )
    pc1 <- round((out.pcoa.log$values$Eigenvalues[1])/sum(out.pcoa.log$values$Eigenvalues),2)
    pc2 <- round((out.pcoa.log$values$Eigenvalues[2])/sum(out.pcoa.log$values$Eigenvalues),2)
    head(eng)
    eng$grp <- factor(eng$grp, levels=c('Real','DM','Semiparametric'))
    
    pcoa = ggplot(eng,aes(x = Axis.1, y = Axis.2,shape = grp, fill = grp)) +
      geom_point(size = 1.5)+
      stat_ellipse(geom = "polygon", alpha = 0.1,type = "norm", color = 'black')+ 
      scale_color_manual(values = cols)+
      scale_shape_manual(values = c(24, 23,21,11))+
      scale_fill_manual(values = cols)+
      labs(fill = '', shape = '') +
      xlab(paste0('PC1, ',(paste0(pc1*100,'%')))) + ylab(paste0('PC2, ',paste0(pc2*100, '%'))) + 
      theme_bw() +
      theme(axis.text = element_text(size = 18, color = "black"),
            axis.title = element_text(size = 18, color = "black"),
            legend.text= element_text(color="black", size = 18),
            legend.title = element_text(color="black", size = 18))
    pcoa
    ggsave(paste0(wd,'/result/CompareDMandSemiparametric/',name,'_PCOA_Overfitting.pdf'),  width = 8, height = 5, dpi = 60, units = 'in')
    
    
    
    
    
    
    ####################### heatmap of abudance 
    h2 = rary2.r; h3 = rary3.r; h4 = rary4.r
    dim(h2);dim(h3);dim(h4)
    h2 = h2[unique(c(rownames(h3),rownames(h4))),]
    r2 = sqrt(h2);r3 = sqrt(h3);r4 = sqrt(h4)
    
    TOPn = T
    if(TOPn){
      top = names(sort(rowMeans(r2), decreasing = T))[1:nrow(r2)]
      r2 = r2[top,];r3 = h3[top,];r4 = r4[top,]
    }
    
    cbind = cbind(r2, r3, r4)
    col.scheme <- brewer.pal(9, 'YlOrRd')
    max(r2);max(r3);max(r4)
    breaks <- c(0, seq(min(cbind[cbind != 0]), 0.2, len=6),seq(0.3,max(cbind),length = 3))
    # breaks <- c(0, seq(min(cbind[cbind != 0]) * 0.99, max(r2, r3, r4), len=6))
    save_pheatmap_png <- function(x, filename, width=1800, height=1700, res = 200) {
      png(filename, width = width, height = height, res = res)
      grid::grid.newpage()
      grid::grid.draw(x$gtable)
      dev.off()
    }
    
    heat2 = pheatmap::pheatmap(r2,legend = TRUE, cluster_rows = F, 
                               cluster_cols = F,clustering_method = 'ward.D2',color = col.scheme,
                               breaks = breaks,
                               annotation_legend = F,show_colnames = F,show_rownames = F,fontsize_row =6,fontsize_col =6,
                               main = paste0('Real'))
    save_pheatmap_png(heat2, paste0(wd,'/result/CompareDMandSemiparametric/',name,'_heat_raw_Overfitting.png'))
    
    heat3 = pheatmap::pheatmap(r3,legend = TRUE, cluster_rows = F, 
                               cluster_cols = F,clustering_method = 'ward.D2', color = col.scheme,
                               breaks = breaks,
                               annotation_legend = F,show_colnames = F,show_rownames = F,fontsize_row =6,fontsize_col =6,
                               main = paste0('DM'))
    
    save_pheatmap_png(heat3, paste0(wd,'/result/CompareDMandSemiparametric/',name,'_heat_DM_Overfitting.png'))
    
    heat4 = pheatmap::pheatmap(r4,legend = TRUE, cluster_rows = F, 
                               cluster_cols = F,clustering_method = 'ward.D2',color = col.scheme,
                               breaks = breaks,
                               annotation_legend = F,show_colnames = F,show_rownames = F,fontsize_row =6,fontsize_col =6,
                               main = paste0('Semiparametric'))
    save_pheatmap_png(heat4, paste0(wd,'/result/CompareDMandSemiparametric/',name,'_heat_SemiSimulation_Overfitting.png'))
    
  }else{
    break
  }
}














