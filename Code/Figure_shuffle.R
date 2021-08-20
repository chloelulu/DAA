#### Below code for Figure11D and supplementary Figure10 in the manuscript

pkg = c('readr','dplyr','tidyr','reshape','RColorBrewer','ggplot2','ggpubr','stringr','scales','tibble','kableExtra','formattable','reactable','htmltools',"htmlwidgets","webshot2","heatmaply","circlize")
suppressPackageStartupMessages(sapply(pkg, require, character = T))
setwd("/Users/m216453/Documents/Mayo_Research/SemiSimulation/Submit/DAA/")
wd <- getwd()
source('~/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Code/Function.R')
sub.methods <- c("Rarefy+Wilcox","TSS+Wilcox",'GMPR+Wilcox', 'GMPR+DESeq2','GMPR+edgeR',
                 'Wrench+MSeq',"RAIDA","corncob",'MaAsLin2',"ANCOM-BC",
                  "DACOMP","LDM","Omnibus",'ZicoSeq', "Aldex2(Wilcox)","GMPR+glm","eBay(Wilcox)")



# files <- list.files('/Users/m216453/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/realdataanalysis/shuffle/',pattern = 'Rdata$')
# f <- NULL
# for(i in 1:length(files)){
#   sub <- strsplit(files[i],split= '_')[[1]]
#   unique_data <- paste0(sub[!(sub %in% tail(sub, n=2))], collapse = '_')
#   iter <- head(tail(sub, n=2),1)
#   f <- c(f, unique_data)
# }
# f1 <- unique(f)

selected <- read_csv("/Users/m216453/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/realdataanalysis/dataset_107res/select.txt", col_names = FALSE)
f1 <- gsub('/Users/m216453/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/realdataanalysis/dataset_107res/','',selected$X1)
f1 = gsub('_res.Rdata','',f1)

SIG <- list();iter_file <- NULL
for(f in f1){
  for(j in 1:50){
    try({
      load(paste0('/Users/m216453/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/realdataanalysis/shuffle/new/',f,'_',j,'_res.Rdata'))
      idx <- c('taxa',colnames(res.sum)[grep('fdr',colnames(res.sum))])
      cat(f, '\n')
      res <- as.data.frame(res.sum[,idx]) %>% column_to_rownames('taxa')
      colnames(res) <- gsub('fdr.','',colnames(res))
      colnames(res)[colnames(res) =='Wilcox'] = 'TSS+Wilcox'
      colnames(res)[colnames(res) =='Rarefy'] = 'Rarefy+Wilcox'
      colnames(res) = gsub('ANCOMBC','ANCOM-BC',colnames(res))
      colnames(res) = gsub('glmquassi','GMPR+glm',colnames(res))
      colnames(res) = gsub('DESeq2.gmpr','GMPR+DESeq2',colnames(res))
      colnames(res) = gsub('DESeq2.Wrench','Wrench+DESeq2',colnames(res))
      colnames(res) = gsub('edgeR.gmpr','GMPR+edgeR',colnames(res))
      colnames(res) = gsub('edgeR.Wrench','Wrench+edgeR',colnames(res))
      colnames(res) = gsub('^MSeq2.Wrench$','Wrench+MSeq',colnames(res))
      colnames(res) = gsub('^MSeq2$','MSeq',colnames(res))
      colnames(res) = gsub('eBayW','eBay(Wilcox)',colnames(res))
      colnames(res) = gsub('eBayt','eBay(t-test)',colnames(res))
      colnames(res) = gsub('BBinomial','corncob',colnames(res))
      colnames(res) = gsub('Aldex2we','Aldex2(t-test)',colnames(res))
      colnames(res) = gsub('^Aldex2$','Aldex2(Wilcox)',colnames(res))
      colnames(res) = gsub('^mbzinb','Omnibus',colnames(res))
      colnames(res) = gsub('^Wilcox.Wrench','Wrench+Wilcox',colnames(res))
      colnames(res) = gsub('^Wilcox.gmpr','GMPR+Wilcox',colnames(res))
      colnames(res) = gsub('Maaslin2','MaAsLin2',colnames(res))
      
      
      # len <- sub.methods[!(sub.methods %in% colnames(res))]
      # if(length(len) > 0){
      #   for(k in 1:length(len)){
      #     res[,len[k]] <- rep(1, nrow(res))
      #   }
      # }
      
      res[is.na(res)] <- 1
      res[grep('NA',res)] <- 1
      len <- sub.methods[(sub.methods %in% colnames(res))]
      
      res <- res[,len]
      
      sig <- apply(res, 2, function(x) sum(x <=0.05))
      iter_file <- c(iter_file, paste0(f,'.',j))
      
      
      sig0 = as.data.frame(sig)
      colnames(sig0) = paste0(f,'.',j)
      sig0 = sig0 %>% rownames_to_column('method')
      
      SIG[[paste0(f,'.',j)]] <- sig0
      # SIG <- rbind(SIG, sig)
    })
  }
}



d1 <- SIG[[1]]
for(i in 2:length(SIG)){
  d1 <- full_join(d1, SIG[[i]])
}

dim(d1)
d1[1:5,1:5]
save(SIG,iter_file,d1, file = '~/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/realdataanalysis/realshuffle_res_new.Rdata')




load('~/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/realdataanalysis/realshuffle_res.Rdata')
datasetinfo <- read.csv('/Users/m216453/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/realdataanalysis/106Realdata.csv')
colnames(datasetinfo)[1] <- 'dataset'
dd <- as.data.frame(t(d1 %>% column_to_rownames('method'))) %>% rownames_to_column('class')
dd$dataset <- gsub('\\..*','',dd$class);length(unique(dd$dataset))
dd <- dd %>% filter(dataset %in% datasetinfo$dataset)
dd1 <- dd %>% filter(dataset %in% datasetinfo$dataset) %>% right_join(datasetinfo %>% dplyr::select(dataset,Taxa.number ))
head(dd1);length(unique(dd1$dataset))
dd.pct <- as.data.frame(apply(dd1[,-c(1,ncol(dd1), ncol(dd1)-1)], 2, function(x) x/dd1[,ncol(dd1)]))
rownames(dd.pct) <- NULL
rownames(dd.pct) <- as.character(dd1$class)
# dd.pct$dataset <- gsub('\\..*','',rownames(dd.pct))
# rownames(dd.pct) <- NULL
head(dd.pct)

m <- apply(dd.pct, 2, function(x) mean(x[!is.na(x)]))
md <- apply(dd.pct, 2, function(x) median(x[!is.na(x)]))
se <- apply(dd.pct, 2, function(x) sd(x[!is.na(x)]) / sqrt(length(x[!is.na(x)])))
ymin <- m - 1.96 * se
ymin <- ifelse(ymin > 0, ymin, 0)
ymax <- m + 1.96 * se
dd.sum <- as.data.frame(cbind(mean= m, median=md, ymin, ymax)) %>% rownames_to_column('methods')
head(dd.sum);dim(dd.sum)


ggplot(dd.sum, aes(x = reorder(methods, mean), y = mean*100, fill = methods)) + geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax=ymax*100, ymin=ymin*100, linetype = NULL),  width=0.7, size = 0.2,position=position_dodge2(.7, preserve = "single")) +
  scale_fill_manual(values = cols) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4,color = 'black', size = 16),
        axis.text.y = element_text(color = 'black', size = 20),
        plot.title = element_text(size=16),
        plot.background = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(color = 'black', size = 20), legend.position = 'none') +
  labs(y = '% of differential taxa', x = '', title = 'Supplementary Figure 10')
ggsave('/Users/m216453/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Result/SupplementaryFigure10.pdf', width =6, height = 6)


## false discover iteration
head(dd)
ff <- dd
ff$iter <- gsub('.*\\.','', ff$class)
ff$dataset <- gsub('\\..*','', ff$class)
head(ff);length(unique(ff$dataset))
ff <- ff %>% dplyr::select(-class)
ff[,-c(ncol(ff):(ncol(ff)-1))] <- apply(ff[,-c(ncol(ff):(ncol(ff)-1))], 2, function(x) ifelse(x>0,1, 0))
ff1 <- aggregate(.~ dataset, ff%>% dplyr::select(-iter), function(x) sum(x[!is.na(x)]==1)/length(x[!is.na(x)]))
setdiff(unique(ff$dataset),ff1$dataset)

ff[is.na(ff)] <- -1
ff1 <- aggregate(.~ dataset, ff%>% dplyr::select(-iter), function(x) sum(x==1)/length(x[x!= -1]))
dim(ff1)
ff1.m <- melt(ff1) %>% na.omit()
head(ff1.m)
ff1.m1 = aggregate(value ~ variable,ff1.m, function(x) mean(x))
ff1.m1 = ff1.m1[order(ff1.m1$value),]
head(ff1.m1)
ggplot(ff1.m1, aes(x = reorder(variable, value), y = value, fill = variable)) + geom_bar(stat = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4,color = 'black', size = 14),
        axis.text.y = element_text(color = 'black', size = 14),
        axis.title.y = element_text(color = 'black', size = 14),
        legend.position = 'none') +
  geom_hline(yintercept = 0.05, color = 'forestgreen',linetype='dashed') + 
  geom_hline(yintercept = 0.1, color = 'yellow',linetype='dashed') + 
  geom_hline(yintercept = 0.2, color = 'red',linetype='dashed') + 
  scale_y_continuous(breaks=c(0.05, 0.1,0.2,0.5,0.75)) + 
  labs(y = 'Observed FDR', x = '', fill = '')
ggsave('~/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Result/realdataanalysis/Figure11D.pdf', width = 4, height = 6)

library(ggpubr)
ggarrange(p0, p1, ncol = 2)




