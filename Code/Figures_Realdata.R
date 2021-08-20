#### Below code for Figure11ABC in the manuscript

pkg = c('readr','dplyr','tidyr','reshape','RColorBrewer','ggplot2','ggpubr','stringr','scales','tibble','kableExtra','formattable','reactable','htmltools',"htmlwidgets","webshot2","heatmaply","circlize","ComplexHeatmap","dendextend")
suppressPackageStartupMessages(sapply(pkg, require, character = T))

source('Code/Function.R')
setwd("/Users/m216453/Documents/Mayo_Research/SemiSimulation/Submit/DAA/")
wd <- getwd()
sub.methods <- c("TSS+Wilcox",'GMPR+Wilcox', 'GMPR+DESeq2','GMPR+edgeR','Wrench+MSeq',
                "RAIDA","corncob",'MaAsLin2',
                "ANCOM-BC","DACOMP","LDM","Omnibus",'ZicoSeq', "Aldex2(Wilcox)","GMPR+glm","eBay(Wilcox)")

sub.methods1 <- c('GMPR+Wilcox', 'GMPR+DESeq2','GMPR+edgeR','Wrench+MSeq',
                 "RAIDA","corncob",'MaAsLin2',
                 "ANCOM-BC","DACOMP","LDM","Omnibus",'ZicoSeq', "Aldex2(Wilcox)","GMPR+glm","eBay(Wilcox)")
selected <- read_csv("Data/realdataanalysis/dataset_107res/select.txt", col_names = FALSE)
files <- gsub('/Users/m216453/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/realdataanalysis/dataset_107res/','',selected$X1)

# files <- list.files('Data/realdataanalysis/dataset_107res/', pattern = 'Rdata$')
for(fdr in c(0.05, 0.1, 0.2)){
  RES <- name <- PCTs <- NULL;RESS <- OVERLAP <- jacrd <- RES.SUM <-failure <-  list();dist.mat <- list()
  for(i in 1:length(files)){
    cat(i,'\n')
    load(paste0(wd,'/Data/realdataanalysis/dataset_107res/',files[i]))
    idx <- colnames(res.sum)[grep('taxa|^fdr',colnames(res.sum))]
    if(!is.null(res.sum)){
      res <- res.sum[,idx] %>% column_to_rownames('taxa')
      colnames(res) <- gsub('fdr.','',colnames(res))
      colnames(res)[colnames(res) =='Wilcox'] = 'TSS+Wilcox'
      colnames(res)[colnames(res) =='Rarefy'] = 'Rarefy+Wilcox'
      colnames(res) = gsub('ANCOMBC','ANCOM-BC',colnames(res))
      colnames(res) = gsub('glmquassi','GMPR+glm',colnames(res))
      colnames(res) = gsub('DESeq2.gmpr','GMPR+DESeq2',colnames(res))
      colnames(res) = gsub('DESeq2.Wrench','Wrench+DESeq2',colnames(res))
      colnames(res) = gsub('edgeR.gmpr','GMPR+edgeR',colnames(res))
      colnames(res) = gsub('edgeR.Wrench','Wrench+edgeR',colnames(res))
      colnames(res) = gsub('MSeq2.Wrench','Wrench+MSeq',colnames(res))
      colnames(res) = gsub('MSeq2','metagenomeSeq',colnames(res))
      colnames(res) = gsub('eBayW','eBay(Wilcox)',colnames(res))
      colnames(res) = gsub('eBayt','eBay(t-test)',colnames(res))
      colnames(res) = gsub('BBinomial','Beta-binomial',colnames(res))
      colnames(res) = gsub('Aldex2we','Aldex2(t-test)',colnames(res))
      colnames(res) = gsub('^Aldex2$','Aldex2(Wilcox)',colnames(res))
      colnames(res) = gsub('^mbzinb','Omnibus',colnames(res))
      colnames(res) = gsub('^Wilcox.Wrench','Wrench+Wilcox',colnames(res))
      colnames(res) = gsub('^Wilcox.gmpr','GMPR+Wilcox',colnames(res))
      colnames(res) = gsub('BBinomial','Beta-binomial',colnames(res))
      colnames(res) = gsub('Maaslin2','MaAsLin2',colnames(res))
      len <- sub.methods[!(sub.methods %in% colnames(res))]
      if(length(len) > 0){
        for(k in 1:length(len)){
          res[,len[k]] <- rep(1, nrow(res))
        }
      }
      cat(files[i],'_',len,'\n')
      
      res <- res %>% dplyr::select(sub.methods)
      dim(res)
      ## fail of the dataset
      failure[[i]] <- which(apply(res, 2, function(x) length(grep('NA',x))==nrow(res)))
      
      RES.SUM[[i]] <- res
      
      res[is.na(res)] <- 1
      res[grep('NA',res)] <- 1
      
      res_no <- res[,sub.methods1];dim(res_no) # exclude TSS+Wilcox
      overlap <- matrix(NA, nrow = ncol(res_no), ncol = ncol(res_no))
      colnames(overlap) <- rownames(overlap) <- colnames(res_no)
      
      
      for(m1 in colnames(res_no)){
        for(m2 in colnames(res_no)){
          x1 <- res_no[,c(m1, m2)]
          overlap[m1,m2] <- ifelse(sum(apply(x1, 1, function(x) sum(x <=fdr)>0)) ==0, 1, sum(apply(x1, 1, function(x) sum(x <=fdr)==2))/sum(apply(x1, 1, function(x) sum(x <=fdr)>0)))
        }
      }
      
      OVERLAP[[i]] <- overlap
      
      ## calculate distance
      iidx <- apply(res, 1, function(x)sum(x <= fdr)) > 0
      res_0 <- res[iidx,,drop =F]
      jacrd[[i]] <- proxy::dist(t(res_0), method = 'Jaccard')
      
      pct <- apply(res, 2, function(x) sum(x<=fdr)/length(x))
      
      PCTs <- rbind(PCTs, pct)
      
      dist.mat[[i]] <- res
      ress <- t(apply(res,1, function(x) {ifelse(x > fdr, 0, sum(x <=fdr))}))
      
      tab0 <- rep(NA, ncol(res));names(tab0) <- 1:ncol(res)
      tab0 <- as.data.frame(tab0) %>% rownames_to_column('freq')
      
      for(j in 1:nrow(tab0)){
        tab <- as.data.frame(table(ress[,j])/nrow(ress))
        colnames(tab)[1:2] <- c('freq',colnames(ress)[j])
        tab0 <- full_join(tab0, tab)
      }
      tab0 <- tab0 %>% dplyr::select(-tab0)
      tab0[is.na(tab0)] <- 0
      
      RESS[[gsub('.Rdata','',files[i])]] <- tab0
      
      comf <- apply(res, 2, function(x) sum(x<=fdr))
      name <- c(name, gsub('.Rdata','',files[i]))
      RES <- rbind(RES,comf)
    }
  }
  
  names <- gsub('.Rdata','',files)
  rownames(RES) <- name
  rownames(PCTs) <- gsub('_res','',name)
  save(RESS,RES, dist.mat, PCTs, OVERLAP, jacrd, RES.SUM, failure, file = paste0(wd,'/Data/realdataanalysis/alldata_',fdr,'_res_norarefy.Rdata'))
}



## Check failing of ZicoSeq in which dataset
txt <- NULL
for(i in 1:length(RES.SUM)){
  df = RES.SUM[[i]][,'ZicoSeq']
  if(length(grep('NA',df))>0){
    txt <- c(txt, files[i])
    cat(files[i],':',length(grep('NA',df)),'\n')
  }
}

###################
## Figure 11b overlap at 0.05, 0.1, 0.2 levels, exclude TSS+Wilcox and Rarefy+Wilcox to decrease bias
output = paste0(wd,'/Result/realdataanalysis/Figure11b/');if(!(dir.exists(output))){dir.create(output)}
sub.methods <- c('GMPR+Wilcox','GMPR+DESeq2','GMPR+edgeR','Wrench+MSeq',"RAIDA",
                  "ANCOM-BC","DACOMP","LDM","Omnibus",'ZicoSeq',
                  "corncob",'MaAsLin2',"Aldex2(Wilcox)","GMPR+glm","eBay(Wilcox)")

for(fdr in c(0.05)){
  load(paste0(wd,'/Data/realdataanalysis/alldata_',fdr,'_res_norarefy.Rdata'))
  ## overlap figure
  for(i in 1:length(OVERLAP)){
    n= sum(is.na(OVERLAP[[i]]))
    if(n >0){
      print(paste0(i,':',n,'\n'))
    }
  }
  
  mean.overlap <- OVERLAP[[1]]
  colnames(mean.overlap) <- gsub("Beta-binomial",'corncob',colnames(mean.overlap))
  rownames(mean.overlap) <- gsub("Beta-binomial",'corncob',rownames(mean.overlap))
  mean.overlap <- mean.overlap[,colnames(mean.overlap) %in% sub.methods]
  mean.overlap <- mean.overlap[rownames(mean.overlap) %in% sub.methods,]
  dim(mean.overlap)
  for(i in 2:length(OVERLAP)){
    mean.overlap1 <- OVERLAP[[i]]
    colnames(mean.overlap1) <- gsub("Beta-binomial",'corncob',colnames(mean.overlap1))
    rownames(mean.overlap1) <- gsub("Beta-binomial",'corncob',rownames(mean.overlap1))
    mean.overlap1 <- mean.overlap1[,colnames(mean.overlap1) %in% sub.methods]
    mean.overlap1 <- mean.overlap1[rownames(mean.overlap1) %in% sub.methods,]
    
    mean.overlap <- mean.overlap + mean.overlap1
  }
  mean.overlap <- mean.overlap/length(OVERLAP)
  colnames(mean.overlap) <- gsub('Wilcox','Wilcox',colnames(mean.overlap))
  colnames(mean.overlap) <- gsub('Wrench\\+metagenomeSeq','Wrench+MSeq',colnames(mean.overlap))
  colnames(mean.overlap) <- gsub('Beta\\-binomial','corconb',colnames(mean.overlap))
  colnames(mean.overlap) <- gsub('GLM\\(quasipoisson\\)','Wrench+glm',colnames(mean.overlap))
  
  rownames(mean.overlap) <- gsub('Wilcox','Wilcox',rownames(mean.overlap))
  rownames(mean.overlap) <- gsub('Wrench\\+metagenomeSeq','Wrench+MSeq',rownames(mean.overlap))
  rownames(mean.overlap) <- gsub('Beta\\-binomial','corconb',rownames(mean.overlap))
  rownames(mean.overlap) <- gsub('GLM\\(quasipoisson\\)','Wrench+glm',rownames(mean.overlap))
  
  sink(paste0(wd,'/Result/realdataanalysis/Figure11b/summary_mean_overlap_fdr_',fdr,'.txt'))
  cat(paste0('-----------When FDR is ', fdr, '.---------\n' ))
  for(k in 1:ncol(mean.overlap)){
    x <- mean.overlap[,k]
    x <- x[!names(x) %in% colnames(mean.overlap)[k]]
    cat(colnames(mean.overlap)[k],': \n')
    print(summary(x))
    cat('\n')
  }
  
  cat('-----------Overall summary------------ \n Max :',max(mean.overlap[lower.tri(mean.overlap)]))
  cat('\n')
  cat('Min :',min(mean.overlap[lower.tri(mean.overlap)]),'\n')
  cat('\n')
  cat('Mean :',mean(mean.overlap[lower.tri(mean.overlap)]),'\n')
  cat('\n')
  cat('Median :',median(mean.overlap[lower.tri(mean.overlap)]),'\n')
  sink()
  
  
  df <- as.matrix(mean.overlap)
  html_file <- paste0(output,'overlap_fdr_',fdr,'.html')
  tb <- heatmaply_cor(
    df,colors = rev(brewer.pal(11, 'RdYlGn')),
    point_size_mat = mean.overlap,
    limits = c(0,1), 
    node_type = "scatter", 
    trace = "none",margins = c(40, 40, 40, 40),
    dist_method = NULL,fontsize_col = 14,fontsize_row = 14,
    heatmap_layers = theme(axis.text=element_text(colour="black", family = 'Arial', angle = 90)),
    hclust_method = "complete")
  saveWidget(widget = tb, file = html_file, selfcontained = TRUE)
  
  img_file <- paste0(output,'Figure11B.png')
  webshot(url = html_file, file = img_file, vwidth = 900, vheight = 800, zoom = 2)
}





sub.methods <- c("TSS+Wilcox",'GMPR+Wilcox', 'GMPR+DESeq2','GMPR+edgeR','Wrench+MSeq',
                  "RAIDA","corncob",'MaAsLin2',
                  "ANCOM-BC","DACOMP","LDM","Omnibus",'ZicoSeq', "Aldex2(Wilcox)","GMPR+glm","eBay(Wilcox)")

load(paste0(wd,'/Data/realdataanalysis/alldata_0.05_res_norarefy.Rdata'))
colnames(RES) <- gsub('Beta-binomial','corncob',colnames(RES))
## importing dataset information table
datainfo <- read.csv(paste0(wd,'/Data/realdataanalysis/SampleTaxaAlpha.csv'), row.names = 1)
rownames(datainfo) <- NULL
colnames(datainfo) <- c("data","Alpha diversity", "Taxa number","Sample size" )

sort(apply(RES, 2, function(x) mean(x)))
sort(apply(RES, 2, function(x) median(x)))


RES2 <- as.data.frame(t(apply(RES, 1,function(x) scale(x, scale = T, center = T))))
rownames(RES2) <- gsub('_res','',rownames(RES))
colnames(RES2) <- colnames(RES)


datainfo1 <- datainfo[datainfo$data %in% rownames(RES2),]
rownames(datainfo1) <- NULL
datainfo1 <- datainfo1 %>% column_to_rownames('data')
datainfo1 <- datainfo1[order(-datainfo1$`Sample size`),,drop =F]

library(dendsort)
RES2 <- RES2[rownames(datainfo1),,drop =F] %>% na.omit() 
RES2 = RES2[,sub.methods]
dim(RES2)
datainfo1 <- datainfo1[rownames(RES2),,drop=F]
# write.csv(datainfo1, file = '/Users/m216453/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/realdataanalysis/106Realdata.csv')
col_fun = colorRamp2(c(seq(min(RES2),0, len = 9), 0,seq(0, max(RES2),len = 9)), c(rev(brewer.pal(9, 'YlGn')), "white", brewer.pal(9, 'YlOrRd')))

col_fun_sample = colorRamp2(c(10,50,100,200,400,500,700,1000,1800), brewer.pal(9, 'Purples'))
# colorRamp2(c(seq(min(datainfo1$`Sample size`),median(datainfo1$`Sample size`)*0.95,len = 4),median(datainfo1$`Sample size`),seq(median(datainfo1$`Sample size`) *1.1,max(datainfo1$`Sample size`), len = 3)), c(brewer.pal(9, 'Reds')[2:5],brewer.pal(9, 'Reds')[6],brewer.pal(9, 'Reds')[7:9]))
col_fun_taxa = colorRamp2(c(seq(min(datainfo1$`Taxa number`),median(datainfo1$`Taxa number`)*0.95,len = 4),
                            median(datainfo1$`Taxa number`),seq(median(datainfo1$`Taxa number`) *1.1,max(datainfo1$`Taxa number`), len = 3)), c(brewer.pal(9, 'Blues')[2:5],brewer.pal(9, 'Blues')[5],brewer.pal(9, 'Blues')[7:9]))
col_fun_alpha = colorRamp2(c(seq(min(datainfo1$`Alpha diversity`),median(datainfo1$`Alpha diversity`)*0.95,len = 4),
                             median(datainfo1$`Alpha diversity`),seq(median(datainfo1$`Alpha diversity`) *1.1,max(datainfo1$`Alpha diversity`), len = 3)), c(brewer.pal(9, 'Purples')[2:5],brewer.pal(9, 'Purples')[5],brewer.pal(9, 'Purples')[7:9]))

ha = rowAnnotation(df=datainfo1[,-1,drop =F], simple_anno_size = unit(0.4, "cm"), 
                   col = list(`Sample size` = col_fun_sample,
                              `Taxa number` = col_fun_taxa))
lgd <- Legend(col_fun = col_fun_sample, title = 'Sample size',at = c(10,50,100,200,400,500,700,1000,1800))
col_dend = (hclust(dist(t(RES2))))
# col_dend = color_branches(col_dend, k =4)
column_ha = HeatmapAnnotation(` ` = anno_boxplot(RES2,axis = F,gp = gpar(fill = cols)))
ht <- Heatmap(as.matrix(RES2), heatmap_legend_param = list(title = 'Standardized No. findings '),
              column_km = 4,
              col = col_fun, row_names_gp = gpar(fontsize =7), show_row_names = F, #clustering_method_columns = "ward.D2",
              row_order = c(1:nrow(RES2)),#cluster_columns = col_dend,
              #cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.1f", RES[i, j]), x, y, gp = gpar(fontsize = 10))},
              border = T,rect_gp = gpar(col = "black", lwd = 0.5), show_row_dend =F, bottom_annotation = column_ha, 
              left_annotation = ha)
ht
pdf(paste0(wd,'/result/realdataanalysis/Figure11a.pdf'),width = 6, height = 13)
draw(ht)#, heatmap_legend_side = "left", annotation_legend_side = "right", annotation_legend_list = lgd
dev.off()
dim(RES2)
dim(RES2)
sort(apply(RES2, 2, function(x) median(x)))
sort(apply(RES2, 2, function(x) mean(x)))

load(paste0(wd,'/data/realdataanalysis/alldata_0.05_res_norarefy.Rdata'))
## #### ## failure of method and dataset
fails <- NULL
for (i in 1:length(failure)){
  if(length(failure[[i]]) > 0){
    fails <- rbind(fails,paste0(names(failure[[i]]),':',files[i]))
  }
}
fails



###############
## Figure11C: covaerage
## RES.SUM is the list contains all res files for selected methods
## for all methods
length(sub.methods)
overlap <- coverage <- list()
for(i in 1:length(RES.SUM)){
  for(fdr in c(0.05, 0.1, 0.2)){
    data <- RES.SUM[[i]]
    colnames(data) <- gsub("Beta-binomial", "corncob",colnames(data))
    data <- data[,sub.methods, drop =F]
    data[is.na(data)] <- 1
    all.comm <- apply(data, 1, function(x) sum(x<=fdr)) > 0# union
    all.comm1 <- apply(data, 1, function(x) sum(x<=fdr)) == length(sub.methods) # intersect 
    data1 <- data[all.comm,,drop =F]
    data2 <- data[all.comm1,,drop =F]
    name <- paste0(i,'_',fdr)
    overlap[[name]] <- nrow(data2)/nrow(data1)
    coverage[[name]] <- nrow(data1)/nrow(data)
  }
}


coverage1 <- melt(coverage)
coverage1$value <- (coverage1$value)*100
colnames(coverage1) <- c('coverage','dataset_fdr')
head(coverage1)
coverage1$dataset <- gsub('\\_.*','',coverage1$dataset_fdr)
coverage1$fdr <- gsub('.*_','',coverage1$dataset_fdr)
head(coverage1)
coverage1 <- coverage1[,-2]
head(coverage1)
aggregate(coverage ~fdr, coverage1, function(x)mean(x))
aggregate(coverage ~fdr, coverage1, function(x)median(x))
summary(coverage1[coverage1$fdr==0.05,][,1])
summary(coverage1[coverage1$fdr==0.1,][,1])
summary(coverage1[coverage1$fdr==0.2,][,1])

ggplot(coverage1, aes(x = fdr, y = coverage, color = fdr)) + 
  theme_bw() +
  stat_boxplot(geom='errorbar', linetype=1, width=0.5)+
  geom_boxplot() +
  geom_jitter(aes(colour = fdr), width = 0.25) +
  scale_color_manual(values = c('forestgreen','orange','red')) + 
  # geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.4, size = 0.5,position=position_dodge2(.7, preserve = "single")) +
  # scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
  labs(y = 'coverage,%', x = 'FDR', color = "", fill = '') +
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
ggsave(paste0(wd,'/result/realdataanalysis/Figure11C.pdf'), width = 6,height = 5)






