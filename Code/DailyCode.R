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



data_summary1 <- function(data, formula){
  ## formula: value ~ variables to be aggregated
  med <- aggregate(as.formula(formula), data, function(x) max(x[!is.na(x)]))
  return(med)
}


reduce_list <- function(data = data){Reduce(
  function(x, y, ...) {
    if (is.null(x)){
      y
    } else if (is.null(y)){
      x
    }
    else{
      merge(x, y, all = TRUE, ...)
    }},
  data
)}



reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}
scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}





count.ct = function(data = x4,j = 2){
  ct = NULL
  for(i in 1:nrow(data)){
    ct=c(ct,str_count(data[i,j], "\\*") %>% sum())
  }
  ct = gsub('^3$',brewer.pal(11,'PiYG')[9],ct)
  ct = gsub('^2$',brewer.pal(8,'Set2')[6],ct)
  ct = gsub('^1$',brewer.pal(8,'RdGy')[2],ct)
  ct = gsub('^0$',brewer.pal(12,'Set3')[9],ct)
  return(ct)
} 

count.score = function(data = x4,j = 11, color1 = brewer.pal(8,'Paired')[1], color2 = 'black', value = 27){
  ct = data[,j]
  ct[ct==value] = color1
  ct[-grep('^\\#',ct)] =color2
  return(ct)
}


cols <- c('ZicoSeq'=brewer.pal(9,'Paired')[5], 'Omnibus' = brewer.pal(9,'PuBuGn')[7], 'mbzinb' = brewer.pal(9,'PuBuGn')[8], 
          'GLM(quasipoisson)' = brewer.pal(9,'PuRd')[3],'LinDA'=brewer.pal(9,'Set1')[8],'LINDA2'=brewer.pal(9,'Set1')[7],
          'RAIDA'=brewer.pal(11,'BrBG')[7],'RioNorm2' = brewer.pal(11,'BrBG')[8],
          'eBay(Wilcoxon)' = brewer.pal(8,'RdPu')[4],'eBay(t-test)' = brewer.pal(8,'RdPu')[6], 'ANCOM-BC' = brewer.pal(8,'Set2')[7],
          'LDM'=brewer.pal(9,'Set1')[1], 'DACOMP'=brewer.pal(8,'Dark2')[6],'Aldex2(Wilcoxon)'=brewer.pal(8,'YlGn')[3],'Aldex2(t-test)'=brewer.pal(8,'YlGn')[2],
          'Beta-binomial'=brewer.pal(11,'BrBG')[2],'DESeq2'=brewer.pal(9,'Greys')[7], 'Wrench+DESeq2'= brewer.pal(9,'Greys')[3], 'GMPR+DESeq2'=brewer.pal(9,'Greys')[5],
          'TSS+Wilcoxon'=brewer.pal(9,'Purples')[7],'Rarefy+Wilcoxon'= brewer.pal(9,'Set1')[5], 'GMPR+Wilcoxon'= brewer.pal(9,'Purples')[5], 'Wrench+Wilcoxon'= brewer.pal(9,'Purples')[3], 
          'TSS+Spearman'=brewer.pal(9,'Purples')[7],'Rarefy+Spearman'=brewer.pal(9,'Set1')[5], 'GMPR+Spearman'= brewer.pal(9,'Purples')[5], 
          'TSS+t-test'=brewer.pal(9,'GnBu')[7],'Rarefy+t-test'= brewer.pal(9,'GnBu')[5], 'GMPR+t-test'= brewer.pal(9,'GnBu')[5], 'Wrench+t-test'= brewer.pal(9,'GnBu')[3], 
          'edgeR'=brewer.pal(9,'Greens')[7], 'Wrench+edgeR'=brewer.pal(9,'Greens')[3], 'GMPR+edgeR'=brewer.pal(9,'Greens')[5], 
          'metagenomeSeq'=brewer.pal(9,'Blues')[7], 'Wrench+metagenomeSeq'=brewer.pal(9,'Blues')[5])
