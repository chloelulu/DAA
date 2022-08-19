## Fig.S6: Comparison of the FDR control and power using different thresholds to select the reference set under (a) balanced and (b) unbalanced settings 
## Run in cluster
# setwd("/research/bsi/projects/staff_analysis/m216453/SemiSim/realanalysis/result03182022/")
# files <- list.files(getwd(), pattern = '*Rdata')
# 
# load(files[1])
# diff <- cbind.data.frame(err.var = out$zicoseq.obj$err.var[which(out$zicoseq.obj$p.adj.fdr<=0.05)])
# ref <- cbind.data.frame(err.var = out$zicoseq.obj$err.var[out$zicoseq.obj$ref.ind])
# err.var <- (rbind(diff %>% mutate(type = 'diff.taxa'),ref %>% mutate(type = 'ref.taxa')) %>% mutate(data = file))
# i = 1
# for(file in files){
#   i = i + 1
#   
#   cat(i,': ',file,'\n')
#   out <- NULL
#   load(file)
#   diff <- cbind.data.frame(err.var = out$zicoseq.obj$err.var[which(out$zicoseq.obj$p.adj.fdr<=0.05)])
#   ref <- cbind.data.frame(err.var = out$zicoseq.obj$err.var[out$zicoseq.obj$ref.ind])
#   err <- (rbind(diff %>% mutate(type = 'diff.taxa'),ref %>% mutate(type = 'ref.taxa')) %>% mutate(data = file))
#   err.var <- rbind(err.var, err)
# }
# 
# head(err.var)
# mean(err.var$err.var == 0)
# colnames(err.var)[1] = 'error variance'
# save(err.var, file = "/research/bsi/projects/staff_analysis/m216453/SemiSim/realanalysis/error.variance.Rdata")
# err.var$err.var <- log(err.var$err.var)

## local plot
load("/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/realdataanalysis/error.variance.Rdata")
selected <- read_csv("/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/realdataanalysis/dataset_107res/select.txt", col_names = FALSE)
files <- gsub('/Users/m216453/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/realdataanalysis/dataset_107res/','',selected$X1)
err.var <- err.var[err.var$data %in% files,]
unique(err.var$data)
err.var$type <- gsub('diff.taxa','differential taxa',err.var$type)
err.var$type <- gsub('ref.taxa','reference taxa',err.var$type)


ggplot(err.var, aes(x = type, y = `error variance`)) + 
  # stat_boxplot( geom='errorbar', linetype=1, width=0.3)+  
  # geom_violin(aes(fill = type)) + 
  geom_boxplot(outlier.size = 0.4,aes(fill = type)) +
  scale_fill_brewer(palette = 'Set1') + 
  scale_y_continuous(trans = sqrt_trans(),
                     breaks = trans_breaks("sqrt", function(x) x^2),
                     labels = trans_format("sqrt", math_format(.x^2))) +
  labs(fill = '', x= '') + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size = 20),
        axis.title = element_text(color="black", size = 20),
        axis.ticks.x = element_blank(),
        title = element_text(size = 20),
        legend.position = 'right',
        plot.margin = margin(1, 1, 1, 1, "cm"),
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 20))  
out <- "/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Result/"
ggsave(file = paste0(out, "SimulationEvaluation03122022/Fig.S6.pdf"), width = 6, height = 5)
  

