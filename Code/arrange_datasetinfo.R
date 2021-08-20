library(readxl);library(dplyr)
Table <- read_excel("~/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Result/SimulationEvaluation/Table.xlsx", 
                    sheet = "Supplementary Table 3") %>% dplyr::select(-c('taxa.number','sample.size'))
head(Table);dim(Table)

data_106 <- read.csv('/Users/m216453/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/realdataanalysis/106Realdata.csv')
colnames(data_106)[1] <- colnames(Table)[1]
head(data_106)



data <- read.csv('/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/datainfo.csv')
colnames(data)[1] <- colnames(Table)[1]
data <- data %>%filter(dataset_name %in% data_106$dataset_name)
head(data);head(Table);dim(data)
head(data)

data1 <- data %>% left_join(Table)
head(data1)
write.csv(data1, '~/Documents/Mayo_Research/SemiSimulation/Submit/DAA/Data/realdataanalysis/DataSetInfo.csv', row.names = F)
