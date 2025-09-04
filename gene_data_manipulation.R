library(dplyr)
library(tidyverse)
library(GEOquery)
library(ggplot2)
dat<-read.csv("C:/Users/Adros/Desktop/Bioinformatics/GSE183947_fpkm.csv")
head(dat)


#get metadata
gse<-getGEO(GEO = 'GSE183947',GSEMatrix = TRUE)
metadata<- pData(phenoData(gse[[1]]))
head(metadata)

metadata.modified<-metadata %>%
  select(1,10,11,17) %>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue: ","",tissue)) %>%
  mutate(metastasis = gsub("metastasis: ","",metastasis))


#reshaping data to longform
data.long<- dat %>%
  rename(gene = X) %>%
  gather(key='samples',value = 'FPKM',-gene)

#join long form data with metadata
data.long<-data.long %>%
  left_join(.,metadata.modified,by = c("samples" = "description"))
#basic analysis
data.long %>%
  filter(gene=='BRCA1'|gene =='BRCA2') %>%
  group_by(gene,tissue) %>%
  summarize(mean_FPKM = mean(FPKM),median_FPKM = median(FPKM)) %>%
  head()

#Generate visualizations with ggplot
ggplot(data.long,aes(x=metastasis,y = samples))+geom_col()
data.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(.,aes(x = samples,y=FPKM,fill = tissue)) +
  geom_col()

data.long %>%
  filter(gene=='BRCA1') %>%
  ggplot(.,aes(x = FPKM,fill= tissue)) + geom_density(alpha = 0.1)


data.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(.,aes(x = metastasis,y = FPKM))+geom_violin()


data.long %>%
  filter(gene == 'BRCA1'| gene== 'BRCA2') %>%
  spread(key = gene,value = FPKM) %>%
  ggplot(.,aes(x = BRCA1,y=BRCA2,color =tissue))+geom_point()+geom_smooth(method = 'lm',se =FALSE)

genes.of.interest<-c('BRCA1','BRCA2','TP53','ALK','MYCN')
data.long %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(.,aes(x = samples,y = gene,fill = FPKM)) +geom_tile()+scale_fill_gradient(low = 'white',high = 'red')