library(GEOquery)
library(affy)
library(tidyverse)
#there are various normalization methods to normalize microarray data, but we will use RMA normalization
#get supplementary file
getGEOSuppFiles("GSE148537", makeDirectory = TRUE)

#Untar files
untar("GSE148537/GSE148537_RAW.tar", exdir = 'data/')

#reading cel files
raw.data<- ReadAffy(celfile.path = "data/")

#performing RMA normalization
normalized_data<-rma(raw.data)

#get expression estimates
normalized.expr<- as.data.frame(exprs(normalized_data))

#map probe id's to gene symbols
gse<- getGEO("GSE148537",GSEMatrix = TRUE)

#fetch feature data to get ID - gene symbol mapping
feature.data <- gse$GSE148537_series_matrix.txt.gz@featureData@data
#subset
feature.data <- feature.data[,c(1,11)]

#merge normalized expression data with feature data
normalized.expr<-normalized.expr %>%
  rownames_to_column(var = 'ID') %>%
  inner_join(.,feature.data, by= 'ID')

