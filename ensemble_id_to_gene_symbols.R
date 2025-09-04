library(annotables)
library(biomaRt)
library(tidyverse)

#input list of ensembl id's
ensembl.ids<-read.delim('ensembl_ids')
#using biomart
listEnsembl()
ensembl<-useEnsembl(biomart = 'genes')
datasets<-listDatasets(ensembl)
ensembl.connection<-useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
attr<-listAttributes(ensembl.connection)
filters<-listFilters(ensembl.connection)
getBM(attributes = c('ensembl_gene_id','external_gene_name'),
      filters = "ensembl_gene_id",
      values = ensembl.ids,
      mart = ensembl.connection)


#using annotables 
grch38 %>%
  filter(ensgene %in% ensembl.ids$column)
