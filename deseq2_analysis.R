#DESeq2 Workflow

#Step 1: Read the data and get it in right format 

counts_data<- read.csv('counts_data.csv')
head(counts_data)

#Read information about the column 
colData<-read.csv('sample_info.csv')

#make sure that row names in coldata matches to column  names in counts data
all(colnames(counts_data) %in%rownames(colData)
)

all(colnames(counts_data)==rownames(colData))

#construct a DESeqdataset object
dds<-DESeqDataSetFromMatrix(countData  = counts_data,
                       colData = colData,
                       design = ~dexamethasone)

#pre-filtering: removing rows with low gene counts
#keeping rows that have atleast 10 reads total
keep<- rowSums(counts(dds))>=10
dds<-dds[keep,]

#set the factor level
dds$dexamethasone<-relevel(dds$dexamethasone, ref="untreated")

#Note: If a dataset has technical replicates we have to collapse them, Never collapse the biological replicates
#Run Deseq2
dds<- DESeq(dds)
res<-results(dds)
res
#baseMean : average of normalized counts taken over all the sample
#log2foldchange: fold change of gene in treatd condition compared to untraeted, +ve means upregulated genes in treated condition 
#lfcSE: Standard error for log2foldchange
#stat: wald test values for genes
#p-value: p-value of test statistic for the gene
#padj: corrected p values of multiple testing
summary(res)
res0.01<- results(dds, alpha = 0.01)
summary(res0.01)

#contrasts
resultsNames(dds)
#e.g. treated_4hrs, treated_8hrs, untreated
results(dds, contrast = c("dexamethasone","treated_4hrs", "untreated")) #if there are different levels

#MA plot
plotMA(res) #here triangles means these genes have highr chanegs
