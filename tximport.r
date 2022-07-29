#install tximport and tximportData from bioconductor

library(tximport)
library(tximportData)
library(readr)
library(jsonlite)
library(edgeR)
#prep file paths
dir <- "/scratch/henrygeo/RNAproject/Results/salmon_quant/"
samples <- read.table("/scratch/henrygeo/RNAproject/Results/samples.txt", header = T)
files <- file.path(dir, samples$sample, "quant.sf")
all(file.exists(files))
names(files)<-samples$sample

#check they all exist


#set up the transcript to gene table
tx2gene <- read_tsv("g2t.txt")
colnames(tx2gene)<- c("GENEID", "TXNAME")
tx2g <- subset(tx2gene, select = c(TXNAME, GENEID))
head(tx2g)

txi <- tximport(files, type = "salmon", tx2gene = tx2g)
names(txi)txi



y<-DGEList(counts = txi$counts, group = samples$sample)
y<-calcNormFactors(y, doWeighting=F)
norm_counts<-cpm(y)

#only transcipts with a summed TPM of >1
subcount <- norm_counts[rowSums(norm_counts)>1,]

#only those with moderate variance (mean is a little above 2 before subset)
highvar <-subcount[apply(subcount,1,var)>=1,]

#scaling for downstream
thvar<-t(highvar)
thvar<-apply(thvar,2,scale)
rownames(thvar)<-samples$sample
#back to normal orientation
scalevar<-t(thvar)

write.csv(scalevar, file = "scaledfiltered.csv")

library(tidyr)
library(ggplot2)
df <-read_csv("scaledfiltered.csv")
names(df)[1]<-"gene"
df<-gather(df, key = "ID", value = "TMM", -gene)
ggplot(df)+geom_point(aes(x=ID,y=TMM))