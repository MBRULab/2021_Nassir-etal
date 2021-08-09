#GSE152418 -PBMC RNAseq data from 17 COVID-19 and 17 healthy subjects. 
library(edgeR)
library(ggplot2)
tab <- read.csv("GSE152418.csv")
summary(tab)
dim(tab)
#Size: 60683    35
tab[1,2]
countdata <- tab[(1:60683),(2:35)]
dim(countdata)
rownames(countdata) <- tab[,1]
rownames(countdata)
colnames(countdata)
myCPM <- cpm(countdata)
head(myCPM)
dgeObj <- DGEList(myCPM)
logcounts <- as.matrix(cpm(dgeObj,log=TRUE))

rownames(logcounts)
dim(logcounts)
row.names(logcounts)
#FCGR3B, FFAR2, CCL3L1, TNFAIP6 
my_gene<-logcounts[c("ENSG00000162747", "ENSG00000126262", "ENSG00000276085", "ENSG00000123610"), ]  
write.csv(t(my_gene),"norm_data_GSE152418.csv")
#Rearranged and renamed the table into 4 diff tables with only 2 columns as below
#fcg.csv
#sample expression
#Covid_FCGR3B	1.44939523
#Covid_FCGR3B	2.240226055
#Control_FCGR3B	2.289875474
#Control_FCGR3B	1.969695405
#ffar.csv
#Covid_FFAR2	4.330644047
#Control_FFAR2	3.459198059
#ccl.csv
#Covid_CCL3L1	5.04868585
#Control_CCL3L1	0.999994229
#tnf.csv
#Covid_TNF	2.762210455
#Control_TNF	2.340863497
#Read each file one by one
expr <- read.csv("gene_name.csv")
#boxplot
g <- ggplot(expr, aes(x=sample, y=expression, fill=sample))
g + geom_boxplot()+geom_point()+scale_fill_manual(values=c("#999999", "indianred1"))+ theme_classic()
#t-test
q<-t.test(expression~sample, data=expr, var.equal = FALSE, 
          alternative = 'less')
q$p.value


