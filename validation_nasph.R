#GSE152075 - consisted of nasopharyngeal RNAseq data from 430 COVID-19 subjects and 54 negative controls 
library(edgeR)
library(ggplot2)
tab <- read.csv("GSE152075.csv")
summary(tab)
dim(tab)
#Size: 35785 485
tab[1,2]
countdata <- tab[(1:35784),(2:485)]
dim(countdata)
rownames(countdata) <- tab[,1]
rownames(countdata)
colnames(countdata)
myCPM <- cpm(countdata)
head(myCPM)
dgeObj <- DGEList(myCPM)
logcounts <- cpm(dgeObj,log=TRUE)

#FCGR3B, FFAR2, CCL3L1, TNFAIP6 
my_gene<-logcounts[c("FCGR3B", "FFAR2", "CCL3L1", "TNFAIP6"), ]  
write.csv(t(my_gene),"norm_data_GSE152075.csv")
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

#75 percentile value to fing ccl3l1+ and ccl3L1- cells
q75<-quantile(logcounts,0.75)
q75

#Use q75 to sort ccl3l1+ and ccl3l1- expression values for FCGR3B. Load the table - plot calculate 
ccl <- read.csv("fcg_ccl.csv")
#boxplot
g <- ggplot(ccl, aes(x=sample, y=expression, fill=sample))
g + geom_boxplot()+geom_point()+scale_fill_manual(values=c("#999999", "indianred1"))+ theme_classic()
#t-test
q<-t.test(expression~sample, data=ccl, var.equal = FALSE, 
          alternative = 'less')
q$p.value

