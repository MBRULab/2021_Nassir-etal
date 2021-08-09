library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(ggplot2)
SeuratDisk::Convert("GSE155249_supplement.h5ad",dest='h5seurat', assay='RNA',overwrite = FALSE, verbose = TRUE)
rd <- LoadH5Seurat("GSE155249_supplement.h5seurat")
rd<-NormalizeData(rd)
nam <- c("FCGR3B", "FFAR2", "CCL3L1", "TNFAIP6")
cntrl <- FetchData(object=rd, vars = nam, cells=names(rd$Patient[rd$Patient==6]))
write.csv(cntrl,'control.csv')
CoV <- FetchData(object=rd, vars = nam,cells=c(names(rd$Patient[rd$Patient==1]),
                                              names(rd$Patient[rd$Patient==2]),
                                              names(rd$Patient[rd$Patient==3]),
                                              names(rd$Patient[rd$Patient==4]),
                                              names(rd$Patient[rd$Patient==5])))
write.csv(CoV,'cov.csv')

#Rearranged and renamed the table "control.csv" and "cov.csv" into 4 diff tables with only 2 columns as below
#fcg.csv
#sample expression
#Covid_FCGR3B	1.44939523
#Control_FCGR3B	2.289875474
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

