#define the format of the result table
tabledata= data.frame( pathwayid= character(), pathwayname= character(), Pvalue = double(), Oddsratio = double(), stringsAsFactors=FALSE)
#load required libraries
library(GeneOverlap)
library(cogena)
library(dplyr)
#load input
allGOnKEGG= gmt2list("allGOnKEGG_June2020.gmt")
brs=read.csv("sev_gene.txt",sep=',')
genelist= as.vector(brs$GeneID)

#loop and pairwase comparison between each pathway and the input genelist using geneoverlap (one sided FET)
for (i in 1:length(allGOnKEGG)) 
{
  pathway= as.data.frame(allGOnKEGG[i])[,1] 
  
  # Restrict analysis to pathways of size 50 to 1000 genes 
  if (length(pathway) <50  | length(pathway) > 1000) 
  {next()}
  

  
  #use geneoverlap to calculate P value and odds ratio
  go.obj <- newGeneOverlap(as.vector(genelist), pathway)
  go.obj <- testGeneOverlap(go.obj)
  # Restrict analysis to pathways intersection > 5
  intersection <- go.obj@intersection
  i_len <- length(intersection)
  if (i_len < 5)
  {next()}
  
  tabledata[i, 3] = go.obj@pval
  tabledata[i, 4] = go.obj@odds.ratio
  #Manipulation of string to get desired pathway name and pathway ID
  tabledata[i, 1] = strsplit(names(allGOnKEGG[i]), "%")[[1]][3] 
  tabledata[i, 2] = strsplit(names(allGOnKEGG[i]), "%")[[1]][1] 
  
}
#remove cells that have NA value (these are the pathways with size greater that 1000 genes and less than 50 genes) and reorder in ascending order of p value

tabledataALL= tabledata[which(tabledata$pathwayid!= "NA" ),]
tabledataALL= tabledataALL %>% arrange(tabledataALL$Pvalue)

#My implementation of FDR control and selecting FDR to be less than 0.01

newp = p.adjust(tabledataALL$Pvalue, method = "fdr", n = length(tabledataALL$Pvalue))
tabledataALL= tabledataALL %>%  mutate(FDR = newp)
tabledataALL = tabledataALL [which(tabledataALL$FDR< 0.01),]

#Save final result table "tabledataALL" to a csv file in the directory
write.csv(tabledataALL, "syed/out_severeGene.csv", row.names=FALSE)

