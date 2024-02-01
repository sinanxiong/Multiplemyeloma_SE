# ---
# title: "RNA-seq differential expression"
# author: "Tuan Zea Tan"
# last update: "Jan 23, 2024"
# ---

library(EBSeq)
print('usage: goEBseq(infile,outfile,conditions)')
print('y = as.factor(rep(c("C1","C2","C3"),each=3))')
print('can use for case where one group has only one sample')
print('can use for paired analysis by making one group one sample')


goEBseq <- function(inFile,outFile,y){
x = read.table(inFile,sep='\t',header=T,row.names=1)
xx = as.matrix(x)

sx = MedianNorm(xx) #normalization size factor


EBOut=EBTest(Data=xx,Conditions=y,sizeFactors=sx, maxround=5)

GeneFC=PostFC(EBOut)
print(GeneFC$Direction)
print(EBOut$ConditionOrder)
#PPMat, 1st column = PPEE, second = PPDE
fx = data.frame(EBOut$PPMat,GeneFC$PostFC,GeneFC$RealFC,unlist(EBOut$C1Mean),unlist(EBOut$C2Mean),check.rows=FALSE)


write.table(fx,file=outFile, sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE)

}