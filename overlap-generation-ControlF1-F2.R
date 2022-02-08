#simple script to see the overlap between F1 and F2 generations and the comparsion between F1 and F2 control
x<-read.csv('~/Desktop/f2-tissues/deseq/generation-DEGs-LRT-jan2022.csv', header=TRUE)
head(x)
nrow(x)
y<-read.csv('~/Desktop/f2-tissues/deseq/DEG/A-tA-F1vsF2-sig.csv', header=TRUE)
head(y)
nrow(y)

xy<-merge(x,y, by='X')
nrow(xy)
#Overlap of 6326 