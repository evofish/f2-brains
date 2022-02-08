setwd('~/Desktop/f2-tissues/deseq/')
library('DESeq2')
library('pheatmap')
library('RColorBrewer')
library('ggplot2')

x<-read.csv('f2_tissues_edit-dec20.txt', row.names="Geneid", header=T, sep='\t')
x$A2_B=NULL
#coldata<-read.table("f2-conditions.txt", header=T,row.names=1)
coldata<-read.table("F2-brain-traits-Reviewer2-Jan2022.txt", header=T,row.names=1) #this is for evaluating the effect of grandparents
head(x)
head(coldata)
coldata <- coldata[-c(2), ]#to remove sample A2 from traits table from the original
coldata

#remove oultier from counts and coldata
x<-as.matrix(x)

#basic count stats
totalCounts=colSums(x) 
totalCounts
min(totalCounts) #29036187
max(totalCounts) #56488723
mean(totalCounts) #41477978

all(rownames(coldata) %in% colnames(x)) #to check that coldata has all the correct names

#Heatmap of all samples to check outliers
#differences between all treatments
dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design = ~condition)
dds
featureData <- data.frame(gene=rownames(x))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
dds <- DESeq(dds)
res <- results(dds)
res
table(res$padj<0.05)
#18591  9344 

head(res)
vsd=getVarianceStabilizedData(dds)
head(vsd)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval", "padj")
vsdpvals=cbind(vsd,vals)
head(vsdpvals)
write.csv(vsdpvals,"F2-brain-LRT_VSDandPVALS_2020.csv", quote=F)

#sample to sample matrix 
sampleDists <- dist(t(vsd))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(vsd))
colnames(sampleDistMatrix) <- paste(colnames(vsd))
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(250)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
#sample A2_B looks like an outlier 

####
#LTR for Temp
dds.tg <- DESeqDataSetFromMatrix(countData=x, colData = coldata, design= ~stage+generation+temp+sex)
dds.tg2 <- DESeq(dds.tg, test="LRT", reduced=~generation+stage+sex)
res.tg <-results(dds.tg2)
write.csv(res.tg, file="Temp-DEG-LRT-table4GO.csv", quote=FALSE) #genes that are different just due to temperature
table(res.tg$padj<0.05)
#FALSE  TRUE
#21107  3691 effect of tem (not considering the effect of grandparents)

#test for grandparents Jan 2022
#temp
dds.tg <- DESeqDataSetFromMatrix(countData=x, colData = coldata, design= ~stage+generation+temp+sex+grandparents)
dds.tg2 <- DESeq(dds.tg, test="LRT", reduced=~generation+stage+sex+grandparents)
res.tg <-results(dds.tg2)
#write.csv(res.tg, file="Family-DEG-LRT-table4GO.csv", quote=FALSE) #genes that are different just due to temperature
table(res.tg$padj<0.05)
#FALSE  TRUE 
#21359   362 effect of temp considering grandparents 

#grandparents
dds.tg <- DESeqDataSetFromMatrix(countData=x, colData = coldata, design= ~stage+generation+temp+sex+grandparents)
dds.tg2 <- DESeq(dds.tg, test="LRT", reduced=~stage+generation+temp+sex)
res.tg <-results(dds.tg2)
write.csv(res.tg, file="grandparent-DEG-LRT-table4GO-jan2022.csv", quote=FALSE) #genes that are different just due to grandparent
table(res.tg$padj<0.05)
#FALSE  TRUE 
#28336   268

dds.tg <- DESeqDataSetFromMatrix(countData=x, colData = coldata, design= ~stage+generation+temp+sex+grandparents)
dds.tg3 <- DESeq(dds.tg, test="LRT", reduced=~temp+stage+sex+grandparents)
res.tg <-results(dds.tg3)
write.csv(res.tg, file="Generation-DEG-LRT-table4GO.csv", quote=FALSE) #genes that are different just due to generation
table(res.tg$padj<0.05)
#FALSE  TRUE 
#14821 13699  #this means that the differences between Generation are huge, which means they should be analyzed separately 
#16918 11060 this is when grandparents are added to the model

sigs<-res.tg[which(res.tg$padj<0.05),] #this command is only to write the table, use sig=which(res$padj<0.05) for heatmap
#write.csv(sigs, file="generation-DEGs-LRT-may2021.csv")
write.csv(sigs, file="generation-DEGs-LRT-jan2022.csv")

#test of just generation without considering other variables 
dds.tg <- DESeqDataSetFromMatrix(countData=x, colData = coldata, design= ~generation)
dds.tg33 <- DESeq(dds.tg)
res.tg <-results(dds.tg33)
#write.csv(res.tg, file="Generation-DEG-LRT-table4GO.csv", quote=FALSE) #genes that are different just due to generation
table(res.tg$padj<0.05)
#FALSE  TRUE 
#13595 15665  #the differences between Generation are still huge, which means they should be analyzed separately 


dds.tg4 <- DESeq(dds.tg, test="LRT", reduced=~temp+generation+sex+grandparents)
res.tg <-results(dds.tg4)
write.csv(res.tg, file="Stage-DEG-LRT-table4GO.csv", quote=FALSE) #genes that are different just due to stage at which stressor is applied
table(res.tg$padj<0.05)
#FALSE  TRUE 
#32240    37 #the stage at which the stressor was applied, Step was included as developmental treatment
#27971     7 once grandparents are included

####Effect of Sex 
dds.tg5 <- DESeq(dds.tg, test="LRT", reduced=~temp+generation+stage+grandparents)
res.tg <-results(dds.tg5)
write.csv(res.tg, file="Sex-DEG-LRT-table4GO.csv", quote=FALSE) #genes that are different just due to stage at which stressor is applied
table(res.tg$padj<0.05)
###False TRUE
##27272    8

############PCA for all: 
library(vegan)
library(rgl)
library(ape) 
library(RColorBrewer)

dds<- DESeqDataSetFromMatrix(countData=x, colData=coldata, design=~condition)
dds.v<-DESeq(dds)
summary(dds.v) # 29681 
res.v<-results(dds.v)
res.v
res.v<-res.v[order(res.v$padj),]
table(res.v$padj<0.05)
#FALSE  TRUE 
#18591  9344

head(res.v)
vsd.v=getVarianceStabilizedData(dds.v)
head(vsd.v)
vals.v=cbind(res.v$pvalue, res.v$padj)
head(vals.v)
colnames(vals.v)=c("pval.v", "padj.v")
vsdpvals.v=cbind(vsd.v,vals.v)
head(vsdpvals.v)
write.csv(vsdpvals.v, "F2-Brains-ConditionsVSDandPVALS_May21.csv", quote=F)

dat<-read.csv("F2-Brains-ConditionsVSDandPVALS_May21.csv", header=T)
head(dat)
colnames(dat)
data=dat1[,2:48]
row.names(data)=dat1$X
head(data)

# extracting experimental conditions 

traits=read.table("F2-brain-treatments-pca.txt",header=T, sep='\t') #this was not used
head(traits)
colorete <-
  with(traits,
       data.frame(condition = levels(condition),
                  color = I(brewer.pal(nlevels(condition), name = 'Paired'))))    
colores<-merge(traits, colorete)
colores

# principal coordinate calculation
dd.pcoa=pcoa(dist(t(data)))
scores=dd.pcoa$vectors
dd.pcoa 

quartz() 
plot(scores[,1], scores[,2],col=colorete$color, pch=19, cex=2, main= "PCoA for All Samples", xlab='PC1=39.69%', ylab='PC2=7.59%')
legend("bottomright", c("A_F1", "B_F1", "C_F1", "A_F2", "B_F2", "C_F2", "sC", "ttB", "ttC"), col=c("#A6CEE3", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#1F78B4"), pch=c(19),cex=0.75, horiz=FALSE)
#this PCA is a mess doesn't say anything 

###### PCA FOR F1 - based on LRT for F1
#for this, the DDS was done just for the 3 categories of the F1, to see the number DEGs 
x<-read.csv('f2_tissues_edit-dec20.txt', row.names="Geneid", header=T, sep='\t')
x$A2_B=NULL
head(x) #this is the main counts file 
f1x<-x[c(6:11,17:22,29:34)]
colnames(f1x)
coldata<-read.table("f2-tissues-conditions-F1.txt", header=T,row.names=1)
coldata
all(rownames(coldata) %in% colnames(f1x)) #check names in both documents 

dds<- DESeqDataSetFromMatrix(countData=f1x, colData=coldata, design=~condition)
dds.f1<-DESeq(dds)
summary(dds.f1) # 33604
res.f1<-results(dds.f1)
res.f1
res.f1<-res.f1[order(res.f1$padj),]
table(res.f1$padj<0.05)
#FALSE  TRUE 
#27229   110 #DEGs between all treatments of the F1 together 

vsd.f1=getVarianceStabilizedData(dds.f1)
head(vsd.f1)
vals.f1=cbind(res.f1$pvalue, res.f1$padj)
head(vals.f1)
colnames(vals.f1)=c("pval.f1", "padj.f1")
vsdpvals.f1=cbind(vsd.f1,vals.f1)
head(vsdpvals.f1)
write.csv(vsdpvals.f1, "F2-VSDandPVALS_F1-May21.csv", quote=F)

dat<-read.csv("F2-VSDandPVALS_F1-May21.csv", header=T)
head(dat )
data=dat[,2:19]
row.names(data)=dat$X
head(data)

# principal coordinate calculation
colorete <-
  with(coldata,
       data.frame(condition = levels(condition),
                  color = I(brewer.pal(nlevels(condition), name = 'Paired'))))    
colores<-merge(coldata, colorete)
colores

dd.pcoa=pcoa(dist(t(data)))
scores=dd.pcoa$vectors
dd.pcoa 
quartz() 
plot(scores[,1], scores[,2],col=colores$color, pch=19, cex=2, main= "PCoA for F1", xlab='PC1=30.19%', ylab='PC2=16.61%')
legend("bottomright", c("A_F1", "B_F1", "C_F1"), col=c('#A6CEE3', '#1F78B4', '#B2DF8A'), pch=c(19),cex=0.75, horiz=FALSE)


###### PCA FOR F2 - based on LRT for F2
#for this, the DDS was done just for the 6 categories of the F2, to see the number DEGs  
x<-read.csv('f2_tissues_edit-dec20.txt', row.names="Geneid", header=T, sep='\t')
x$A2_B=NULL
colnames(x) #this is the main counts file 

f2x<-x[c(1:5,12:16,23:28,35:47)]
colnames(f2x)
coldata<-read.table("f2-conditions-pca-f2.txt", header=T,row.names=1)
coldata
all(rownames(coldata) %in% colnames(f2x)) #check names in both documents 

dds<- DESeqDataSetFromMatrix(countData=f2x, colData=coldata, design=~condition)
dds.f2<-DESeq(dds)
summary(dds.f2) # 33604
res.f2<-results(dds.f2)
res.f2
res.f2<-res.f2[order(res.f2$padj),]
table(res.f2$padj<0.05)
#FALSE  TRUE 
#23091   188 #DEGs between all treatments of the F2 together 

vsd.f2=getVarianceStabilizedData(dds.f2)
head(vsd.f2)
vals.f2=cbind(res.f2$pvalue, res.f2$padj)
head(vals.f2)
colnames(vals.f2)=c("pval.f2", "padj.f2")
vsdpvals.f2=cbind(vsd.f2,vals.f2)
head(vsdpvals.f2)
write.csv(vsdpvals.f2, "F2-VSDandPVALS_F2-May21.csv", quote=F)

dat<-read.csv("F2-VSDandPVALS_F2-May21.csv", header=T)
colnames(dat)
data=dat[,2:30]
row.names(data)=dat$X
head(data)

# principal coordinate calculation
colorete <-
     with(coldata,
         data.frame(condition = levels(condition),
                    color = I(brewer.pal(nlevels(condition), name = 'Paired'))))    
colores<-merge(coldata, colorete)

dd.pcoa=pcoa(dist(t(data)))
scores=dd.pcoa$vectors
dd.pcoa 
quartz() 
plot(scores[,1], scores[,2],col= colores$color , pch=19, cex=2, main= "PCoA for F2", xlab='PC1=39.69%', ylab='PC2=7.58%')
legend("bottomright", c("A_F2","B_F2", "C_F2", "sC_F2","ttB_F2", "ttC_F2"), col=c('#33A02C', '#A6CEE3', '#1F78B4', '#B2DF8A', '#FB9A99', '#E31A1C'), pch=c(19),cex=0.75, horiz=FALSE)
#legend was set manually based on "colores" frame

#### Pairewise contrasts done with the original DESEQ2 Command based on individual temperature "conditions" 
###A vs tA - F1 vs F2
res<-results(dds, contrast=c('condition', 'A', 'tA')) #here is where the two contrasting conditions get defined
write.csv(res, file="A-tA-F1vsF2-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.05) 
#FALSE  TRUE 
#19140  7550  
sigs<-res[which(res$padj<0.05),] #this command is only to write the table, use sig=which(res$padj<0.05) for heatmap
write.csv(sigs, file='A-tA-F1vsF2-sig.csv')

###A vs B - F1
res<-results(dds, contrast=c('condition', 'A', 'B')) #here is where the two contrasting conditions get defined
write.csv(res, file="A-B-F1-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.05) 
#FALSE  TRUE 
#21585   122    
sigs<-res[which(res$padj<0.05),] #this command is only to write the table, use sig=which(res$padj<0.05) for heatmap
write.csv(sigs, file='AvsB-F1-sig.csv')

###A vs C - F1
res<-results(dds, contrast=c('condition', 'A', 'C')) #here is where the two contrasting conditions get defined
write.csv(res, file="A-C-F1-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.05) 
#FALSE  TRUE 
#25819   249    
sigs<-res[which(res$padj<0.05),] #this command is only to write the table, use sig=which(res$padj<0.05) for heatmap
write.csv(sigs, file='AvsC-F1-sig.csv')

###B vs C - F1
res<-results(dds, contrast=c('condition', 'B', 'C')) #here is where the two contrasting conditions get defined
write.csv(res, file="B-C-F1-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.05) 
#FALSE  TRUE 
#26049    19    
sigs<-res[which(res$padj<0.05),] #this command is only to write the table, use sig=which(res$padj<0.05) for heatmap
write.csv(sigs, file='BvsC-F1-sig.csv')

####
##For F2
# tA vs dB - F2
res<-results(dds, contrast=c('condition', 'tA', 'dB')) #here is where the two contrasting conditions get defined
write.csv(res, file="tA-dB-F2-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.05) 
#FALSE  TRUE 
#20431    27   
sigs<-res[which(res$padj<0.05),] #this command is only to write the table, use sig=which(res$padj<0.05) for heatmap
write.csv(sigs, file='tAvsdB-F1-sig.csv')

# tA vs dC - F2
res<-results(dds, contrast=c('condition', 'tA', 'dC')) #here is where the two contrasting conditions get defined
write.csv(res, file="tA-dC-F2-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.05) 
#FALSE  TRUE 
#24014   187  
sigs<-res[which(res$padj<0.05),] #this command is only to write the table, use sig=which(res$padj<0.05) for heatmap
write.csv(sigs, file='tAvsdC-F1-sig.csv')

# dB vs dC - F2
res<-results(dds, contrast=c('condition', 'dB', 'dC')) #here is where the two contrasting conditions get defined
write.csv(res, file="dB-dC-F2-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.05) 
#FALSE  TRUE 
#27301    10     
sigs<-res[which(res$padj<0.05),] #this command is only to write the table, use sig=which(res$padj<0.05) for heatmap
write.csv(sigs, file='dBvsdC-F2-sig.csv')

# tA vs tB - F2
res<-results(dds, contrast=c('condition', 'tA', 'tB')) #here is where the two contrasting conditions get defined
write.csv(res, file="tA-tB-F2-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.05) 
#FALSE  TRUE 
#26640    50    
sigs<-res[which(res$padj<0.05),] #this command is only to write the table, use sig=which(res$padj<0.05) for heatmap
write.csv(sigs, file='tAvstB-F2-sig.csv')

#tA vs tC - F2
res<-results(dds, contrast=c('condition', 'tA', 'tC')) #here is where the two contrasting conditions get defined
write.csv(res, file="tA-tC-F2-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.05) 
#FALSE  TRUE 
#27199   112    
sigs<-res[which(res$padj<0.05),] #this command is only to write the table, use sig=which(res$padj<0.05) for heatmap
write.csv(sigs, file='tAvstC-F2-sig.csv')

#tC vs sC - F2
res<-results(dds, contrast=c('condition', 'tC', 'sC')) #here is where the two contrasting conditions get defined
write.csv(res, file="tC-sC-F2-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.05) 
#FALSE  TRUE 
#26688     2  
sigs<-res[which(res$padj<0.05),] #this command is only to write the table, use sig=which(res$padj<0.05) for heatmap
write.csv(sigs, file='tCvs-sC-F2-sig.csv')

#dB vs tB - F2
res<-results(dds, contrast=c('condition', 'dB', 'tB')) #here is where the two contrasting conditions get defined
write.csv(res, file="dB-tB-F2-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.05) 
#FALSE  TRUE 
#26030    38    
sigs<-res[which(res$padj<0.05),] #this command is only to write the table, use sig=which(res$padj<0.05) for heatmap
write.csv(sigs, file='dBvstB-F2-sig.csv')

#dC vs tC - F2
res<-results(dds, contrast=c('condition', 'dC', 'tC')) #here is where the two contrasting conditions get defined
write.csv(res, file="dC-tC-F2-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.05) 
#FALSE  TRUE 
#26037    31   
sigs<-res[which(res$padj<0.05),] #this command is only to write the table, use sig=which(res$padj<0.05) for heatmap
write.csv(sigs, file='dCvstC-F2-sig.csv')

#dC vs sC - F2
res<-results(dds, contrast=c('condition', 'dC', 'sC')) #here is where the two contrasting conditions get defined
write.csv(res, file="dC-sC-F2-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.05) 
#FALSE  TRUE 
#32317     2     
sigs<-res[which(res$padj<0.05),] #this command is only to write the table, use sig=which(res$padj<0.05) for heatmap
write.csv(sigs, file='dCvs-sC-F2-sig.csv')

######
#tB vs sC - F2
res<-results(dds, contrast=c('condition', 'tB', 'sC')) #here is where the two contrasting conditions get defined
write.csv(res, file="tB-sC-F2-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.05) 
#FALSE  TRUE 
#32317     149     
sigs<-res[which(res$padj<0.05),] #this command is only to write the table, use sig=which(res$padj<0.05) for heatmap
write.csv(sigs, file='tBvs-sC-F2-sig.csv')
#######

#dB vs sC - F2
res<-results(dds, contrast=c('condition', 'dB', 'sC')) #here is where the two contrasting conditions get defined
write.csv(res, file="dB-sC-F2-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.05) 
#FALSE  TRUE 
#25427    17     
sigs<-res[which(res$padj<0.05),] #this command is only to write the table, use sig=which(res$padj<0.05) for heatmap
write.csv(sigs, file='dBvs-sC-F2-sig.csv')

#tA vs sC - F2
res<-results(dds, contrast=c('condition', 'tA', 'sC')) #here is where the two contrasting conditions get defined
write.csv(res, file="tA-sC-F2-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.05) 
#FALSE  TRUE 
#26642    48     
sigs<-res[which(res$padj<0.05),] #this command is only to write the table, use sig=which(res$padj<0.05) for heatmap
write.csv(sigs, file='tA-vs-sC-F2-sig.csv')

#tB vs tC - F2
res<-results(dds, contrast=c('condition', 'tB', 'tC')) #here is where the two contrasting conditions get defined
write.csv(res, file="tB-tC-F2-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.05) 
#FALSE  TRUE 
#24461   361    
sigs<-res[which(res$padj<0.05),] #this command is only to write the table, use sig=which(res$padj<0.05) for heatmap
write.csv(sigs, file='tBvstC-F2-sig.csv')

######## 
##loop for annotation of DEGs of the pairwise comparisons
setwd('~/Desktop/f2-tissues/deseq/DEG')
w<-read.table(file='../../annotation/apoly-annot-columns-jul18.txt', header=T, sep='\t')
w <- w[-c(5), ]
head(w)

y<-read.csv(file='A-tA-F1vsF2-sig.csv', header=T)
out <- merge(y, w, by.y="x", all.x=TRUE)
write.table(out, file='A-tA-merged-F1F2.txt',sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

files<-list.files(pattern="*-sig.csv", full.names=TRUE, recursive=FALSE)
files

#make sure that the headers are the same

lapply(files, function(x) {
  y <- read.csv(x, header=TRUE, sep=',') # load file
  # apply function
  out <- merge(y, w, by.y="x",all.x=TRUE)
  # write to file
  write.table(out, file=paste(x,"merged.txt",sep="\t"), quote=FALSE, row.names=FALSE, col.names=TRUE)
}) #this generates the merged files, the name of the files is not great needs to be edited further
