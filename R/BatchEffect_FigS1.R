### plot FigS1 for 8 datasets by Qianyi on July 2020
#8 datasets:
#  4 time points: 8 dpi, 13 dpi, 16 dpi, 21 dpi
#  2 types: WT vs SCARKO (sertoli-cell Androgen-receptor knockout)

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(cowplot)
library(RColorBrewer)
redblue100<-rgb(read.table('../redblue100.txt',sep='\t',row.names=1,header=T))
gg_color_hue <- function(n) {
hues = seq(15, 375, length = n + 1)
hcl(h = hues, l = 65, c = 100)[1:n]
}
myBrewerPalette<-c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#f781bf")


### load Seurat object for all merged data corrected by batch effect using CCA
load(file="all.cells.cca_umap_dim1.3_res.0.3_no5_merg2.3_rename_101619_FINALOBJ_merg.Leydig_031320.Robj")
ls()
dge=all.cells.cca_2
dgeall=dge


### plot S1A: heatmap for rank correlation of cluster centroids across the 8 datasets with 4 time points and 2 genotypes
dpi8WT<-SubsetData(dge, subset.name = 'group', accept.value = 'dpi8WT', subset.raw = T)
dpi8SCAR<-SubsetData(dge, subset.name = 'group', accept.value = 'dpi8SCAR', subset.raw = T)
dpi13WT<-SubsetData(dge, subset.name = 'group', accept.value = 'dpi13WT', subset.raw = T)
dpi13SCAR<-SubsetData(dge, subset.name = 'group', accept.value = 'dpi13SCAR', subset.raw = T)
dpi16WT<-SubsetData(dge, subset.name = 'group', accept.value = 'dpi16WT', subset.raw = T)
dpi16SCAR<-SubsetData(dge, subset.name = 'group', accept.value = 'dpi16SCAR', subset.raw = T)
dpp21WT<-SubsetData(dge, subset.name = 'group', accept.value = 'dpp21WT', subset.raw = T)
dpp21SCAR<-SubsetData(dge, subset.name = 'group', accept.value = 'dpp21SCAR', subset.raw = T)

dpi8WT_centroids<-log(AverageExpression(dpi8WT)+1)
dpi8SCAR_centroids<-log(AverageExpression(dpi8SCAR)+1)
dpi13WT_centroids<-log(AverageExpression(dpi13WT)+1)
dpi13SCAR_centroids<-log(AverageExpression(dpi13SCAR)+1)
dpi16WT_centroids<-log(AverageExpression(dpi16WT)+1)
dpi16SCAR_centroids<-log(AverageExpression(dpi16SCAR)+1)
dpp21WT_centroids<-log(AverageExpression(dpp21WT)+1)
dpp21SCAR_centroids<-log(AverageExpression(dpp21SCAR)+1)

colnames(dpi8WT_centroids) <- paste("dpi8_WT", colnames(dpi8WT_centroids), sep = "_")
colnames(dpi8SCAR_centroids) <- paste("dpi8_SCARKO", colnames(dpi8SCAR_centroids), sep = "_")
colnames(dpi13WT_centroids) <- paste("dpi13_WT", colnames(dpi13WT_centroids), sep = "_")
colnames(dpi13SCAR_centroids) <- paste("dpi13_SCARKO", colnames(dpi13SCAR_centroids), sep = "_")
colnames(dpi16WT_centroids) <- paste("dpi16_WT", colnames(dpi16WT_centroids), sep = "_")
colnames(dpi16SCAR_centroids) <- paste("dpi16_SCARKO", colnames(dpi16SCAR_centroids), sep = "_")
colnames(dpp21WT_centroids) <- paste("dpp21_WT", colnames(dpp21WT_centroids), sep = "_")
colnames(dpp21SCAR_centroids) <- paste("dpp21_SCARKO", colnames(dpp21SCAR_centroids), sep = "_")

dall<-cbind(dpi8WT_centroids,dpi8SCAR_centroids,dpi13WT_centroids,dpi13SCAR_centroids,dpi16WT_centroids,dpi16SCAR_centroids,dpp21WT_centroids,dpp21SCAR_centroids)

all_cor <- cor(dall,  method = 'spearman')
write.table(all_cor,file="Global_batch_centroids_rho.txt",row.names=T,col.names=T,quote=F,sep="\t")

df = data.frame(strsplit(colnames(all_cor),split = '_'))
df <- t(df)
rownames(df)<-NULL
colnames(df) <- c('Time', 'Genotype', 'CellType')

levels=paste(df[,1],df[,2],sep="_")
colsep.use=cumsum(table(levels)[unique(levels)]))
row.lab=df[,3]
col.lab=rep("",length(levels))
col.lab[round(cumsum(table(levels)[unique(levels)])-table(levels)[unique(levels)]/2)+0]=gsub("_","",unique(levels))

n=length(levels)
ntime=length(unique(df[,1]))
ncluster=length(unique(df[,3]))
sidecol=matrix(0,3,length(levels))
timecol=c("white","gray75","gray50","black")
sidecol[1,]=rep(timecol,each=2*ncluster)
sidecol[2,]=rep(rep(rev(gg_color_hue(2)),each=ncluster),ntime)
sidecol[3,]=rep(myBrewerPalette[1:sum(ncluster)],ntime*2)
clab=cbind(sidecol[3,],sidecol[2,],sidecol[1,])
rlab=sidecol
rownames(rlab)=colnames(df)
colnames(clab)=rev(colnames(df))

col.use=redblue100
pdf(file=paste("Global_Batch_Centroid_RankedCorrelation.pdf",sep=""),height=6,width=6.5)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(all_cor,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=.8,cexRow=.8,ColSideColorsSize = 1.8,RowSideColorsSize = 1.8,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(5,7))
dev.off()
# save as FigS1A
