### plot Fig4A-C, FigS1 and FigS2 and generated Table S2 by Qianyi on July 2020
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
myBrewerPalette1=gg_color_hue(4)
myBrewerPalette<-c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#f781bf")


### load Seurat object for all merged data corrected by batch effect using CCA
load(file="all.cells.cca_umap_dim1.3_res.0.3_no5_merg2.3_rename_101619_FINALOBJ_merg.Leydig_031320.Robj")
ls()
dge=all.cells.cca_2
dgeall=dge

table(dge@meta.data$group)
summary(dge@meta.data$nUMI)
summary(dge@meta.data$nGene)
table(dge@meta.data$group,dge@ident)[unique(dge@meta.data$group),]



### plot Fig 4A: visualization of 7 major cell types in global view 
tiff(paste0("Global_clusters.tiff"),res=300,height=1500,width=2000)
DimPlot(dgeall,pt.size=0.4,reduction.use="umap",no.legend=FALSE,cols.use=myBrewerPalette,do.return=TRUE,do.label=FALSE)
dev.off()
### check if tSNE view seperate germ cells from leydig better
tiff(paste0("Global_clusters_tSNE.tiff"),res=300,height=1500,width=2000)
TSNEPlot(dgeall,pt.size=0.4,no.legend=FALSE,colors.use=myBrewerPalette,do.return=TRUE,do.label=FALSE)
dev.off()


### plot Fig 4B: separated two genotypes (WT and SCARKO) with the other genotypes colored in light grey in the background
dge=SetAllIdent(dge,"stim")
tiff(paste0("Global_SCARKOvsWT_tSNE.tiff"),res=300,height=1500,width=1800)
TSNEPlot(dge,pt.size=0.4,no.legend=FALSE,do.return=TRUE,do.label=FALSE)
dev.off()
# Fig 4B: make background grey lighter in order to distinguish 
library(scales)
dim=dge@dr$umap@cell.embeddings;pos="bottomleft"
dim=dge@dr$tsne@cell.embeddings;pos="topleft"
xlim=range(dim[,1])
ylim=range(dim[,2])
sets=levels(dge@ident)
cols=gg_color_hue(length(sets))
size=length(sets)
  pdf(paste0("Global_SCARKOvsWT_indiv.pdf"),height=3.5,width=3.5*size)
  par(mfrow=c(1,size),mar=c(2.5,2.5,0.5,0.5),mgp=c(1.2, 0.5, 0))
for(seti in 1:length(sets)){
set=sets[seti]
ident=as.character(dge@ident)
names(ident)=names(dge@ident)
ident[which(!(ident %in% set))] <- "Others"
ident=factor(ident,levels=c("Others",set))
ident=sort(ident)
tmp=dim[names(ident),]
plot(tmp,pch=16,cex=0.3,col=c("grey80",alpha(cols[seti],0.8))[ident],xlim=xlim,ylim=ylim)
legend(pos,pch=16,set,col=cols[seti])
}
dev.off()


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

myBrewerPalette<-c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#f781bf")

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

library(gplots)
library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

pdf(file=paste("Global_Batch_Centroid_RankedCorrelation.pdf",sep=""),height=6,width=6.5)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(all_cor,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=.8,cexRow=.8,ColSideColorsSize = 1.8,RowSideColorsSize = 1.8,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(5,7))
dev.off()


### plot Fig 4C: Visualize standardized expression of markers for each of the 7 cell types across 7 cell type centroids in heatmap 
markers=FindAllMarkers(dge,only.pos=TRUE,min.pct=0.2,min.diff.pct=0.2,logfc.threshold=log(2),test.use="bimod",do.print = TRUE)
write.table(markers, file="markers_all.cells.cca.txt")
# save as Table S2 Sheet 1: markers for 7 major cell types
centroid.std=read.table(file="AllCellsCCA_centroid_allgenes_std_032320.txt",header=T,row.names=1,stringsAsFactors=F,sep="\t")
table(dge@ident)
table(markers$cluster)
data.use=centroid.std[markers$gene,]
write.table(data.use,"AllCellsCCA_centroid_std_markers.txt",row.names=T,col.names=T,quote=F,sep="\t")
# save as TableS2 3rd sheet
genes=markers$gene
data.use=centroid.std[markers$gene,]
colnames(centroid.std)=unique(markers$cluster)
levels=colnames(centroid.std)
colsep.use=cumsum(table(levels)[levels])
col.lab=rep("",length(levels))
col.lab=levels
row.lab=NULL
ncluster=length(levels)
sidecol=matrix(0,2,length(levels))
sidecol[1,]=rep(rep(c("white","white"),each=12),3)[1:sum(ncluster)]
sidecol[2,]=myBrewerPalette[1:sum(ncluster)]
clab=cbind(sidecol[2,],sidecol[1,])
rlab=sidecol
rownames(rlab)=c("","Cell Type")
colnames(clab)=c("Cell Type","")

redblue100<-rgb(read.table(paste0(home,"data_DGE/redblue100.txt"),sep='\t',row.names=1,header=T))
col.use=redblue100

library(gplots)
library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

jpeg(file=paste0("Fig4C_centroid_std_markersall.jpeg"),res=300,height=2600,width=1600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=0.3,ColSideColorsSize = 2,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()
jpeg(file=paste0("Fig4C_centroid_std_markersall2.jpeg"),res=300,height=1800,width=1600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=0.3,ColSideColorsSize = 2,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()


### plot Fig S2: visualize known markers used to classify each cell type
knownmarkers=c("Stra8","Dazl","Mki67","Piwil1","Zbtb16","Gfra1","Id4","Sohlh1","Sohlh2",
"Syce1","Setx","Sycp2","Spag6","Sycp3",
"Amhr2","Sox9","Rhox8","Clu","Wt1","Dmrt1","Inha","Tubb3",
"Tcf21","Dcn","Pdgfra","Nr2f2","Ly6a",
"Lipg","Prss35","Tpm1","Smoc2",
"Col1a2","Lama2","S100b","Sct",
"Cyp11a1","Cyp17a1","Star","Lhcgr","Hsd3b1",
"Cd68","Csf1r","Lyz1","Adgre1"  )
length(knownmarkers) # 44
knownmarkers[which(!(knownmarkers %in% rownames(dge@data)))] 
knownmarkers=knownmarkers[which(knownmarkers %in% rownames(dge@data))]
length(knownmarkers) # 44

pdf("knownmarkers_Violin.pdf",height=12,width=24)
VlnPlot(dge,knownmarkers,cols.use=myBrewerPalette,nCol=8,point.size.use=-1)
dev.off()
# save as Fig S2

### Finding local differentially-expressed markers between ImmLeydig and Leydig
gene=FindMarkers(dge,"Immature Leydig", "Adult Leydig",test.use="bimod",logfc.threshold = log(2),min.pct=0.1,min.diff.pct=0.2)
print(c(length(which(gene$avg_logFC>0)),length(which(gene$avg_logFC<0))))
# 5 1
gene=FindMarkers(dge,"Immature Leydig", "Adult Leydig",test.use="bimod",logfc.threshold = log(2),min.pct=0.1)
print(c(length(which(gene$avg_logFC>0)),length(which(gene$avg_logFC<0))))
# 9 3
gene=FindMarkers(dge,"Immature Leydig", "Adult Leydig",test.use="bimod",logfc.threshold = log(1.5),min.pct=0.1)
print(c(length(which(gene$avg_logFC>0)),length(which(gene$avg_logFC<0))))
# 42 19 
write.table(gene,paste0("ImmVsAdultLeydig_de_minpct0.1_log1.5_07.2020.txt"),col.names=T,row.names=T,quote=F,sep="\t")
markers=gene
markers$gene=rownames(markers)
markersp=markers[markers$avg_logFC>0,]
markersn=markers[markers$avg_logFC<0,]
markers=rbind(markersp,markersn)
pdf("markersImmVsAdultLeydig_Violin.pdf",height=1.5*round(sqrt(length(markers$gene))),width=1.7*ceiling(sqrt(length(markers$gene))))
VlnPlot(dge,markers$gene,cols.use=myBrewerPalette,nCol=ceiling(sqrt(length(markers$gene))),point.size.use=-1)
dev.off()
markers[markers$avg_logFC>0,] %>% top_n(8, avg_logFC)  -> topp
markers[markers$avg_logFC<0,] %>% top_n(8, -avg_logFC)  -> topn
markers.use=c(topp$gene,topn$gene)
size=sqrt(length(markers.use))
jpeg("markersImmVsAdultLeydig_Feature.jpeg",res=300,height=500*round(size),width=500*ceiling(size))
FeaturePlot(object = dge, reduction.use="umap",features.plot = markers.use,col=c("lightblue","red"),nCol=ceiling(size))
dev.off()

### check cell size factor between Immature Leydig and Adult Leydig
pdf(file=paste0("CellSizeFactors_ViolinPlot.pdf"),height=3.5,width=7)
	plotlist=list()
plotlist[[1]]=VlnPlot(dgeall, features.plot="nGene", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
plotlist[[2]]=VlnPlot(dgeall, features.plot="nUMI", y.log=T, nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
multiplot(plotlist,cols = 2)
dev.off()


