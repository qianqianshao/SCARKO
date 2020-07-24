### Differentially-expressed markers between SCARKO and WT for germ cell types by Qianyi on 5/24/2020
# directly merged germ cells from 8 datasets, did focused germ cell subset clustering -> 7 germ cell types

#1 species: mouse
#8 datasets:
#  4 time points: 8 dpi, 13 dpi, 16 dpi, 21 dpi
#  2 types: WT vs SCARKO (sertoli-cell Androgen-receptor knockout)


library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(cowplot)
library(RColorBrewer)
myBrewerPalette=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2")[c(4,8,1)],brewer.pal(8,"Set2")[c(4,8,1)])
redblue100<-rgb(read.table('../redblue100.txt',sep='\t',row.names=1,header=T))

load(file="GermCells_DirectMergeAll_May2020.Robj") 
### add metadata column of WT vs SCARKO
typelabel=paste(dge@meta.data$type,as.character(ident),sep="_")
names(typelabel)=names(dge@ident)
levels=NULL
for(label in levels(ident)){
for(type in c("SCAR","WT")){
  levels=c(levels,paste(type,label,sep="_"))
}
}
typelabel=factor(typelabel,levels=levels)
dge=AddMetaData(dge,typelabel,"typelabel")

### visualize cell size factors for 7 germ cell types, 4 time points and/or 2 genotypes separately
plotlist=list()
pdf(file=paste0("GC_SCARKOvsWT_PerCellAttributes_ViolinPlot0.pdf"),height=14,width=14)
dge=SetAllIdent(dge,id="typelabel")
dge@ident=typelabel
plotlist[[1]]=VlnPlot(dge, features.plot="nGene", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
plotlist[[2]]=VlnPlot(dge, features.plot="nUMI", y.log=T, nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
plotlist[[3]]=VlnPlot(dge, features.plot="percent.mito", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
multiplot(plotlist,cols = 1)
dev.off()

pdf(file=paste0("GC_SCARKOvsWT_PerCellAttributes_ViolinPlot1.pdf"),height=12,width=7.2)
dge=SetAllIdent(dge,id="celltype")
dge@ident=ident
plotlist[[1]]=VlnPlot(dge, features.plot="nGene", nCol = 1,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
plotlist[[2]]=VlnPlot(dge, features.plot="nUMI", y.log=T, nCol = 1,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
plotlist[[3]]=VlnPlot(dge, features.plot="percent.mito", nCol = 1,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
multiplot(plotlist,cols = 1)
dev.off()

pdf(file=paste0("GC_SCARKOvsWT_PerCellAttributes_ViolinPlot2.pdf"),height=8,width=7.5)
  plotlist=list()
dge=SetAllIdent(dge,id="orig.ident")
dge@ident=factor(dge@ident,levels=c("dpi8SCAR","dpi8WT","dpi13SCAR","dpi13WT","dpi16SCAR","dpi16WT","dpp21SCAR","dpp21WT"))
plotlist[[1]]=VlnPlot(dge, features.plot="nGene", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
plotlist[[2]]=VlnPlot(dge, features.plot="nUMI", y.log=T, nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
plotlist[[3]]=VlnPlot(dge, features.plot="percent.mito", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
multiplot(plotlist,cols = 1)
dev.off()

pdf(file=paste0("GC_SCARKOvsWT_PerCellAttributes_ViolinPlot3.pdf"),height=8,width=3)
dge=SetAllIdent(dge,id="type")
plotlist[[1]]=VlnPlot(dge, features.plot="nGene", nCol = 1,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
plotlist[[2]]=VlnPlot(dge, features.plot="nUMI", y.log=T, nCol = 1,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
plotlist[[3]]=VlnPlot(dge, features.plot="percent.mito", nCol = 1,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
multiplot(plotlist,cols = 1)
dev.off()


### Differential expression between SCARKO and WT for each of 5 germ cell types
### p-value and fold change for all genes
dge=SetAllIdent(dge,id="celltype")
dge@ident=ident
dgeall=dge


markerslist=list()
for(i in 1:5){
  label=levels(ident)[i]
  dge=SubsetData(dgeall,ident.use=label)
  dge=SetAllIdent(dge,id="type")
### Differential expression for all genes and volcano plot
  gene=FindMarkers(dge,ident.1="SCAR", only.pos=FALSE,test.use="bimod",min.pct=0.1,min.diff.pct=0.2,logfc.threshold=log(2))
write.table(gene,paste0("GC_SCARKOvsWT_C",i,"_de_minpct0.1_mindiff0.2_log2_05.2020.txt"),col.names=T,row.names=T,quote=F,sep="\t")
print(c(length(which(gene$avg_logFC>0)),length(which(gene$avg_logFC<0))))
gene1=rownames(gene[which(gene$avg_logFC>0),])
gene2=rownames(gene[which(gene$avg_logFC<0),])
genes1=c(gene1,gene2)
  gene=FindMarkers(dge,ident.1="SCAR", only.pos=FALSE,test.use="bimod",min.pct=0,min.diff.pct=0,logfc.threshold=log(1.5))
write.table(gene,paste0("GC_SCARKOvsWT_C",i,"_de_log1.5_05.2020.txt"),col.names=T,row.names=T,quote=F,sep="\t")
print(c(length(which(gene$avg_logFC>0)),length(which(gene$avg_logFC<0))))
gene1=rownames(gene[which(gene$avg_logFC>0),])
gene2=rownames(gene[which(gene$avg_logFC<0),])
genes2=c(gene1,gene2)
genes2=genes2[which(!(genes2 %in% genes1))]
  markers1=FindMarkers(dge,ident.1="SCAR", only.pos=FALSE,test.use="bimod",min.pct=-Inf,min.diff.pct=-Inf,logfc.threshold=-Inf,min.cells.gene = -Inf,min.cells.group = -Inf)
a=markers1
jpeg(file=paste0("GC_SCARKOvsWT_C",i,"_Volcano.jpeg"),res=300,height=1600,width=1600)
par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(a$avg_logFC,-log10(a$p_val),pch=16,cex=0.8,col=rgb(0,0,0,0.3),xlab="log Fold Change",ylab="-log10(p-value)")
if(length(genes1)!=0){
points(a[genes1,]$avg_logFC,-log10(a[genes1,]$p_val),cex=0.8,pch=16,col=rgb(1,0,0,0.8))
}
if(length(genes2)!=0){
points(a[genes2,]$avg_logFC,-log10(a[genes2,]$p_val),cex=0.8,pch=16,col=rgb(1,0,1,0.8))
}
plot((a$pct.1-a$pct.2)*100,-log10(a$p_val),pch=16,cex=0.8,col=rgb(0,0,0,0.3),xlab="%Cells Diff",ylab="-log10(p-value)")
if(length(genes1)!=0){
points((a[genes1,]$pct.1-a[genes1,]$pct.2)*100,-log10(a[genes1,]$p_val),pch=16,cex=0.8,col=rgb(1,0,0,0.8))
}
if(length(genes2)!=0){
points((a[genes2,]$pct.1-a[genes2,]$pct.2)*100,-log10(a[genes2,]$p_val),pch=16,cex=0.8,col=rgb(1,0,1,0.8))
}
plot(a$avg_logFC,(a$pct.1-a$pct.2)*100,pch=16,cex=0.8,col=rgb(0,0,0,0.3),xlab="log Fold Change",ylab="%Cells Diff")
if(length(genes1)!=0){
points(a[genes1,]$avg_logFC,(a[genes1,]$pct.1-a[genes1,]$pct.2)*100,pch=16,cex=0.8,col=rgb(1,0,0,0.8))
}
if(length(genes2)!=0){
points(a[genes2,]$avg_logFC,(a[genes2,]$pct.1-a[genes2,]$pct.2)*100,pch=16,cex=0.8,col=rgb(1,0,1,0.8))
}
dev.off()
### calculate average expression for each gene
centroid=log(AverageExpression(dge)+1)
print(which(!(rownames(centroid) %in% rownames(markers1)))) # integer(0)
centroid=data.frame(gene=rownames(centroid),centroid)
markers2=data.frame(gene=rownames(markers1),markers1)
markers=merge(markers2,centroid,by="gene",all=TRUE)
rownames(markers)=markers$gene
markers$diff=markers$pct.1-markers$pct.2
markers1=markers
write.table(markers1,paste0("GC_SCARKOvsWT_C",i,"_allgenes_de_centroid_05.2020.txt"),col.names=T,row.names=T,quote=F,sep="\t")
markerslist[[i]]=markers1
}

# saved the markers for C3 as Table 1.
