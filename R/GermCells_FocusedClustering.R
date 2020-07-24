### Merge and do clustering for only germ cells on 6.24.2019 by Qianyi
# directly merge germ cells from 8 datasets, do focused germ cell subset clustering 

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

### load raw data for germ cells from each of the 8 datasets
label=c("dpi8WT","dpi8SCARKO","dpi13WT","dpi13SCARKO","dpi16WT","dpi16SCARKO","dpp21WT","dpp21SCARKO")
rawgc=list()
for(i in 1:length(label)){
  tmp=read.table(paste0(label[i],"_GCs_062419.txt"),header=T,row.names=1) # data sent to me by Hailey at 11:30 AM
  nCellperGene <- rowSums(tmp>0)
  genes.use=names(nCellperGene[which(nCellperGene!=0)])
  rawgc[[i]]=tmp[genes.use,]
  print(c(dim(tmp),dim(rawgc[[i]]))))
}
label=c("dpi8WT","dpi8SCAR","dpi13WT","dpi13SCAR","dpi16WT","dpi16SCAR","dpp21WT","dpp21SCAR")
time=c("dpi8","dpi13","dpi16","dpp21")
type=c("WT","SCAR")
### merge all germ cells from each individual dataset
set=list()
for(i in 1:length(label)){
  set[[i]]=data.frame(GENE=rownames(rawgc[[i]]),rawgc[[i]])
}
dgedata=Reduce(function(x,y) merge(x,y,all=TRUE), set)
dgedata[is.na(dgedata)] <- 0
row.names(dgedata)=dgedata[,1]
dgedata=dgedata[,-1]
dim(dgedata)                   #[1] 32943  3351
nCellperGene <- rowSums(dgedata>0)
length(which(nCellperGene==0))     #0
nCellperGene[which(nCellperGene==0)] #numeric(0)


### create new Seurat object
dge <- CreateSeuratObject(raw.data = dgedata, project="GC_DirectMergeAll_Jun2019", min.cells=1, min.genes=1)
mito.genes <- grep(pattern = "^MT-|^mt-", x = rownames(x = dge@data), value = TRUE)
percent.mito <- Matrix::colSums(dge@raw.data[mito.genes, ])/Matrix::colSums(dge@raw.data)
dge <- AddMetaData(object = dge, metadata = percent.mito, col.name = "percent.mito")
### Normalize data
dge <- NormalizeData(object = dge)
### Highly variable genes
pdf("GC_DirectMergeAll_HVG0.2.pdf",height=7.5,width=11)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
dge <- FindVariableGenes(object = dge, x.low.cutoff = 0.2, x.high.cutoff = 10, y.cutoff = 0.2, do.plot = TRUE, col.use=rgb(0,0,0,0.8))
dev.off()
print(length(dge@var.genes)) # 2988
### Scale data 
dge <- ScaleData(object = dge)
dge                    # 32943 genes across 3351 samples          
### Add a column of time points in meta data sheet
ident<-gsub("WT","",gsub("SCAR","",as.character(dge@ident)))
names(ident)=names(dge@ident)
ident=factor(ident,levels=time)
table(ident)
# dpi8 dpi13 dpi16 dpp21 
#  216   518  1242  1375
dge=AddMetaData(dge,ident,"time")
### Add a column of wild type in meta data sheet
ident<-gsub("dpi8","",gsub("dpi13","",gsub("dpi16","",gsub("dpp21","",as.character(dge@ident))) ))
names(ident)=names(dge@ident)
ident=factor(ident,levels=type)
table(ident)
#  WT SCAR 
#2389  962
dge=AddMetaData(dge,ident,"type")

save(dge,file="GermCells_DirectMergeAll_Jun2019.Robj")

### PCA
print(Sys.time())    # [1] "2019-06-24 12:02:46 EDT"
dge <- RunPCA(dge, pc.genes = dge@var.genes,pcs.compute = 40,  do.print = TRUE, pcs.print = 5, genes.print = 5)
dge <- ProjectPCA(object = dge, do.print = FALSE)

### Scree Plot for PCA
numPCs=8;i=1 # HVG
pdf("dge_PCA_Variablel_variation.pdf",height=4.5,width=4.5)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(dge@dr$pca@sdev,type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
abline(v=numPCs[i]+0.5,col="red")
dev.off()

### Louvain-Jaccard clustering, tSNE, UMAP using top PCs
dge <- FindClusters(dge, reduction.type = "pca", dims.use = 1:numPCs[i], resolution = seq(0.1,3,by=0.1), save.SNN = T)
dge <- RunTSNE(dge, dims.use = 1:numPCs[i], do.fast = T)
dge <- RunUMAP(dge, reduction.use = "pca", dims.use = 1:numPCs[i])
save(dge,file="GermCells_DirectMergeAll_Jun2019.Robj")
Sys.time() # [1] "2019-06-24 12:09:50 EDT"


### check the number of clusters
print(c( length(unique(dge@meta.data$res.0.1)),length(unique(dge@meta.data$res.0.2)),length(unique(dge@meta.data$res.0.3)),length(unique(dge@meta.data$res.0.4)),length(unique(dge@meta.data$res.0.5)),length(unique(dge@meta.data$res.0.6)),length(unique(dge@meta.data$res.0.7)),length(unique(dge@meta.data$res.0.8)),length(unique(dge@meta.data$res.0.9)),length(unique(dge@meta.data$res.1)),length(unique(dge@meta.data$res.1.1)),length(unique(dge@meta.data$res.1.2)) ))
### decide to use 13 clusters, double-check the number of clusters
res="res.0.6";j=1;resi=1;
   dge=SetAllIdent(dge,id=res)
   print(length(unique(dge@ident)))
   TSNEPlot(dge)

## order cell by cluster ID and randomly shuffle cells within each batch
levels=levels(dge@ident)
ident=factor(dge@ident,levels=levels)

### randomly shuffling cells within each cluster
cells=sort(ident)
cells.use=NULL
for(i in 1:length(levels)){
   set.seed(i)
   tmp=cells[which(cells == levels[i])]
   if(length(tmp)>0){
      tmpname=sample(names(tmp),length(tmp),replace=FALSE)
      cells.use=c(cells.use,tmpname)
   }
}
cells.ident=as.factor(cells)
names(cells.ident)=cells.use
levels(cells.ident)=levels

### for each cluster, calculate average normalized expression of each gene
tmpdge=data.frame(t(as.matrix(dge@data[,cells.use])))
# make sure same order for cells.ident and dge before combining
which(names(cells.ident)!=colnames(dge@data)[cells.use])  # nothing
mouseclustersall=data.frame(ident=cells.ident,t(as.matrix(dge@data[,cells.use])))
genecountsall=matrix(,dim(dge@data)[1],length(unique(mouseclustersall$ident)))
rownames(genecountsall)=rownames(dge@data)
colnames(genecountsall)=unique(mouseclustersall$ident)
for(i in unique(mouseclustersall$ident)){
    genecountsall[,i]=apply(mouseclustersall[mouseclustersall$ident==i,-1],2,function(x) ExpMean(as.numeric(x))) # log(mean(exp(x) - 1) + 1)
    print(i)
}

### Reordering cluster centroid using dissimilarity matrix
library(seriation)
n=ncluster=length(levels)
nbatch=1 # nbatch=length(dataset)
bb=1
tmp=genecountsall[,levels]
tmp=tmp
colnames(tmp)=gsub(".*_","",colnames(tmp))
da <- dist(t(as.matrix(tmp)), method = "euclidean")
# note: dist calculate distance between each row
length(da) # 91
da
###### plot with seriation
pdf(file=paste0("MergedAllGermCells_Centroid_norm_Seriation_",res,".pdf"))
 dissplot(da, method="OLO",options = list(main = paste("Dissimilarity with seriation OLO"),col=bluered(10)))
dev.off()

### get order of seriation
do=seriate(da,method="OLO")
print(get_order(do))
levelss=get_order(do)
levelss=rev(levelss)
levels=levelss
print(levels-1)

### Reordered clusters for all cells
cells.use=colnames(dge@data)
# random shuffling cells within ordered clusters
ident=factor(dge@ident,levels=levels-1)

cells=sort(ident)
cells.use=NULL
for(i in 1:length(levels)){
   set.seed(i)
   tmp=cells[which(cells == levels[i]-1)]
   if(length(tmp)>0){
      tmpname=sample(names(tmp),length(tmp),replace=FALSE)
      cells.use=c(cells.use,tmpname)
   }
}
cells.ident.ordered=factor(as.numeric(cells),ordered=TRUE)
names(cells.ident.ordered)=cells.use

### save ordered cluster ID in dge object
print(which(levels(cells)!=levelss[[j]]-1)) # integer(0)

ordered=paste0(res,"order")

dge=AddMetaData(dge,cells.ident.ordered,ordered)
dge@meta.data[,ordered]=factor(dge@meta.data[,ordered])
dge=SetAllIdent(dge,ordered)
# save the dge file

save(dge,file="GermCells_DirectMergeAll_Jun2019.Robj")


## double-check if I need to reverse the cluster ID orders
### use Acrv1 as marker for the last cluster (round spermatid)
knownmarkers=c("Zbtb16","Gfra1","Tcf3","Morc1","Id4","Kit","Sohlh1","Stra8","Prdm9","Zcwpw1","Acrv1","Prm1")
length(knownmarkers)# [1] 12
pdf("clusters_ordered1_Prm1.pdf",height=10,width=10)
VlnPlot(dge,knownmarkers,cols.use=myBrewerPalette,nCol=3,point.size.use=-1)
dev.off()
### visualize known markers in t-SNE heatmap
pdf("knownmarkers.pdf",height=9,width=12)
FeaturePlot(object = dge, reduction.use="pca",features.plot = knownmarkers,cols.use = c("lightblue", 
    "red"),nCol=4)
FeaturePlot(object = dge, reduction.use="umap",features.plot = knownmarkers,cols.use = c("lightblue", 
    "red"),nCol=4)
dev.off()


###### assign cell types based on markers for each ordered clusters in each dataset
markers=FindAllMarkers(dge,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,"MergedAllGermCells_res.0.6order_markersall_mindiff0.2_logfc2fold_2019.txt",col.names=T,row.names=T,quote=F,sep="\t")
table(dge@ident)
table(markers$cluster)

#ClusterID  1   2   3   4   5   6   7   8   9  10  11  12  13 
#nCell     166 405 444 356 375  75 355  63 150 200 362 219 181 
#nMarker   290 134 106 107  37 163  55  46  15 234 348 261 778

markers %>% group_by(cluster) %>% top_n(5, avg_logFC)  -> top5
data.frame(top5)
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)  -> top2
pdf("markerstop.pdf",height=12,width=21)
FeaturePlot(object = dge, reduction.use="pca",features.plot = top2$gene, min.cutoff = "q9", cols.use = c("lightblue", 
    "red"), pt.size = .8,nCol=7)
FeaturePlot(object = dge, reduction.use="umap",features.plot = top2$gene, min.cutoff = "q9", cols.use = c("lightblue", 
    "red"), pt.size = .8,nCol=7)
dev.off()

### later on decided to remove cluster 1 (N=166 cells) and merge the rest 12 ordered clusters to 7 germ cell types as below 
table(dgerename@ident,dge@ident[names(dgerename@ident)])                      
#                          1   2   3   4   5   6   7   8   9  10  11  12  13
#  Undifferentiated SPGs   0 405   0   0   0   0   0   0   0   0   0   0   0
#  Differentiating SPGs    0   0 444   0   0   0   0   0   0   0   0   0   0
#  Leptotene/Zygotene      0   0   0 356 375   0   0   0   0   0   0   0   0
#  Early-Mid Pachytene     0   0   0   0   0  75 355  63   0   0   0   0   0
#  Late Pachytene          0   0   0   0   0   0   0   0 150 200 362   0   0
#  Diplotene               0   0   0   0   0   0   0   0   0   0   0 219   0
#  Early Rounds            0   0   0   0   0   0   0   0   0   0   0   0 181

### PCA and tSNE plot
dge=SetAllIdent(dge,id="orig.ident")
pdf("PCA_tSNE_UMAP_MergedAllGermCells_Batch.pdf",width=10.5,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=F)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
dge=SetAllIdent(dge,id="time")
dge@ident=factor(dge@ident,levels=time)
pdf("PCA_tSNE_UMAP_MergedAllGermCells_Time.pdf",width=10,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=F)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
dge=SetAllIdent(dge,id="type")
dge@ident=factor(dge@ident,levels=type)
pdf("PCA_tSNE_UMAP_MergedAllGermCells_Type.pdf",width=10,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=F)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
dge=SetAllIdent(dge,id="res.0.6order")
pdf("PCA_tSNE_UMAP_MergedAllGermCells_res.0.6order.pdf",width=9.5,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=T,colors.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
