#Run CCA on ALL cells (both somatic and germ cells) by Hailey
#There was an issue before with duplicate cell names, to circumvent this, name each cell with the time point and genotype. 
#Dates ran: 9/3/19, 9/4/19, 9/11/19

#8 datasets:
#  4 time points: 8 dpi, 13 dpi, 16 dpi, 21 dpi
#  2 types: WT vs SCARKO (sertoli-cell Androgen-receptor knockout)

library(Seurat)
myBrewerPalette<-c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#f781bf")

#### Step 1: load each dataset. 

dpi8WT<-read.table("/home/hlarose/learning/117507/out_2500cells_gene_exon_tagged_cleaned5000cells.dge.txt", sep = "\t", header = T, stringsAsFactors=F, row.names=1)

dpi8SCAR<-read.table("/home/hlarose/learning/117508/out_3000cells_gene_exon_tagged_cleaned5000cells.dge.txt", sep = "\t", header = T, stringsAsFactors=F, row.names=1)

dpi13WT<-read.table("/home/hlarose/learning/117697/out_2500cells_gene_exon_tagged_cleaned5000cells.dge.txt", sep = "\t", header = T, row.names = 1)

dpi13SCAR<-read.table("/home/hlarose/learning/117698/out_5000cells_gene_exon_tagged_cleaned5000cells.dge.txt", sep = "\t", header = T, row.names = 1)

dpi16WT<-read.table("/home/hlarose/learning/115237/out_2500cells_gene_exon_tagged_cleaned5000cells.dge.txt", sep = "\t", header = T, row.names = 1)

dpi16SCAR<-read.table("/home/hlarose/learning/115238/out_4500cells_gene_exon_tagged_cleaned5000cells.dge.txt", sep = "\t", header = T, row.names = 1)

dpp21WT<-read.table("/home/hlarose/learning/115235/out_2500cells_gene_exon_tagged_cleaned5000cells.dge.txt", sep = "\t", header = T, row.names = 1)

dpp21SCAR<-read.table("/home/hlarose/learning/115236/out_4000cells_gene_exon_tagged_cleaned5000cells.dge.txt", sep = "\t", header = T, row.names = 1)

### Change cell names

colnames(dpi8WT) <- paste("dpi8WT", colnames(dpi8WT), sep = "_")
colnames(dpi8SCAR) <- paste("dpi8SCAR", colnames(dpi8SCAR), sep = "_")
colnames(dpi13WT) <- paste("dpi13WT", colnames(dpi13WT), sep = "_")
colnames(dpi13SCAR) <- paste("dpi13SCAR", colnames(dpi13SCAR), sep = "_")
colnames(dpi16WT) <- paste("dpi16WT", colnames(dpi16WT), sep = "_")
colnames(dpi16SCAR) <- paste("dpi16SCAR", colnames(dpi16SCAR), sep = "_")
colnames(dpp21WT) <- paste("dpp21WT", colnames(dpp21WT), sep = "_")
colnames(dpp21SCAR) <- paste("dpp21SCAR", colnames(dpp21SCAR), sep = "_")

########## Filter each data set and then set up the Seurat object for each in prep for CCA alignment. 

@@### 8dpi WT - mouse 1

nGeneperCell <- colSums(dpi8WT>0)
dpi8WT.tmp=dpi8WT[,nGeneperCell>500]
mito.genes <- grep("^mt-", rownames(dpi8WT.tmp), value = T) 
percent.mito <- colSums(dpi8WT.tmp[mito.genes, ])/colSums(dpi8WT.tmp)
length(which(percent.mito>0.1))    #[1] 127
dpi8WT.tmp=dpi8WT[,nGeneperCell>500][,percent.mito<0.1]

### Filter for genes: 

nCellperGene <- rowSums(dpi8WT>0)
nUMIperGene <- rowSums(dpi8WT)
dgedata1=dpi8WT.tmp[which(nCellperGene>3 & nUMIperGene>3),]

### Set up Seurat object.

dge1 <- CreateSeuratObject(raw.data = dgedata1, project = "8dpiWT", min.cells = 3, min.genes = 500)
dge1 <- NormalizeData(dge1)
dge1 <- ScaleData(dge1, display.progress = F)
dge1@meta.data$stim <- "WT"

@@### 8dpi SCARKO - mouse 2

nGeneperCell <- colSums(dpi8SCAR>0)
dpi8SCAR.tmp=dpi8SCAR[,nGeneperCell>500]
mito.genes <- grep("^mt-", rownames(dpi8SCAR.tmp), value = T) 
percent.mito <- colSums(dpi8SCAR.tmp[mito.genes, ])/colSums(dpi8SCAR.tmp)
length(which(percent.mito>0.1))    #[1] 10
dpi8SCAR.tmp=dpi8SCAR[,nGeneperCell>500][,percent.mito<0.1]

### Filter for genes: 

nCellperGene <- rowSums(dpi8SCAR>0)
nUMIperGene <- rowSums(dpi8SCAR)
dgedata2=dpi8SCAR.tmp[which(nCellperGene>3 & nUMIperGene>3),]

### Set up Seurat object. 

dge2 <- CreateSeuratObject(raw.data = dgedata2, project = "8dpiSCARKO", min.cells = 3, min.genes = 500)
dge2 <- NormalizeData(dge2)
dge2 <- ScaleData(dge2, display.progress = F)
dge2@meta.data$stim <- "SCARKO"

@@### 13dpi WT - mouse 3

nGeneperCell <- colSums(dpi13WT>0)
dpi13WT.tmp=dpi13WT[,nGeneperCell>500]
mito.genes <- grep("^mt-", rownames(dpi13WT.tmp), value = T) 
percent.mito <- colSums(dpi13WT.tmp[mito.genes, ])/colSums(dpi13WT.tmp)
length(which(percent.mito>0.1))    #[1] 31
dpi13WT.tmp=dpi13WT[,nGeneperCell>500][,percent.mito<0.1]

### Filter for genes: 

nCellperGene <- rowSums(dpi13WT>0)
nUMIperGene <- rowSums(dpi13WT)
dgedata3=dpi13WT.tmp[which(nCellperGene>3 & nUMIperGene>3),]

### Set up Seurat object.

dge3 <- CreateSeuratObject(raw.data = dgedata3, project = "13dpi_WT", min.cells = 3, min.genes = 500)
dge3 <- NormalizeData(dge3)
dge3 <- ScaleData(dge3, display.progress = F)
dge3@meta.data$stim <- "WT"

@@### 13dpi SCARKO - mouse 4

nGeneperCell <- colSums(dpi13SCAR>0)
dpi13SCAR.tmp=dpi13SCAR[,nGeneperCell>500]
mito.genes <- grep("^mt-", rownames(dpi13SCAR.tmp), value = T) 
percent.mito <- colSums(dpi13SCAR.tmp[mito.genes, ])/colSums(dpi13SCAR.tmp)
length(which(percent.mito>0.1))    
#[1] 89
dpi13SCAR.tmp=dpi13SCAR[,nGeneperCell>500][,percent.mito<0.1]

### Filter for genes: 

nCellperGene <- rowSums(dpi13SCAR>0)
nUMIperGene <- rowSums(dpi13SCAR)
dgedata4=dpi13SCAR.tmp[which(nCellperGene>3 & nUMIperGene>3),]

### Set up Seurat object.

dge4 <- CreateSeuratObject(raw.data = dgedata4, project = "13dpi_SCARKO", min.cells = 3, min.genes = 500)
dge4 <- NormalizeData(dge4)
dge4 <- ScaleData(dge4, display.progress = F)
dge4@meta.data$stim <- "SCARKO"

@@### 16dpi WT - mouse 5

nGeneperCell <- colSums(dpi16WT>0)
dpi16WT.tmp=dpi16WT[,nGeneperCell>500]
mito.genes <- grep("^mt-", rownames(dpi16WT.tmp), value = T) 
percent.mito <- colSums(dpi16WT.tmp[mito.genes, ])/colSums(dpi16WT.tmp)
length(which(percent.mito>0.1))    #[1] 27
dpi16WT.tmp=dpi16WT[,nGeneperCell>500][,percent.mito<0.1]

### Filter for genes: 

nCellperGene <- rowSums(dpi16WT>0)
nUMIperGene <- rowSums(dpi16WT)
dgedata5=dpi16WT.tmp[which(nCellperGene>3 & nUMIperGene>3),]

### Set up Seurat object.

dge5 <- CreateSeuratObject(raw.data = dgedata5, project = "16dpi_WT", min.cells = 3, min.genes = 500)
dge5 <- NormalizeData(dge5)
dge5 <- ScaleData(dge5, display.progress = F)
dge5@meta.data$stim <- "WT"

@@### 16dpi SCARKO - mouse 6

nGeneperCell <- colSums(dpi16SCAR>0)
dpi16SCAR.tmp=dpi16SCAR[,nGeneperCell>500]
mito.genes <- grep("^mt-", rownames(dpi16SCAR.tmp), value = T) 
percent.mito <- colSums(dpi16SCAR.tmp[mito.genes, ])/colSums(dpi16SCAR.tmp)
length(which(percent.mito>0.1))    #[1] 27
dpi16SCAR.tmp=dpi16SCAR[,nGeneperCell>500][,percent.mito<0.1]

### Filter for genes: 

nCellperGene <- rowSums(dpi16SCAR>0)
nUMIperGene <- rowSums(dpi16SCAR)
dgedata6=dpi16SCAR.tmp[which(nCellperGene>3 & nUMIperGene>3),]

### Set up Seurat object.

dge6 <- CreateSeuratObject(raw.data = dgedata6, project = "16dpi_SCARKO", min.cells = 3, min.genes = 500)
dge6 <- NormalizeData(dge6)
dge6 <- ScaleData(dge6, display.progress = F)
dge6@meta.data$stim <- "SCARKO"

@@### 21dpp WT - mouse 7

nGeneperCell <- colSums(dpp21WT>0)
dpp21WT.tmp=dpp21WT[,nGeneperCell>500]
mito.genes <- grep("^mt-", rownames(dpp21WT.tmp), value = T) 
percent.mito <- colSums(dpp21WT.tmp[mito.genes, ])/colSums(dpp21WT.tmp)
length(which(percent.mito>0.1))    #[1] 64
dpp21WT.tmp=dpp21WT[,nGeneperCell>500][,percent.mito<0.1]

### Filter for genes: 

nCellperGene <- rowSums(dpp21WT>0)
nUMIperGene <- rowSums(dpp21WT)
dgedata7=dpp21WT.tmp[which(nCellperGene>3 & nUMIperGene>3),]

### Set up Seurat object.

dge7 <- CreateSeuratObject(raw.data = dgedata7, project = "21dpp_WT", min.cells = 3, min.genes = 500)
dge7 <- NormalizeData(dge7)
dge7 <- ScaleData(dge7, display.progress = F)
dge7@meta.data$stim <- "WT"

### 21dpp SCARKO - mouse 8

nGeneperCell <- colSums(dpp21SCAR>0)
dpp21SCAR.tmp=dpp21SCAR[,nGeneperCell>500]
mito.genes <- grep("^mt-", rownames(dpp21SCAR.tmp), value = T) 
percent.mito <- colSums(dpp21SCAR.tmp[mito.genes, ])/colSums(dpp21SCAR.tmp)
length(which(percent.mito>0.1))    #[1] 16
dpp21SCAR.tmp=dpp21SCAR[,nGeneperCell>500][,percent.mito<0.1]

### Filter for genes: 

nCellperGene <- rowSums(dpp21SCAR>0)
nUMIperGene <- rowSums(dpp21SCAR)
dgedata8=dpp21SCAR.tmp[which(nCellperGene>3 & nUMIperGene>3),]

### Set up Seurat object.

dge8 <- CreateSeuratObject(raw.data = dgedata8, project = "21dpp_SCARKO", min.cells = 3, min.genes = 500)
dge8 <- NormalizeData(dge8)
dge8 <- ScaleData(dge8, display.progress = F)
dge8@meta.data$stim <- "SCARKO"

### Add metadata information - should have done this 62019, not in original

dge1@meta.data$group <- "dpi8WT"
dge2@meta.data$group <- "dpi8SCAR"
dge3@meta.data$group <- "dpi13WT"
dge4@meta.data$group <- "dpi13SCAR"
dge5@meta.data$group <- "dpi16WT"
dge6@meta.data$group <- "dpi16SCAR"
dge7@meta.data$group <- "dpp21WT"
dge8@meta.data$group <- "dpp21SCAR"

#Find Variable genes for each sample. 
dge1<-FindVariableGenes(dge1,do.plot=F)
dge2<-FindVariableGenes(dge2,do.plot=F)
dge3<-FindVariableGenes(dge3,do.plot=F)
dge4<-FindVariableGenes(dge4,do.plot=F)
dge5<-FindVariableGenes(dge5,do.plot=F)
dge6<-FindVariableGenes(dge6,do.plot=F)
dge7<-FindVariableGenes(dge7,do.plot=F)
dge8<-FindVariableGenes(dge8,do.plot=F)

# Gene selection for input to CCA
# Here we take the union of the top 2,000 genes with the highest dispersion (var/mean) from both datasets.

g.1 <- head(rownames(dge1@hvg.info), 2000)
g.2 <- head(rownames(dge2@hvg.info), 2000)
g.3 <- head(rownames(dge3@hvg.info), 2000)
g.4 <- head(rownames(dge4@hvg.info), 2000)
g.5 <- head(rownames(dge5@hvg.info), 2000)
g.6 <- head(rownames(dge6@hvg.info), 2000)
g.7 <- head(rownames(dge7@hvg.info), 2000)
g.8 <- head(rownames(dge8@hvg.info), 2000)


genes.use <- unique(c(g.1, g.2, g.3, g.4, g.5, g.6, g.7, g.8))
genes.use <- intersect(genes.use, rownames(dge1@scale.data))
genes.use <- intersect(genes.use, rownames(dge2@scale.data))
genes.use <- intersect(genes.use, rownames(dge3@scale.data))
genes.use <- intersect(genes.use, rownames(dge4@scale.data))
genes.use <- intersect(genes.use, rownames(dge5@scale.data))
genes.use <- intersect(genes.use, rownames(dge6@scale.data))
genes.use <- intersect(genes.use, rownames(dge7@scale.data))
genes.use <- intersect(genes.use, rownames(dge8@scale.data))

### RunMuliCCA for all 8 samples. 

samples.list<-list(dge1, dge2, dge3, dge4, dge5, dge6, dge7, dge8)

all.cca <- RunMultiCCA(object.list = samples.list, genes.use = genes.use, num.ccs = 10)

### Run biocor to determine dimensions to use...

pdf(file = 'all.cells.cca_bicor_090319.pdf')
MetageneBicorPlot(all.cca, grouping.var = "stim", dims.eval = 1:8, display.progress = FALSE)
dev.off()

### Align samples! Pick dims based on the bicor plot above.. this didn't look great, so chose 3. 

all.cca <- AlignSubspace(all.cca, reduction.type = "cca", grouping.var = "stim", dims.align = 1:3)

all.cca <- RunTSNE(all.cca, reduction.use = "cca.aligned", dims.use = 1:3, do.fast = T, check_duplicates = FALSE) #Unclear why there are duplicates... cells should have unique names
all.cca <- RunPCA(all.cca, reduction.use = "cca.aligned", dims.use = 1:3, do.fast = T)
all.cca <- FindClusters(all.cca, reduction.type = "cca.aligned", resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), dims.use = 1:3)

print(c( length(unique(all.cells.cca@meta.data$res.0.1)),length(unique(all.cells.cca@meta.data$res.0.2)),length(unique(all.cells.cca@meta.data$res.0.3)),length(unique(all.cells.cca@meta.data$res.0.4)),length(unique(all.cells.cca@meta.data$res.0.5)),length(unique(all.cells.cca@meta.data$res.0.6)),length(unique(all.cells.cca@meta.data$res.0.7)),length(unique(all.cells.cca@meta.data$res.0.8)),length(unique(all.cells.cca@meta.data$res.0.9)),length(unique(all.cells.cca@meta.data$res.1))))

# [1]  5  7 10 14 15 16 17 19 20 20

#Lets try res.0.2 and res.0.3 and see what the TSNE plots and marker tables look like. 

all.cells.cca<-all.cca
all.cells.cca<-SetAllIdent(all.cells.cca, id='res.0.2')
pdf(file="all.cells.cca_tsne_dim1.3_res.0.2_090319.pdf")
TSNEPlot(all.cells.cca, do.label = T, do.return = T, pt.size = 0.5, colors.use= myBrewerPalette)
dev.off()

all.cells.cca<-SetAllIdent(all.cells.cca, id='res.0.3')
pdf(file="all.cells.cca_tsne_dim1.3_res.0.3_090319.pdf")
TSNEPlot(all.cells.cca, do.label = T, do.return = T, pt.size = 0.5, colors.use= myBrewerPalette)
dev.off()

# save Robj - take this to desktop if not too large (ha)
save(all.cells.cca, file='all.cells.cca_dim1.3_res.0.2_090319.Robj')

# Selected res.0.3 because res.0.2 was under-clustered. Call markers!

dge<-SetAllIdent(all.cells.cca, id='res.0.3')

#Took back to laptop to run UMAP, and its clear that cluster 5 is a mess... doublets? Remove! 

dge<-SubsetData(dge, ident.remove='5', subset.raw=T)

### Find markers for the clusters, so I know what they are. 

markersall=FindAllMarkers(dge,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)

write.table(markersall, file="all.cells.cca_dim1.3_res.0.2_no5_090319_markers.txt")

#### Run cluster cluster correlation to get the order of the clusters. 
##### order cell by cluster ID and randomly shuffle cells within each batch

dge=SetAllIdent(dge, id="res.0.3") #MODIFY
ident=dget@ident
levels=levels(ident)
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

##### Calculate centroid of each cluster

tmpdge=data.frame(t(as.matrix(dget@data[,cells.use])))

## make sure same order for cells.ident and dge before combining

which(names(cells.ident)!=colnames(dget@data)[cells.use])  # nothing
mouseclustersall=data.frame(ident=cells.ident,t(as.matrix(dget@data[,cells.use])))

##### for each cluster, calculate average normalized expression of each gene

genecountsall=matrix(,dim(dget@data)[1],length(unique(mouseclustersall$ident)))

rownames(genecountsall)=rownames(dget@data)

colnames(genecountsall)=unique(mouseclustersall$ident)

for(i in unique(mouseclustersall$ident)){
  genecountsall[,i]=apply(mouseclustersall[mouseclustersall$ident==i,-1],2,function(x) ExpMean(as.numeric(x))) # log(mean(exp(x) - 1) + 1)
  print(i)
}

### Save the R image do you don't have to load/create everything again...this can take up a lot of memory, don't do this unless you're in 'final' mode. 

# save.image(file = "117507_11918.RData")

##### Ordering the cluster centroids using dissimilarity matrix #####

library(seriation)

tmp=genecountsall[,levels]

da <- dist(t(as.matrix(tmp)), method = "euclidean")

pdf(file=paste("all.cells.cca_no5_merg2.3_renamed_090419_Centroid_norm_Seriation_Dissimilarity.pdf")) 

### plot original matrix

dissplot(da, method = NA,options = list(main = paste("All Cells CCA")))

### plot with seriation using method "ARSA"

dissplot(da, method="ARSA", options = list(main = paste("Dissimilarity with seriation ARSA")))

### plot with seriation using method "OLO"

dissplot(da, method="OLO",options = list(main = paste("Dissimilarity with seriation OLO")))

dev.off()

### Rename idents based off the chart above... 

dge1<-dge
current.cluster.ids <- c(0, 1, 2, 3, 4, 6, 7,8,9)
new.cluster.ids <- c("C4", "C5", "C2", "C2", "C3", "C1", "C7", "C6", "C8")
dge1@ident <- plyr::mapvalues(x = dge@ident, from = current.cluster.ids, to = new.cluster.ids)

dge2<-dge1
current.cluster.ids <- c("C4", "C5", "C2", "C2", "C3", "C1", "C7", "C6", "C8")
new.cluster.ids <- c("4", "5", "2", "2", "3", "1", "7", "6", "8")
dge2@ident <- plyr::mapvalues(x = dge1@ident, from = current.cluster.ids, to = new.cluster.ids)

dge1<-dge
current.cluster.ids <- c(0, 1, 2, 3, 4, 6, 7,8,9)
new.cluster.ids <- c("Immature Leydig", "Adult Leydig", "Interstitial Progenitor", "Interstitial Progenitor", "Stem Leydig", "Macrophage", "Spermatogonia", "Sertoli", "Spermatocytes")
dget@ident <- plyr::mapvalues(x = all.cells.cca@ident, from = current.cluster.ids, to = new.cluster.ids)

dge1<-dge
current.cluster.ids <- c(0, 1, 2, 3, 4, 6, 7,8,9)
new.cluster.ids <- c(0, 1, 2, 2, 4, 6, 7,8,9)
dge1@ident <- plyr::mapvalues(x = dge@ident, from = current.cluster.ids, to = new.cluster.ids)


### save ordered cluster ID in dge object

save(all.cells.cca, file='all.cells.cca_dim1.3_res.0.2_ordered_090319.Robj')

@#%$^&*()*&^%$^&*()_(*&^%$#@$%^&*()_) Fix this damn plot so it is in order!

### Find markers for the clusters, so I know what they are. 

markersall=FindAllMarkers(dge,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)

write.table(markersall, file="test_all.cells.cca_dim1.3_res.0.2_ordered_090319_markers.txt")

### Make a heat map. 

#Add this function
expMean <- function(x) {
  return(log(x = mean(x = exp(x = x) - 1) + 1))} 

# Calculate cluster centroid vals
centroid=matrix(,dim(dget@data)[1],length(unique(dget@ident)))
rownames(centroid)=rownames(dget@data)
colnames(centroid)=unique(dget@ident)
foss Cell Types
centroid.std=(centroid-apply(centroid,1,mean))/apply(centroid,1,sd)
write.table(centroid.std,"AllCellsCCA_centroid_allgenes_std_090419.txt",row.names=T,col.names=T,quote=Fr(i in levels(dget@ident)){
    centroid[,i]=apply(dget@data[,which(dget@ident==i)],1,function(x) expMean(as.numeric(x))) 
    print(i)}

### Standardize genes Acro,sep="\t") #Can load this again later. 

# set up variables for the heatmap. 
levels=levels(dget@ident)
 colsep.use=cumsum(table(gsub("_.*","",levels))[levels])
 col.lab=rep("",length(levels))
 col.lab=gsub(".*_","",levels)
 library(RColorBrewer)
 myBrewerPalette <- c(brewer.pal(12,"Paired"))[c(12,11,10,6,9,5,7,1:4)]
 ncluster=length(levels)
 sidecol=matrix(0,2,length(levels))
 sidecol[1,]=rep(rep(c("white","white"),each=12),3)[1:sum(ncluster)]
 sidecol[2,]=myBrewerPalette[1:sum(ncluster)]
 clab=cbind(sidecol[2,],sidecol[1,])
 rlab=sidecol
 rownames(rlab)=c("","Cell Type")
 colnames(clab)=c("Cell Type","")
 redblue100<- scale_color_gradient(low = "blue", high = "red")
 col.use=redblue100
markersall=FindAllMarkers(dget,only.pos=TRUE,min.pct=0.2,min.diff.pct=0.2,logfc.threshold=log(2),test.use="bimod",do.print = TRUE)
data.use=centroid.std[markersall$gene,]
row.lab=rownames(data.use)
 
write.table(markersall, file="all.cells.cca_dim1.3_res.0.3_no5_merged2.3_090419_markers.txt")

#plot the heatmap
 jpeg((file="workplease.jpeg"),res=300,height=2600,width=1600)
 par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",colsep = colsep.use,sep.color="black",sepwidth=c(0.001,0.001),ColSideColors=clab,labCol=col.lab,labRow=row.lab,ColSideColorsSize = 2,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()

## Plot in TSNE 
pdf(file="all.cells.cca_tsne_dim1.3_res.0.2_090319.pdf")
TSNEPlot(dge, do.label = T, do.return = T, pt.size = 0.5, colors.use= myBrewerPalette)
dev.off()

##Plot each sample over the TSNE. 
all.cca<-dge
p="dpi8WT" 
  all.cca=SetAllIdent(all.cca,id="meta.data$group")
  ident=as.character(all.cca@meta.data$group)
  names(ident)=names(all.cca@meta.data$group)
  ident[grep(p,ident)]<-p
  ident[which(!grepl(p,ident))]<-"Others"
  ident=factor(ident,levels=c(p,"Others"),order=T)
  all.cca@meta.data$group=ident
p1<-DimPlot(all.cca,group.by='group', cols.use=c("grey","#a6cee3"), plot.order='dpi8WT', reduction.use='umap')

all.cca<-dge
p="dpi8SCAR" 
  all.cca=SetAllIdent(all.cca,id="meta.data$group")
  ident=as.character(all.cca@meta.data$group)
  names(ident)=names(all.cca@meta.data$group)
  ident[grep(p,ident)]<-p
  ident[which(!grepl(p,ident))]<-"Others"
  ident=factor(ident,levels=c(p,"Others"),order=T)
  all.cca@meta.data$group=ident
p2<-DimPlot(all.cca,group.by='group', cols.use=c("grey","#1f78b4"), plot.order='dpi8SCAR', reduction.use='umap')

all.cca<-dge
p="dpi13WT" 
  all.cca=SetAllIdent(all.cca,id="meta.data$group")
  ident=as.character(all.cca@meta.data$group)
  names(ident)=names(all.cca@meta.data$group)
  ident[grep(p,ident)]<-p
  ident[which(!grepl(p,ident))]<-"Others"
  ident=factor(ident,levels=c(p,"Others"),order=T)
  all.cca@meta.data$group=ident
p3<-DimPlot(all.cca,group.by='group', cols.use=c("grey","#b2df8a"), plot.order='dpi13WT', reduction.use='umap')

all.cca<-dge
p="dpi13SCAR" 
  all.cca=SetAllIdent(all.cca,id="meta.data$group")
  ident=as.character(all.cca@meta.data$group)
  names(ident)=names(all.cca@meta.data$group)
  ident[grep(p,ident)]<-p
  ident[which(!grepl(p,ident))]<-"Others"
  ident=factor(ident,levels=c(p,"Others"),order=T)
  all.cca@meta.data$group=ident
p4<-DimPlot(all.cca,group.by='group', cols.use=c("grey","#33a02c"), plot.order='dpi13SCAR', reduction.use='umap')

all.cca<-dge
p="dpi13WT" 
  all.cca=SetAllIdent(all.cca,id="meta.data$group")
  ident=as.character(all.cca@meta.data$group)
  names(ident)=names(all.cca@meta.data$group)
  ident[grep(p,ident)]<-p
  ident[which(!grepl(p,ident))]<-"Others"
  ident=factor(ident,levels=c(p,"Others"),order=T)
  all.cca@meta.data$group=ident
p5<-DimPlot(all.cca,group.by='group', cols.use=c("grey","#fb9a99"), plot.order='dpi13WT', reduction.use='umap')

all.cca<-dge
p="dpi16SCAR" 
  all.cca=SetAllIdent(all.cca,id="meta.data$group")
  ident=as.character(all.cca@meta.data$group)
  names(ident)=names(all.cca@meta.data$group)
  ident[grep(p,ident)]<-p
  ident[which(!grepl(p,ident))]<-"Others"
  ident=factor(ident,levels=c(p,"Others"),order=T)
  all.cca@meta.data$group=ident
p6<-DimPlot(all.cca,group.by='group', cols.use=c("grey","#e31a1c"), plot.order='dpi16SCAR', reduction.use='umap')

all.cca<-dge
p="dpp21WT" 
  all.cca=SetAllIdent(all.cca,id="meta.data$group")
  ident=as.character(all.cca@meta.data$group)
  names(ident)=names(all.cca@meta.data$group)
  ident[grep(p,ident)]<-p
  ident[which(!grepl(p,ident))]<-"Others"
  ident=factor(ident,levels=c(p,"Others"),order=T)
  all.cca@meta.data$group=ident
p7<-DimPlot(all.cca,group.by='group', cols.use=c("grey","#fdbf6f"), plot.order='dpp21WT', reduction.use='umap')

all.cca<-dge
p="dpp21SCAR" 
  all.cca=SetAllIdent(all.cca,id="meta.data$group")
  ident=as.character(all.cca@meta.data$group)
  names(ident)=names(all.cca@meta.data$group)
  ident[grep(p,ident)]<-p
  ident[which(!grepl(p,ident))]<-"Others"
  ident=factor(ident,levels=c(p,"Others"),order=T)
  all.cca@meta.data$group=ident
p8<-DimPlot(all.cca,group.by='group', cols.use=c("grey","#ff7f00"), plot.order='dpp21SCAR', reduction.use='umap')

plot_grid(p1,p2,p3,p4,p5,p6,p7,p8)
pdf((file='all.cells.cca.grey.background_landscape_090319.pdf'), width=20, height=10)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, ncol=4, nrow=2)
dev.off()

## Figures

dge<-all.cells.cca
rm(all.cells.cca)
mito.genes <- grep("^mt-", rownames(dge@data), value = T)
length(mito.genes) #24
[1] 26
percent.mito <- colSums(expm1(dge@data[mito.genes, ]))/colSums(expm1(dge@data))
dge <- AddMetaData(dge, percent.mito, "percent.mito")

# add %ChrX and %ChrY genes to the dge datainfo
x=read.table(("mouseChrXgenes"),stringsAsFactors=FALSE)[,1]
y=read.table(("mouseChrYgenes"),stringsAsFactors=FALSE)[,1]
length(x) # 2605
length(y) # 1570

chrx.genes <- x[which(x %in% rownames(dge@data))]
chry.genes <- y[which(y %in% rownames(dge@data))]
length(chrx.genes)  # 1520
length(chry.genes)  # 215
percent.x <- colSums(expm1(dge@data[chrx.genes, ]))/colSums(expm1(dge@data))
percent.y <- colSums(expm1(dge@data[chry.genes, ]))/colSums(expm1(dge@data))
dge <- AddMetaData(dge, percent.x, "percent.x")
dge <- AddMetaData(dge, percent.y, "percent.y")

v1<-VlnPlot(dge, "nGene", nCol = 1,cols.use=myBrewerPalette, point.size.use = F)
v2<-VlnPlot(dge, "nUMI", nCol = 1,cols.use=myBrewerPalette, point.size.use = F)
v3<-VlnPlot(dge, "percent.mito", nCol = 1,cols.use=myBrewerPalette, point.size.use = F)
v4<-VlnPlot(dge, "percent.x", nCol = 1,cols.use=myBrewerPalette, point.size.use = F)
v5<-VlnPlot(dge, "percent.y", nCol = 1,cols.use=myBrewerPalette, point.size.use = F)
plot_grid(v1, v2, v3, v4, v5)
pdf(file='all.cells.cca_NGeneUMImtXY.pdf')
grid.arrange(v1,v2,v3,v4,v5)
dev.off()

############# Pull out those Sertoli Cells...

##Plot the ICs over the germ cells (gray)

dget<-all.cells.cca
current.cluster.ids <- c(0, 1, 2, 4, 6, 7,8,9)
new.cluster.ids <- c("Immature Leydig", "Adult Leydig", "Interstitial Progenitor", "Stem Leydig", "Macrophage", "Spermatogonia", "Sertoli", "Spermatocytes")
dget@ident <- plyr::mapvalues(x = all.cells.cca@ident, from = current.cluster.ids, to = new.cluster.ids)

## All the cells:
pdf("all.cells.cca_res0.2_dim1.3_no5_merg2.3_renamed.pdf")
DimPlot(dget, reduction.use = 'umap', do.label = T, cols.use = myBrewerPalette)
dev.off()

##Over grey germ cells 
pdf("all.cells.cca_res0.2_dim1.3_no5_merg2.3_renamed_greyGCs_090919.pdf")
DimPlot(dget, reduction.use = 'umap', do.label = T, cols.use = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "grey", "#fdbf6f"))
dev.off()

dget<-all.cells.cca
current.cluster.ids <- c(0, 1, 2, 4, 6, 7,8,9)
new.cluster.ids <- c("Immature Leydig", "Adult Leydig", "Interstitial Progenitor", "Stem Leydig", "Macrophage", "Germ Cells", "Sertoli", "Germ Cells")
dget@ident <- plyr::mapvalues(x = all.cells.cca@ident, from = current.cluster.ids, to = new.cluster.ids)

#VlnPlot using Sox9 and Rhox8 to identify sertoli cells. 

pdf("all.cells.cca_res0.2_dim1.3_no5_merg2.3_renamed_greyGCs_090919.pdf", height=3)
VlnPlot(dget, features.plot = c("Sox9", "Rhox8"), cols.use = myBrewerPalette, point.size.use = -1, x.lab.rot = 90, size.x.use = 8, size.y.use = 8, size.title.use = 12)
dev.off()

scs<-SubsetData(dget, ident.use='Sertoli', subset.raw = T)

#Make a pdf of just the sertoli cells, in tsne space (looks better than umap)

pdf(file="sertoli_cells_tsne_fromall.cells.cca_bystim_090919.pdf")
DimPlot(scs, reduction.use='tsne', group.by='stim')
dev.off()

# Make a pdf of the contribution of each sample. 

scs<-SetAllIdent(scs, id= 'group')
table(scs@ident)

# dpi8SCAR    dpi8WT dpi13SCAR   dpi13WT dpi16SCAR   dpi16WT dpp21SCAR   dpp21WT 
#      557       205       190       112        82         9       251         4 

pdf("sertoli_contribution.per.timepoint_tsne_090919.pdf")
DimPlot(scs, reduction.use='tsne', do.label=T)
dev.off()

#Plot the Sox9 ratio.

pdf("sertoli_sox9_wtvscarko_090919.pdf", height=3)
VlnPlot(scs, features.plot = c("Sox9", "Rhox8"), point.size.use = -1, x.lab.rot = 90, size.x.use = 8, size.y.use = 8, size.title.use = 12)
dev.off()

#Call markers between the genotypes. 

markers=FindMarkers(scs,only.pos=TRUE,ident.1='WT', ident.2 = 'SCARKO', test.use="bimod",logfc.threshold = log(1.5),min.diff.pct=0.18,do.print = TRUE)
write.table(markers,"SCs_fromall.cells.cca_res0.3_markersWTvSCAR_logfc1.5fold_09092019.txt",col.names=T,row.names=T,quote=F,sep="\t")

scs_final<-save(scs, file="scs_final.all.cells.cca_res0.2_dim1.3_no5_merg2.3_091119.Robj")

@@@@@@@@@@@@@@@ Leydig cells: 

setwd("/Volumes/HL/all.cells.cca_Sept19/090419/")
load("all.cells.cca_umap_dim1.3_res.0.3_no5_merg2.3_rename_090319.Robj")
table(all.cells.cca@ident)

leydigs<-SubsetData(dget, ident=c("Immature Leydig", "Adult Leydig", "Stem Leydig"))

%%%%%%%%%%% Returning to analysis and making figures: 

setwd("/Volumes/HL/all.cells.cca_Sept19/090419/")
load("all.cells.cca_umap_dim1.3_res.0.3_no5_merg2.3_rename_090319.Robj")
table(dget@ident)

#Remake the nUMI, nGENE, mt%, etc: Picked up 9/11/19

mito.genes <- grep("^mt-", rownames(dget@data), value = T)
length(mito.genes)
[1] 26
percent.mito <- colSums(expm1(dget@data[mito.genes, ]))/colSums(expm1(dget@data))
dget <- AddMetaData(dget, percent.mito, "percent.mito")

# add %ChrX and %ChrY genes to the dge datainfo
x=read.table(("mouseChrXgenes"),stringsAsFactors=FALSE)[,1]
y=read.table(("mouseChrYgenes"),stringsAsFactors=FALSE)[,1]
length(x) # 2605
length(y) # 1570

# add %ChrX and %ChrY genes to the dge datainfo
x=read.table(("mouseChrXgenes"),stringsAsFactors=FALSE)[,1]
y=read.table(("mouseChrYgenes"),stringsAsFactors=FALSE)[,1]
length(x) # 2605
length(y) # 1570

chrx.genes <- x[which(x %in% rownames(dget@data))]
chry.genes <- y[which(y %in% rownames(dget@data))]
length(chrx.genes)  # 1520
length(chry.genes)  # 215
percent.x <- colSums(expm1(dget@data[chrx.genes, ]))/colSums(expm1(dget@data))
percent.y <- colSums(expm1(dget@data[chry.genes, ]))/colSums(expm1(dget@data))
dget <- AddMetaData(dget, percent.x, "percent.x")
dget <- AddMetaData(dget, percent.y, "percent.y")

v1<-VlnPlot(dget, "nGene", nCol = 1,cols.use=myBrewerPalette, point.size.use = F, size.x.use = 8, x.lab.rot = 90, size.title.use = 16, size.y.use = 8)
v2<-VlnPlot(dget, "nUMI", nCol = 1,cols.use=myBrewerPalette, point.size.use = F, size.x.use = 8, x.lab.rot = 90, size.title.use = 16, size.y.use = 8)
v3<-VlnPlot(dget, "percent.mito", nCol = 1,cols.use=myBrewerPalette, point.size.use = F, size.x.use = 8, x.lab.rot = 90, size.title.use = 16, size.y.use = 8)
v4<-VlnPlot(dget, "percent.x", nCol = 1,cols.use=myBrewerPalette, point.size.use = F, size.x.use = 8, x.lab.rot = 90, size.title.use = 16, size.y.use = 8)
v5<-VlnPlot(dget, "percent.y", nCol = 1,cols.use=myBrewerPalette, point.size.use = F, size.x.use = 8, x.lab.rot = 90, size.title.use = 16, size.y.use = 8)
plot_grid(v1, v2, v3, v4, v5)
pdf(file='all.cells.cca_NGeneUMImtXY_renamed_merged_res0.3_091119.pdf')
grid.arrange(v1,v2,v3,v4,v5)
dev.off()

### Plot SCs against everything else gray. 

all.cca<-dget
p="Sertoli" 
  ident=as.character(all.cca@ident)
  names(ident)=names(all.cca@ident)
  ident[grep(p,ident)]<-p
  ident[which(!grepl(p,ident))]<-"Others"
  ident=factor(ident,levels=c(p,"Others"),order=T)
  all.cca@ident=ident
pdf("Sertoli.all.cells.cca_greybackground_101619.pdf")
DimPlot(all.cca, cols.use=c("grey","#fdbf6f"), plot.order='Sertoli', reduction.use='umap')
dev.off()

rm(all.cca)

### Make a heat map with new labels - note this was done on the server.  

load("all.cells.cca_umap_dim1.3_res.0.3_no5_merg2.3_rename_090319.Robj")

new_order<-c("Spermatogonia", "Spermatocytes", "Sertoli", "Interstitial Progenitor", "Progenitor Leydig", "Immature Leydig", "Adult Leydig", "Macrophage")

dget@ident <- factor(dget@ident, levels = new_order)

dge1<-dget

#Add this function
expMean <- function(x) {
  return(log(x = mean(x = exp(x = x) - 1) + 1))} 

# Calculate cluster centroid vals
centroid=matrix(,dim(dge1@data)[1],length(unique(dge1@ident)))
rownames(centroid)=rownames(dge1@data)
colnames(centroid)=new_order
for(i in levels(dge1@ident)){
    centroid[,i]=apply(dge1@data[,which(dge1@ident==i)],1,function(x) expMean(as.numeric(x))) 
    print(i)}

### Standardize genes Across Cell Types
centroid.std=(centroid-apply(centroid,1,mean))/apply(centroid,1,sd)
write.table(centroid.std,"AllCellsCCA_centroid_allgenes_std_091219.txt",row.names=T,col.names=T,quote=F,sep="\t") #Can load this again later. 

# set up variables for the heatmap. 
levels=levels(dge1@ident)
 colsep.use=cumsum(table(gsub("_.*","",levels))[levels])
 col.lab=rep("",length(levels))
 col.lab=gsub(".*_","",levels)
 ncluster=length(levels)
 sidecol=matrix(0,2,length(levels))
 sidecol[1,]=rep(rep(c("white","white"),each=12),3)[1:sum(ncluster)]
 sidecol[2,]=myBrewerPalette[1:sum(ncluster)]
 clab=cbind(sidecol[2,],sidecol[1,])
 rlab=sidecol
 rownames(rlab)=c("","Cell Type")
 colnames(clab)=c("Cell Type","")
 redblue100<- scale_color_gradient(low = "blue", high = "red")
 col.use=redblue100
markersall=FindAllMarkers(dge1,only.pos=TRUE,min.pct=0.2,min.diff.pct=0.2,logfc.threshold=log(2),test.use="bimod",do.print = TRUE)
data.use=centroid.std[markersall$gene,]
 row.lab=rownames(data.use)

write.table(markersall, file="all.cells.cca_dim1.3_res.0.3_no5_merged2.3_renamed_101619_markers.txt")

#plot the heatmap
jpeg((file="GermCells_all.cells.cca_Oct162019_final_heatmap.jpeg"),res=300,height=2600,width=1600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use, Colv=F, dendrogram='none', cluster.by.col=F)


######### Lets go back to the germ cells. 

setwd("~/Documents/Hammoud_Lab/learning/raw.data_forGCs_fromCCAof.indv.timepts")
load("GermCells_DirectMergeAll_Jun2019_Final.Robj")

load("GermCells_DirectMergeAll_Jun2019_Final.Robj")
current.cluster.ids <- c(2, 3, 4.5, 6.7, 8, 9, 10, 11,12,13)
new.cluster.ids <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
dge_merg@ident <- plyr::mapvalues(x = dge_merg@ident, from = current.cluster.ids, to = new.cluster.ids)
table(dge_merg@ident)

#1   2   3   4   5   6   7   8   9  10 
#405 444 731 430  63 150 200 362 219 181 

save(all_gcs_merged_final, file="GermCells_DirectMergeAll_Sept2019.Robj")

markersall<-FindAllMarkers(dge_merg,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2), min.diff.pct=0.2,do.print = TRUE)

write.table(markersall, file="GermCells_DirectMergeAll_Sept2019_markers.txt")

### Make a heat map with new labels - 

dge1<-dge_merg

#Add this function
expMean <- function(x) {
  return(log(x = mean(x = exp(x = x) - 1) + 1))} 

# Calculate cluster centroid vals
centroid=matrix(,dim(dge1@data)[1],length(unique(dge1@ident)))
rownames(centroid)=rownames(dge1@data)
colnames(centroid)=new.cluster.ids ##((From above - issue was that the original command was putting the clusters in the wrong order!!!))
for(i in levels(dge1@ident)){
    centroid[,i]=apply(dge1@data[,which(dge1@ident==i)],1,function(x) expMean(as.numeric(x))) 
    print(i)}

### Standardize genes Across Cell Types
centroid.std=(centroid-apply(centroid,1,mean))/apply(centroid,1,sd)


# set up variables for the heatmap. 
levels=levels(dge1@ident)
 colsep.use=cumsum(table(gsub("_.*","",levels))[levels])
 col.lab=rep("",length(levels))
 col.lab=gsub(".*_","",levels)
 ncluster=length(levels)
 sidecol=matrix(0,2,length(levels))
 sidecol[1,]=rep(rep(c("white","white"),each=12),3)[1:sum(ncluster)]
 sidecol[2,]=myBrewerPalette[1:sum(ncluster)]
 clab=cbind(sidecol[2,],sidecol[1,])
 rlab=sidecol
 rownames(rlab)=c("","Cell Type")
 colnames(clab)=c("Cell Type","")
 redblue100<- scale_color_gradient(low = "blue", high = "red")
 col.use=redblue100
markersall=FindAllMarkers(dge1,only.pos=TRUE,min.pct=0.2,min.diff.pct=0.2,logfc.threshold=log(2),test.use="bimod",do.print = TRUE)
data.use=centroid.std[markersall$gene,]
 row.lab=rownames(data.use)

write.table(markersall, file="all.cells.cca_dim1.3_res.0.3_no5_merged2.3_renamed_091219_markers.txt")

#plot the heatmap
 jpeg((file="GermCells_DirectMergeAll_Sept2019_final_heatmap.jpeg"),res=300,height=2600,width=1600)
 par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
 heatmap.3(data.use, Colv=FALSE)
 dev.off()

### 10/16/19

dget<-RenameIdent(dget, old.ident.name = 'Stem Leydig', new.ident.name = 'Progenitor Leydig')

all.cells.cca_final<-dget

save(all.cells.cca_final, file="all.cells.cca_umap_dim1.3_res.0.3_no5_merg2.3_rename_101619_FINALOBJ.Robj")

new_order<-c("Spermatogonia", "Spermatocytes", "Sertoli", "Interstitial Progenitor", "Progenitor Leydig", "Immature Leydig", "Adult Leydig", "Macrophage")

## Note this is the order from the heatmap. Ran heatmap 'blind' then set order based on blocks of genes.

all.cells.cca_final@ident <- factor(all.cells.cca_final@ident, levels = new_order)

v1<-VlnPlot(dget, "nGene", nCol = 1,cols.use=myBrewerPalette, point.size.use = F, size.x.use = 8, x.lab.rot = 90, size.title.use = 16, size.y.use = 8)

v2<-VlnPlot(dget, "nUMI", nCol = 1,cols.use=myBrewerPalette, point.size.use = F, size.x.use = 8, x.lab.rot = 90, size.title.use = 16, size.y.use = 8)

pdf(file='all.cells.cca_NGeneUMI_101619.pdf')
grid.arrange(v1,v2)
dev.off()

## Heatmap

load("all.cells.cca_umap_dim1.3_res.0.3_no5_merg2.3_rename_101619_FINALOBJ.Robj")
new_order<-c("Spermatogonia", "Spermatocytes", "Sertoli", "Interstitial Progenitor", "Progenitor Leydig", "Immature Leydig", "Adult Leydig", "Macrophage")
dget@ident <- factor(dget@ident, levels = new_order)
dge1<-dget
expMean <- function(x) {
+   return(log(x = mean(x = exp(x = x) - 1) + 1))} 
centroid=matrix(,dim(dge1@data)[1],length(unique(dge1@ident)))
rownames(centroid)=rownames(dge1@data)
colnames(centroid)=new_order
for(i in levels(dge1@ident)){
+     centroid[,i]=apply(dge1@data[,which(dge1@ident==i)],1,function(x) expMean(as.numeric(x))) 
+     print(i)}
centroid.std=(centroid-apply(centroid,1,mean))/apply(centroid,1,sd)
levels=levels(dge1@ident)
 colsep.use=cumsum(table(gsub("_.*","",levels))[levels])
 col.lab=rep("",length(levels))
 col.lab=gsub(".*_","",levels)
 ncluster=length(levels)
 sidecol=matrix(0,2,length(levels))
 sidecol[1,]=rep(rep(c("white","white"),each=12),3)[1:sum(ncluster)]
 sidecol[2,]=myBrewerPalette[1:sum(ncluster)]
 clab=cbind(sidecol[2,],sidecol[1,])
 rlab=sidecol
 rownames(rlab)=c("","Cell Type")
 colnames(clab)=c("Cell Type","")
 redblue100<- scale_color_gradient(low = "blue", high = "red")
 col.use=redblue100
markersall=FindAllMarkers(dge1,only.pos=TRUE,min.pct=0.2,min.diff.pct=0.2,logfc.threshold=log(2),test.use="bimod",do.print = TRUE)

data.use=centroid.std[markersall$gene,]
 row.lab=rownames(data.use)
jpeg((file="GermCells_all.cells.cca_Oct162019_final_heatmap.jpeg"),res=300,height=2600,width=1600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use, Colv=F, dendrogram='none', cluster.by.col=F)
dev.off()


## All the cells:
pdf("all.cells.cca_res0.2_dim1.3_no5_merg2.3_renamed_FINALOBJ_101619.pdf")
DimPlot(all.cells.cca_final, reduction.use = 'umap', do.label = T, cols.use = myBrewerPalette)
dev.off()

#Subset only the somatic genes:

somatic<-SubsetData(all.cells.cca_final, ident.use = c('Sertoli', 'Interstitial Progenitor', 'Progenitor Leydig', 'Immature Leydig', 'Adult Leydig', 'Macrophage'))
rm(all.cells.cca_final)

##Somatic cell only plots:

pdf("all.cells.cca_res0.2_dim1.3_umap_only_somatic_101619.pdf")
DimPlot(somatic, reduction.use = 'umap', do.label = T, cols.use = c("#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00"))
dev.off()

pdf("all.cells.cca_res0.2_dim1.3_tsne_only_somatic_101619.pdf")
TSNEPlot(somatic, do.label = T, colors.use = c("#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00"))
dev.off()

pdf("all.cells.cca_res0.2_dim1.3_tsne_only_sertoli_oversomatic_101619.pdf")
TSNEPlot(somatic, do.label = T, colors.use = c("#b2df8a", "grey", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00"))
dev.off()

all.cca<-somatic
p="Sertoli" 
  ident=as.character(all.cca@ident)
  names(ident)=names(all.cca@ident)
  ident[grep(p,ident)]<-p
  ident[which(!grepl(p,ident))]<-"Interstitial Cells"
  ident=factor(ident,levels=c(p,"Interstitial Cells"),order=T)
  all.cca@ident=ident
pdf("Sertoli.all.cells.cca_greybackground_101619.pdf")
DimPlot(all.cca, cols.use=c("grey","#b2df8a"), plot.order='Sertoli', reduction.use='tsne')
dev.off()

rm(all.cca)


setwd("/Volumes/HL/Final R objects for MAH/")
load("all.cells.cca_umap_dim1.3_res.0.3_no5_merg2.3_rename_101619_FINALOBJ.Robj")
somatic<-SubsetData(all.cells.cca_final, ident.use = c('Sertoli', 'Interstitial Progenitor', 'Progenitor Leydig', 'Immature Leydig', 'Adult Leydig', 'Macrophage'))
rm(all.cells.cca_final)

p1<-DimPlot(somatic, reduction.use = ‘umap’, cols.use= myBrewerPalette)
remove.me<-FeatureLocator(p1, data.plot = somatic@dr$tsne@cell.embeddings)
cellnames<-somatic@raw.data@Dimnames[[2]]
keepcells<-setdiff(cellnames,remove.me)

somatic2<-SubsetData(somatic, accept.value = keepcells)
DimPlot(somatic2, reduction.use = 'umap')

save(somatic2,"all.cells.cca_umap_dim1.3_res.0.3_no5_merg2.3_rename_101619_FINALOBJ_FILTEREDSOMATIC.Robj")

scs<-SubsetData(somatic2,ident.use = 'Sertoli', subset.raw = T)
scs<-SetAllIdent(scs,id='stim')
markers<-FindAllMarkers(scs,logfc.threshold = log(2),min.diff.pct = .2,only.pos = T)

markers<-FindAllMarkers(scs,logfc.threshold = log(1.5),min.diff.pct = .2,only.pos = T)

write.table(centroids,file="sertoli_cell_cluster_centroids_from_allcca_101619_FILTEREDFINALOBJ.txt")
write.table(centroids,file="sertoli_cell_cluster_centroids_from_allcca_120619_FILTEREDFINALOBJ.txt")

while read file; do grep -w "$file" sertoli_cell_cluster_centroids_from_allcca_120619_FILTEREDFINALOBJ.txt; done < SC_markersall_Oct2019_WTvSCAR_geneIDs.txt > scs_WTvSCAR_upreggenes_centroidexp_120619.txt

save(scs, file="scs_final.all.cells.cca_res0.2_dim1.3_no5_merg2.3_FINALOBJ_FILTERED_101619.Robj")

markers<-FindAllMarkers(scs,logfc.threshold = log(1.5),min.diff.pct = .2,only.pos = T)
centroids<-log(AverageExpression(scs)+1)


##Somatic cell only plots:

pdf("all.cells.cca_res0.2_dim1.3_umap_only_somatic_filtered_101619.pdf")
DimPlot(somatic2, reduction.use = 'umap', do.label = T, cols.use = myBrewerPalette[3:12])
dev.off()

pdf("all.cells.cca_res0.2_dim1.3_tsne_only_somatic_filtered_101619.pdf")
TSNEPlot(somatic2, do.label = T, colors.use = myBrewerPalette[3:12])
dev.off()

all.cca<-somatic2
p="Sertoli" 
  ident=as.character(all.cca@ident)
  names(ident)=names(all.cca@ident)
  ident[grep(p,ident)]<-p
  ident[which(!grepl(p,ident))]<-"Interstitial Cells"
  ident=factor(ident,levels=c(p,"Interstitial Cells"),order=T)
  all.cca@ident=ident
pdf("all.cells.cca_res0.2_dim1.3_tsne_only_sertoli_oversomatic_filtered_101619.pdf")
DimPlot(all.cca, cols.use=c("grey","#b2df8a"), plot.order='Sertoli', reduction.use='tsne')
dev.off()

v1<-VlnPlot(somatic2, "nGene", nCol = 1,cols.use=myBrewerPalette, point.size.use = F, size.x.use = 8, x.lab.rot = 90, size.title.use = 16, size.y.use = 8)

v2<-VlnPlot(somatic2, "nUMI", nCol = 1,cols.use=myBrewerPalette, point.size.use = F, size.x.use = 8, x.lab.rot = 90, size.title.use = 16, size.y.use = 8)

pdf(file='all.cells.cca_somatic_NGeneUMI_101619_filtered.pdf')
grid.arrange(v1,v2)
dev.off()

pdf("all.cells.cca_rhox8sox9_somatic_101619.pdf", height=3)
VlnPlot(somatic2, features.plot = c("Sox9", "Rhox8"), point.size.use = -1, x.lab.rot = 90, size.x.use = 8, size.y.use = 8, size.title.use = 12, cols.use = myBrewerPalette[3:12])
dev.off()

##### 02/25/20

Recluster somatic cells and call markers, again... 

setwd("/Volumes/HL/Final R objects for MAH/")
load("all.cells.cca_umap_dim1.3_res.0.3_no5_merg2.3_rename_101619_FINALOBJ_FILTEREDSOMATIC.Robj")

somatic2 <- FindClusters(somatic2, reduction.type = "cca.aligned", resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), dims.use = 1:3, force.recalc = T)

print(c( length(unique(somatic2@meta.data$res.0.1)),length(unique(somatic2@meta.data$res.0.2)),length(unique(somatic2@meta.data$res.0.3)),length(unique(somatic2@meta.data$res.0.4)),length(unique(somatic2@meta.data$res.0.5)),length(unique(somatic2@meta.data$res.0.6)),length(unique(somatic2@meta.data$res.0.7)),length(unique(somatic2@meta.data$res.0.8)),length(unique(somatic2@meta.data$res.0.9)),length(unique(somatic2@meta.data$res.1))))

 [1]  4  5  8 11 12 12 13 17 15 19
 
somatic2_rc <- FindClusters(somatic2, reduction.type = "cca.aligned", resolution = c(0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29), dims.use = 1:3, force.recalc = T)


print(c(length(unique(somatic2_rc@meta.data$res.0.2)),length(unique(somatic2_rc@meta.data$res.0.21)),length(unique(somatic2_rc@meta.data$res.0.22)),length(unique(somatic2_rc@meta.data$res.0.23)),length(unique(somatic2_rc@meta.data$res.0.24)),length(unique(somatic2_rc@meta.data$res.0.25)),length(unique(somatic2_rc@meta.data$res.0.26)),length(unique(somatic2_rc@meta.data$res.0.27)),length(unique(somatic2_rc@meta.data$res.0.28)), length(unique(somatic2_rc@meta.data$res.0.29))))

[1] 6 6 6 7 7 7 7 7 7 7

somatic2_rc <- FindClusters(somatic2, reduction.type = "cca.aligned", resolution = c(0.21, 0.22, 0.23), dims.use = 1:3, force.recalc = T)

print(c(length(unique(somatic2_rc@meta.data$res.0.21)),length(unique(somatic2_rc@meta.data$res.0.22)),length(unique(somatic2_rc@meta.data$res.0.23))))

$$$%%%% How does this compare to the somatic2 original cluster? 

somatic2_rc<-SetAllIdent(somatic2_rc, id='res.0.22')
pdf(file="all.cells.cca_tsne_dim1.3_res.0.22_090319.pdf")
TSNEPlot(all.cells.cca, do.label = T, do.return = T, pt.size = 0.5, colors.use= myBrewerPalette)
dev.off()

markersall=FindAllMarkers(somatic2_rc,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)

write.table(markersall, file="all.cells.cca_umap_dim1.3_res.0.3_no5_merg2.3_rename_101619_FINALOBJ_FILTEREDSOMATIC_RECLUSTERED_022820_MARKERS.txt")


save(somatic2_rc,"all.cells.cca_umap_dim1.3_res.0.3_no5_merg2.3_rename_101619_FINALOBJ_FILTEREDSOMATIC_reclustered_022920.Robj")


@@@@@@@@@@@ ALL CELLS CCA EDITS MARCH 2020  

setwd("/Volumes/HL/Final R objects for MAH/")
load("all.cells.cca_umap_dim1.3_res.0.3_no5_merg2.3_rename_101619_FINALOBJ.Robj")
new_order<-c("Spermatogonia", "Spermatocytes", "Sertoli", "Interstitial Progenitor", "Progenitor Leydig", "Immature Leydig", "Adult Leydig", "Macrophage")
all.cells.cca_final@ident <- factor(all.cells.cca_final@ident, levels = new_order)
save(all.cells.cca_final, file="all.cells.cca_umap_dim1.3_res.0.3_no5_merg2.3_rename_101619_FINALOBJ_merg.Leydig_030120.Robj")

# Recalled markers.... 

#Rename/Combined the Progenitor Leydig cluster to Immature Leydig cluster. 

load("all.cells.cca_umap_dim1.3_res.0.3_no5_merg2.3_rename_101619_FINALOBJ_merg.Leydig_030120.Robj")
all.cells.cca_2<-RenameIdent(all.cells.cca_2, old.ident.name = 'Progenitor Leydig', new.ident.name = 'Immature Leydig')
new_order<-c("Spermatogonia", "Spermatocytes", "Sertoli", "Interstitial Progenitor", "Immature Leydig", "Adult Leydig", "Macrophage")
all.cells.cca_2@ident <- factor(all.cells.cca_2@ident, levels = new_order)
save(all.cells.cca_2, file="all.cells.cca_umap_dim1.3_res.0.3_no5_merg2.3_rename_101619_FINALOBJ_merg.Leydig_031320.Robj")

#### NOTE: rename interstitial progenitor to interstitial fibroblast (tcf21+) before making figures. 

new<-RenameIdent(all.cells.cca_2, old.ident.name = 'Interstitial Progenitor', new.ident.name = 'Interstitial Fibroblast (Tcf21+)')

markersall<-FindAllMarkers(new,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)

write.table(markersall, file="all.cells.cca_umap_dim1.3_res.0.3_no5_merg2.3_rename_101619_FINALOBJ_merg.Leydig_030920_MARKERS.txt")
