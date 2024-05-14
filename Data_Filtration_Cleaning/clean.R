#### CLEAN
#### Code for publication unifying all 'clean' scripts
#### Replicable method for taking filtered objects to fully cleaned objects




## Note: For cleaning, we grouped by 10x chemistry and sample similarity (cellular diversity),
##  and attempted minimally-invasive grouped clustering, including merge, ComBat, or CCA integration,
##  in that order. Once sufficient sample alignment was reached, the grouped objects were cleaned of
##  partial cells and doublets by manual cluster removal. Samples with no comparable counterparts were
##  cluster-cleaned alone.

## Final object groupings and alignment methods:
# 'BC1P3','BC1P6' = StartingBC (v2) - Integration
# 'FB13','FB14' = StartingFB (v2) - ComBat
# 'BAL' (v3) = MacStart - Alone
# 'RLMVEC' (v3) = ECStart - Alone
# 'BEF1','BEF2','BCL5' = v2EngLungs - Integration
# 'BCEC2','BEF3','BEF12','BEF14','BEF15' = v3EngLungs - Integration
# 'BEFM1','BEFM2','BEFM4','BEFM5','BEFM6' = v3BEFMLungs - Integration




#### Cleaning objects by combining and removing trash clusters
# Load filtered objects
load('./bef.metadata.filter1.2022-06-09.Robj')
load('./bef.seurat.objs.filter1.2022-06-09.Robj')




#### ComBat batch correction - FB13 & FB14
# Setup sample combinations
sample.names <- c('FB13','FB14')
dataset <- c('StartingFB')
dir_name <- paste(sample.names,collapse='_')

# Merge specified objects
objects <- new.seurat.data[sample.names]
object <- merge(x = objects[[1]], y = objects[2:length(objects)])
object$Dataset <- dataset

# Implementing ComBat
count_matrix <- as.matrix(GetAssayData(object = object, slot = "counts"))
adjusted <- ComBat_seq(count_matrix, batch=object$Sample, group=NULL)

object[['ComBat']] <- CreateAssayObject(counts = adjusted)
DefaultAssay(object) <- 'ComBat'

# Process & plot objects together
object <- NormalizeData(object)
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, npcs = 50, verbose = F)
png(paste0(dir_name,'_combat_elbow.png'), width = 800, height = 500)
ElbowPlot(object,ndims = 50)+labs(title = dir_name)
dev.off()
pdf(paste0(dir_name,'_combat_pcheatmap1.pdf'), width = 12, height = 16)
DimHeatmap(object, dims = 1:15, cells = 200, balanced = T)
dev.off()
pdf(paste0(dir_name,'_combat_pcheatmap2.pdf'), width = 12, height = 16)
DimHeatmap(object, dims = 16:30, cells = 200, balanced = T)
dev.off()

pcs <- 25
res <- 0.5
tag <- paste('.pcs',pcs,'res',res,sep='.')
object <- FindNeighbors(object, dims = 1:pcs)
object <- FindClusters(object, resolution = res)
object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
p1 <- UMAPPlot(object,label = T) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) +
        theme(plot.title = element_text(size = 24),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoLegend() + NoAxes()
p2 <- UMAPPlot(object,group.by = 'Sample') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))

# histogram of proportion sample per cluster
nclusters <- length(levels(object$seurat_clusters))
object$Sample <- factor(object$Sample)
object$seurat_clusters <- factor(object$seurat_clusters)
data = table(object$Sample,object$seurat_clusters)
fill_col <- hue_pal()(length(objects))
p3 <- barplot(data,main=dir_name,xlab='Cluster',ylab='nCells',col=rep(fill_col,nclusters),beside=T)
        legend("topright",sample.names,fill=fill_col) 

# QC checks
p4 <- FeaturePlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','perc.spliced','Mki67','Top2a'),ncol=2,label = T,repel=T)
p5 <- VlnPlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','perc.spliced','Mki67','Top2a'),ncol=2,pt.size=0.1)

# Backtrack to original Filtration plots
# ViolinPlots whole sample, jitter colored by cluster
p <- VlnPlot(object, features = c("nCount_RNA"),pt.size=0,log=T,group.by='Dataset',cols=1) +
    geom_jitter(aes(color = factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    theme(legend.position = 'none',axis.title.x = element_blank())
q <- VlnPlot(object, features = c("nFeature_RNA"),pt.size=0,log=T,group.by='Dataset',cols=1) +
    geom_jitter(aes(color = factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    theme(legend.position = 'none',axis.title.x = element_blank())
r <- VlnPlot(object, features = c("percent.mt"),pt.size=0,group.by='Dataset',cols=1,y.max=(max(object$percent.mt)+0.5)) +
    geom_jitter(aes(color = factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    theme(legend.position = 'none',axis.title.x = element_blank())
p6 <- plot_grid(p,q,r,ncol=3)

# Lineage checks
p7 <- FeaturePlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'), label = T,repel=T)
p8 <- VlnPlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'),ncol=2,pt.size=0.1)

# Scatter Plot colored by cluster
list <- subset(bef.metadata,Sample %in% sample.names)
p9 <- ggplot() +
    geom_point(data=list, aes(x=nCount_RNA, y=nFeature_RNA, color=factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    labs(title = dir_name) +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
        plot.subtitle = element_text(color='red'),
        plot.caption = element_text(color='red'),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size = 12,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Assemble plots into a cowplot
plots <- list(p1,p2,p3,p4,p5,p6,p7,p8,p9)
ht <- ceiling(length(plots)/3)*600
png(paste0(dir_name,tag,'_combat_allplots.png'), width = 2000,height=ht)
plot_grid(plotlist=plots,ncol=3)
dev.off()

save(object,file = paste0(dir_name,".combat.",Sys.Date(),".Robj"))




#### Seurat CCA Integration - BC1P3/BC1P6, BEF1/BEF2/BCL5,
##      BCEC2/BEF3/BEF12/BEF14/BEF15, BEFM1/BEFM2/BEFM4/BEFM5/BEFM6
# Setup sample combinations
sample.names <- c('BC1P3','BC1P6')
dataset <- c('StartingBC')
dir_name <- paste(sample.names,collapse='_')

# Normalize and identify variable features for each object independently
objects <- new.seurat.data[sample.names]
objects <- lapply(X = objects, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x,nfeatures = 2000)
})

# Select variable features across datasets for integration
ngenes <- 2000
dims <- 30
tag <- paste('.nfts',ngenes,'dims',dims,sep='.')
features <- SelectIntegrationFeatures(object.list = objects, nfeatures = ngenes)

# Perform integration
object.anchors <- FindIntegrationAnchors(object.list=objects, anchor.features=features,
        reduction='cca',dims=1:dims)
object <- IntegrateData(anchorset = object.anchors, dims=1:dims)

# Plot
DefaultAssay(object) <- "integrated"
object$Dataset <- dataset
object <- ScaleData(object)
object <- RunPCA(object, npcs = 100, verbose = F)
png(paste(dir_name,tag,'_integrated_elbow.png',sep=""), width = 800, height = 500)
ElbowPlot(object,ndims = 50)+labs(title = dir_name)
dev.off()
pdf(paste(dir_name,tag,'_integrated_pcheatmap1.pdf',sep=""), width = 12, height = 16)
DimHeatmap(object, dims = 1:15, cells = 200, balanced = T)
dev.off()
pdf(paste(dir_name,tag,'_integrated_pcheatmap2.pdf',sep=""), width = 12, height = 16)
DimHeatmap(object, dims = 16:30, cells = 200, balanced = T)
dev.off()

pcs <- 15
res <- 0.5
tag2 <- paste(tag,'pcs',pcs,'res',res,sep='.')
object <- FindNeighbors(object, dims = 1:pcs)
object <- FindClusters(object, resolution = res)
object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
p1 <- UMAPPlot(object,label = T) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) +
        theme(plot.title = element_text(size = 24),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoLegend() + NoAxes()
p2 <- UMAPPlot(object,group.by = 'Sample') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste('nfeatures = ',ngenes,', & dims = ',dims,sep='')) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red'))

# histogram of proportion sample per cluster
nclusters <- length(levels(object$seurat_clusters))
object$Sample <- factor(object$Sample)
object$seurat_clusters <- factor(object$seurat_clusters)
data = table(object$Sample,object$seurat_clusters)
fill_col <- hue_pal()(length(objects))
p3 <- barplot(data,main=dir_name,xlab='Cluster',ylab='nCells',col=rep(fill_col,nclusters),beside=T)
        legend("topright",sample.names,fill=fill_col)

# QC checks
DefaultAssay(object) <- "RNA"
p4 <- FeaturePlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','perc.spliced','Mki67','Top2a'),ncol=2,label = T,repel=T)
p5 <- VlnPlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','perc.spliced','Mki67','Top2a'),ncol=2,pt.size=0.1)

# Backtrack to original Filtration plots
# ViolinPlots whole sample, jitter colored by cluster
p <- VlnPlot(object, features = c("nCount_RNA"),pt.size=0,log=T,group.by='Dataset',cols=1) +
    geom_jitter(aes(color = factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    theme(legend.position = 'none',axis.title.x = element_blank())
q <- VlnPlot(object, features = c("nFeature_RNA"),pt.size=0,log=T,group.by='Dataset',cols=1) +
    geom_jitter(aes(color = factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    theme(legend.position = 'none',axis.title.x = element_blank())
r <- VlnPlot(object, features = c("percent.mt"),pt.size=0,group.by='Dataset',cols=1,y.max=(max(object$percent.mt)+0.5)) +
    geom_jitter(aes(color = factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    theme(legend.position = 'none',axis.title.x = element_blank())
p6 <- plot_grid(p,q,r,ncol=3)

# Lineage checks
p7 <- FeaturePlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'), label = T,repel=T)
p8 <- VlnPlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'),ncol=2,pt.size=0.1)

# Scatter Plot colored by cluster
list <- subset(bef.metadata,Sample %in% sample.names)
p9 <- ggplot() +
    geom_point(data=list, aes(x=nCount_RNA, y=nFeature_RNA, color=factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    labs(title = dir_name) +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
        plot.subtitle = element_text(color='red'),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size = 12,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Assemble plots into a cowplot
plots <- list(p1,p2,p3,p4,p5,p6,p7,p8,p9)
ht <- ceiling(length(plots)/3)*600
png(paste0(dir_name,tag2,'_integrated_allplots.png'), width = 2000,height=ht)
plot_grid(plotlist=plots,ncol=3)
dev.off()




#### Clean individual objects - MacAlv & RLMVEC
# Setup sample to clean
sample.name <- c('RLMVEC')
dataset <- c('StartingEC')

# Pull object & dimensional reduction
x <- which(sample.name == names(new.seurat.data))
object <- new.seurat.data[[x]]
object$Dataset <- dataset
object <- NormalizeData(object)
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, npcs = 50, verbose = F)
png(paste0(sample.name,'_elbow.png'), width = 800, height = 500)
ElbowPlot(object,ndims = 50)+labs(title = sample.name)
dev.off()
pdf(paste0(sample.name,'_pcheatmap1.pdf'), width = 12, height = 16)
DimHeatmap(object, dims = 1:15, cells = 200, balanced = T)
dev.off()
pdf(paste0(sample.name,'_pcheatmap2.pdf'), width = 12, height = 16)
DimHeatmap(object, dims = 16:30, cells = 200, balanced = T)
dev.off()

pcs <- 15
res <- 0.5
tag <- paste('.pcs',pcs,'res',res,sep='.')
object <- FindNeighbors(object, dims = 1:pcs)
object <- FindClusters(object, resolution = res)
object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
p1 <- UMAPPlot(object,label = T) + labs(title = sample.name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) +
        theme(plot.title = element_text(size = 24),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoLegend() + NoAxes()
p2 <- UMAPPlot(object,group.by = 'Sample') + labs(title = sample.name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))

# histogram of proportion sample per cluster
nclusters <- length(levels(object$seurat_clusters))
object$seurat_clusters <- factor(object$seurat_clusters)
data = table(object$Sample,object$seurat_clusters)
p3 <- barplot(data,main=sample.name,xlab='Cluster',ylab='nCells')

# QC checks
p4 <- FeaturePlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','perc.spliced','Mki67','Top2a'),ncol=2,label = T,repel=T)
p5 <- VlnPlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','perc.spliced','Mki67','Top2a'),ncol=2,pt.size=0.1)

# Backtrack to original Filtration plots
# ViolinPlots not broken apart, jitter points colored by cluster
p <- VlnPlot(object, features = c("nCount_RNA"),pt.size=0,log=T,group.by='Dataset',cols=1) +
    geom_jitter(aes(color = factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    theme(legend.position = 'none',axis.title.x = element_blank())
q <- VlnPlot(object, features = c("nFeature_RNA"),pt.size=0,log=T,group.by='Dataset',cols=1) +
    geom_jitter(aes(color = factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    theme(legend.position = 'none',axis.title.x = element_blank())
r <- VlnPlot(object, features = c("percent.mt"),pt.size=0,group.by='Dataset',cols=1,y.max=(max(object$percent.mt)+0.5)) +
    geom_jitter(aes(color = factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    theme(legend.position = 'none',axis.title.x = element_blank())
p6 <- plot_grid(p,q,r,ncol=3)

# Lineage checks
p7 <- FeaturePlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'), label = T,repel=T)
p8 <- VlnPlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'),ncol=2,pt.size=0.1)

# Scatter Plot colored by cluster
list <- subset(bef.metadata,Sample %in% sample.name)
p9 <- ggplot() +
    geom_point(data=list, aes(x=nCount_RNA, y=nFeature_RNA, color=factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    labs(title = sample.name) +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
        plot.subtitle = element_text(color='red'),
        plot.caption = element_text(color='red'),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size = 12,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())


## Assemble plots into a cowplot
plots <- list(p1,p2,p3,p4,p5,p6,p7,p8,p9)
ht <- ceiling(length(plots)/3)*600
png(paste0(sample.name,tag,'_allplots.png'), width = 2000,height=ht)
plot_grid(plotlist=plots,ncol=3)
dev.off()




#### Generate supplemental figure to visualize object cleaning
# Setup sample combinations
sample.names <- c('BEFM1','BEFM2','BEFM4','BEFM5','BEFM6')
dataset <- c('v3BEFMLungs')
cleantype <- c('Integration')
dir_name <- paste(sample.names,collapse='_')
objects <- new.seurat.data[sample.names]

# Alone -
object <- objects[[1]]
object <- NormalizeData(object)
object <- FindVariableFeatures(object)
tag <- c('')

# ComBat -
load('./FB13_FB14/FB13_FB14.combat.2022-06-13.Robj')
DefaultAssay(object) <- 'ComBat'
tag <- c('')

# Integration -
objects <- lapply(X = objects, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x)
})
features <- SelectIntegrationFeatures(object.list = objects)
object.anchors <- FindIntegrationAnchors(object.list=objects, anchor.features=features,reduction='cca')
object <- IntegrateData(anchorset = object.anchors)
DefaultAssay(object) <- "integrated"
tag <- c('.nfts.2000.dims.30')

# all -
object$Dataset <- dataset
object$Cleantype <- cleantype
object <- ScaleData(object)
object <- RunPCA(object, npcs = 100, verbose = F)
png(paste(dir_name,'_elbow.png',sep=""), width = 800, height = 500)
ElbowPlot(object,ndims = 50)+labs(title = dir_name)
dev.off()
pdf(paste(dir_name,'_pcheatmap1.pdf',sep=""), width = 12, height = 16)
DimHeatmap(object, dims = 1:15, cells = 200, balanced = T)
dev.off()
pdf(paste(dir_name,'_pcheatmap2.pdf',sep=""), width = 12, height = 16)
DimHeatmap(object, dims = 16:30, cells = 200, balanced = T)
dev.off()

# Relevant QC metrics for identifying trash clusters
pcs <- 50
res <- 0.85
tag2 <- paste(tag,'pcs',pcs,'res',res,sep='.')
object <- FindNeighbors(object, dims = 1:pcs)
object <- FindClusters(object, resolution = res)
object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
p1 <- UMAPPlot(object,label = T,label.size = 6) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) +
        theme(plot.title = element_text(size = 24),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoLegend()
p2 <- UMAPPlot(object,group.by = 'Sample') + labs(title = dir_name) +
        theme(plot.title = element_text(size = 24,hjust = 0)) +
        # Alone
        #labs(subtitle = c("Clean Alone")) + theme(plot.subtitle = element_text(size = 20))
        # ComBat
        #labs(subtitle = c("ComBat")) + theme(plot.subtitle = element_text(size = 20))
        # Integration
        labs(subtitle = c("nfeatures = 2000\ndims = 30")) + theme(plot.subtitle = element_text(size = 20))

data <- prop.table(table(object$seurat_clusters, object$Sample),1)*100
data <- as.data.frame(data)
p3 <- ggplot(data, aes(x = Var1, y = Freq, fill = Var2)) + NoLegend() +
  geom_col(position = "dodge", width = 0.5) + scale_fill_manual(values = hue_pal()(length(objects))) +
  labs(title = c('Cluster Composition by Sample')) + xlab("Cluster") + ylab("% cluster per sample") +
  theme(plot.title = element_text(size = 24),
        legend.title = element_blank(),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
p3 <- plot_grid(NULL,p3,NULL,nrow=3,rel_heights = c(1,2,1))

toprow <- plot_grid(p1,p2,p3,ncol=3,labels = c('A','B','C'),label_size = 40)

DefaultAssay(object) <- "RNA"
p4 <- FeaturePlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','Mki67'),ncol=2,label = T,repel=T)
p5 <- VlnPlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','Mki67'),ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())
midrow <- plot_grid(p4,p5,ncol=2,labels = c('D','E'),rel_widths = c(1,2),label_size = 40)

p6 <- FeaturePlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'), label = T,repel=T)
p7 <- VlnPlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'),ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())
midrow2 <- plot_grid(p6,p7,ncol=2,labels = c('F','G'),rel_widths = c(1,2),label_size = 40)

# Identify clusters to be removed
trash <- c(1,4,8,12,16,19)
contaminant <- c()  # cells of unexpected class
new.labels <- levels(object$seurat_clusters)
for (i in 1:length(trash)){
    x <- which(new.labels == trash[i])
    new.labels[x] <- c('Trash')
}
for (i in 1:length(contaminant)){
    x <- which(new.labels == contaminant[i])
    new.labels[x] <- c('Contaminant')
}
names(new.labels) <- levels(object$seurat_clusters)
object <- RenameIdents(object, new.labels)
p8 <- UMAPPlot(object,cols = c('Trash' = 'red','Contaminant' = 'dodgerblue')) + NoAxes() +
        labs(title = c('Clusters to Remove'),caption = paste("Remove:",paste(trash,collapse=', ')))  +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.caption = element_text(size = 18,color='red'))
p8 <- plot_grid(NULL,p8,nrow=2,rel_heights = c(1,9))

# Cluster marker characterization
Idents(object) <- object$seurat_clusters
myplots <- list()
for (i in 1:length(trash)){
    x <- trash[i]
    cluster.markers <- FindMarkers(object, ident.1 = trash[i],only.pos = T)
    cluster.markers$ratio <- cluster.markers$pct.1/cluster.markers$pct.2
    cluster.markers <- cluster.markers[order(-cluster.markers$ratio),]
    marks <- cluster.markers[1:8, ]
    marks[, 1:6] <- round(marks[, 1:6], digits = 4)
    marks$cluster <- rep(trash[i],nrow(marks))
    myplots[[i]] <- as.data.frame(marks)
}
length(myplots)

# plot marker lists
p9 <- tableGrob(myplots[[1]],theme = ttheme_minimal(base_size = 16))
p10 <- tableGrob(myplots[[2]],theme = ttheme_minimal(base_size = 16))

p11 <- tableGrob(myplots[[3]],theme = ttheme_minimal(base_size = 16))
p12 <- tableGrob(myplots[[4]],theme = ttheme_minimal(base_size = 16))
p13 <- tableGrob(myplots[[5]],theme = ttheme_minimal(base_size = 16))
p14 <- tableGrob(myplots[[6]],theme = ttheme_minimal(base_size = 16))
listplot <- plot_grid(p9,p10,p11,p12,p13,p14,nrow=3,labels = c('I','J','K','L','M','N'),label_size = 40)

lastrow <- plot_grid(p8,listplot,ncol=2,labels = c('H',''),label_size = 40,rel_widths = c(1,2))

# Cowplot all panels
plots <- list(toprow,midrow,midrow2,lastrow)
png(paste0(dir_name,tag2,'_suppfig.png'), width = 2000,height=2750)
plot_grid(plotlist=plots,ncol=1)
dev.off()




#### Perform first round of cleaning by cluster removal
## Recreate clean-able object (Alone, ComBat, Integration)
# Setup sample combinations
sample.names <- c('BEFM1','BEFM2','BEFM4','BEFM5','BEFM6')
dataset <- c('v3BEFMLungs')
cleantype <- c('Integration')
chemistry <- c('v3')
dir_name <- paste(sample.names,collapse='_')
objects <- new.seurat.data[sample.names]

# Alone -
object <- objects[[1]]
object <- NormalizeData(object)
object <- FindVariableFeatures(object)
tag <- c('')

# ComBat -
load('./FB13_FB14.combat.2022-06-13.Robj')
DefaultAssay(object) <- 'ComBat'
tag <- c('')

# Integration -
objects <- lapply(X = objects, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x)
})
features <- SelectIntegrationFeatures(object.list = objects)
object.anchors <- FindIntegrationAnchors(object.list=objects, anchor.features=features,reduction='cca')
object <- IntegrateData(anchorset = object.anchors)
DefaultAssay(object) <- "integrated"
tag <- c('.nfts.2000.dims.30')

# all -
object$Dataset <- dataset
object$Cleantype <- cleantype
object$Chemistry <- chemistry
object <- ScaleData(object)
object <- RunPCA(object, npcs = 100, verbose = F)

pcs <- 50
res <- 0.85
tag2 <- paste(tag,'pcs',pcs,'res',res,sep='.')
object <- FindNeighbors(object, dims = 1:pcs)
object <- FindClusters(object, resolution = res)
object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)

# Identify clusters to be removed
trash <- c(1,4,8,12,16,19)
contaminant <- c()
new.labels <- levels(object$seurat_clusters)
for (i in 1:length(trash)){
    x <- which(new.labels == trash[i])
    new.labels[x] <- c('Trash')
}
for (i in 1:length(contaminant)){
    x <- which(new.labels == contaminant[i])
    new.labels[x] <- c('Contaminant')
}
names(new.labels) <- levels(object$seurat_clusters)
object <- RenameIdents(object, new.labels)
object$rnd1_contaminants <- Idents(object)
object

# Remove trash, preserve contaminants, re-cluster
object <- subset(object, idents = c('Trash'), invert = T)

# Alone -
DefaultAssay(object) <- 'RNA'
object <- NormalizeData(object)
object <- FindVariableFeatures(object)
tag <- c('')

# ComBat - re-run ComBat without trash clusters
count_matrix <- as.matrix(GetAssayData(object = object, slot = "counts"))
adjusted <- ComBat_seq(count_matrix, batch=object$Sample, group=NULL)  # takes a long time...
object[['ComBat']] <- CreateAssayObject(counts = adjusted)
DefaultAssay(object) <- 'ComBat'
object <- NormalizeData(object)
object <- FindVariableFeatures(object)
tag <- c('')

# Integration - re-run integration without trash clusters
objects <- SplitObject(object, split.by = "Sample")
features <- SelectIntegrationFeatures(object.list = objects)
object.anchors <- FindIntegrationAnchors(object.list=objects, anchor.features=features,reduction='cca')
object <- IntegrateData(anchorset = object.anchors)
DefaultAssay(object) <- "integrated"
tag <- c('.nfts.2000.dims.30')

# all -
object <- ScaleData(object)
object <- RunPCA(object, npcs = 100, verbose = F)

pcs <- 50
res <- 0.85
tag2 <- paste(tag,'pcs',pcs,'res',res,sep='.')
object <- FindNeighbors(object, dims = 1:pcs)
object <- FindClusters(object, resolution = res)
object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
p1 <- UMAPPlot(object,label = T,label.size = 6) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) +
        theme(plot.title = element_text(size = 24),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoLegend()
p2 <- UMAPPlot(object,group.by = 'Sample') + labs(title = dir_name) +
        theme(plot.title = element_text(size = 24,hjust = 0)) +
        # Alone
        #labs(subtitle = c("Clean Alone")) + theme(plot.subtitle = element_text(size = 20))
        # ComBat
        #labs(subtitle = c("ComBat")) + theme(plot.subtitle = element_text(size = 20))
        # Integration
        labs(subtitle = c("nfeatures = 2000\ndims = 30")) + theme(plot.subtitle = element_text(size = 20))

data <- prop.table(table(object$seurat_clusters, object$Sample),1)*100
data <- as.data.frame(data)
p3 <- ggplot(data, aes(x = Var1, y = Freq, fill = Var2)) + NoLegend() +
  geom_col(position = "dodge", width = 0.5) + scale_fill_manual(values = hue_pal()(length(objects))) +
  labs(title = c('Cluster Composition by Sample')) + xlab("Cluster") + ylab("% cluster per sample") +
  theme(plot.title = element_text(size = 24),
        legend.title = element_blank(),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
p3 <- plot_grid(NULL,p3,NULL,nrow=3,rel_heights = c(1,2,1))

toprow <- plot_grid(p1,p2,p3,ncol=3,labels = c('A','B','C'),label_size = 40)

DefaultAssay(object) <- "RNA"
p4 <- FeaturePlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','Mki67'),ncol=2,label = T,repel=T)
p5 <- VlnPlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','Mki67'),ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())
midrow <- plot_grid(p4,p5,NULL,ncol=3,labels = c('D','E'),rel_widths = c(2,2,1),label_size = 40)

p6 <- FeaturePlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'), label = T,repel=T)
p7 <- VlnPlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'),ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())
midrow2 <- plot_grid(p6,p7,NULL,ncol=3,labels = c('F','G'),rel_widths = c(2,2,1),label_size = 40)

# Visualize remaining contaminants
p8 <- UMAPPlot(object,group.by = 'rnd1_contaminants',cols = c('Contaminant' = 'dodgerblue')) + NoAxes() +
        labs(title = c('Remaining Contaminant Populations')) + theme(plot.title = element_text(size = 24,hjust = 0))
p8 <- plot_grid(NULL,p8,nrow=2,rel_heights = c(1,9))

lastrow <- plot_grid(p8,NULL,NULL,ncol=3,labels = c('H'),label_size = 40)

# Cowplot all panels
plots <- list(toprow,midrow,midrow2,lastrow)
png(paste0(dir_name,tag2,'_clean_allplots.png'), width = 2000,height=2750)
plot_grid(plotlist=plots,ncol=1)
dev.off()

# Save combined & cleaned object
save(object,file = paste(dir_name,cleantype,"clean1",Sys.Date(),"Robj",sep="."))




#### Split, finalize and save each individual object
## Generate final cleaning supplement
## Separate to individual, "clean" objects, and confirm/characterize
# Pull combined sample objects
sample.names <- c('BEFM1','BEFM2','BEFM4','BEFM5','BEFM6')
dir_name <- paste(sample.names,collapse='_') ; dir_name
load(list.files(path = ".", pattern = "clean1"))

# Split and characterize individually
# Epi: ['Krt5','Hopx'] | 'Krt5','Sftpc','Hopx','Aqp5' ; 'Sftpc','Sftpb','Defb4','Lyz2'
# Endo: ['Lyve1','Prx'] | 'Lyve1','Prx','Plvap','Apln'
# Mes: ['Dcn','Acta2'] | 'Dcn','Acta2','Itga8','Msln'
# Imm: 'Mrc1','Naaa','Cd14','Cd3e'
# BEF: 'Krt5','Hopx','Lyve1','Prx','Dcn','Acta2' | BEFM: 'Krt5','Hopx','Lyve1','Dcn','Mrc1','Naaa'
objects <- SplitObject(objects, split.by = "Sample")

pcs <- 40
res <- 0.5
tag <- paste('.pcs',pcs,'res',res,sep='.')
goi <- c('Krt5','Hopx','Lyve1','Dcn','Mrc1','Naaa')
plots <- list() ; new.data <- list()
for (i in 1:length(objects)){
    object <- objects[[i]]
    sample.name <- names(objects[i])
    message(sample.name)
    object <- NormalizeData(object)
    object <- FindVariableFeatures(object)
    object <- ScaleData(object)
    object <- RunPCA(object, npcs = 100, verbose = F)
    object <- FindNeighbors(object, dims = 1:pcs)
    object <- FindClusters(object, resolution = res)
    object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)

    p1 <- UMAPPlot(object,label = T,label.size = 6) + labs(title = sample.name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) +
            theme(plot.title = element_text(size = 24),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoLegend()
    p2 <- FeaturePlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','Mki67'),ncol=2,label = T,repel=T)
    p3 <- VlnPlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','Mki67'),ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())
    p4 <- FeaturePlot(object, goi, label = T,repel=T)
    p5 <- FeaturePlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'), label = T,repel=T)
    p6 <- VlnPlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'),ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())
    plots[[i]] <- list(p1,p2,p3,p4,p5,p6)
    new.data[[i]] <- object
}
names(new.data) <- sample.names

for (i in 1:length(objects)){
    object <- new.data[[i]]
    sample.name <- names(objects[i])
    plotlist <- plots[[i]]
    message(sample.name)
    setwd(paste0("../",sample.name))

    png(paste0(sample.name,tag,'_cleancheck.png'), width = 2000,height=1375)  # 2750
    print(plot_grid(plotlist=plotlist,ncol=3))
    dev.off()
}

# iterate through each object
i <- 3
object <- new.data[[i]]
sample.name <- names(new.data[i])
# Identify clusters to be removed & label clusters suspicious (sus)
trash <- c(0,6)
sus <- c()
new.labels <- levels(object$seurat_clusters)
for (i in 1:length(sus)){
    x <- which(new.labels == sus[i])
    new.labels[x] <- c('sus')
}
for (i in 1:length(trash)){
    x <- which(new.labels == trash[i])
    new.labels[x] <- c('Trash')
}
names(new.labels) <- levels(object$seurat_clusters)
object <- RenameIdents(object, new.labels)
object$rnd2_suspicious <- Idents(object)
p1 <- UMAPPlot(object,group.by = 'rnd2_suspicious',cols = c('sus' = 'springgreen4','Trash' = 'red')) +
        labs(title = c('Remaining Contaminant Populations')) + theme(plot.title = element_text(size = 24,hjust = 0))

# Remove trash, preserve contaminants, re-cluster
object <- subset(object, idents = c('Trash'), invert = T)

tag <- paste('.pcs',pcs,'res',res,sep='.')
plots <- list()
object <- NormalizeData(object)
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, npcs = 100, verbose = F)
object <- FindNeighbors(object, dims = 1:pcs)
object <- FindClusters(object, resolution = res)
object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)

p2 <- UMAPPlot(object,label = T,label.size = 6) + labs(title = paste(sample.name,'- Clean'),subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) +
        theme(plot.title = element_text(size = 24),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoLegend()
p3 <- FeaturePlot(object, goi, label = T,repel=T)
p4 <- FeaturePlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','Mki67'),ncol=2,label = T,repel=T)
p5 <- VlnPlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','Mki67'),ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())
p6 <- FeaturePlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'), label = T,repel=T)
p7 <- VlnPlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'),ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())

# Class objects
new.labels <- levels(object$seurat_clusters)
epi <- c(1,3,5,7)
endo <- c(2,4,9)
mes <- c(6)
imm <- c(0,8)
for (i in 1:length(levels(object$seurat_clusters))){
    for (j in 1:length(epi)){
        x <- which(new.labels == epi[j])
        new.labels[x] <- c('Epithelium')
    }
    for (j in 1:length(endo)){
        x <- which(new.labels == endo[j])
        new.labels[x] <- c('Endothelium')
    }
    for (j in 1:length(mes)){
        x <- which(new.labels == mes[j])
        new.labels[x] <- c('Mesenchyme')
    }
    for (j in 1:length(imm)){
        x <- which(new.labels == imm[j])
        new.labels[x] <- c('Immune')
    }
}
names(new.labels) <- levels(object$seurat_clusters)
object <- RenameIdents(object, new.labels)
object$class <- Idents(object)
cols <- hue_pal()(4)
class_colors <- c('Epithelium'=cols[2],'Endothelium'=cols[3],'Mesenchyme'=cols[1],'Immune'=cols[4])
p8 <- UMAPPlot(object,group.by = 'class',cols = class_colors) +
        labs(title = c('Cell Class')) + theme(plot.title = element_text(size = 24,hjust = 0))

toprow <- plot_grid(p1,p2,p3,ncol=3,labels = c('A','B','C'),label_size = 40)
midrow <- plot_grid(p4,p5,NULL,ncol=3,labels = c('D','E'),label_size = 40)
lastrow <- plot_grid(p6,p7,p8,ncol=3,labels = c('F','G','H'),label_size = 40)
plots <- list(toprow,midrow,lastrow)

png(paste0(sample.name,tag,'_cleancheck2.png'), width = 2000,height=2063)
print(plot_grid(plotlist=plots,nrow=3))
dev.off()

# Save clean individual objects
save(object,file = paste(sample.name,"clean",Sys.Date(),"Robj",sep="."))
