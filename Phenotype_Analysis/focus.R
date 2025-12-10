#### FOCUS
#### Code for publication unifying all 'focus' scripts
#### Replicable methods for focused analyses of classed objects,
####    including archetype definition, pseudotime via slingshot, & RNA velocity via velocyto-R




## Note: Engineered cell archetypes and cell types were identified through iterative clustering of integrated class
##  objects, combined with targeted transcriptomic characterizations. Iterative analyses that contributed to archetype
##  definition are not presented here, but rather the final steps in assigning archetype labels, which should yield
##  metadata that is consistent with the final processed objects uploaded to GEO.




#### EPITHELIUM ####

### Archetype and cell type definition of integrated object clustering
## Re-integrate epithelial object including starting cells and engineered populations
# Load object containing putative labels from previous analyses
dir_name <- c('epi')
load('./epi.seurat.objs.int.2022-12-15.Robj')

object$Sample <- factor(object$Sample,levels=c('BC1P3','BC1P6','BCL5','BCEC2','BEF1',
    'BEF2','BEF3','BEF12','BEF14','BEF15','BEFM1','BEFM2','BEFM4','BEFM5','BEFM6'))
object$Dataset2 <- factor(object$Dataset2,levels=c('Start','Mono','Co','Tri_E','Tri_L',
    'Quad_E','Quad_L'))

# Re-embed and -cluster on existing integration
DefaultAssay(object) <- "integrated"
object <- ScaleData(object)
object <- RunPCA(object, npcs = 100, verbose = F)

pcs <- 37
res <- 0.6
tag <- paste('pcs',pcs,'res',res,sep='.')
DefaultAssay(object) <- "integrated"
object <- FindNeighbors(object, dims = 1:pcs)
object <- FindClusters(object, resolution = res)
object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
p1 <- UMAPPlot(object,group.by = 'seurat_clusters',label = T) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoAxes()
p2 <- UMAPPlot(object,group.by = 'Sample') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
p3 <- UMAPPlot(object,group.by = 'putative_labels') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
Idents(object) <- object$Sample
toprow <- plot_grid(p1,p2,p3,ncol=3)

data <- prop.table(table(object$seurat_clusters,object$Sample),1)*100
data <- as.data.frame(data)
p5 <- ggplot(data, aes(x = Var1, y = Freq, fill = Var2)) + NoLegend() +
  geom_col(position = "dodge") + scale_fill_manual(values = hue_pal()(length(unique(object$Sample)))) +
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
midrow1 <- plot_grid(p5)

DefaultAssay(object) <- "RNA"
Idents(object) <- object$seurat_clusters
p6 <- FeaturePlot(object, c('Krt5','Hopx','Aqp5','Sftpc'),ncol=2,label = F)  # Epi
p7 <- FeaturePlot(object, c('Nkx2-1','Sftpb','Defb4','Lyz2'),ncol=2,label = F)  # ATII
p8 <- FeaturePlot(object, c('Ager','Aqp5','Hopx','Pdpn'),ncol=2,label = F)  # ATI
p9 <- FeaturePlot(object, c('Krt5','Krt8','Igfbp2','Sox2'),ncol=2,label = F)  # Basal
p10 <- FeaturePlot(object, c('Scgb1a1','Rarres1','Muc20','Clca1'),ncol=2,label = F)  # Secretory
p11 <- FeaturePlot(object, c('Dclk1','Espn','Sox9','Pou2f3'),ncol=2,label = F)  # Tuft1
p12 <- FeaturePlot(object, c('Top2a','Trpm5','Sox9','Lgr5'),ncol=2,label = F)  # Tuft2
p13 <- FeaturePlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','Mki67'),ncol=2,label = F)  # QC
p14 <- FeaturePlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'),ncol=2,label = F)  # Lineage
midrow2 <- plot_grid(p6,p7,p8,ncol=3)
midrow3 <- plot_grid(p9,p10,p11,ncol=3)
lastrow <- plot_grid(p12,p13,p14,ncol=3)

# Cowplot all panels
plots <- list(toprow,midrow1,midrow2,midrow3,lastrow)
png(paste(dir_name,tag,'_integrated7_allplots.png',sep=""), width = 2000,height=2750)
plot_grid(plotlist=plots,ncol=1,rel_heights = c(2,1,2,2,2))
dev.off()


## Remake putative labels
table(object$seurat_clusters,object$putative_labels)
Idents(object) <- object$putative_labels
cycling <- subset(object,idents='Cycling')
basc <- subset(object,idents='BASC')
Idents(object) <- object$seurat_clusters

# create new column in metadata
object$sub_cluster <- as.character(object$seurat_clusters)
object$sub_cluster[Cells(cycling)] <- paste(Idents(cycling))
object$sub_cluster[Cells(basc)] <- paste(Idents(basc))

table(object$sub_cluster,object$putative_labels)

p1 <- UMAPPlot(object,group.by = 'sub_cluster') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20)) + NoAxes()
p2 <- UMAPPlot(object,group.by = 'putative_labels') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
toprow <- plot_grid(p1,p2,ncol=2)
p3 <- FeaturePlot(object, c('Krt5','Hopx'),ncol=2,order=T)
p4 <- VlnPlot(object,c('Krt5','Hopx'),group.by='sub_cluster',ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())

png(paste(dir_name,tag,'_integrated7_subcluster.png',sep=""), width = 1000,height=1500)
plot_grid(plotlist=list(toprow,p3,p4),ncol=1)
dev.off()

Idents(object) <- object$sub_cluster
object <- RenameIdents(object,
                '0'='Hopx',
                '1'='Krt5',
                '2'='Trans',
                '3'='Hopx',
                '4'='Trans',
                '5'='Krt5',
                '6'='Trans',
                '7'='Trans',
                '8'='Krt5',
                '9'='Trans',
                '10'='Krt5',
                '11'='Krt5',
                '12'='Krt5',
                '13'='Secretory')
object$putative_labels <- Idents(object)

table(object$putative_labels)
object$putative_labels <- factor(object$putative_labels,levels=c('Krt5','Hopx',
    'Trans','Cycling','BASC','Secretory'))

save(object,file = paste("epi.seurat.objs.int.",Sys.Date(),".Robj",sep=""))



### Pseudotime via Slingshot on final epithelial object
## Slingshot v2.2.1
# Load and prepare object and color palettes 
dir_name <- c('epi')
load('./epi.refined.2023-03-06.Robj')
object <- epi ; rm(epi)
DefaultAssay(object) <- "RNA"
# pull and add metadata of interest
load('./epi.clean.annotations_2023-04-03.Robj')
tosave <- c('Final1','Final2','Dataset3')
epi.metadata <- epi@meta.data
epi.metadata <- epi.metadata[c(tosave)]
object <- AddMetaData(object,metadata = epi.metadata,col.name=names(epi.metadata))
rm(epi,epi.metadata)
object$putative_labels <- object$stash
object$putative_labels <- factor(object$putative_labels,levels=c('Krt5','Hopx',
    'Trans','Cycling','BASC'))
object$Final1 <- factor(object$Final1,levels=c('Basal-like','ATI-like',
    'Transitional','Cycling','BASC'))
load('../all.color_palettes_1.R')
sample_colors <- color_pals$sample_colors
names(sample_colors)[6] <- "BAL"
dataset_colors <- color_pals$dataset_colors
class_colors <- color_pals$class_colors
class <- dir_name
pcs <- c(1,4,5,7,9,10,11,30,35)
res <- 0.25
if (class == 'imm'){
        types <- levels(object$putative_labels)  # n=10
        pt1 <- brewer.pal(6,"Dark2")
        pt2 <- brewer.pal(8,"Set2")[-7]
        cols <- c('#7570B3','#B3B3B3','#66C2A5','#8DA0CB','#FC8D62','#FFD92F','#E7298A',
                '#E78AC3','#A6D854','#66A61E')
        type_colors <- cols
        names(type_colors) <- types       
} else if (class != 'imm'){
        types <- levels(object$Final1)  # n<10
        cols <- brewer.pal(length(levels(object$Final1)),'Set2')
        type_colors <- cols
        names(type_colors) <- types
}

# Run Slingshot pipeline
DefaultAssay(object) <- "integrated"  # recommended
slobject <- as.SingleCellExperiment(object,assay = c('integrated','RNA'))

slobject <- slingshot(slobject,reducedDim = reducedDim(slobject,'PCA')[,pcs],
    approx_points = 40,start.clus = 'Cycling',clusterLabels = 'putative_labels')

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(slobject$slingPseudotime_1, breaks=100)]

png(paste(dir_name,'sling_int3_umap_1.png',sep='_'), width = 800, height = 800)
plot(reducedDims(slobject)$UMAP, col = plotcol, pch=16, asp = 1)
dev.off()

save(slobject,file = paste("epi.seurat.objs.sling",Sys.Date(),"Robj",sep="."))

# Pull and save Slingshot pseudotime metadata
sling.metadata <- slingPseudotime(slobject)

sling.metadata2 <- slingCurveWeights(slobject)

object <- AddMetaData(object,metadata = sling.metadata,col.name=colnames(sling.metadata))

png(paste(dir_name,'sling_int3_fp_1.png',sep='_'), width = 800, height = 400)
FeaturePlot(object,c('Lineage1','Lineage2'),ncol=2,label = T,repel=T,order=T)
dev.off()

sling.metadata <- slingPseudotime(slobject)
colnames(sling.metadata) <- paste0(colnames(sling.metadata),'_ps')

sling.metadata2 <- slingCurveWeights(slobject)
colnames(sling.metadata2) <- paste0(colnames(sling.metadata2),'_cw')

sling.metadata <- data.frame(sling.metadata,sling.metadata2)

save(sling.metadata,file = paste("sling.metadata.",Sys.Date(),".Robj",sep=""))

object <- AddMetaData(object,metadata = sling.metadata,col.name=names(sling.metadata))




### RNA Velocity via Velocyto on final epithelial object
## Velocyto-R v0.6 & Pagoda2
# Load and prepare object and color palettes (see above)

# Run Velocyto-R pipeline
library(velocyto.R)
library('pagoda2')

# Setup with Pagoda2
DefaultAssay(object) <- "integrated"
Idents(object) <- object$Final1
object <- subset(object, downsample = 2000)
object <- FindVariableFeatures(object)
length(object@assays$integrated@var.features)
object <- ScaleData(object)
nrow(GetAssayData(object,slot="scale.data"))

r <- Pagoda2$new(as(GetAssayData(object,slot="scale.data"),"dgCMatrix"),modelType="raw")
r$reductions[['PCA']] <- object@reductions$pca
r$reductions[['umap']] <- object@reductions$umap
r$embeddings[['PCA']][['umap']] <- Embeddings(object[['umap']])

cluster.label <- object$Final1
cell.colors <- sccore:::fac2col(cluster.label,
        level.colors=type_colors[levels(object$Final1)])

png(paste(dir_name,'_pagumap_3.png',sep=''), width = 800, height = 800)
r$plotEmbedding(type='PCA',embeddingType='umap',colors=cell.colors,mark.clusters=T,
        min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,
        main='cell clusters')
dev.off()

# Calculating Velocity
emat <- object$spliced; nmat <- object$unspliced;
emat <- emat[,rownames(r$counts)]; nmat <- nmat[,rownames(r$counts)];

cluster.label <- object$Final1
cell.colors <- sccore:::fac2col(cluster.label,
        level.colors=type_colors[levels(object$Final1)])

cell.dist <- as.dist(1-armaCor(t(r$embeddings[['PCA']][['umap']])))

emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.08)
length(intersect(rownames(emat),rownames(nmat)))
# [1] 912

fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,
        cell.dist=cell.dist,fit.quantile=fit.quantile)
emb <- r$embeddings$PCA$umap

x <- show.velocity.on.embedding.cor(emb,rvel.cd,n=200,n.cores=1,scale='sqrt',
        cell.colors=ac(cell.colors,alpha=0.5),cex=0.8,arrow.scale=8,show.grid.flow=T,
        min.grid.cell.mass=5,grid.n=20,arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)

res <- 400
png(paste(dir_name,'_pagvel_3.png',sep=''), width=650*res/72, height=650*res/72,res=res)
show.velocity.on.embedding.cor(emb,rvel.cd,n=800,n.cores=1,scale='sqrt',cc=x$cc,
        cell.colors=ac(cell.colors),cex=1,arrow.scale=15,show.grid.flow=T,
        min.grid.cell.mass=30,grid.n=20,arrow.lwd=4,do.par=F,cell.border.alpha = 0.1,
        xlab = c('UMAP_1'),ylab= c('UMAP_2'))
dev.off()

# Velocyto on BASC subset of epithelium
# Subset based on embedding
cut <- Embeddings(object[['umap']])[which(Embeddings(object[['umap']])[,'UMAP_1'] > -3),]
cut <- cut[which(cut[,'UMAP_1'] < 0.5),]
cut <- cut[which(cut[,'UMAP_2'] > -4),]
cut <- cut[which(cut[,'UMAP_2'] < 0),]

sub <- subset(object,cells=rownames(cut))  # 628 cells

# Setup with Pagoda2
DefaultAssay(sub) <- "integrated"
Idents(sub) <- sub$Final1
sub <- FindVariableFeatures(sub)
length(sub@assays$integrated@var.features)
sub <- ScaleData(sub)
nrow(GetAssayData(sub,slot="scale.data"))

r <- Pagoda2$new(as(GetAssayData(sub,slot="scale.data"),"dgCMatrix"),modelType="raw")
r$reductions[['PCA']] <- sub@reductions$pca
r$reductions[['umap']] <- sub@reductions$umap
r$embeddings[['PCA']][['umap']] <- Embeddings(sub[['umap']])

cluster.label <- sub$Final1
cell.colors <- sccore:::fac2col(cluster.label,
        level.colors=type_colors[levels(sub$Final1)])

png(paste(dir_name,'_pagumap_sub_1.png',sep=''), width = 800, height = 800)
r$plotEmbedding(type='PCA',embeddingType='umap',colors=cell.colors,mark.clusters=T,
        min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,
        main='cell clusters')
dev.off()

# Calculating Velocity
emat <- sub$spliced; nmat <- sub$unspliced;
emat <- emat[,rownames(r$counts)]; nmat <- nmat[,rownames(r$counts)];

cluster.label <- sub$Final1
cell.colors <- sccore:::fac2col(cluster.label,
        level.colors=type_colors[levels(sub$Final1)])

cell.dist <- as.dist(1-armaCor(t(r$embeddings[['PCA']][['umap']])))

emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.08)
length(intersect(rownames(emat),rownames(nmat)))
# [1] 667

fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,
        cell.dist=cell.dist,fit.quantile=fit.quantile)
emb <- r$embeddings$PCA$umap

x <- show.velocity.on.embedding.cor(emb,rvel.cd,n=200,n.cores=1,scale='sqrt',
        cell.colors=ac(cell.colors,alpha=0.5),cex=0.8,arrow.scale=8,show.grid.flow=T,
        min.grid.cell.mass=5,grid.n=20,arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)

res <- 400
png(paste(dir_name,'_pagvel_sub_1.png',sep=''), width=400*res/72, height=400*res/72,res=res)
show.velocity.on.embedding.cor(emb,rvel.cd,n=25,n.cores=1,scale='sqrt',cc=x$cc,
        cell.colors=ac(cell.colors),cex=1,arrow.scale=5,show.grid.flow=T,
        min.grid.cell.mass=5,grid.n=18,arrow.lwd=4,do.par=F,cell.border.alpha = 0.1,
        min.arrow.size = 0.25,xlab = c('UMAP_1'),ylab= c('UMAP_2'))
dev.off()




#### ENDOTHELIUM ####

### Archetype and cell type definition of integrated object clustering
## Re-integrate endothelial object including starting cells and engineered populations
# Load object containing putative labels from previous analyses
dir_name <- c('endo')
load('./endo.seurat.objs.int.2022-10-28.Robj')

object$putative_labels <- factor(object$putative_labels,levels=c('Interesting','Cycling2','Cycling3','Lymphatic','Homeostatic'))
object$Dataset2 <- Idents(object)
object$Dataset2 <- factor(object$Dataset2,levels=c('Start','Co','Tri_E','Tri_L','Quad_E','Quad_L'))

# Re-embed and -cluster on existing integration
DefaultAssay(object) <- "integrated"
object <- ScaleData(object)
object <- RunPCA(object, npcs = 100, verbose = F)
png(paste(dir_name,'_integrated4_elbow.png',sep=""), width = 800, height = 500)
ElbowPlot(object,ndims = 50)+labs(title = dir_name)
dev.off()
pdf(paste(dir_name,'_integrated4_pcheatmap1.pdf',sep=""), width = 12, height = 16)
DimHeatmap(object, dims = 1:15, cells = 200, balanced = T)
dev.off()
pdf(paste(dir_name,'_integrated4_pcheatmap2.pdf',sep=""), width = 12, height = 16)
DimHeatmap(object, dims = 16:30, cells = 200, balanced = T)
dev.off()
pdf(paste(dir_name,'_integrated4_pcheatmap3.pdf',sep=""), width = 12, height = 16)
DimHeatmap(object, dims = 31:45, cells = 200, balanced = T)
dev.off()

pcs <- 10
res <- 0.3
tag <- paste('pcs',pcs,'res',res,sep='.')
DefaultAssay(object) <- "integrated"
object <- FindNeighbors(object, dims = 1:pcs)
object <- FindClusters(object, resolution = res)
object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
p1 <- UMAPPlot(object,group.by = 'seurat_clusters',label = T) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoAxes()
p2 <- UMAPPlot(object,group.by = 'Sample') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
p3 <- UMAPPlot(object,group.by = 'putative_labels') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
Idents(object) <- object$Sample
toprow <- plot_grid(p1,p2,p3,ncol=3)

data <- prop.table(table(object$seurat_clusters,object$Sample),1)*100
data <- as.data.frame(data)
p5 <- ggplot(data, aes(x = Var1, y = Freq, fill = Var2)) + NoLegend() +
  geom_col(position = "dodge") + scale_fill_manual(values = hue_pal()(length(unique(object$Sample)))) +
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
midrow1 <- plot_grid(p5)

DefaultAssay(object) <- "RNA"
Idents(object) <- object$seurat_clusters
p6 <- FeaturePlot(object, c('Lyve1','Prx','Plvap','Apln'),ncol=2,label = F)  # Endothelium
p7 <- FeaturePlot(object, c('Plvap','Aplnr','Cdh5','Top2a'),ncol=2,label = F)  # gCap1
p8 <- FeaturePlot(object, c('Kit','Wnt5b','Nostrin','Vegfa'),ncol=2,label = F)  # gCap2
p9 <- FeaturePlot(object, c('Apln','Prx','Ca4','Ednrb'),ncol=2,label = F)  # aCap
p10 <- FeaturePlot(object, c('Efnb2','Vwf','Lyve1','Gja5'),ncol=2,label = F)  # Arterial
p11 <- FeaturePlot(object, c('Vcam1','Vwf','Lyve1','Selp'),ncol=2,label = F)  # Venous
p12 <- FeaturePlot(object, c('Prox1','Mmrn1','nCount_RNA','nFeature_RNA'),ncol=4,label = F)  # Lymphatic
midrow2 <- plot_grid(p6,p7,p8,ncol=3)
midrow3 <- plot_grid(p9,p10,p11,ncol=3)
lastrow <- plot_grid(p12,NULL,ncol=2,rel_widths=c(2,1))

# Cowplot all panels
plots <- list(toprow,midrow1,midrow2,midrow3,lastrow)
png(paste(dir_name,tag,'_integrated4_allplots.png',sep=""), width = 2000,height=2750)
plot_grid(plotlist=plots,ncol=1,rel_heights = c(2,1,2,2,1))
dev.off()


## Remake putative labels
table(object$seurat_clusters,object$putative_labels)
Idents(object) <- object$seurat_clusters
object <- RenameIdents(object,
                '0'='Interesting',
                '1'='Cycling2',
                '2'='Lymphatic',
                '3'='Cycling3',
                '4'='Homeostatic')
object$putative_labels <- Idents(object)

table(object$putative_labels)
object$putative_labels <- factor(object$putative_labels,levels=c('Interesting',
    'Cycling2','Cycling3','Lymphatic','Homeostatic'))

table(object$Sample)
object$Sample <- factor(object$Sample,levels=c('RLMVEC','BCEC2','BEF1','BEF2','BEF3',
    'BEF12','BEF14','BEF15','BEFM1','BEFM2','BEFM4','BEFM5','BEFM6'))

# Get new marker list for endothelial putative_labels
Idents(object) <- object$putative_labels
markers <- FindAllMarkers(object, only.pos = T)
markers$ratio <- markers$pct.1/markers$pct.2
markers <- markers[order(-markers$ratio),]

save(markers,file=paste(dir_name,'seurat.objs.int.marks.',Sys.Date(),".Robj",sep=""))

save(object,file = paste("endo.seurat.objs.int.",Sys.Date(),".Robj",sep=""))



### Pseudotime via Slingshot on final endothelial object
## Slingshot v2.2.1
# Load and prepare object and color palettes 
dir_name <- c('endo')
load('./endo.seurat.objs.int.2023-01-25.Robj')
load('../all.color_palettes_1.R')
sample_colors <- color_pals$sample_colors
names(sample_colors)[6] <- "BAL"
dataset_colors <- color_pals$dataset_colors
class_colors <- color_pals$class_colors
class <- dir_name
object$putative_labels <- factor(object$putative_labels,levels=c('Interesting','Cycling2','Lymphatic','Cycling3','Homeostatic'))
pcs <- 10
res <- 0.3
if (class == 'imm'){
        types <- levels(object$putative_labels)  # n=10
        pt1 <- brewer.pal(6,"Dark2")
        pt2 <- brewer.pal(8,"Set2")[-7]
        cols <- c('#7570B3','#B3B3B3','#66C2A5','#8DA0CB','#FC8D62','#FFD92F','#E7298A',
                '#E78AC3','#A6D854','#66A61E')
        type_colors <- cols
        names(type_colors) <- types       
} else if (class != 'imm'){
        types <- levels(object$putative_labels)  # n<10
        cols <- brewer.pal(length(levels(object$putative_labels)),'Set2')
        type_colors <- cols
        names(type_colors) <- types
}

# Run Slingshot pipeline
Idents(object) <- object$Final1
DefaultAssay(object) <- "integrated"
slobject <- as.SingleCellExperiment(object,assay = c('integrated','RNA'))

slobject <- slingshot(slobject,reducedDim = reducedDim(slobject,'PCA')[,1:10],
    approx_points = 40,start.clus = 'Progenitor',clusterLabels = 'Final1')

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(slobject$slingPseudotime_1, breaks=100)]

png(paste(dir_name,'sling_int_umap_4.png',sep='_'), width = 800, height = 800)
plot(reducedDims(slobject)$UMAP, col = plotcol, pch=16, asp = 1)
dev.off()

save(slobject,file = paste(dir_name,"seurat.objs.sling",Sys.Date(),"Robj",sep="."))

# Featureplots of pseudotime by lineage, split by Sample
if (ncol(object) < 10000){pt.size=1.5} else if (ncol(object) >= 10000){pt.size=NULL}
p1 <- UMAPPlot(object,group.by = 'Final1',cols = type_colors,pt.size=pt.size) + 
        labs(title = class,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) + 
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red'))
p2 <- UMAPPlot(object,group.by = 'Sample',cols = sample_colors[levels(object$Sample)],pt.size=pt.size) + 
        labs(title = class,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = c(' ')) + 
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red'))
p3 <- FeaturePlot(object,'Lineage1_ps',cols = colorRampPalette(brewer.pal(11,'Spectral')[-6])(100),pt.size=pt.size,label = T,repel=T,order=T)
toprow <- plot_grid(p1,p2,p3,ncol=3)
p4 <- FeaturePlot(object,'Lineage1_ps',split.by='Sample',combine=F,cols = colorRampPalette(brewer.pal(11,'Spectral')[-6])(100),
        keep.scale='all',pt.size=pt.size,label = T,repel=T,order=T)
lastrow <- plot_grid(plotlist = p4, ncol=3)

png(paste(dir_name,'sling_int_lin_fp_1.png',sep='_'), width = 2000,height=2750)
plot_grid(toprow,lastrow,nrow=2,rel_heights=c(1,5))
dev.off()

# Pull and save Slingshot pseudotime metadata
sling.metadata <- slingPseudotime(slobject)
colnames(sling.metadata) <- paste0(colnames(sling.metadata),'_ps')

sling.metadata2 <- slingCurveWeights(slobject)
colnames(sling.metadata2) <- paste0(colnames(sling.metadata2),'_cw')

sling.metadata <- data.frame(sling.metadata,sling.metadata2)

save(sling.metadata,file = paste("sling.metadata.",Sys.Date(),".Robj",sep=""))



### RNA Velocity via Velocyto on final epithelial object
## Velocyto-R v0.6 & Pagoda2
# Load and prepare object and color palettes (see above)

# Run Velocyto-R pipeline
library(velocyto.R)
library('pagoda2')

# Setup with Pagoda2
DefaultAssay(object) <- "integrated"
object <- FindVariableFeatures(object)
length(object@assays$integrated@var.features)
object <- ScaleData(object)
nrow(GetAssayData(object,slot="scale.data"))

r <- Pagoda2$new(as(GetAssayData(object,slot="scale.data"),"dgCMatrix"),modelType="raw")
r$reductions[['PCA']] <- object@reductions$pca
r$reductions[['umap']] <- object@reductions$umap
r$embeddings[['PCA']][['umap']] <- Embeddings(object[['umap']])

cluster.label <- object$Final1
cell.colors <- sccore:::fac2col(cluster.label,
        level.colors=type_colors[levels(object$Final1)])

png(paste(dir_name,'_pagumap_2.png',sep=''), width = 800, height = 800)
r$plotEmbedding(type='PCA',embeddingType='umap',colors=cell.colors,mark.clusters=T,
        min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,
        main='cell clusters')
dev.off()

# Calculating Velocity
emat <- object$spliced; nmat <- object$unspliced;
emat <- emat[,rownames(r$counts)]; nmat <- nmat[,rownames(r$counts)];

cluster.label <- object$Final1
cell.colors <- sccore:::fac2col(cluster.label,
        level.colors=type_colors[levels(object$Final1)])

cell.dist <- as.dist(1-armaCor(t(r$embeddings[['PCA']][['umap']])))

emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 1)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.5)
length(intersect(rownames(emat),rownames(nmat)))
# [1] 796

fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,
        cell.dist=cell.dist,fit.quantile=fit.quantile)
emb <- r$embeddings$PCA$umap

x <- show.velocity.on.embedding.cor(emb,rvel.cd,n=1000,n.cores=1,scale='sqrt',
        cell.colors=ac(cell.colors,alpha=0.5),cex=0.8,arrow.scale=8,show.grid.flow=T,
        min.grid.cell.mass=5,grid.n=20,arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)

png(paste(dir_name,'_pagvel_4.png',sep=''), width = 900, height = 800)
show.velocity.on.embedding.cor(emb,rvel.cd,n=800,n.cores=1,scale='sqrt',cc=x$cc,
        cell.colors=ac(cell.colors),cex=1,arrow.scale=11,show.grid.flow=T,
        min.grid.cell.mass=50,grid.n=15,arrow.lwd=4,do.par=F,cell.border.alpha = 0.1,
        xlab = c('UMAP_1'),ylab= c('UMAP_2'))
dev.off()

res <- 400
png(paste(dir_name,'_pagvel_4.png',sep=''), width=650*res/72, height=650*res/72,res=res)
show.velocity.on.embedding.cor(emb,rvel.cd,n=800,n.cores=1,scale='sqrt',cc=x$cc,
        cell.colors=ac(cell.colors),cex=1,arrow.scale=11,show.grid.flow=T,
        min.grid.cell.mass=50,grid.n=15,arrow.lwd=4,do.par=F,cell.border.alpha = 0.1,
        xlab = c('UMAP_1'),ylab= c('UMAP_2'))
dev.off()



### Plot velocyto colored by slingshot
load('./sling.metadata.2023-04-04.Robj')
object <- AddMetaData(object,metadata = sling.metadata,col.name=names(sling.metadata))

cluster.label <- sort(object$Lineage1_ps)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(ncol(object))
cell.colors <- sccore:::fac2col(cluster.label,shuffle=FALSE,
        level.colors=colors)

png(paste(dir_name,'_pagvel_3.png',sep=''), width = 800, height = 800)
show.velocity.on.embedding.cor(emb,rvel.cd,n=1000,n.cores=1,scale='sqrt',cc=x$cc,
        cell.colors=ac(cell.colors,alpha=0.5),cex=1,arrow.scale=12,show.grid.flow=T,
        min.grid.cell.mass=50,grid.n=15,arrow.lwd=4,do.par=F,cell.border.alpha = 0.1)
dev.off()




#### MESENCHYME ####

### Archetype and cell type definition of integrated object clustering
## Re-integrate mesenchymal object including starting cells and engineered populations
# Load object containing putative labels from previous analyses
dir_name <- c('mes')
load('./mes.seurat.objs.int.2022-10-14.Robj')
load('./mes.metadata.2022-11-18.Robj')
object <- AddMetaData(object,metadata = mes.metadata,col.name=names(mes.metadata))

object$Sample <- factor(object$Sample,levels=c('FB13','FB14','BEF1','BEF2','BEF3',
    'BEF12','BEF14','BEF15','BEFM1','BEFM2','BEFM4','BEFM5','BEFM6'))
object$Dataset2 <- factor(object$Dataset2,levels=c('Start','Tri_E','Tri_L','Quad_E',
    'Quad_L'))

# Re-embed and -cluster on existing integration
DefaultAssay(object) <- "integrated"
object <- ScaleData(object)
object <- RunPCA(object, npcs = 100, verbose = F)

pcs <- 35
res <- 0.35
tag <- paste('pcs',pcs,'res',res,sep='.')
DefaultAssay(object) <- "integrated"
object <- FindNeighbors(object, dims = 1:pcs)
object <- FindClusters(object, resolution = res)
object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
p1 <- UMAPPlot(object,group.by = 'seurat_clusters',label = T) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoAxes()
p2 <- UMAPPlot(object,group.by = 'Sample') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
p3 <- UMAPPlot(object,group.by = 'putative_labels') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
Idents(object) <- object$Sample
toprow <- plot_grid(p1,p2,p3,ncol=3)

data <- prop.table(table(object$seurat_clusters,object$Sample),1)*100
data <- as.data.frame(data)
p5 <- ggplot(data, aes(x = Var1, y = Freq, fill = Var2)) + NoLegend() +
  geom_col(position = "dodge") + scale_fill_manual(values = hue_pal()(length(unique(object$Sample)))) +
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
midrow1 <- plot_grid(p5)

DefaultAssay(object) <- "RNA"
Idents(object) <- object$seurat_clusters
p6 <- FeaturePlot(object, c('Dcn','Acta2','Itga8','Wnt11'),ncol=2,label = F)  # Mesenchyme
p7 <- FeaturePlot(object, c('Col14a1','Dcn','Pi16','Col13a1'),ncol=2,label = F)  # Pi16+ Fibs
p8 <- FeaturePlot(object, c('Msln','Rspo1','Dclk1','Krt8'),ncol=2,label = F)  # Mesothelium
p9 <- FeaturePlot(object, c('Col13a1','Itga8','Hsd11b1','Col14a1'),ncol=2,label = F)  # Itga8+ Fibs
p10 <- FeaturePlot(object, c('Wnt11','Dcn','Acta2','Myh11'),ncol=2,label = F)  # Myofibroblasts
p11 <- FeaturePlot(object, c('Lgr5','Notum','Col13a1','Col14a1'),ncol=2,label = F)  # Notum+ Fibs
p12 <- FeaturePlot(object, c('Acta2','Myh11','Acta1','Itga8'),ncol=2,label = F)  # Smooth Muscle
p13 <- FeaturePlot(object, c('Gucy1a1','Gucy1a2','Gucy1b1','Acta2'),ncol=2,label = F)  # Pericytes
midrow2 <- plot_grid(p6,p7,p8,ncol=3)
midrow3 <- plot_grid(p9,p10,p11,ncol=3)
lastrow <- plot_grid(p12,p13,NULL,ncol=3)

# Cowplot all panels
plots <- list(toprow,midrow1,midrow2,midrow3,lastrow)
png(paste(dir_name,tag,'_integrated4_allplots.png',sep=""), width = 2000,height=2750)
plot_grid(plotlist=plots,ncol=1,rel_heights = c(2,1,2,2,2))
dev.off()


## Remake putative labels
table(object$seurat_clusters,object$putative_labels)
Idents(object) <- object$putative_labels
peri <- subset(object,idents='Pericytes')
meso <- subset(object,idents='Mesothelium')
homeo <- subset(object,idents='Homeostatic')

Idents(object) <- object$seurat_clusters
# create new column in metadata
object$sub_cluster <- as.character(object$seurat_clusters)
object$sub_cluster[Cells(peri)] <- paste(Idents(peri))
object$sub_cluster[Cells(meso)] <- paste(Idents(meso))
object$sub_cluster[Cells(homeo)] <- paste(Idents(homeo))

table(object$sub_cluster,object$putative_labels)

p1 <- UMAPPlot(object,group.by = 'sub_cluster') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20)) + NoAxes()
p2 <- UMAPPlot(object,group.by = 'putative_labels') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
toprow <- plot_grid(p1,p2,ncol=2)
png(paste(dir_name,tag,'_integrated7_subcluster.png',sep=""), width = 1000,height=500)
plot_grid(plotlist=list(toprow),ncol=1)
dev.off()

Idents(object) <- object$sub_cluster
object <- RenameIdents(object,
                '0'='General',
                '1'='General',
                '2'='General',
                '3'='General',
                '4'='General',
                '5'='Homeostatic',
                '6'='Cycling',
                '7'='Cycling',
                '8'='General',
                '9'='Mesothelium',
                '10'='Pericytes')
object$putative_labels <- Idents(object)

table(object$putative_labels)
object$putative_labels <- factor(object$putative_labels,levels=c('General','Cycling',
    'Homeostatic','Mesothelium','Pericytes'))


## Get new marker list for mesenchymal putative_labels
Idents(object) <- object$putative_labels
markers <- FindAllMarkers(object, only.pos = T)
markers$ratio <- markers$pct.1/markers$pct.2
markers <- markers[order(-markers$ratio),]

save(markers,file=paste(dir_name,'seurat.objs.int.marks',Sys.Date(),"Robj",sep="."))

save(object,file = paste(dir_name,"seurat.objs.int",Sys.Date(),"Robj",sep="."))



### Pseudotime via Slingshot on final mesenchymal object
## Slingshot v2.2.1
# Load and prepare object and color palettes 
dir_name <- c('mes')
load('./mes.seurat.objs.int.2023-01-27.Robj')
load('../all.color_palettes_1.R')
sample_colors <- color_pals$sample_colors
names(sample_colors)[6] <- "BAL"
dataset_colors <- color_pals$dataset_colors
class_colors <- color_pals$class_colors
class <- dir_name
object$putative_labels <- factor(object$putative_labels,levels=c('General','Cycling','Homeostatic','Mesothelium','Pericytes'))
pcs <- 35
res <- 0.35
if (class == 'imm'){
        types <- levels(object$putative_labels)  # n=10
        pt1 <- brewer.pal(6,"Dark2")
        pt2 <- brewer.pal(8,"Set2")[-7]
        cols <- c('#7570B3','#B3B3B3','#66C2A5','#8DA0CB','#FC8D62','#FFD92F','#E7298A',
                '#E78AC3','#A6D854','#66A61E')
        type_colors <- cols
        names(type_colors) <- types       
} else if (class != 'imm'){
        types <- levels(object$putative_labels)  # n<10
        cols <- brewer.pal(length(levels(object$putative_labels)),'Set2')
        type_colors <- cols
        names(type_colors) <- types
}

# Run Slingshot pipeline
DefaultAssay(object) <- "integrated"
slobject <- as.SingleCellExperiment(object,assay = c('integrated','RNA'))

slobject <- slingshot(slobject,reducedDim = reducedDim(slobject,'PCA')[,1:35],
    approx_points = 40,start.clus = 'Cycling',clusterLabels = 'putative_labels')

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(slobject$slingPseudotime_1, breaks=100)]

png(paste(dir_name,'sling_int_umap_1.png',sep='_'), width = 800, height = 800)
plot(reducedDims(slobject)$UMAP, col = plotcol, pch=16, asp = 1)
dev.off()

save(slobject,file = paste(dir_name,"seurat.objs.sling",Sys.Date(),"Robj",sep="."))
load('./mes.seurat.objs.sling.2023-02-14.Robj')

# Pull and save Slingshot pseudotime metadata
sling.metadata <- slingPseudotime(slobject)
sling.metadata2 <- slingCurveWeights(slobject)
object <- AddMetaData(object,metadata = sling.metadata,col.name=colnames(sling.metadata))

png(paste(dir_name,'sling_int_fp_1.png',sep='_'), width = 800, height = 800)
FeaturePlot(object,c('Lineage1','Lineage2','Lineage3'),ncol=2,label = T,repel=T,order=T)
dev.off()

sling.metadata <- slingPseudotime(slobject)
colnames(sling.metadata) <- paste0(colnames(sling.metadata),'_ps')

sling.metadata2 <- slingCurveWeights(slobject)
colnames(sling.metadata2) <- paste0(colnames(sling.metadata2),'_cw')

sling.metadata <- data.frame(sling.metadata,sling.metadata2)

save(sling.metadata,file = paste("sling.metadata.",Sys.Date(),".Robj",sep=""))

object <- AddMetaData(object,metadata = sling.metadata,col.name=names(sling.metadata))

# Featureplots of pseudotime by lineage
if (ncol(object) < 10000){pt.size=1.5} else if (ncol(object) >= 10000){pt.size=NULL}
p1 <- UMAPPlot(object,group.by = 'putative_labels',cols = type_colors,pt.size=pt.size) + 
        labs(title = class,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) + 
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red'))
p2 <- UMAPPlot(object,group.by = 'Sample',cols = sample_colors[levels(object$Sample)],pt.size=pt.size) + 
        labs(title = class,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = c(' ')) + 
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red'))
toprow <- plot_grid(p1,p2,NULL,ncol=3)
p3 <- FeaturePlot(object,'Lineage1_ps',cols = colorRampPalette(brewer.pal(11,'Spectral')[-6])(100),pt.size=pt.size,label = T,repel=T,order=T)
p4 <- FeaturePlot(object,'Lineage2_ps',cols = colorRampPalette(brewer.pal(11,'Spectral')[-6])(100),pt.size=pt.size,label = T,repel=T,order=T)
p5 <- FeaturePlot(object,'Lineage3_ps',cols = colorRampPalette(brewer.pal(11,'Spectral')[-6])(100),pt.size=pt.size,label = T,repel=T,order=T)
lastrow <- plot_grid(p3,p4,p5,ncol=3)

# Cowplot all panels
png(paste(dir_name,'sling_int_allplots_1.png',sep='_'), width = 2000,height=1250)
plot_grid(toprow,lastrow,nrow=2)
dev.off()

# Featureplots of pseudotime by lineage, split by Sample
if (ncol(object) < 10000){pt.size=1.5} else if (ncol(object) >= 10000){pt.size=NULL}
p1 <- UMAPPlot(object,group.by = 'putative_labels',cols = type_colors,pt.size=pt.size) + 
        labs(title = class,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) + 
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red'))
p2 <- UMAPPlot(object,group.by = 'Sample',cols = sample_colors[levels(object$Sample)],pt.size=pt.size) + 
        labs(title = class,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = c(' ')) + 
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red'))
p3 <- FeaturePlot(object,'Lineage1_ps',cols = colorRampPalette(brewer.pal(11,'Spectral')[-6])(100),pt.size=pt.size,label = T,repel=T,order=T)
toprow <- plot_grid(p1,p2,p3,ncol=3)
p4 <- FeaturePlot(object,'Lineage1_ps',split.by='Sample',combine=F,cols = colorRampPalette(brewer.pal(11,'Spectral')[-6])(100),
        keep.scale='all',pt.size=pt.size,label = T,repel=T,order=T)
lastrow <- plot_grid(plotlist = p4, ncol=3)

png(paste(dir_name,'sling_int_lin1_fp_1.png',sep='_'), width = 2000,height=2750)
plot_grid(toprow,lastrow,nrow=2,rel_heights=c(1,5))
dev.off()

# Ridge plots of pseudotime by lineage
object$putative_labels <- factor(object$putative_labels,levels=c('Pericytes','Homeostatic','Mesothelium','General','Cycling'))
type_colors <- type_colors[as.character(levels(object$putative_labels))]
p1 <- RidgePlot(object,'Lineage1_ps',group.by='putative_labels',cols = type_colors) + 
    theme(plot.title = element_text(size = 30),
        legend.text = element_text(size = 22),
        legend.key.size = unit(2, "lines"),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
p2 <- RidgePlot(object,'Lineage1_ps',group.by='Sample',cols = sample_colors[levels(object$Sample)]) + 
    theme(plot.title = element_text(size = 30),
        legend.text = element_text(size = 22),
        legend.key.size = unit(2, "lines"),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
p3 <- RidgePlot(object,'Lineage1_ps',group.by='Dataset2',cols = dataset_colors[levels(object$Dataset2)]) + 
    theme(plot.title = element_text(size = 30),
        legend.text = element_text(size = 22),
        legend.key.size = unit(2, "lines"),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
png(paste(dir_name,'sling_int_lin1_allplots_1.png',sep='_'), width = 2000,height=2000)
plot_grid(p1,p2,p3,nrow=3)
dev.off()

object$putative_labels <- factor(object$putative_labels,levels=c('Pericytes','Homeostatic','Mesothelium','General','Cycling'))
type_colors <- type_colors[as.character(levels(object$putative_labels))]
p1 <- RidgePlot(object,'Lineage2_ps',group.by='putative_labels',cols = type_colors) + 
    theme(plot.title = element_text(size = 30),
        legend.text = element_text(size = 22),
        legend.key.size = unit(2, "lines"),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
p2 <- RidgePlot(object,'Lineage2_ps',group.by='Sample',cols = sample_colors[levels(object$Sample)]) + 
    theme(plot.title = element_text(size = 30),
        legend.text = element_text(size = 22),
        legend.key.size = unit(2, "lines"),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
p3 <- RidgePlot(object,'Lineage2_ps',group.by='Dataset2',cols = dataset_colors[levels(object$Dataset2)]) + 
    theme(plot.title = element_text(size = 30),
        legend.text = element_text(size = 22),
        legend.key.size = unit(2, "lines"),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
png(paste(dir_name,'sling_int_lin2_allplots_1.png',sep='_'), width = 2000,height=2000)
plot_grid(p1,p2,p3,nrow=3)
dev.off()

object$putative_labels <- factor(object$putative_labels,levels=c('Pericytes','Homeostatic','Mesothelium','General','Cycling'))
type_colors <- type_colors[as.character(levels(object$putative_labels))]
p1 <- RidgePlot(object,'Lineage3_ps',group.by='putative_labels',cols = type_colors) + 
    theme(plot.title = element_text(size = 30),
        legend.text = element_text(size = 22),
        legend.key.size = unit(2, "lines"),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
p2 <- RidgePlot(object,'Lineage3_ps',group.by='Sample',cols = sample_colors[levels(object$Sample)]) + 
    theme(plot.title = element_text(size = 30),
        legend.text = element_text(size = 22),
        legend.key.size = unit(2, "lines"),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
p3 <- RidgePlot(object,'Lineage3_ps',group.by='Dataset2',cols = dataset_colors[levels(object$Dataset2)]) + 
    theme(plot.title = element_text(size = 30),
        legend.text = element_text(size = 22),
        legend.key.size = unit(2, "lines"),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
png(paste(dir_name,'sling_int_lin3_allplots_1.png',sep='_'), width = 2000,height=2000)
plot_grid(p1,p2,p3,nrow=3)
dev.off()

# Summary plot
if (ncol(object) < 10000){pt.size=1.5} else if (ncol(object) >= 10000){pt.size=NULL}
p1 <- UMAPPlot(object,group.by = 'putative_labels',cols = type_colors,pt.size=pt.size) + 
        labs(title = class,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) + 
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red'))
p2 <- UMAPPlot(object,group.by = 'Sample',cols = sample_colors[levels(object$Sample)],pt.size=pt.size) + 
        labs(title = class,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = c(' ')) + 
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red'))
p3 <- FeaturePlot(object,'Lineage3_ps',cols = colorRampPalette(brewer.pal(11,'Spectral')[-6])(100),pt.size=pt.size,label = T,repel=T,order=T)
toprow <- plot_grid(p1,p2,p3,ncol=3)

object$putative_labels <- factor(object$putative_labels,levels=c('Pericytes','Homeostatic','Mesothelium','General','Cycling'))
type_colors <- type_colors[as.character(levels(object$putative_labels))]
p4 <- RidgePlot(object,'Lineage3_ps',group.by='putative_labels',cols = type_colors) + 
    theme(plot.title = element_text(size = 30),
        legend.text = element_text(size = 22),
        legend.key.size = unit(2, "lines"),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
p5 <- RidgePlot(object,'Lineage3_ps',group.by='Sample',cols = sample_colors[levels(object$Sample)]) + 
    theme(plot.title = element_text(size = 30),
        legend.text = element_text(size = 22),
        legend.key.size = unit(2, "lines"),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
p6 <- RidgePlot(object,'Lineage3_ps',group.by='Dataset2',cols = dataset_colors[levels(object$Dataset2)]) + 
    theme(plot.title = element_text(size = 30),
        legend.text = element_text(size = 22),
        legend.key.size = unit(2, "lines"),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
midrow <- plot_grid(p4,p5,p6,nrow=3)

data <- prop.table(table(object$Sample,object$Lineage3_bi),1)*100
data <- as.data.frame(data)
data_sub <- data[data$Var2 %in% 1, ]
p7 <- ggplot(data_sub, aes(x = Var1, y = Freq, fill = Var1)) + NoLegend() +
  geom_col(position = "dodge") + scale_fill_manual(values = sample_colors[levels(object$Sample)]) +
  labs(title = paste('Lineage 3 Composition per Sample')) + xlab("") + ylab("% lineage per sample") +
  theme(plot.title = element_text(size = 24),
        legend.title = element_blank(),
        axis.text.x=element_text(size=18,angle = 90,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

data2 <- prop.table(table(object$Dataset2,object$Lineage3_bi),1)*100
data2 <- as.data.frame(data2)
data2_sub <- data2[data2$Var2 %in% 1, ]
p8 <- ggplot(data2_sub, aes(x = Var1, y = Freq, fill = Var1)) + NoLegend() +
  geom_col(position = "dodge") + scale_fill_manual(values = dataset_colors[levels(object$Dataset2)]) +
  labs(title = paste('Lineage 3 Composition per Dataset')) + xlab("") + ylab("% lineage per dataset") +
  theme(plot.title = element_text(size = 24),
        legend.title = element_blank(),
        axis.text.x=element_text(size=18,angle = 90,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
lastrow <- plot_grid(p7,p8,ncol=2)

png(paste(dir_name,'sling_int_lin3_allplots_2.png',sep='_'), width = 2000,height=2750)
plot_grid(toprow,midrow,lastrow,nrow=3,rel_heights=c(1,2,1))
dev.off()



### RNA Velocity via Velocyto on final epithelial object
## Velocyto-R v0.6 & Pagoda2
# Load and prepare object and color palettes (see above and the following)
load('./mes.clean.annotations_2023-04-03.Robj')
tosave <- c('Final1','Final2','Dataset3')
mes.metadata <- mes@meta.data
mes.metadata <- mes.metadata[c(tosave)]
object <- AddMetaData(object,metadata = mes.metadata,col.name=names(mes.metadata))
rm(mes,mes.metadata)
object$Final1 <- factor(object$Final1,levels=c('Alveolar','Cycling',
    'Aberrant','Remodeling','Pericytes'))
type_colors <- c('Alveolar'='#66C2A5','Cycling'='#FC8D62','Remodeling'='#E78AC3',
    'Aberrant'='#8DA0CB','Pericytes'='#A6D854')

# Run Velocyto-R pipeline
library(velocyto.R)
library('pagoda2')

# Setup with Pagoda2
DefaultAssay(object) <- "integrated"
Idents(object) <- object$Final1
object <- subset(object, downsample = 5000)
object <- FindVariableFeatures(object)
length(object@assays$integrated@var.features)
object <- ScaleData(object)
nrow(GetAssayData(object,slot="scale.data"))

r <- Pagoda2$new(as(GetAssayData(object,slot="scale.data"),"dgCMatrix"),modelType="raw")
r$reductions[['PCA']] <- object@reductions$pca
r$reductions[['umap']] <- object@reductions$umap
r$embeddings[['PCA']][['umap']] <- Embeddings(object[['umap']])

cluster.label <- object$Final1
cell.colors <- sccore:::fac2col(cluster.label,
        level.colors=type_colors[levels(object$Final1)])

png(paste(dir_name,'_pagumap_2.png',sep=''), width = 800, height = 800)
r$plotEmbedding(type='PCA',embeddingType='umap',colors=cell.colors,mark.clusters=T,
        min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,
        main='cell clusters')
dev.off()

# Calculating Velocity
emat <- object$spliced; nmat <- object$unspliced;
emat <- emat[,rownames(r$counts)]; nmat <- nmat[,rownames(r$counts)];

cluster.label <- object$Final1
cell.colors <- sccore:::fac2col(cluster.label,
        level.colors=type_colors[levels(object$Final1)])

cell.dist <- as.dist(1-armaCor(t(r$embeddings[['PCA']][['umap']])))

emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.9)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.3)
length(intersect(rownames(emat),rownames(nmat)))
# [1] 885

fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,
        cell.dist=cell.dist,fit.quantile=fit.quantile)
emb <- r$embeddings$PCA$umap

x <- show.velocity.on.embedding.cor(emb,rvel.cd,n=1000,n.cores=1,scale='sqrt',
        cell.colors=ac(cell.colors,alpha=0.5),cex=0.8,arrow.scale=8,show.grid.flow=T,
        min.grid.cell.mass=5,grid.n=20,arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)

res <- 400
png(paste(dir_name,'_pagvel_4.png',sep=''), width=650*res/72, height=650*res/72,res=res)
show.velocity.on.embedding.cor(emb,rvel.cd,n=800,n.cores=1,scale='sqrt',cc=x$cc,
        cell.colors=ac(cell.colors),cex=1,arrow.scale=9,show.grid.flow=T,
        min.grid.cell.mass=40,grid.n=15,arrow.lwd=4,do.par=F,cell.border.alpha = 0.1,
        xlab = c('UMAP_1'),ylab= c('UMAP_2'))
dev.off()

save(x,file=paste(dir_name,'_velocity_',Sys.Date(),".Robj",sep=""))

# Velocyto on Remodeling/Pericyte subset of mesenchyme
# Subset based on embedding
cut <- Embeddings(object[['umap']])[which(Embeddings(object[['umap']])[,'UMAP_1'] > 1),]
cut <- cut[which(cut[,'UMAP_1'] < 5),] 
cut <- cut[which(cut[,'UMAP_2'] > 1),]
cut <- cut[which(cut[,'UMAP_2'] < 5),]

sub <- subset(object,cells=rownames(cut))  # 966 cells

# Setup with Pagoda2
DefaultAssay(sub) <- "integrated"
Idents(sub) <- sub$Final1
sub <- FindVariableFeatures(sub)
length(sub@assays$integrated@var.features)
sub <- ScaleData(sub)
nrow(GetAssayData(sub,slot="scale.data"))

r <- Pagoda2$new(as(GetAssayData(sub,slot="scale.data"),"dgCMatrix"),modelType="raw")
r$reductions[['PCA']] <- sub@reductions$pca
r$reductions[['umap']] <- sub@reductions$umap
r$embeddings[['PCA']][['umap']] <- Embeddings(sub[['umap']])

cluster.label <- sub$Final1
cell.colors <- sccore:::fac2col(cluster.label,
        level.colors=type_colors[levels(sub$Final1)])

png(paste(dir_name,'_pagumap_sub_2.png',sep=''), width = 800, height = 800)
r$plotEmbedding(type='PCA',embeddingType='umap',colors=cell.colors,mark.clusters=T,
        min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,
        main='cell clusters')
dev.off()

# Calculating Velocity
emat <- sub$spliced; nmat <- sub$unspliced;
emat <- emat[,rownames(r$counts)]; nmat <- nmat[,rownames(r$counts)];

cluster.label <- sub$Final1
cell.colors <- sccore:::fac2col(cluster.label,
        level.colors=type_colors[levels(sub$Final1)])

cell.dist <- as.dist(1-armaCor(t(r$embeddings[['PCA']][['umap']])))

emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.9)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.3)
length(intersect(rownames(emat),rownames(nmat)))
# [1] 952

fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,
        cell.dist=cell.dist,fit.quantile=fit.quantile)
emb <- r$embeddings$PCA$umap

x <- show.velocity.on.embedding.cor(emb,rvel.cd,n=200,n.cores=1,scale='sqrt',
        cell.colors=ac(cell.colors,alpha=0.5),cex=0.8,arrow.scale=8,show.grid.flow=T,
        min.grid.cell.mass=5,grid.n=20,arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)

res <- 400
png(paste(dir_name,'_pagvel_sub_2.png',sep=''), width=400*res/72, height=400*res/72,res=res)
show.velocity.on.embedding.cor(emb,rvel.cd,n=35,n.cores=1,scale='sqrt',cc=x$cc,
        cell.colors=ac(cell.colors),cex=1,arrow.scale=5,show.grid.flow=T,
        min.grid.cell.mass=5,grid.n=15,arrow.lwd=4,do.par=F,cell.border.alpha = 0.1,
        min.arrow.size = 0.1,xlab = c('UMAP_1'),ylab= c('UMAP_2'))
dev.off()



### Plot velocyto colored by slingshot
load('./sling.metadata.2023-02-14.Robj')
object <- AddMetaData(object,metadata = sling.metadata,col.name=names(sling.metadata))
lins <- data.frame(object$Lineage1_ps,object$Lineage2_ps,object$Lineage3_ps)
object$Lineage_avg <- rowMeans(lins,na.rm=TRUE)

cluster.label <- sort(object$Lineage_avg)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(ncol(object))
cell.colors <- sccore:::fac2col(cluster.label,shuffle=FALSE,
        level.colors=colors)

png(paste(dir_name,'_pagvel_2.png',sep=''), width = 800, height = 800)
show.velocity.on.embedding.cor(emb,rvel.cd,n=1000,n.cores=1,scale='sqrt',cc=x$cc,
        cell.colors=ac(cell.colors,alpha=0.5),cex=1,arrow.scale=12,show.grid.flow=T,
        min.grid.cell.mass=50,grid.n=15,arrow.lwd=4,do.par=F,cell.border.alpha = 0.1)
dev.off()




#### IMMUNE ####

### Archetype and cell type definition of integrated object clustering
## Re-integrate immune object including starting cells and engineered populations
# Load object containing putative labels from previous analyses
dir_name <- c('imm')
load('./imm.seurat.objs.int.2023-01-18.Robj')

Idents(object) <- object$Sample
object <- RenameIdents(object,'MacAlv'='BAL')
object$Sample <- Idents(object)
object$Sample <- factor(object$Sample,levels=c('FB13','FB14','BAL','BEF12','BEFM1',
    'BEFM2','BEFM4','BEFM5','BEFM6'))
object$Dataset2 <- factor(object$Dataset2,levels=c('Start','Tri_L','Quad_E','Quad_L'))

# Re-embed and -cluster on existing integration
DefaultAssay(object) <- "integrated"
object <- ScaleData(object)
object <- RunPCA(object, npcs = 100, verbose = F)

pcs <- 10
res <- 0.5
tag <- paste('pcs',pcs,'res',res,sep='.')
DefaultAssay(object) <- "integrated"
object <- FindNeighbors(object, dims = 1:pcs)
object <- FindClusters(object, resolution = res)
object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
p1 <- UMAPPlot(object,group.by = 'seurat_clusters',label = T,pt.size=2) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoAxes()
p2 <- UMAPPlot(object,group.by = 'Sample',pt.size=2) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
p3 <- UMAPPlot(object,group.by = 'putative_labels',pt.size=2) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
Idents(object) <- object$Sample
toprow <- plot_grid(p1,p2,p3,ncol=3)

data <- prop.table(table(object$seurat_clusters,object$Sample),1)*100
data <- as.data.frame(data)
p5 <- ggplot(data, aes(x = Var1, y = Freq, fill = Var2)) + NoLegend() +
  geom_col(position = "dodge") + scale_fill_manual(values = hue_pal()(length(unique(object$Sample)))) +
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
midrow1 <- plot_grid(p5)

DefaultAssay(object) <- "RNA"
Idents(object) <- object$seurat_clusters
p6 <- FeaturePlot(object, c('Mrc1','Tac1','Naaa','Top2a'),ncol=2,label = F,order=T)  # Alveolar_Mac
p7 <- FeaturePlot(object, c('C1qc','C1qb','Cd163','Stab1'),ncol=2,label = F,order=T)  # Interstitial_Mac
p8 <- FeaturePlot(object, c('Defb14','Capn3','Eno3','Slamf7'),ncol=2,label = F,order=T)  # Mono_NC_Act
p9 <- FeaturePlot(object, c('Rnase2','Hrh3','F5','Flt3'),ncol=2,label = F,order=T)  # Mono_C_Act
p10 <- FeaturePlot(object, c('Ms4a1','Cd19','Cd79b','Ebf1'),ncol=2,label = F,order=T)  # B
p11 <- FeaturePlot(object, c('Cd3d','Cd3e','Cd8b','Cd27'),ncol=2,label = F,order=T)  # T
p12 <- FeaturePlot(object, c('S100a8','S100a9','Nkg7','Gzma'),ncol=2,label = F,order=T)  # Neutro/NK
p13 <- FeaturePlot(object, c('Il17rb','Il2ra','Gata3','Ltb4r'),ncol=2,label = F,order=T)  # ILC
p14 <- FeaturePlot(object, c('Siglech','Smim5','Kmo','Flt3'),ncol=2,label = F,order=T)  # DC
p15 <- FeaturePlot(object, c('Tph1','Jchain','Asz1','Dnase2b'),ncol=2,label = F,order=T)  # Mast/Plasma/Eos/Baso
midrow2 <- plot_grid(p6,p7,p8,ncol=3)
midrow3 <- plot_grid(p9,p10,p11,ncol=3)
lastrow <- plot_grid(p12,p13,p14,p15,ncol=4)

# Cowplot all panels
plots <- list(toprow,midrow1,midrow2,midrow3,lastrow)
png(paste(dir_name,tag,'_integrated8_allplots.png',sep=""), width = 2000,height=2750)
plot_grid(plotlist=plots,ncol=1,rel_heights = c(2,1,2,2,2))
dev.off()


## Remake putative labels
table(object$seurat_clusters,object$putative_labels)
p1 <- UMAPPlot(object,group.by = 'seurat_clusters') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20)) + NoAxes()
p2 <- UMAPPlot(object,group.by = 'putative_labels') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
toprow <- plot_grid(p1,p2,ncol=2)
png(paste(dir_name,tag,'_integrated8_subcluster.png',sep=""), width = 1000,height=500)
plot_grid(plotlist=list(toprow),ncol=1)
dev.off()

Idents(object) <- object$putative_labels
ilc <- subset(object,idents='ILC')
mac_int <- subset(object,idents='Mac_Int')
mo_act <- subset(object,idents='Mo_Act')
mo_nc <- subset(object,idents='Mo_NC')
neutro <- subset(object,idents='Neutro')
plasma <- subset(object,idents='Plasma')
mac_alv_pro <- subset(object,idents='Mac_Alv_Pro')

Idents(object) <- object$seurat_clusters
# create new column in metadata
object$sub_cluster <- as.character(object$seurat_clusters)
object$sub_cluster[Cells(ilc)] <- paste(Idents(ilc))
object$sub_cluster[Cells(mac_int)] <- paste(Idents(mac_int))
object$sub_cluster[Cells(mo_act)] <- paste(Idents(mo_act))
object$sub_cluster[Cells(mo_nc)] <- paste(Idents(mo_nc))
object$sub_cluster[Cells(neutro)] <- paste(Idents(neutro))
object$sub_cluster[Cells(plasma)] <- paste(Idents(plasma))
object$sub_cluster[Cells(mac_alv_pro)] <- paste(Idents(mac_alv_pro))

table(object$sub_cluster,object$putative_labels)
p1 <- UMAPPlot(object,group.by = 'sub_cluster') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20)) + NoAxes()
p2 <- UMAPPlot(object,group.by = 'putative_labels') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
toprow <- plot_grid(p1,p2,ncol=2)
png(paste(dir_name,tag,'_integrated8_subcluster2.png',sep=""), width = 1000,height=500)
plot_grid(plotlist=list(toprow),ncol=1)
dev.off()

Idents(object) <- object$sub_cluster
object <- RenameIdents(object,
                '0'='Mac_Alv',
                '1'='Mac_Alv',
                '2'='Mac_Alv',
                '3'='Mac_Alv',
                '4'='Mac_Alv',
                '5'='Mac_Alv',
                '6'='Mac_Alv',
                '8'='T',
                '9'='B')
object$putative_labels <- Idents(object)

table(object$putative_labels)
object$putative_labels <- factor(object$putative_labels,levels=sort(levels(object$putative_labels)))


## Get new marker list for immune putative_labels
Idents(object) <- object$putative_labels
markers <- FindAllMarkers(object, only.pos = T)
markers$ratio <- markers$pct.1/markers$pct.2
markers <- markers[order(-markers$ratio),]

save(markers,file=paste(dir_name,'seurat.objs.int.marks',Sys.Date(),"Robj",sep="."))

save(object,file = paste(dir_name,"seurat.objs.int",Sys.Date(),"Robj",sep="."))



### Macrophage Analysis
## Pull and re-integrate macrophage populations from immune object
# Load object containing putative labels from previous analyses and subset macrophages
dir_name <- c('imm')
load('./imm.seurat.objs.int.2023-01-27.Robj')

Idents(object) <- object$putative_labels
mac <- subset(object,idents=c('Mac_Alv','Mac_Int','Mac_Alv_Pro'))  # 5199 cells

# Re-embed and -cluster on existing integration
DefaultAssay(mac) <- "integrated"
mac <- ScaleData(mac)
mac <- RunPCA(mac, npcs = 100, verbose = F)

pcs <- 15
res <- 0.4
tag <- paste('pcs',pcs,'res',res,sep='.')
DefaultAssay(mac) <- "integrated"
mac <- FindNeighbors(mac, dims = 1:pcs)
mac <- FindClusters(mac, resolution = res)
mac <- RunUMAP(mac, reduction = "pca", dims = 1:pcs)
p1 <- UMAPPlot(mac,group.by = 'seurat_clusters',label = T,pt.size=2) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(mac))) +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoAxes()
p2 <- UMAPPlot(mac,group.by = 'Sample',pt.size=2) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
p3 <- UMAPPlot(mac,group.by = 'putative_labels',pt.size=2) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
Idents(mac) <- mac$Sample
toprow <- plot_grid(p1,p2,p3,ncol=3)

data1 <- prop.table(table(mac$Sample,mac$seurat_clusters),1)*100
data1 <- as.data.frame(data1)
p5 <- ggplot(data1, aes(x = Var1, y = Freq, fill = Var2)) + NoLegend() +
  geom_col(position = "dodge") + scale_fill_manual(values = hue_pal()(length(unique(mac$seurat_clusters)))) +
  labs(title = c('Sample Composition by Cluster')) + xlab("Sample") + ylab("% sample per cluster") +
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
        
data2 <- prop.table(table(mac$seurat_clusters,mac$Sample),1)*100
data2 <- as.data.frame(data2)
p6 <- ggplot(data2, aes(x = Var1, y = Freq, fill = Var2)) + NoLegend() +
  geom_col(position = "dodge") + scale_fill_manual(values = hue_pal()(length(unique(mac$Sample)))) +
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
midrow1 <- plot_grid(p5,p6,ncol=2)

DefaultAssay(mac) <- "RNA"
Idents(mac) <- mac$seurat_clusters
p6 <- FeaturePlot(mac, c('Mrc1','Cx3cr1','Siglec1','Arg1'),ncol=2,order=T)  # Mononuclear Phagocytes
p6 <- FeaturePlot(mac, c('Naaa','Siglec10','Gpd1','Nrcam'),ncol=2,order=T)  # Alveolar Macrophage
p7 <- FeaturePlot(mac, c('C1qc','C1qb','Cd163','Stab1'),ncol=2,order=T)  # Interstitial Macrophage
p8 <- FeaturePlot(mac, c('G0s2','Eno3','Itgax','Cd79b'),ncol=2,order=T)  # Nonclassical Monocyte
p9 <- FeaturePlot(mac, c('Sell','Rnase2','F5','Vcan'),ncol=2,order=T)  # Classical Monocyte
p10 <- FeaturePlot(mac, c('Slamf7','Jaml','Cd200','Cd226'),ncol=2,order=T)  # Activated Monocyte
p11 <- FeaturePlot(mac, c('Defb14','Capn3','Eno3','Top2a'),ncol=2,order=T)  # Mono_NC_Act
p12 <- FeaturePlot(mac, c('Rnase2','Hrh3','F5','Flt3'),ncol=2,order=T)  # Mono_C_Act
p13 <- FeaturePlot(mac, c('Fcgr1a','Ccl20','Cxcl10','Cxcl11'),ncol=2,order=T)  # M1
p14 <- FeaturePlot(mac, c('Ccl24','Cxcr4','Tgfb1','Mmp9'),ncol=2,order=T)  # M2
midrow2 <- plot_grid(p6,p7,p8,ncol=3)
midrow3 <- plot_grid(p9,p10,p11,ncol=3)
lastrow <- plot_grid(p12,p13,p14,ncol=3)

# Cowplot all panels
plots <- list(toprow,midrow1,midrow2,midrow3,lastrow)
png(paste(dir_name,'mac',tag,'integrated_allplots.png',sep="_"), width = 2000,height=2750)
plot_grid(plotlist=plots,ncol=1,rel_heights = c(2,1,2,2,2))
dev.off()


#### LEFT OFF HERE ####


### Evaluate clusters - by seurat_clusters
p1 <- UMAPPlot(mac,group.by = 'seurat_clusters',label = T) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(mac))) +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoAxes()
p2 <- UMAPPlot(mac,group.by = 'Sample') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
p3 <- UMAPPlot(mac,group.by = 'putative_labels') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
toprow <- plot_grid(p1,p2,p3,NULL,ncol=4)

clusters <- levels(mac$seurat_clusters)
goi <- c(); marks <- c(); genes <- c()
myplots <- vector('list',length(clusters))
for (i in 1:length(clusters)) {
    cluster <- clusters[[i]]
    message(cluster)

    goi <- markers[which(markers$cluster == cluster),]
    goi <- goi[1:15, ]  # 15
    row.names(goi) <- goi$gene
    marks <- round(goi[, 2:5], digits = 4)
    marks$ratio <- round(goi$ratio, digits = 4)
    marks <- as.data.frame(marks)
    p <- tableGrob(marks,theme = ttheme_minimal(base_size = 16))

    #marks <- marks[order(marks$avg_log2FC),]
    genes <- row.names(marks)
    q <- FeaturePlot(mac,c(genes[1:4]),ncol=2,order=T)
    myplots[[i]] <- plot_grid(p,q,ncol=2)
}
lastrow <- plot_grid(plotlist=myplots,ncol=2,labels=paste('Cluster',clusters),label_size = 40)

# Cowplot all panels
plots <- list(toprow,lastrow)
png(paste(dir_name,'mac',tag,'integrated_clusters_1.png',sep="_"), width = 2000,height=2750)
plot_grid(plotlist=plots,ncol=1,rel_heights = c(1,5))
dev.off()





### Investigate markers - by Dataset3 (Start/Eng)
Idents(mac) <- mac$Dataset2
mac <- RenameIdents(mac,
                'Tri_L'='Eng',
                'Quad_E'='Eng',
                'Quad_L'='Eng')
mac$Dataset3 <- Idents(mac)
mac$Dataset3 <- factor(mac$Dataset3,levels=c('Start','Eng'))

Idents(mac) <- mac$Dataset3
markers <- FindAllMarkers(mac, only.pos = T)
markers$ratio <- markers$pct.1/markers$pct.2
markers <- markers[order(-markers$ratio),]
save(markers,file=paste(dir_name,'mac','seurat.objs.int.marks',Sys.Date(),"Robj",sep="."))
write.table(markers,file = paste(dir_name,'mac','seurat.objs.int.marks',Sys.Date(),"txt",sep="."),sep="\t")
load('./imm.mac.seurat.objs.int.marks..2023-01-30.Robj')  # differentiate macrophages by Start/Eng

markers[which(markers$gene == 'Il10'),]


## Evaluate clusters - by Dataset3 (Start/Eng)
p1 <- UMAPPlot(mac,group.by = 'Dataset3',label = T) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(mac))) +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoAxes()
p2 <- UMAPPlot(mac,group.by = 'Sample') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
p3 <- UMAPPlot(mac,group.by = 'putative_labels') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
toprow <- plot_grid(p1,p2,p3,NULL,ncol=4)

clusters <- levels(mac$Dataset3)
goi <- c(); marks <- c(); genes <- c()
myplots <- vector('list',length(clusters))
for (i in 1:length(clusters)) {
    cluster <- clusters[[i]]
    message(cluster)

    goi <- markers[which(markers$cluster == cluster),]
    row.names(goi) <- goi$gene
    marks <- round(goi[, 2:5], digits = 4)
    marks$ratio <- round(goi$ratio, digits = 4)
    marks <- as.data.frame(marks)
    marks_r <- marks[1:40, ]  # 15
    p <- tableGrob(marks_r,theme = ttheme_minimal(base_size = 16))
    marks_a <- marks[order(-marks$avg_log2FC),]
    marks_a <- marks_a[1:40, ]  # 15
    q <- tableGrob(marks_a,theme = ttheme_minimal(base_size = 16))
    pq <- plot_grid(p,q,nrow=2)

    #marks <- marks[order(marks$avg_log2FC),]
    genes <- row.names(marks)
    r <- FeaturePlot(mac,c(genes[1:20]),ncol=2,order=T)
    myplots[[i]] <- plot_grid(pq,r,ncol=2)
}
lastrow <- plot_grid(plotlist=myplots,ncol=2,labels=paste('Cluster',clusters),label_size = 40)

# Cowplot all panels
plots <- list(toprow,lastrow)
png(paste(dir_name,'mac',tag,'integrated_clusters_3.png',sep="_"), width = 2000,height=2750)
plot_grid(plotlist=plots,ncol=1,rel_heights = c(1,5))
dev.off()




### Investigate M1/M2 Mac_Alv heterogeneity - by Start/Eng
DefaultAssay(mac) <- "RNA"
## Shaykhiev 2009 (human) - lung-specific
m1.genes <- c('Fcgr3a','Fcgr2a','Fcgr1a','Il15ra','C3ar1','Cd69','Cd80','Cd86','Tlr2','Tlr4','Icam1','Il1b',
    'Il6','Il12b','Il18','Il23a','Il32','Tnf','Tnfsf10','Tnfaip6','Cxcl1','Cxcl9','Cxcl10','Cxcl11','Ccl4',
    'Ccl5','Ccl20','Gbp1','Gbp2','Gbp3','Gbp4','Gbp5','Acod1','Irf1','Irf7','Socs3','Nos2','Pde4b','Apol3','Cfb')
m2.genes <- c('Msr1','Mrc1','Mrc2','Cxcr4','Ccr5','Clec7a','Lpar6','Stab1','Cd9','Cd36','Cd163','Mertk',
    'Adora3','Fcer2','Il4r','Il10','Il1rn','Tgfb1','Ccl17','Ccl18','Ccl22','Ccl23','Ccl24','Rgs1','Gas7',
    'Arg1','Mmp2','Mmp7','Mmp9','Hs3st1','Hs3st2','Col6a2','Fn1')
# Score macrophage activation
genes.use <- intersect(m1.genes,rownames(mac))
mac <- AddModuleScore(mac, features = list(genes.use),name = 'm1.score')
genes.use <- intersect(m2.genes,rownames(mac))
mac <- AddModuleScore(mac, features = list(genes.use),name = 'm2.score')

Idents(mac) <- mac$seurat_clusters
p1 <- FeaturePlot(mac,'m1.score1',cols=c('lightgrey','red'),label = T,order=T)
p2 <- FeaturePlot(mac,'m2.score1',cols=c('lightgrey','red'),label = T,order=T)
toprow <- plot_grid(p1,p2,ncol=2)
p3 <- VlnPlot(mac,c('m1.score1','m2.score1'),group.by='Dataset3',ncol=2,pt.size=0.1)

png(paste(dir_name,'mac',tag,'integrated_m1m2_1.png',sep="_"), width = 1000,height=1000)
plot_grid(plotlist=list(toprow,p3),nrow=2)
dev.off()

# Save Shaykhiev M1/M2 metadata
tosave <- c('m1.score1','m2.score1')
mac.metadata <- mac@meta.data
mac.metadata <- mac.metadata[c(tosave)]
save(mac.metadata,file = paste("mac.metadata.",Sys.Date(),".Robj",sep=""))

#load('./mac.metadata.2023-02-03.Robj')  # calculated incorrectly!
load('./mac.metadata.2023-03-09.Robj')  # calculated correctly w list()
object <- AddMetaData(object,metadata = mac.metadata,col.name=names(mac.metadata))




## Weagel 2015
m1.genes <- c('Il12a','Il12b','Il23a','Tnf','Irf3','Jun','Nfkb1','Nfkb2','Cxcl10',
    'Nos2','Stat1','Cxcl9','Irf1','Socs1')
m2a.genes <- c('Ccl17','Arg1','Irf4','Il10','Socs3','Tgfb1','Tgfb2')
m2b.genes <- c('Il10','Ccl1','Il1a','Il1b','Il6','Tnf')
m2c.genes <- c('Il10','Tgfb1','Tgfb2','Cd163','Mrc1')
m2d.genes <- c('Vegfa','Vegfb','Il10','Il12a','Il12b','Tnf','Tgfb1','Tgfb2')
m2.genes <- unique(c(m2a.genes,m2b.genes,m2c.genes,m2d.genes))
# Score macrophage activation
mac3 <- mac
genes.use <- intersect(m1.genes,rownames(mac3))
mac3 <- AddModuleScore(mac3, features = list(genes.use),name = 'm1.score')
genes.use <- intersect(m2.genes,rownames(mac3))
mac3 <- AddModuleScore(mac3, features = list(genes.use),name = 'm2.score')
genes.use <- intersect(m2a.genes,rownames(mac3))
mac3 <- AddModuleScore(mac3, features = list(genes.use),name = 'm2a.score')
genes.use <- intersect(m2b.genes,rownames(mac3))
mac3 <- AddModuleScore(mac3, features = list(genes.use),name = 'm2b.score')
genes.use <- intersect(m2c.genes,rownames(mac3))
mac3 <- AddModuleScore(mac3, features = list(genes.use),name = 'm2c.score')
genes.use <- intersect(m2d.genes,rownames(mac3))
mac3 <- AddModuleScore(mac3, features = list(genes.use),name = 'm2d.score')

Idents(mac3) <- mac3$seurat_clusters
p1 <- FeaturePlot(mac3,'m1.score1',cols=c('lightgrey','red'),label = T,order=T)
p2 <- FeaturePlot(mac3,'m2.score1',cols=c('lightgrey','red'),label = T,order=T)
toprow <- plot_grid(p1,p2,NULL,ncol=3)
p3 <- VlnPlot(mac3,c('m1.score1','m2.score1'),group.by='Dataset3',ncol=2,pt.size=0.1)
midrow1 <- plot_grid(p3,NULL,ncol=2,rel_widths=c(2,1))

p4 <- FeaturePlot(mac3,'m2d.score1',cols=c('lightgrey','red'),label = T,order=T)
p5 <- FeaturePlot(mac3,'m2b.score1',cols=c('lightgrey','red'),label = T,order=T)
p6 <- FeaturePlot(mac3,'m2c.score1',cols=c('lightgrey','red'),label = T,order=T)
midrow2 <- plot_grid(p4,p5,p6,ncol=3)
p7 <- VlnPlot(mac3,c('m2d.score1','m2b.score1','m2c.score1'),group.by='Dataset3',ncol=3,pt.size=0.1)
lastrow <- plot_grid(p7)

plots <- list(toprow,midrow1,midrow2,lastrow)
png(paste(dir_name,'mac',tag,'integrated_m1m2_4.png',sep="_"), width = 1000,height=1500)
plot_grid(plotlist=plots,ncol=1)
dev.off()


## Murray 2014 - Not great
m1.genes <- c('Stat1','Socs1','Nfkbiz','Irf5','Irf1','Tnf','Il6','Il27','Il23a',
    'Il12a','Il1b','Il12b','Cxcl10','Cxcl8','Ccl5','Cxcl9','Cxcl11','Marco',
    'Mmp9','Arg1','Nos2','Ido1','Kynu','Ptx3','Gbp1','Ccr7','Cd40')
m2a.genes <- c('Stat6','Irf4','Socs2','Socs1','Gata3','Ccl17','Ccl24','Ccl22','Ccl4',
    'Ccl13','Ccl17','Ccl18','Mrc1','Stab1','Fn1','Tgfb1','Mmp1','Mmp12','Tg','F13a1',
    'Arg1','Retnla','Chi3l3','Alox15','Tgm2','Adora3','Il17rb','Cd200r1','Cd200r1l')
m2b.genes <- c('Il10','Il6','Cxcl13','Ccl1','Ccl20','Nos2')
m2c.genes <- c('Stat3','Nfil3','Sbno2','Socs3','Il10','Il4r')
m2c2.genes <- c(m2c.genes,'Id3','Rgs1','Smad2','Cd163','Stab1','Marco','F13a1','Tgfbr2',
    'Alox5ap','Il17rb','Adora3')
m2.genes <- unique(c(m2a.genes,m2b.genes,m2c2.genes))
# Score macrophage activation
mac4 <- mac
genes.use <- intersect(m1.genes,rownames(mac4))
mac4 <- AddModuleScore(mac4, features = list(genes.use),name = 'm1.score')
genes.use <- intersect(m2.genes,rownames(mac4))
mac4 <- AddModuleScore(mac4, features = list(genes.use),name = 'm2.score')
genes.use <- intersect(m2a.genes,rownames(mac4))
mac4 <- AddModuleScore(mac4, features = list(genes.use),name = 'm2a.score')
genes.use <- intersect(m2b.genes,rownames(mac4))
mac4 <- AddModuleScore(mac4, features = list(genes.use),name = 'm2b.score')
genes.use <- intersect(m2c.genes,rownames(mac4))
mac4 <- AddModuleScore(mac4, features = list(genes.use),name = 'm2c.score')
genes.use <- intersect(m2c2.genes,rownames(mac4))
mac4 <- AddModuleScore(mac4, features = list(genes.use),name = 'm2c2.score')

Idents(mac4) <- mac4$seurat_clusters
p1 <- FeaturePlot(mac4,'m1.score1',cols=c('lightgrey','red'),label = T,order=T)
p2 <- FeaturePlot(mac4,'m2.score1',cols=c('lightgrey','red'),label = T,order=T)
p3 <- FeaturePlot(mac4,'m2a.score1',cols=c('lightgrey','red'),label = T,order=T)
toprow <- plot_grid(p1,p2,p3,ncol=3)
p4 <- VlnPlot(mac4,c('m1.score1','m2.score1','m2a.score1'),group.by='Dataset3',ncol=3,pt.size=0.1)
midrow1 <- plot_grid(p4)

p5 <- FeaturePlot(mac4,'m2b.score1',cols=c('lightgrey','red'),label = T,order=T)
p6 <- FeaturePlot(mac4,'m2c.score1',cols=c('lightgrey','red'),label = T,order=T)
p7 <- FeaturePlot(mac4,'m2c2.score1',cols=c('lightgrey','red'),label = T,order=T)
midrow2 <- plot_grid(p5,p6,p7,ncol=3)
p8 <- VlnPlot(mac4,c('m2b.score1','m2c.score1','m2c2.score1'),group.by='Dataset3',ncol=3,pt.size=0.1)
lastrow <- plot_grid(p8)

plots <- list(toprow,midrow1,midrow2,lastrow)
png(paste(dir_name,'mac',tag,'integrated_m1m2_5.png',sep="_"), width = 1000,height=1500)
plot_grid(plotlist=plots,ncol=1)
dev.off()







# focus_imm_7.R
#### Focused Analysis of Engineered Immune
## Parse Macrophage Heterogeneity


# focus_imm_8.R
#### Focused Analysis of Engineered Immune
## Pseudotime by Slingshot & Velocyto


# focus_imm_9.R
#### Focused Analysis of Engineered Immune
## Re-do Velocyto




#### ADDITIONAL ####

# focus_rus_1.R
#### Focused Analysis guided by Ruslan
## May 2024, Regional ECM expression in engineered epithelium & ECM gene list


# focus_rus_2.R
#### Focused Analysis guided by Ruslan
## May 2024, Generating native v2 object to evaluate regional ECM expression in native


# focus_rus_3.R
#### Focused Analysis guided by Ruslan
## May 2024, Regional ECM expression in native & engineered epithelium


# focus_rus_4.R
#### Focused Analysis guided by Ruslan
## June 2024, Regional ECM expression in native & engineered epithelium WO starting peBC


# focus_rus_5.R
#### Focused Analysis guided by Ruslan
## June 2024, ECM Figure


# focus_rus_6.R
#### Focused Analysis guided by Ruslan
## June 2024, Functional Delegation
