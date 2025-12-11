#### EXPLORE
#### Code for publication unifying cleaned objects into integrated class objects
#### Replicable method for exploring and integrating cleaned objects by class




## Note: We tried generating unified class objects by merge, reference mapping, or CCA integration, in that order.
##  In all classes, useful sample alignment was only achieved by integration. Integration was iterated across 
##  multiple sample combinations, including with samples that are not included in this manuscript. All scripts have
##  been edited to run only with published samples. Only code that contributed to final processed objects is shown.
##  A sample of initial exploratory analyses are also included.




#### Assembling relevant datasets and lists of objects
## Generate list of all individual cleaned objects
allobjects <- c("BC1P3","BC1P6","RLMVEC","FB13","FB14","BAL","BCL5","BCEC2","BEF1","BEF2","BEF3",
    "BEF12","BEF14","BEF15","BEFM1","BEFM2","BEFM4","BEFM5","BEFM6")
seurat.data <- list()
for (i in 1:length(allobjects)){
    sample.name <- allobjects[i]
    message(sample.name)
    setwd(paste0("../",sample.name))
    load(list.files(path = ".", pattern = "clean.2022-07-23"))
    seurat.data[[i]] <- object
}
names(seurat.data) <- allobjects

save(seurat.data,file = paste0("bef.seurat.objs.",Sys.Date(),".Robj"))


## Object lists by cell class - Epithelium, Endothelium, Mesenchyme, Immune
epi <- list() ; epi.names <- c()
endo <- list() ; endo.names <- c()
mes <- list() ; mes.names <- c()
imm <- list() ; imm.names <- c()
for (i in 1:length(seurat.data)){
    object <- seurat.data[[i]]
    sample.name <- names(seurat.data[i])
    message(sample.name)
    Idents(object) <- object$class
    if ("Epithelium" %in% object$class){
        epi[i] <- subset(object,idents = c('Epithelium'))
        epi.names[i] <- sample.name
    }
    if ("Endothelium" %in% object$class){
        endo[i] <- subset(object,idents = c('Endothelium'))
        endo.names[i] <- sample.name
    }
    if ("Mesenchyme" %in% object$class){
        mes[i] <- subset(object,idents = c('Mesenchyme'))
        mes.names[i] <- sample.name
    }
    if ("Immune" %in% object$class){
        imm[i] <- subset(object,idents = c('Immune'))
        imm.names[i] <- sample.name
    }
}
names(epi) <- epi.names
names(endo) <- endo.names
names(mes) <- mes.names
names(imm) <- imm.names

# Remove empty slots from list
epi <- epi[lengths(epi) != 0]
endo <- endo[lengths(endo) != 0]
mes <- mes[lengths(mes) != 0]
imm <- imm[lengths(imm) != 0]

save(epi,file = paste0("epi.seurat.objs.",Sys.Date(),".Robj"))
save(endo,file = paste0("endo.seurat.objs.",Sys.Date(),".Robj"))
save(mes,file = paste0("mes.seurat.objs.",Sys.Date(),".Robj"))
save(imm,file = paste0("imm.seurat.objs.",Sys.Date(),".Robj"))


## Table of metadata info for each object
df <- data.frame()
Name <- c() ; Chemistry <- c() ; Cells <- c() ; Epithelium <- c() ; Endothelium <- c() ; Mesenchyme <- c() ; Immune <- c()
Mean_Counts <- c() ; Mean_Features <- c()
for (i in 1:length(seurat.data)){
    object <- seurat.data[[i]]
    sample.name <- names(seurat.data[i])
    message(sample.name)

    Name[i] <- sample.name
    Chemistry[i] <- unique(object$Chemistry)
    Mean_Counts[i] <- round(mean(object$nCount_RNA))
    Mean_Features[i] <- round(mean(object$nFeature_RNA))
    Cells[i] <- ncol(object)
    Epithelium[i] <- length(which(object$class == 'Epithelium'))
    Endothelium[i] <- length(which(object$class == 'Endothelium'))
    Mesenchyme[i] <- length(which(object$class == 'Mesenchyme'))
    Immune[i] <- length(which(object$class == 'Immune'))

}

df <- cbind(Name,Chemistry,Mean_Counts,Mean_Features,Cells,Epithelium,Endothelium,Mesenchyme,Immune)

p <- tableGrob(df,theme = ttheme_minimal(base_size = 20))

png(paste("bef.seurat.objs_clean_table.png"), width = 1200,height=800)
plot(p)
dev.off()




#### EPITHELIUM ####

### Initial exploration of epi objects
## Initial embeddings of each object
# Load all epithelium
load("./epi.seurat.objs.2022-07-27.Robj")

# Cowplot of every epithelial object UMAP
cols <- hue_pal()(4)
class_colors <- c('Epithelium'=cols[2],'Endothelium'=cols[3],'Mesenchyme'=cols[1],'Immune'=cols[4])

# Order objects
allobjects <- c("BC1P3","BC1P6","FB13","FB14","MacAlv","BCL5","BCEC2","BEF1","BEF2","BEF3",
    "BEF12","BEF14","BEF15","BEFM1","BEFM2","BEFM4","BEFM5","BEFM6","Native")
epi <- epi[allobjects]

# Plot without re-embedding
myplots <- vector('list',length(epi))
for (i in 1:length(epi)){
    object <- epi[[i]]
    sample.name <- names(epi[i])
    message(sample.name)
    p <- UMAPPlot(object,group.by = 'class',cols = class_colors) +
            labs(title = paste(sample.name,'- Epithelium'), subtitle=paste('nCells =',ncol(object))) +
            theme(plot.title = element_text(size = 24,hjust = 0))
    myplots[[i]] <- p
}

png(paste("epi_all_umaps2.png"), width = 3000,height=2000)
print(plot_grid(plotlist=myplots,ncol=6))
dev.off()

# Plot with re-embedding (skip bc too small: FB13, FB14, MacAlv)
epi2 <- epi[c("BC1P3","BC1P6","BCL5","BCEC2","BEF1","BEF2","BEF3","BEF12","BEF14","BEF15","BEFM1",
    "BEFM2","BEFM4","BEFM5","BEFM6","Native")]
pcs <- 10
myplots <- vector('list',length(epi2))
new.epi <- list() ; final.names <- c()
for (i in 1:length(epi2)) {
    object <- epi2[[i]]
    sample.name <- names(epi2[i])
    message(sample.name)
    object <- NormalizeData(object)
    object <- FindVariableFeatures(object)
    object <- ScaleData(object)
    object <- RunPCA(object, npcs = 50, verbose = F)
    object <- FindNeighbors(object, dims = 1:pcs)
    object <- FindClusters(object, resolution = 0.5)
    object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
    p <- UMAPPlot(object,label = T) + labs(title = paste(sample.name,'- Epithelium'),subtitle = paste("PC's =",pcs,"\nRes = 0.5",'\nCells =',ncol(object))) +
        NoLegend() + theme(plot.title = element_text(size = 24))
    myplots[[i]] <- p
    new.epi[[i]] <- object
    final.names[[i]] <- sample.name
}
names(new.epi) <- final.names
epi3 <- c(new.epi,epi[c('FB13','FB14','MacAlv')])
epi3 <- epi3[order(names(epi))]

png(paste("epi_all_umaps_22.png"), width = 3000,height=2000)
print(plot_grid(plotlist=myplots,ncol=6))
dev.off()


## Explore top native markers via FeaturePlots
natepi <- epi$Native
Idents(natepi) <- natepi$CellType_Final

markers <- FindAllMarkers(natepi, only.pos = T)
markers$ratio <- markers$pct.1/markers$pct.2
markers <- markers[order(-markers$ratio),]

# Epithelial markers of interest:
# Epi: ['Krt5','Hopx','Aqp5','Sftpc']
# ATII: 'Sftpc','Sftpb','Defb4','Lyz2'
# ATI: 'Ager','Aqp5','Hopx','Pdpn'
# Basal: 'Krt5','Krt8','Igfbp2','Sox2'
# Ciliated: 'Ak9','Pifo','Ccdc153'
# Secretory: 'Scgb1a1','Rarres1','Muc20','Clca1'
# Tuft1: 'Dclk1','Espn','Sox9','Pou2f3'
# Tuft2: 'Dclk1','Trpm5','Sox9','Lgr5'
# BASC: 'Sftpc','Scgb1a1','Sox9','Lgr5'

# Makes a sheet of individual gene fps across all objects
genes <- c('Krt5','Hopx','Aqp5','Sftpc')
myplots <- vector('list',length(epi3))
for (i in 1:length(genes)) {
    gene <- genes[i]
    message(gene)
    for (j in 1:length(epi3)) {
        object <- epi3[[j]]
        sample.name <- names(epi3[j])
        message(sample.name)
        p <- FeaturePlot(object,features=gene,label = F) &
            labs(title = paste(sample.name,'-',gene))
        myplots[[j]] <- p
    }
    png(paste("epi_all_fp_",gene,".png",sep=''),width = 2000,height=3000)
    print(plot_grid(plotlist=myplots,ncol=4))
    dev.off()
}

# Makes a sheet of a set of gene fps across all objects
type <- c('ATI')
genes <- c('Ager','Aqp5','Hopx','Pdpn')
myplots1 <- vector('list',length(genes))
myplots2 <- vector('list',length(epi3))
for (i in 1:length(epi3)) {
    object <- epi3[[i]]
    sample.name <- names(epi3[i])
    message(sample.name)
    for (j in 1:length(genes)) {
        gene <- genes[j]
        message(gene)
        p <- FeaturePlot(object,features=gene,label = F) +
            labs(title = paste(sample.name,'-',gene))
        myplots1[[j]] <- p
    }
    myplots2[[i]] <- plot_grid(plotlist=myplots1,ncol=2)
}

plot <- plot_grid(plotlist=myplots2,ncol=4)
title <- ggdraw() + draw_label(type,fontface = 'bold',size = 24,x = 0,hjust = 0)
png(paste("epi_all_fp_",type,".png",sep=''),width = 2000,height=3000)
print(plot_grid(title,plot,ncol=1,rel_heights = c(0.1, 10)))
dev.off()


## Explore top native markers via Heatmaps
load("/nobackup1/greaneya/DGEs_30K_Rank_Cutoff/Native/native.marks.2022-08-02.Robj")

# AverageExpression Heatmaps
# Select top genes per native epi type
table(new.epi$Native$CellType_Final)
native.epi.names <- c('ATII','ATI','Tuft','Ciliated','Secretory','ATII-ATI','BASC')
marks.short <- markers %>% group_by(cluster) %>% top_n(50, avg_log2FC) %>%
    select(cluster, gene, avg_log2FC, p_val,ratio) %>%
    arrange(cluster, desc(avg_log2FC))

# Loop to generate heatmaps of each native cell type with engineered (wo negative native ctrls)
epi4 <- new.epi[c("BC1P3","BC1P6","BCEC2","BCL5","BEF1","BEF12","BEF14","BEF15","BEF2","BEF3",
    "BEFM1","BEFM2","BEFM4","BEFM5","BEFM6")]  # minus native, tiny samples
myplots <- vector('list', length(native.epi.names))
unityNormalize <- function(x){(x-min(x))/(max(x)-min(x))}
for (i in 1:length(native.epi.names)) {
    hm.short <- list()
    type <- native.epi.names[i]
    message(type)
    identities <- c(type,names(epi4))
    ctrl <- subset(new.epi$Native,subset = CellType_Final == c(type))
    hm.object <- c(ctrl,epi4)
    names(hm.object) <- identities
    topgenes <- marks.short[which(marks.short$cluster==type), ]
    genes.use <- topgenes$gene
    for (j in 1:length(hm.object)) {
        object <- hm.object[[j]]
        genes.int <- intersect(topgenes$gene,rownames(object))
        genes.use <- intersect(genes.use,genes.int)
    }
    genes.use <- head(genes.use,50)
    for (k in 1:length(hm.object)) {
        object <- hm.object[[k]]
        Idents(object) <- object$Sample
        hm.short[[k]] <- subset(object,downsample = 120)
    }
    hm.merge <- merge(x=hm.short[[1]],y=hm.short[2:length(hm.short)])
    hm.merge <- ScaleData(hm.merge,features = genes.use)
    hm.avg <- AverageExpression(hm.merge,assays='RNA',slot='scale.data',features=genes.use)
    hm.avg <- t(scale(t(hm.avg$RNA)))
    hm.avg <- t(apply(as.matrix(hm.avg), 1, unityNormalize))
    colnames(hm.avg)[which(colnames(hm.avg)=='Native')] <- type
    hm.avg <- hm.avg[,c(identities)]
    message(dim(hm.avg))

    colors.inferno <- colorRamp2(breaks = c(seq(min(hm.avg),max(hm.avg),length.out=60)), inferno(n=60), space = "RGB")
    myplots[[i]] <- local({
        i <- i
        p1 <- Heatmap(as.matrix(hm.avg),
                                col= colors.inferno,
                                row_names_gp = gpar(fontsize = 18),
                                show_column_dend = T,
                                clustering_method_columns = 'complete',
                                cluster_rows=T,
                                cluster_columns=T,
                                cluster_column_slices=F,
                                show_column_names=T,
                                column_names_side = "top",
                                column_dend_side = "bottom",
                                column_names_gp = gpar(fontsize = 30),
                                show_row_names = T,
                                show_heatmap_legend = F)
    })
    png(paste('nativeepi_heatmap_',type,'_2.png',sep=""), width = 800, height = 1000)  #400,1000
    draw(myplots[[i]],padding = unit(c(2, 2, 18, 2), "mm"))  #12
    dev.off()
}


## Embedding-agnostic marker comparisons
# Merge without embedding
objects <- epi[c("BC1P3","BC1P6","BCL5","BCEC2","BEF1","BEF2","BEF3","BEF12","BEF14","BEF15","BEFM1","BEFM2",
    "BEFM4","BEFM5","BEFM6")]
dir_name <- c('epi')

object <- merge(x = objects[[1]], y = objects[2:length(objects)])
object <- NormalizeData(object)
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, npcs = 100, verbose = F)

object$Sample <- factor(object$Sample,levels = c("BC1P3","BC1P6","BCL5","BCEC2","BEF1","BEF2","BEF3","BEF12",
    "BEF14","BEF15","BEFM1","BEFM2","BEFM4","BEFM5","BEFM6"))

Idents(object) <- object$Sample
object <- RenameIdents(object,
        "BC1P3"='Start',
        "BC1P6"='Start',
        "BCL5"='Mono',
        "BCEC2"='Co',
        "BEF1"='Tri',
        "BEF2"='Tri',
        "BEF3"='Tri',
        "BEF12"='Tri',
        "BEF14"='Tri',
        "BEF15"='Tri',
        "BEFM1"='Quad',
        "BEFM2"='Quad',
        "BEFM4"='Quad',
        "BEFM5"='Quad',
        "BEFM6"='Quad')
object$Dataset2 <- Idents(object)

DefaultAssay(object) <- "RNA"
p1 <- VlnPlot(object,c('Krt5','Hopx','Aqp5','Sftpc'),group.by='Sample',ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())  # Epi
p2 <- VlnPlot(object,c('Sftpc','Sftpb','Defb4','Lyz2'),group.by='Sample',ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())  # ATII
p3 <- VlnPlot(object,c('Ager','Aqp5','Hopx','Pdpn'),group.by='Sample',ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())  # ATI
p4 <- VlnPlot(object,c('Krt5','Krt8','Igfbp2','Sox2'),group.by='Sample',ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())  # Basal
p5 <- VlnPlot(object,c('Ak9','Pifo','Ccdc153','Top2a'),group.by='Sample',ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())  # Ciliated
p6 <- VlnPlot(object,c('Scgb1a1','Rarres1','Muc20','Clca1'),group.by='Sample',ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())  # Secretory
p7 <- VlnPlot(object,c('Dclk1','Espn','Sox9','Pou2f3'),group.by='Sample',ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())  # Tuft1
p8 <- VlnPlot(object,c('Dclk1','Trpm5','Sox9','Lgr5'),group.by='Sample',ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())  # Tuft2
p9 <- VlnPlot(object,c('Sftpc','Scgb1a1','Sox9','Lgr5'),group.by='Sample',ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())  # BASC
toprow <- plot_grid(p1,p2,p3,ncol=3)
midrow2 <- plot_grid(p4,p5,p6,ncol=3)
midrow4 <- plot_grid(p7,p8,p9,ncol=3)

p10 <- VlnPlot(object,c('Krt5','Hopx','Aqp5','Sftpc'),group.by='Dataset2',ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())  # Epi
p11 <- VlnPlot(object,c('Sftpc','Sftpb','Defb4','Lyz2'),group.by='Dataset2',ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())  # ATII
p12 <- VlnPlot(object,c('Ager','Aqp5','Hopx','Pdpn'),group.by='Dataset2',ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())  # ATI
p13 <- VlnPlot(object,c('Krt5','Krt8','Igfbp2','Sox2'),group.by='Dataset2',ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())  # Basal
p14 <- VlnPlot(object,c('Ak9','Pifo','Ccdc153','Top2a'),group.by='Dataset2',ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())  # Ciliated
p15 <- VlnPlot(object,c('Scgb1a1','Rarres1','Muc20','Clca1'),group.by='Dataset2',ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())  # Secretory
p16 <- VlnPlot(object,c('Dclk1','Espn','Sox9','Pou2f3'),group.by='Dataset2',ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())  # Tuft1
p17 <- VlnPlot(object,c('Dclk1','Trpm5','Sox9','Lgr5'),group.by='Dataset2',ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())  # Tuft2
p18 <- VlnPlot(object,c('Sftpc','Scgb1a1','Sox9','Lgr5'),group.by='Dataset2',ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())  # BASC
midrow1 <- plot_grid(p10,p11,p12,ncol=3)
midrow3 <- plot_grid(p13,p14,p15,ncol=3)
lastrow <- plot_grid(p16,p17,p18,ncol=3)

# Cowplot all panels
plots <- list(toprow,midrow1,midrow2,midrow3,midrow4,lastrow)
png(paste(dir_name,'_merge_allplots.png',sep=""), width = 2000,height=2750)
plot_grid(plotlist=plots,ncol=1)
dev.off()

save(object,file = paste("epi.seurat.objs.mer.",Sys.Date(),".Robj",sep=""))



### Generating integrated object of all Epithelium
# Normal Seurat CCA Integration of all engineered objects, without native or tiny samples
objects <- epi[c("BC1P3","BC1P6","BCL5","BCEC2","BEF1","BEF2","BEF3","BEF12","BEF14","BEF15",
    "BEFM1","BEFM2","BEFM4","BEFM5","BEFM6")]
dir_name <- c('epi')
objects <- lapply(X = objects, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x,nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = objects, nfeatures = 2000)
object.anchors <- FindIntegrationAnchors(object.list=objects, anchor.features=features,
        reduction='cca',dims=1:30)
object <- IntegrateData(anchorset = object.anchors, dims=1:30, k.weight = 57)

# Perform integration
DefaultAssay(object) <- "integrated"
object <- ScaleData(object)
object <- RunPCA(object, npcs = 100, verbose = F)
png(paste(dir_name,'_integrated2_elbow.png',sep=""), width = 800, height = 500)
ElbowPlot(object,ndims = 50)+labs(title = dir_name)
dev.off()
pdf(paste(dir_name,'_integrated2_pcheatmap1.pdf',sep=""), width = 12, height = 16)
DimHeatmap(object, dims = 1:15, cells = 200, balanced = T)
dev.off()
pdf(paste(dir_name,'_integrated2_pcheatmap2.pdf',sep=""), width = 12, height = 16)
DimHeatmap(object, dims = 16:30, cells = 200, balanced = T)
dev.off()
pdf(paste(dir_name,'_integrated2_pcheatmap3.pdf',sep=""), width = 12, height = 16)
DimHeatmap(object, dims = 31:45, cells = 200, balanced = T)
dev.off()

pcs <- 20  # iterate
res <- 0.3  # iterate
tag <- paste('pcs',pcs,'res',res,sep='.')
DefaultAssay(object) <- "integrated"
object <- FindNeighbors(object, dims = 1:pcs)
object <- FindClusters(object, resolution = res)
object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
p1 <- UMAPPlot(object,group.by = 'seurat_clusters',label = T) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoAxes()
p2 <- UMAPPlot(object,group.by = 'Sample') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
Idents(object) <- object$Sample
toprow <- plot_grid(p1,p2,NULL,ncol=3)

data <- prop.table(table(object$seurat_clusters,object$Sample),1)*100
data <- as.data.frame(data)
p5 <- ggplot(data, aes(x = Var1, y = Freq, fill = Var2)) + NoLegend() +
  geom_col(position = "dodge", width = 1) + scale_fill_manual(values = hue_pal()(length(unique(object$Sample)))) +
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
p7 <- FeaturePlot(object, c('Sftpc','Sftpb','Defb4','Lyz2'),ncol=2,label = F)  # ATII
p8 <- FeaturePlot(object, c('Ager','Aqp5','Hopx','Pdpn'),ncol=2,label = F)  # ATI
p9 <- FeaturePlot(object, c('Krt5','Krt8','Igfbp2','Sox2'),ncol=2,label = F)  # Basal
p10 <- FeaturePlot(object, c('Ak9','Pifo','Ccdc153','Top2a'),ncol=2,label = F)  # Ciliated
p11 <- FeaturePlot(object, c('Scgb1a1','Rarres1','Muc20','Clca1'),ncol=2,label = F)  # Secretory
p12 <- FeaturePlot(object, c('Dclk1','Espn','Sox9','Pou2f3'),ncol=2,label = F)  # Tuft1
p13 <- FeaturePlot(object, c('Dclk1','Trpm5','Sox9','Lgr5'),ncol=2,label = F)  # Tuft2
p14 <- FeaturePlot(object, c('Sftpc','Scgb1a1','Sox9','Lgr5'),ncol=2,label = F)  # BASC
midrow2 <- plot_grid(p6,p7,p8,ncol=3)
midrow3 <- plot_grid(p9,p10,p11,ncol=3)
lastrow <- plot_grid(p12,p13,p14,ncol=3)

# Cowplot all panels
plots <- list(toprow,midrow1,midrow2,midrow3,lastrow)
png(paste(dir_name,tag,'_integrated2_allplots.png',sep=""), width = 2000,height=2750)
plot_grid(plotlist=plots,ncol=1,rel_heights = c(2,1,2,2,2))
dev.off()


## Evaluate clustering within integrated object
markers <- FindAllMarkers(object, only.pos = T)
markers$ratio <- markers$pct.1/markers$pct.2
markers <- markers[order(-markers$ratio),]

save(markers,file=paste(dir_name,'.seurat.objs.int.marks.',Sys.Date(),".Robj",sep=""))

# Plot top markers
p1 <- UMAPPlot(object,group.by = 'seurat_clusters',label = T) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoAxes()
p2 <- UMAPPlot(object,group.by = 'Sample') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
toprow <- plot_grid(p1,p2,NULL,NULL,ncol=4)

clusters <- levels(object$seurat_clusters)
#clusters <- clusters[1:10]
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
    q <- FeaturePlot(object,c(genes[1:4]),ncol=2,label = F)
    myplots[[i]] <- plot_grid(p,q,ncol=2)
}
lastrow <- plot_grid(plotlist=myplots,ncol=2,labels=paste('Cluster',clusters),label_size = 40)

# Cowplot all panels
plots <- list(toprow,lastrow)
png(paste(dir_name,tag,'_integrated2_clusters_1.png',sep=""), width = 2000,height=2750)
plot_grid(plotlist=plots,ncol=1,rel_heights = c(1,5))  # 2.5,3.5
dev.off()

# Extra plots - QC & lineage
png(paste(dir_name,tag,'_integrated3_fp_qc.png',sep=""), width = 800, height = 1200)
FeaturePlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','perc.spliced','Mki67','Top2a'),ncol=2,label = T,repel=T)
dev.off()
png(paste(dir_name,tag,'_integrated3_vln_qc.png',sep=""), width = 800, height = 1200)
VlnPlot(object,c('percent.mt','nCount_RNA','nFeature_RNA','perc.spliced','Mki67','Top2a'),group.by='seurat_clusters',ncol=2,pt.size=0.1)
dev.off()

png(paste(dir_name,tag,'_integrated3_fp_lin.png',sep=""), width = 800, height = 800)
FeaturePlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'),ncol=2,label = T,repel=T)
dev.off()
png(paste(dir_name,tag,'_integrated3_vln_lin.png',sep=""), width = 800, height = 800)
VlnPlot(object,c('Cdh5','Col1a1','Epcam','Ptprc'),group.by='seurat_clusters',ncol=2,pt.size=0.1)
dev.off()

save(object,file = paste("epi.seurat.objs.int.",Sys.Date(),".Robj",sep=""))




#### ENDOTHELIUM ####

### Initial exploration of endo objects
## Initial embeddings of each object
# Load all endothelium
load("./endo.seurat.objs.2022-07-27.Robj")

# Cowplot of every endothelial object UMAP
cols <- hue_pal()(4)
class_colors <- c('Epithelium'=cols[2],'Endothelium'=cols[3],'Mesenchyme'=cols[1],'Immune'=cols[4])

# Order objects
allobjects <- c("RLMVEC","Native","BCEC2","BEF1","BEF2","BEF3","BEF12","BEF14","BEF15",
    "BEFM1","BEFM2","BEFM4","BEFM5","BEFM6")
endo <- endo[allobjects]

# Plot without re-embedding
myplots <- vector('list',length(endo))
for (i in 1:length(endo)){
    object <- endo[[i]]
    sample.name <- names(endo[i])
    message(sample.name)
    p <- UMAPPlot(object,group.by = 'class',cols = class_colors) +
            labs(title = paste(sample.name,'- Endothelium'), subtitle=paste('nCells =',ncol(object))) +
            theme(plot.title = element_text(size = 24,hjust = 0))
    myplots[[i]] <- p
}

png(paste("endo_all_umaps2.png"), width = 2000,height=2000)
print(plot_grid(plotlist=myplots,ncol=4))
dev.off()

# Plot with re-embedding
pcs <- 10
myplots <- vector('list',length(endo))
new.endo <- list() ; final.names <- c()
for (i in 1:length(endo)) {
    object <- endo[[i]]
    sample.name <- names(endo[i])
    message(sample.name)
    object <- NormalizeData(object)
    object <- FindVariableFeatures(object)
    object <- ScaleData(object)
    object <- RunPCA(object, npcs = 50, verbose = F)
    object <- FindNeighbors(object, dims = 1:pcs)
    object <- FindClusters(object, resolution = 0.5)
    object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
    p <- UMAPPlot(object,label = T) + labs(title = paste(sample.name,'- Endothelium'),subtitle = paste("PC's =",pcs,"\nRes = 0.5",'\nCells =',ncol(object))) + NoLegend() +
        theme(plot.title = element_text(size = 24))
    myplots[[i]] <- p
    new.endo[[i]] <- object
    final.names[[i]] <- sample.name
}
names(new.endo) <- final.names

png(paste("endo_all_umaps_22.png"), width = 2000,height=2000)
print(plot_grid(plotlist=myplots,ncol=4))
dev.off()


## Explore top native markers via FeaturePlots
natendo <- endo$Native
Idents(natendo) <- natendo$CellType_Final

markers <- FindAllMarkers(natendo, only.pos = T)
markers$ratio <- markers$pct.1/markers$pct.2
markers <- markers[order(-markers$ratio),]

## Endothelial markers of interest
# Endothelium: 'Lyve1','Prx','Plvap','Apln'
# gCap1: 'Plvap','Aplnr','Cdh5','Lyve1'
# gCap2: 'Kit','Wnt5b','Nostrin','Vegfa'
# aCap: 'Apln','Prx','Ca4','Ednrb'
# Arterial: 'Efnb2','Vwf','Lyve1','Gja5'
# Venous: 'Vcam1','Vwf','Lyve1','Selp'
# Lymphatic: 'Prox1','Mmrn1'

# Makes a sheet of individual gene fps across all objects
endo2 <- new.endo
genes <- c('Prox1','Mmrn1')
myplots <- vector('list',length(endo2))
for (i in 1:length(genes)) {
    gene <- genes[i]
    message(gene)
    for (j in 1:length(endo2)) {
        object <- endo2[[j]]
        sample.name <- names(endo2[j])
        message(sample.name)
        p <- FeaturePlot(object,features=gene,label = F) &
            labs(title = paste(sample.name,'-',gene))
        myplots[[j]] <- p
    }
    png(paste("endo_all_fp_",gene,".png",sep=''),width = 2000,height=1000)
    print(plot_grid(plotlist=myplots,ncol=4))
    dev.off()
}

# Makes a sheet of a set of gene fps across all objects
type <- c('Lymphatic')
genes <- c('Prox1','Mmrn1')
myplots1 <- vector('list',length(genes))
myplots2 <- vector('list',length(endo2))
for (i in 1:length(endo2)) {
    object <- endo2[[i]]
    sample.name <- names(endo2[i])
    message(sample.name)
    for (j in 1:length(genes)) {
        gene <- genes[j]
        message(gene)
        p <- FeaturePlot(object,features=gene,label = F) +
            labs(title = paste(sample.name,'-',gene))
            #theme(plot.title = element_text(size = 20,hjust = 0,vjust = 0))
        myplots1[[j]] <- p
    }
    myplots2[[i]] <- plot_grid(plotlist=myplots1,ncol=2)
}

plot <- plot_grid(plotlist=myplots2,ncol=4)
title <- ggdraw() + draw_label(type,fontface = 'bold',size = 24,x = 0,hjust = 0)
png(paste("endo_all_fp_",type,".png",sep=''),width = 2000,height=1000)
print(plot_grid(title,plot,ncol=1,rel_heights = c(0.1, 10)))
dev.off()


## Explore top native markers via Heatmaps
load("/nobackup1/greaneya/DGEs_30K_Rank_Cutoff/Native/native.marks.2022-08-02.Robj")

# AverageExpression Heatmaps
# Select top genes per native endo type
table(new.endo$Native$CellType_Final)
native.endo.names <- c('gCap','aCap','Arterial','Venous','Lymphatic')
marks.short <- markers %>% group_by(cluster) %>% top_n(50, avg_log2FC) %>%
    select(cluster, gene, avg_log2FC, p_val,ratio) %>%
    arrange(cluster, desc(avg_log2FC))

# Loop to generate heatmaps of each native cell type with engineered (wo negative native ctrls)
endo4 <- new.endo[c("BCEC2","BEF1","BEF12","BEF14","BEF15","BEF2","BEF3","BEFM1","BEFM2","BEFM4","BEFM5",
    "BEFM6","RLMVEC")]  # minus Native
myplots <- vector('list', length(native.endo.names))
unityNormalize <- function(x){(x-min(x))/(max(x)-min(x))}
for (i in 1:length(native.endo.names)) {
    hm.short <- list()
    type <- native.endo.names[i]
    message(type)
    identities <- c(type,names(endo4))
    ctrl <- subset(new.endo$Native,subset = CellType_Final == c(type))
    hm.object <- c(ctrl,endo4)
    names(hm.object) <- identities
    topgenes <- marks.short[which(marks.short$cluster==type), ]
    genes.use <- topgenes$gene
    for (j in 1:length(hm.object)) {
        object <- hm.object[[j]]
        genes.int <- intersect(topgenes$gene,rownames(object))
        genes.use <- intersect(genes.use,genes.int)
    }
    genes.use <- head(genes.use,50)
    for (k in 1:length(hm.object)) {
        object <- hm.object[[k]]
        Idents(object) <- object$Sample
        hm.short[[k]] <- subset(object,downsample = 120)
    }
    hm.merge <- merge(x=hm.short[[1]],y=hm.short[2:length(hm.short)])
    hm.merge <- ScaleData(hm.merge,features = genes.use)
    hm.avg <- AverageExpression(hm.merge,assays='RNA',slot='scale.data',features=genes.use)
    hm.avg <- t(scale(t(hm.avg$RNA)))
    hm.avg <- t(apply(as.matrix(hm.avg), 1, unityNormalize))
    colnames(hm.avg)[which(colnames(hm.avg)=='Native')] <- type
    hm.avg <- hm.avg[,c(identities)]
    message(dim(hm.avg))

    colors.inferno <- colorRamp2(breaks = c(seq(min(hm.avg),max(hm.avg),length.out=60)), inferno(n=60), space = "RGB")
    myplots[[i]] <- local({
        i <- i
        p1 <- Heatmap(as.matrix(hm.avg),
                                col= colors.inferno,
                                row_names_gp = gpar(fontsize = 18),
                                show_column_dend = T,
                                clustering_method_columns = 'complete',
                                cluster_rows=T,
                                cluster_columns=T,
                                cluster_column_slices=F,
                                show_column_names=T,
                                column_names_side = "top",
                                column_dend_side = "bottom",
                                column_names_gp = gpar(fontsize = 30),
                                show_row_names = T,
                                show_heatmap_legend = F)
    })
    png(paste('nativeendo_heatmap_',type,'_2.png',sep=""), width = 800, height = 1000)  #400,1000
    draw(myplots[[i]],padding = unit(c(2, 2, 18, 2), "mm"))  #12
    dev.off()
}



### Generating integrated object of all Endothelium
# Normal Seurat CCA Integration of all engineered objects, without native
objects <- endo[c("RLMVEC","BCEC2","BEF1","BEF2","BEF3","BEF12","BEF14","BEF15","BEFM1",
    "BEFM2","BEFM4","BEFM5","BEFM6")]
dir_name <- c('endo')
objects <- lapply(X = objects, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x,nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = objects, nfeatures = 2000)
object.anchors <- FindIntegrationAnchors(object.list=objects, anchor.features=features,
        reduction='cca',dims=1:30)
object <- IntegrateData(anchorset = object.anchors, dims=1:30, k.weight = 54)

# Perform integration
DefaultAssay(object) <- "integrated"
object <- ScaleData(object)
object <- RunPCA(object, npcs = 100, verbose = F)
png(paste(dir_name,'_integrated2_elbow.png',sep=""), width = 800, height = 500)
ElbowPlot(object,ndims = 50)+labs(title = dir_name)
dev.off()
pdf(paste(dir_name,'_integrated2_pcheatmap1.pdf',sep=""), width = 12, height = 16)
DimHeatmap(object, dims = 1:15, cells = 200, balanced = T)
dev.off()
pdf(paste(dir_name,'_integrated2_pcheatmap2.pdf',sep=""), width = 12, height = 16)
DimHeatmap(object, dims = 16:30, cells = 200, balanced = T)
dev.off()
pdf(paste(dir_name,'_integrated2_pcheatmap3.pdf',sep=""), width = 12, height = 16)
DimHeatmap(object, dims = 31:45, cells = 200, balanced = T)
dev.off()

pcs <- 10  # iterate
res <- 0.4  # iterate
tag <- paste('pcs',pcs,'res',res,sep='.')
DefaultAssay(object) <- "integrated"
object <- FindNeighbors(object, dims = 1:pcs)
object <- FindClusters(object, resolution = res)
object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
p1 <- UMAPPlot(object,group.by = 'seurat_clusters',label = T) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoAxes()
p2 <- UMAPPlot(object,group.by = 'Sample') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
Idents(object) <- object$Sample
toprow <- plot_grid(p1,p2,NULL,ncol=3)

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
p7 <- FeaturePlot(object, c('Plvap','Aplnr','Cdh5','Lyve1'),ncol=2,label = F)  # gCap1
p8 <- FeaturePlot(object, c('Kit','Wnt5b','Nostrin','Vegfa'),ncol=2,label = F)  # gCap2
p9 <- FeaturePlot(object, c('Apln','Prx','Ca4','Ednrb'),ncol=2,label = F)  # aCap
p10 <- FeaturePlot(object, c('Efnb2','Vwf','Lyve1','Gja5'),ncol=2,label = F)  # Arterial
p11 <- FeaturePlot(object, c('Vcam1','Vwf','Lyve1','Selp'),ncol=2,label = F)  # Venous
p12 <- FeaturePlot(object, c('Prox1','Mmrn1'),ncol=2,label = F)  # Lymphatic
midrow2 <- plot_grid(p6,p7,p8,ncol=3)
midrow3 <- plot_grid(p9,p10,p11,ncol=3)
lastrow <- plot_grid(p12,NULL,NULL,ncol=3)

# Cowplot all panels
plots <- list(toprow,midrow1,midrow2,midrow3,lastrow)
png(paste(dir_name,tag,'_integrated2_allplots.png',sep=""), width = 2000,height=2750)
plot_grid(plotlist=plots,ncol=1,rel_heights = c(2,1,2,2,1))
dev.off()


## Evaluate clustering within integrated object
markers <- FindAllMarkers(object, only.pos = T)
markers$ratio <- markers$pct.1/markers$pct.2
markers <- markers[order(-markers$ratio),]

save(markers,file=paste(dir_name,'seurat.objs.int.marks.',Sys.Date(),".Robj",sep=""))

# Plot top markers
p1 <- UMAPPlot(object,group.by = 'seurat_clusters',label = T) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoAxes()
p2 <- UMAPPlot(object,group.by = 'Sample') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
toprow <- plot_grid(p1,p2,NULL,NULL,ncol=4)

clusters <- levels(object$seurat_clusters)
#clusters <- clusters[1:10]
goi <- c(); marks <- c(); genes <- c()
myplots <- vector('list',length(clusters))
for (i in 1:length(clusters)) {
    cluster <- clusters[[i]]
    message(cluster)

    goi <- markers[which(markers$cluster == cluster),]
    goi <- goi[1:20, ]  # 15
    row.names(goi) <- goi$gene
    marks <- round(goi[, 2:5], digits = 4)
    marks$ratio <- round(goi$ratio, digits = 4)
    marks <- as.data.frame(marks)
    p <- tableGrob(marks,theme = ttheme_minimal(base_size = 16))

    #marks <- marks[order(marks$avg_log2FC),]
    genes <- row.names(marks)
    q <- FeaturePlot(object,c(genes[1:4]),ncol=2,label = F)
    myplots[[i]] <- plot_grid(p,q,ncol=2)
}
lastrow <- plot_grid(plotlist=myplots,ncol=2,labels=paste('Cluster',clusters),label_size = 40)

# Cowplot all panels
plots <- list(toprow,lastrow)
png(paste(dir_name,tag,'_integrated2_clusters_3.png',sep=""), width = 2000,height=2750)
plot_grid(plotlist=plots,ncol=1,rel_heights = c(1,5))  # 2.5,3.5
dev.off()

# Extra plots - QC & lineage
png(paste(dir_name,tag,'_integrated2_fp_qc.png',sep=""), width = 800, height = 1200)
FeaturePlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','perc.spliced','Mki67','Top2a'),ncol=2,label = T,repel=T)
dev.off()
png(paste(dir_name,tag,'_integrated2_vln_qc.png',sep=""), width = 800, height = 1200)
VlnPlot(object,c('percent.mt','nCount_RNA','nFeature_RNA','perc.spliced','Mki67','Top2a'),group.by='seurat_clusters',ncol=2,pt.size=0.1)
dev.off()

png(paste(dir_name,tag,'_integrated2_fp_lin.png',sep=""), width = 800, height = 800)
FeaturePlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'),ncol=2,label = T,repel=T)
dev.off()
png(paste(dir_name,tag,'_integrated2_vln_lin.png',sep=""), width = 800, height = 800)
VlnPlot(object,c('Cdh5','Col1a1','Epcam','Ptprc'),group.by='seurat_clusters',ncol=2,pt.size=0.1)
dev.off()

save(object,file = paste("endo.seurat.objs.int.",Sys.Date(),".Robj",sep=""))




#### MESENCHYME ####

### Initial exploration of mes objects
## Initial embeddings of each object
# Load all mesenchyme
load("./mes.seurat.objs.2022-07-27.Robj")

# Cowplot of every mesenchyme object UMAP
cols <- hue_pal()(4)
class_colors <- c('Epithelium'=cols[2],'Endothelium'=cols[3],'Mesenchyme'=cols[1],'Immune'=cols[4])

# Order objects
allobjects <- c("FB13","FB14","Native","BEF1","BEF2","BEF3","BEF12","BEF14","BEF15",
    "BEFM1","BEFM2","BEFM4","BEFM5","BEFM6")
mes <- mes[allobjects]

# Plot without re-embedding
myplots <- vector('list',length(mes))
for (i in 1:length(mes)){
    object <- mes[[i]]
    sample.name <- names(mes[i])
    message(sample.name)
    p <- UMAPPlot(object,group.by = 'class',cols = class_colors) +
            labs(title = paste(sample.name,'- Mesenchyme'), subtitle=paste('nCells =',ncol(object))) +
            theme(plot.title = element_text(size = 24,hjust = 0))
    myplots[[i]] <- p
}

png(paste("mes_all_umaps2.png"), width = 2500,height=2000)
print(plot_grid(plotlist=myplots,ncol=5))
dev.off()

# Plot with re-embedding
pcs <- 10
myplots <- vector('list',length(mes))
new.mes <- list() ; final.names <- c()
for (i in 1:length(mes)) {
    object <- mes[[i]]
    sample.name <- names(mes[i])
    message(sample.name)
    object <- NormalizeData(object)
    object <- FindVariableFeatures(object)
    object <- ScaleData(object)
    object <- RunPCA(object, npcs = 50, verbose = F)
    object <- FindNeighbors(object, dims = 1:pcs)
    object <- FindClusters(object, resolution = 0.5)
    object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
    p <- UMAPPlot(object,label = T) + labs(title = paste(sample.name,'- Mesenchyme'),subtitle = paste("PC's =",pcs,"\nRes = 0.5",'\nCells =',ncol(object))) + NoLegend() +
        theme(plot.title = element_text(size = 24))
    myplots[[i]] <- p
    new.mes[[i]] <- object
    final.names[[i]] <- sample.name
}
names(new.mes) <- final.names

png(paste("mes_all_umaps_22.png"), width = 2500,height=2000)
print(plot_grid(plotlist=myplots,ncol=5))
dev.off()


## Explore top native markers via FeaturePlots
natmes <- mes$Native
Idents(natmes) <- natmes$CellType_Final

markers <- FindAllMarkers(natmes, only.pos = T)
markers$ratio <- markers$pct.1/markers$pct.2
markers <- markers[order(-markers$ratio),]

## Mesenchymal markers of interest
# Mesenchyme: 'Dcn','Acta2','Itga8','Wnt11'
# Pi16+ Fibs: 'Col14a1','Dcn','Pi16','Col13a1'
# Mesothelium: 'Msln','Rspo1','Dclk1','Krt8'
# Itga8+ Fibs: 'Col13a1','Itga8','Hsd11b1','Col14a1'
# Myofibroblasts: 'Wnt11','Dcn','Acta2','Myh11'
# Notum+ Fibs: 'Lgr5','Notum','Col13a1','Col14a1'
# Smooth Muscle: 'Acta2','Myh11','Acta1','Itga8'
# Pericytes: 'Gucy1a1','Gucy1a2','Gucy1b1','Acta2'

# Makes a sheet of individual gene fps across all objects
mes2 <- new.mes
genes <- c('Acta2','Myh11','Acta1','Itga8')
myplots <- vector('list',length(mes2))
for (i in 1:length(genes)) {
    gene <- genes[i]
    message(gene)
    for (j in 1:length(mes2)) {
        object <- mes2[[j]]
        sample.name <- names(mes2[j])
        message(sample.name)
        p <- FeaturePlot(object,features=gene,label = F) &
            labs(title = paste(sample.name,'-',gene))
        myplots[[j]] <- p
    }
    png(paste("mes_all_fp_",gene,".png",sep=''),width = 2000,height=2000)
    print(plot_grid(plotlist=myplots,ncol=4))
    dev.off()
}

# Makes a sheet of a set of gene fps across all objects
type <- c('Smooth Muscle')
genes <- c('Acta2','Myh11','Acta1','Itga8')
myplots1 <- vector('list',length(genes))
myplots2 <- vector('list',length(mes2))
for (i in 1:length(mes2)) {
    object <- mes2[[i]]
    sample.name <- names(mes2[i])
    message(sample.name)
    for (j in 1:length(genes)) {
        gene <- genes[j]
        message(gene)
        p <- FeaturePlot(object,features=gene,label = F) +
            labs(title = paste(sample.name,'-',gene))
        myplots1[[j]] <- p
    }
    myplots2[[i]] <- plot_grid(plotlist=myplots1,ncol=2)
}

plot <- plot_grid(plotlist=myplots2,ncol=4)
title <- ggdraw() + draw_label(type,fontface = 'bold',size = 24,x = 0,hjust = 0)
png(paste("mes_all_fp_",type,".png",sep=''),width = 2000,height=2000)
print(plot_grid(title,plot,ncol=1,rel_heights = c(0.1, 10)))
dev.off()


## Explore top native markers via Heatmaps
load("/nobackup1/greaneya/DGEs_30K_Rank_Cutoff/Native/native.marks.2022-08-02.Robj")

# AverageExpression Heatmaps
# Select top genes per native mes type
table(new.mes$Native$CellType_Final)
native.mes.names <- c('Pi16+_Fibroblasts','Mesothelium','Itga8+_Fibroblasts','Myofibroblasts','Notum+_Fibroblasts','Smooth_Muscle','Pericytes')
marks.short <- markers %>% group_by(cluster) %>% top_n(50, avg_log2FC) %>%
    select(cluster, gene, avg_log2FC, p_val,ratio) %>%
    arrange(cluster, desc(avg_log2FC))

# Loop to generate heatmaps of each native cell type with engineered (wo negative native ctrls)
mes4 <- new.mes[c("BEF1","BEF2","BEF3","BEF12","BEF14","BEF15","BEFM1","BEFM2","BEFM4","BEFM5",
    "BEFM6","FB13","FB14")]  # minus Native
myplots <- vector('list', length(native.mes.names))
unityNormalize <- function(x){(x-min(x))/(max(x)-min(x))}
for (i in 1:length(native.mes.names)) {
    hm.short <- list()
    type <- native.mes.names[i]
    message(type)
    identities <- c(type,names(mes4))
    ctrl <- subset(new.mes$Native,subset = CellType_Final == c(type))
    hm.object <- c(ctrl,mes4)
    names(hm.object) <- identities
    topgenes <- marks.short[which(marks.short$cluster==type), ]
    genes.use <- topgenes$gene
    for (j in 1:length(hm.object)) {
        object <- hm.object[[j]]
        genes.int <- intersect(topgenes$gene,rownames(object))
        genes.use <- intersect(genes.use,genes.int)
    }
    genes.use <- head(genes.use,50)
    for (k in 1:length(hm.object)) {
        object <- hm.object[[k]]
        Idents(object) <- object$Sample
        hm.short[[k]] <- subset(object,downsample = 120)
    }
    hm.merge <- merge(x=hm.short[[1]],y=hm.short[2:length(hm.short)])
    hm.merge <- ScaleData(hm.merge,features = genes.use)
    hm.avg <- AverageExpression(hm.merge,assays='RNA',slot='scale.data',features=genes.use)
    hm.avg <- t(scale(t(hm.avg$RNA)))
    hm.avg <- t(apply(as.matrix(hm.avg), 1, unityNormalize))
    colnames(hm.avg)[which(colnames(hm.avg)=='Native')] <- type
    hm.avg <- hm.avg[,c(identities)]
    message(dim(hm.avg))

    colors.inferno <- colorRamp2(breaks = c(seq(min(hm.avg),max(hm.avg),length.out=60)), inferno(n=60), space = "RGB")
    myplots[[i]] <- local({
        i <- i
        p1 <- Heatmap(as.matrix(hm.avg),
                                col= colors.inferno,
                                row_names_gp = gpar(fontsize = 18),
                                show_column_dend = T,
                                clustering_method_columns = 'complete',
                                cluster_rows=T,
                                cluster_columns=T,
                                cluster_column_slices=F,
                                show_column_names=T,
                                column_names_side = "top",
                                column_dend_side = "bottom",
                                column_names_gp = gpar(fontsize = 30),
                                show_row_names = T,
                                show_heatmap_legend = F)
    })
    png(paste('nativemes_heatmap_',type,'_2.png',sep=""), width = 800, height = 1000)  #400,1000
    draw(myplots[[i]],padding = unit(c(2, 2, 18, 2), "mm"))  #12
    dev.off()
}



### Generating integrated object of all Mesenchyme
# Normal Seurat CCA Integration of all engineered objects, without native
objects <- mes[c("BEF1","BEF2","BEF3","BEF12","BEF14","BEF15","BEFM1","BEFM2","BEFM4","BEFM5","BEFM6")]
dir_name <- c('mes')
objects <- lapply(X = objects, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x,nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = objects, nfeatures = 2000)
object.anchors <- FindIntegrationAnchors(object.list=objects, anchor.features=features,
        reduction='cca',dims=1:30)
object <- IntegrateData(anchorset = object.anchors, dims=1:30, k.weight = 88)

# Perform integration
DefaultAssay(object) <- "integrated"
object <- ScaleData(object)
object <- RunPCA(object, npcs = 100, verbose = F)
png(paste(dir_name,'_integrated3_elbow.png',sep=""), width = 800, height = 500)
ElbowPlot(object,ndims = 50)+labs(title = dir_name)
dev.off()
pdf(paste(dir_name,'_integrated3_pcheatmap1.pdf',sep=""), width = 12, height = 16)
DimHeatmap(object, dims = 1:15, cells = 200, balanced = T)
dev.off()
pdf(paste(dir_name,'_integrated3_pcheatmap2.pdf',sep=""), width = 12, height = 16)
DimHeatmap(object, dims = 16:30, cells = 200, balanced = T)
dev.off()
pdf(paste(dir_name,'_integrated3_pcheatmap3.pdf',sep=""), width = 12, height = 16)
DimHeatmap(object, dims = 31:45, cells = 200, balanced = T)
dev.off()

pcs <- 20  # iterate
res <- 0.4  # iterate
tag <- paste('pcs',pcs,'res',res,sep='.')
DefaultAssay(object) <- "integrated"
object <- FindNeighbors(object, dims = 1:pcs)
object <- FindClusters(object, resolution = res)
object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
p1 <- UMAPPlot(object,group.by = 'seurat_clusters',label = T) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoAxes()
p2 <- UMAPPlot(object,group.by = 'Sample') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
Idents(object) <- object$Sample
toprow <- plot_grid(p1,p2,NULL,ncol=3)

data <- prop.table(table(object$seurat_clusters,object$Sample),1)*100
data <- as.data.frame(data)
p5 <- ggplot(data, aes(x = Var1, y = Freq, fill = Var2)) + NoLegend() +
  geom_col(position = "dodge", width = 1) + scale_fill_manual(values = hue_pal()(length(unique(object$Sample)))) +
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
png(paste(dir_name,tag,'_integrated3_allplots.png',sep=""), width = 2000,height=2750)
plot_grid(plotlist=plots,ncol=1,rel_heights = c(2,1,2,2,2))
dev.off()


## Evaluate clustering within integrated object
markers <- FindAllMarkers(object, only.pos = T)
markers$ratio <- markers$pct.1/markers$pct.2
markers <- markers[order(-markers$ratio),]

save(markers,file=paste(dir_name,'seurat.objs.int.marks.',Sys.Date(),".Robj",sep=""))

# Plot top markers
p1 <- UMAPPlot(object,group.by = 'seurat_clusters',label = T) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoAxes()
p2 <- UMAPPlot(object,group.by = 'Sample') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
toprow <- plot_grid(p1,p2,NULL,NULL,ncol=4)

clusters <- levels(object$seurat_clusters)
#clusters <- clusters[1:10]
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
    q <- FeaturePlot(object,c(genes[1:4]),ncol=2,label = F)
    myplots[[i]] <- plot_grid(p,q,ncol=2)
}
lastrow <- plot_grid(plotlist=myplots,ncol=2,labels=paste('Cluster',clusters),label_size = 40)

# Cowplot all panels
plots <- list(toprow,lastrow)
png(paste(dir_name,tag,'_integrated3_clusters_1.png',sep=""), width = 2000,height=2750)
plot_grid(plotlist=plots,ncol=1,rel_heights = c(1,5))  # 2.5,3.5
dev.off()

# Extra plots - QC & lineage
png(paste(dir_name,tag,'_integrated3_fp_qc.png',sep=""), width = 800, height = 1200)
FeaturePlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','perc.spliced','Mki67','Top2a'),ncol=2,label = T,repel=T)
dev.off()
png(paste(dir_name,tag,'_integrated3_vln_qc.png',sep=""), width = 800, height = 1200)
VlnPlot(object,c('percent.mt','nCount_RNA','nFeature_RNA','perc.spliced','Mki67','Top2a'),group.by='seurat_clusters',ncol=2,pt.size=0.1)
dev.off()

png(paste(dir_name,tag,'_integrated3_fp_lin.png',sep=""), width = 800, height = 800)
FeaturePlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'),ncol=2,label = T,repel=T)
dev.off()
png(paste(dir_name,tag,'_integrated3_vln_lin.png',sep=""), width = 800, height = 800)
VlnPlot(object,c('Cdh5','Col1a1','Epcam','Ptprc'),group.by='seurat_clusters',ncol=2,pt.size=0.1)
dev.off()

save(object,file = paste("mes.seurat.objs.int.",Sys.Date(),".Robj",sep=""))




#### IMMUNE ####

### Initial exploration of imm objects
## Initial embeddings of each object
# Load all immune
load("./imm.seurat.objs.2022-07-27.Robj")

# Cowplot of every immune object UMAP
cols <- hue_pal()(4)
class_colors <- c('Epithelium'=cols[2],'Endothelium'=cols[3],'Mesenchyme'=cols[1],'Immune'=cols[4])

# Order objects
allobjects <- c("MacAlv","FB13","FB14","Native","BEF12","BEFM1","BEFM2","BEFM4","BEFM5","BEFM6")
imm <- imm[allobjects]

# Plot without re-embedding
myplots <- vector('list',length(imm))
for (i in 1:length(imm)){
    object <- imm[[i]]
    sample.name <- names(imm[i])
    message(sample.name)
    p <- UMAPPlot(object,group.by = 'class',cols = class_colors) +
            labs(title = paste(sample.name,'- Immune'), subtitle=paste('nCells =',ncol(object))) +
            theme(plot.title = element_text(size = 24,hjust = 0))
    myplots[[i]] <- p
}

png(paste("imm_all_umaps2.png"), width = 2000,height=2000)
print(plot_grid(plotlist=myplots,ncol=4))
dev.off()

# Plot with re-embedding (skip bc too small: BEF12)
imm2 <- imm[c("MacAlv","FB13","FB14","Native","BEFM1","BEFM2","BEFM4","BEFM5","BEFM6")]
pcs <- 10
myplots <- vector('list',length(imm2))
new.imm <- list() ; final.names <- c()
for (i in 1:length(imm2)) {
    object <- imm2[[i]]
    sample.name <- names(imm2[i])
    message(sample.name)
    object <- NormalizeData(object)
    object <- FindVariableFeatures(object)
    object <- ScaleData(object)
    object <- RunPCA(object, npcs = 50, verbose = F)
    object <- FindNeighbors(object, dims = 1:pcs)
    object <- FindClusters(object, resolution = 0.5)
    object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
    p <- UMAPPlot(object,label = T) + labs(title = paste(sample.name,'- Immune'),subtitle = paste("PC's =",pcs,"\nRes = 0.5",'\nCells =',ncol(object))) + NoLegend() +
        theme(plot.title = element_text(size = 24))
    myplots[[i]] <- p
    new.imm[[i]] <- object
    final.names[[i]] <- sample.name
}
names(new.imm) <- final.names

png(paste("imm_all_umaps_22.png"), width = 2000,height=2000)
print(plot_grid(plotlist=myplots,ncol=4))
dev.off()


## Explore top native markers via FeaturePlots
natimm <- imm$Native
Idents(natimm) <- natimm$CellType_Final

markers <- FindAllMarkers(natimm, only.pos = T)
markers$ratio <- markers$pct.1/markers$pct.2
markers <- markers[order(-markers$ratio),]

## Immune markers of interest
# Alveolar_Mac: 'Mrc1','Arg1','Naaa','Top2a'
# Interstitial_Mac: 'C1qc','C1qb','Cd163','Arg1'
# Mono_NC_C: 'Itgax','Cd79b','Sell','Rnase2' ['F5','Vcan']
# Mono_Act: 'Slamf7','Jaml','Cd200','Cd226'
# B: 'Ms4a1','Cd19','Cd79b','Ebf1'
# T: 'Cd3d','Cd3e','Cd8a','Cd40lg'
# Neutrophils: 'S100a8','S100a9','Mgam','Cxcr2'
# NK: 'Nkg7','Gzma','Klrd1','Ncr1'
# ILC: 'Il2ra','Gata3','Il17rb','Il1rl1'
# DC: 'Cd209','Siglech','Kmo','Irf8'
# Mast_Plasma: 'Tph1','Tpsb2','Jchain'
# Eosinophil: 'Ccl24','Adgre1','Asz1','Selp'
# Basophil: 'Il3ra','Cd200r1','Dnase2b','Cdh1'

# Makes a sheet of individual gene fps across all objects
imm3 <- c(new.imm,imm[c('BEF12')])
imm3 <- imm3[order(names(imm))]
genes <- c('Slamf7','Jaml','Cd200','Cd226')
myplots <- vector('list',length(imm3))
for (i in 1:length(genes)) {
    gene <- genes[i]
    message(gene)
    for (j in 1:length(imm3)) {
        object <- imm3[[j]]
        sample.name <- names(imm3[j])
        message(sample.name)
        p <- FeaturePlot(object,features=gene,label = F) &
            labs(title = paste(sample.name,'-',gene))
        myplots[[j]] <- p
    }
    png(paste("imm_all_fp_",gene,".png",sep=''),width = 2000,height=2000)
    print(plot_grid(plotlist=myplots,ncol=4))
    dev.off()
}

# Makes a sheet of a set of gene fps across all objects
type <- c('Mono_Act')
genes <- c('Slamf7','Jaml','Cd200','Cd226')
myplots1 <- vector('list',length(genes))
myplots2 <- vector('list',length(imm3))
for (i in 1:length(imm3)) {
    object <- imm3[[i]]
    sample.name <- names(imm3[i])
    message(sample.name)
    for (j in 1:length(genes)) {
        gene <- genes[j]
        message(gene)
        p <- FeaturePlot(object,features=gene,label = F) +
            labs(title = paste(sample.name,'-',gene))
        myplots1[[j]] <- p
    }
    myplots2[[i]] <- plot_grid(plotlist=myplots1,ncol=2)
}

plot <- plot_grid(plotlist=myplots2,ncol=4)
title <- ggdraw() + draw_label(type,fontface = 'bold',size = 24,x = 0,hjust = 0)
png(paste("imm_all_fp_",type,".png",sep=''),width = 2000,height=2000)
print(plot_grid(title,plot,ncol=1,rel_heights = c(0.1, 10)))
dev.off()


## Explore top native markers via Heatmaps
load("/nobackup1/greaneya/DGEs_30K_Rank_Cutoff/Native/native.marks.2022-08-02.Robj")

# AverageExpression Heatmaps
# Select top genes per native imm type
table(new.imm$Native$CellType_Final)
native.imm.names <- c('B','Mac Alv','T','Mono 1','NK','Mono 2','Neutro','Mac Inter','Mono 3',
    'Cell cycle','ILC','pDC','Killer T','Mast','Ccl24+')
marks.short <- markers %>% group_by(cluster) %>% top_n(50, avg_log2FC) %>%
    select(cluster, gene, avg_log2FC, p_val,ratio) %>%
    arrange(cluster, desc(avg_log2FC))

# Loop to generate heatmaps of each native cell type with engineered (wo negative native ctrls)
imm4 <- new.imm[c("BEFM1","BEFM2","BEFM4","BEFM5","BEFM6","FB13","FB14","MacAlv",
    "TXP3_L","TXP3_R","TXP4_L","TXP4_R")]  # minus Native
myplots <- vector('list', length(native.imm.names))
unityNormalize <- function(x){(x-min(x))/(max(x)-min(x))}
for (i in 1:length(native.imm.names)) {
    hm.short <- list()
    type <- native.imm.names[i]
    message(type)
    identities <- c(type,names(imm4))
    ctrl <- subset(new.imm$Native,subset = CellType_Final == c(type))
    hm.object <- c(ctrl,imm4)
    names(hm.object) <- identities
    topgenes <- marks.short[which(marks.short$cluster==type), ]
    genes.use <- topgenes$gene
    for (j in 1:length(hm.object)) {
        object <- hm.object[[j]]
        genes.int <- intersect(topgenes$gene,rownames(object))
        genes.use <- intersect(genes.use,genes.int)
    }
    genes.use <- head(genes.use,50)
    for (k in 1:length(hm.object)) {
        object <- hm.object[[k]]
        Idents(object) <- object$Sample
        hm.short[[k]] <- subset(object,downsample = 120)
    }
    hm.merge <- merge(x=hm.short[[1]],y=hm.short[2:length(hm.short)])
    hm.merge <- ScaleData(hm.merge,features = genes.use)
    hm.avg <- AverageExpression(hm.merge,assays='RNA',slot='scale.data',features=genes.use)
    hm.avg <- t(scale(t(hm.avg$RNA)))
    hm.avg <- t(apply(as.matrix(hm.avg), 1, unityNormalize))
    colnames(hm.avg)[which(colnames(hm.avg)=='Native')] <- type
    hm.avg <- hm.avg[,c(identities)]
    message(dim(hm.avg))

    colors.inferno <- colorRamp2(breaks = c(seq(min(hm.avg),max(hm.avg),length.out=60)), inferno(n=60), space = "RGB")
    myplots[[i]] <- local({
        i <- i
        p1 <- Heatmap(as.matrix(hm.avg),
                                col= colors.inferno,
                                row_names_gp = gpar(fontsize = 18),
                                show_column_dend = T,
                                clustering_method_columns = 'complete',
                                cluster_rows=T,
                                cluster_columns=T,
                                cluster_column_slices=F,
                                show_column_names=T,
                                column_names_side = "top",
                                column_dend_side = "bottom",
                                column_names_gp = gpar(fontsize = 30),
                                show_row_names = T,
                                show_heatmap_legend = F)
    })
    png(paste('nativeimm_heatmap_',type,'_2.png',sep=""), width = 800, height = 1000)  #400,1000
    draw(myplots[[i]],padding = unit(c(2, 2, 18, 2), "mm"))  #12
    dev.off()
}



### Generating integrated object of all Immune
# Normal Seurat CCA Integration of all engineered objects, without native or tiny samples
objects <- imm[c("MacAlv","BEFM1","BEFM2","BEFM4","BEFM5","BEFM6")]
dir_name <- c('imm')
objects <- lapply(X = objects, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x,nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = objects, nfeatures = 2000)
object.anchors <- FindIntegrationAnchors(object.list=objects, anchor.features=features,
        reduction='cca',dims=1:30)
object <- IntegrateData(anchorset = object.anchors, dims=1:30, k.weight = 55)

# Perform integration
DefaultAssay(object) <- "integrated"
object <- ScaleData(object)
object <- RunPCA(object, npcs = 100, verbose = F)
png(paste(dir_name,'_integrated2_elbow.png',sep=""), width = 800, height = 500)
ElbowPlot(object,ndims = 50)+labs(title = dir_name)
dev.off()
pdf(paste(dir_name,'_integrated2_pcheatmap1.pdf',sep=""), width = 12, height = 16)
DimHeatmap(object, dims = 1:15, cells = 200, balanced = T)
dev.off()
pdf(paste(dir_name,'_integrated2_pcheatmap2.pdf',sep=""), width = 12, height = 16)
DimHeatmap(object, dims = 16:30, cells = 200, balanced = T)
dev.off()
pdf(paste(dir_name,'_integrated2_pcheatmap3.pdf',sep=""), width = 12, height = 16)
DimHeatmap(object, dims = 31:45, cells = 200, balanced = T)
dev.off()

pcs <- 30  # iterate
res <- 0.5  # iterate
tag <- paste('pcs',pcs,'res',res,sep='.')
DefaultAssay(object) <- "integrated"
object <- FindNeighbors(object, dims = 1:pcs)
object <- FindClusters(object, resolution = res)
object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
p1 <- UMAPPlot(object,group.by = 'seurat_clusters',label = T) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoAxes()
p2 <- UMAPPlot(object,group.by = 'Sample') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
Idents(object) <- object$Sample
toprow <- plot_grid(p1,p2,NULL,ncol=3)

data <- prop.table(table(object$seurat_clusters,object$Sample),1)*100
data <- as.data.frame(data)
p5 <- ggplot(data, aes(x = Var1, y = Freq, fill = Var2)) + NoLegend() +
  geom_col(position = "dodge", width = 1) + scale_fill_manual(values = hue_pal()(length(unique(object$Sample)))) +
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
p6 <- FeaturePlot(object, c('Mrc1','Tac1','Naaa','Top2a'),ncol=2,label = F)  # Alveolar_Mac
p7 <- FeaturePlot(object, c('C1qc','C1qb','Cd163','Stab1'),ncol=2,label = F)  # Interstitial_Mac
p8 <- FeaturePlot(object, c('Defb14','Capn3','Eno3','Slamf7'),ncol=2,label = F)  # Mono_NC_Act
p9 <- FeaturePlot(object, c('Rnase2','Hrh3','F5','Flt3'),ncol=2,label = F)  # Mono_C_Act
p10 <- FeaturePlot(object, c('Ms4a1','Cd19','Cd79b','Ebf1'),ncol=2,label = F)  # B
p11 <- FeaturePlot(object, c('Cd3d','Cd3e','Cd8b','Cd27'),ncol=2,label = F)  # T
p12 <- FeaturePlot(object, c('S100a8','S100a9','Nkg7','Gzma'),ncol=2,label = F)  # Neutro/NK
p13 <- FeaturePlot(object, c('Il17rb','Il2ra','Gata3','Ltb4r'),ncol=2,label = F)  # ILC
p14 <- FeaturePlot(object, c('Siglech','Smim5','Kmo','Flt3'),ncol=2,label = F)  # DC
p15 <- FeaturePlot(object, c('Tph1','Jchain','Asz1','Dnase2b'),ncol=2,label = F)  # Mast/Plasma/Eos/Baso
midrow2 <- plot_grid(p6,p7,p8,ncol=3)
midrow3 <- plot_grid(p9,p10,p11,ncol=3)
lastrow <- plot_grid(p12,p13,p14,p15,ncol=4)

# Cowplot all panels
plots <- list(toprow,midrow1,midrow2,midrow3,lastrow)
png(paste(dir_name,tag,'_integrated2_allplots.png',sep=""), width = 2000,height=2750)
plot_grid(plotlist=plots,ncol=1,rel_heights = c(2,1,2,2,2))
dev.off()

save(object,file = paste("imm.seurat.objs.int.",Sys.Date(),".Robj",sep=""))











