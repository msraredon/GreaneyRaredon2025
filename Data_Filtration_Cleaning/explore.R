#### EXPLORE
#### Code for publication unifying all 'explore' scripts
#### Replicable method for integrating cleaned objects by class




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

# explore_endo_1.R
#### Initial exploration of all Endothelium
## Re-clustering individual objects & Gene Expression Heatmaps


# explore_endo_2.R
#### Initial exploration of all Endothelium
## CCA Integration & Merge


# explore_endo_3.R
#### Initial exploration of all Endothelium
## Integration by Reference Mapping


# explore_endo_4.R
#### Exploration of all Endothelium
## Integration by Reference Mapping wo TXP_R's


# explore_endo_5.R
#### Exploration of Engineered Endothelium
## Re-do CCA Integration







#### MESENCHYME ####

# explore_mes_1.R
#### Initial exploration of all Mesenchyme
## Re-clustering individual objects & Gene Expression Heatmaps


# explore_mes_2.R
#### Initial exploration of all Mesenchyme
## CCA Integration & Merge


# explore_mes_3.R
#### Initial exploration of all Mesenchyme
## Integration by Reference Mapping


# explore_mes_4.R
#### Exploration of all Mesenchyme
## Integration by Reference Mapping


# explore_mes_5.R
#### Exploration of Engineered Mesenchyme
## Re-do CCA Integration




#### IMMUNE ####

# explore_imm_1.R
#### Initial exploration of all Immune
## Re-clustering individual objects & Gene Expression Heatmaps


# explore_imm_2.R
#### Initial exploration of all Immune
## Update Cell Labels & CCA Integration


# explore_imm_3.R
#### Initial exploration of all Immune
## Integration by Reference Mapping


# explore_imm_4.R
#### Exploration of all Immune
## Integration by Reference Mapping wo TXP_R's


# explore_imm_5.R
#### Exploration of Engineered Immune
## Re-do CCA Integration









