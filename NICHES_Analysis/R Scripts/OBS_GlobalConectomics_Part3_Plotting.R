# Require packages
require(Seurat)
require(NICHES)
require(dplyr)
require(cowplot)

# Set WD
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics")

# Load color palette (does not include native colors yet)
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Allie's Object Filtration/Final Objects/all.color_palettes_2.R")
color_pals

# load system to cell data
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/system.to.cell.BEFM.project.2023-08-03.Robj")
system.to.cell$Tissue <- factor(system.to.cell$Tissue,
                                levels = c('Pseudo.Start','BCL5','BCEC2',
                                           'BEF1','BEF2','BEF3',
                                           'BEF12','BEF14','BEF15',
                                           'BEFM1','BEFM2',
                                           'BEFM4','BEFM5','BEFM6',
                                           'fRat','mRat','P0-f','P0-m'))

# Set idents to class and subset out each class on its own
Idents(system.to.cell) <- system.to.cell$class
epi <- subset(system.to.cell,subset=class=='Epithelium')
end <- subset(system.to.cell,subset=class=='Endothelium')
mes <- subset(system.to.cell,subset=class=='Mesenchyme')
imm <- subset(system.to.cell,subset=class=='Immune')
epi <- ScaleData(epi)
end <- ScaleData(end)
mes <- ScaleData(mes)
imm <- ScaleData(imm)

# load previously computed markers lists
load('class.markers.total.2023-08-15.Robj')

##### explorations of significant markers, legacy # exclude certain mechanisms from analysis ####
# total.mechs <- c(class.markers.total$Epithelium$Dataset3$gene,
#                  class.markers.total$Epithelium$Eng_Dataset2$gene,
#                  class.markers.total$Epithelium$Native$gene)
# to.exclude <- grep('Itgb|Itga|Cav1|Sdc1|Cd44|Lrp1|Lrp2',total.mechs)
# mechs.included <- total.mechs[-to.exclude]
# 
# # determine a data-derived gene list
# start.moi <- class.markers.total$Epithelium$Dataset3[class.markers.total$Epithelium$Dataset3$gene %in% mechs.included,] %>% subset(cluster=='Start') %>% top_n(20,power)
# eng.moi <- class.markers.total$Epithelium$Eng_Dataset2[class.markers.total$Epithelium$Eng_Dataset2$gene %in% mechs.included,] %>% group_by(cluster) %>% top_n(10,power)
# nat.moi <- class.markers.total$Epithelium$Native[class.markers.total$Epithelium$Native$gene %in% mechs.included,] %>% group_by(cluster) %>% top_n(10,power)
# epi.moi <- c(start.moi$gene,
#              eng.moi$gene,
#              nat.moi$gene)
# # determine a dataset3 gene list
# dataset3.moi <- class.markers.total$Epithelium$Dataset3[class.markers.total$Epithelium$Dataset3$gene %in% mechs.included,] %>% group_by(cluster )%>% top_n(50,power)
# epi.moi <- unique(dataset3.moi$gene)
# 
# # determine a native cell type gene list
# native.moi <- class.markers.total$Epithelium$Native[class.markers.total$Epithelium$Native$gene %in% mechs.included,] %>% group_by(cluster)%>% top_n(20,power)
# epi.moi <- unique(native.moi$gene)
# 
# # or make a custom moi
# #epi.moi <- c('Dll4—Notch3','Hgf—Met','Wnt4—Fzd2','Tgfb1—Itgb6','Cthrc1—Fzd6')
# 
# # make a heatmap
# Idents(epi) <- epi$Tissue
# DoHeatmap(epi,features = epi.moi,cells = WhichCells(epi,downsample = 100),group.by = 'Tissue',slot = 'scale.data')
# 
# # what if we look at the data in native space
# i=1
# temp <- subset(system.to.cell,subset=class==class.list[i])
# VlnPlot(temp,'nFeature_SystemToCell',group.by = 'Tissue',log=T)
# temp <- subset(temp,nFeature_SystemToCell>50)
# native <- subset(temp,subset = Dataset3=='Native')
# start <- subset(temp,subset = Dataset3=='Start')
# eng <- subset(temp,subset = Dataset3=='Eng')
# native <- ScaleData(native)
# start <- ScaleData(start)
# eng <- ScaleData(eng)
# temp <- ScaleData(temp)
# temp <- FindVariableFeatures(temp)
# temp <- RunPCA(temp)
# ElbowPlot(temp,ndims=50)
# temp <- RunUMAP(temp,dims = 1:50)
# DimPlot(temp)
# DimPlot(temp,group.by = 'Dataset3')
# DimPlot(temp,group.by = 'Final1')
# FeaturePlot(temp,'Tgfb1—Itgb6')
#####

# Epithelium integration byt Dataset3
temp <- subset(epi,subset = nFeature_SystemToCell>100)
table(temp$Dataset3)
table(temp$Final1)
ifnb.list <- SplitObject(temp, split.by = "Dataset3")
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = ifnb.list)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
epi.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(epi.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
epi.combined <- ScaleData(epi.combined, verbose = FALSE)
epi.combined <- RunPCA(epi.combined, npcs = 100, verbose = FALSE)
ElbowPlot(epi.combined,ndims=100)
# Look carefully at PCs
pdf(file='epi.combined.PCs.pdf',width=10,height=8)
ElbowPlot(epi.combined,ndims = 100)
PCHeatmap(epi.combined,cells=200,balanced=T,dims=1:9)
PCHeatmap(epi.combined,cells=200,balanced=T,dims=10:18)
PCHeatmap(epi.combined,cells=200,balanced=T,dims=19:27)
PCHeatmap(epi.combined,cells=200,balanced=T,dims=28:36)
PCHeatmap(epi.combined,cells=200,balanced=T,dims=37:45)
PCHeatmap(epi.combined,cells=200,balanced=T,dims=46:54)
PCHeatmap(epi.combined,cells=200,balanced=T,dims=55:63)
PCHeatmap(epi.combined,cells=200,balanced=T,dims=64:72)
PCHeatmap(epi.combined,cells=200,balanced=T,dims=73:81)
PCHeatmap(epi.combined,cells=200,balanced=T,dims=82:90)
PCHeatmap(epi.combined,cells=200,balanced=T,dims=91:99)
dev.off()

epi.combined <- RunUMAP(epi.combined, reduction = "pca", dims = 23:53) # trying a new approach to get rid of non-indepednet pcs?
epi.combined <- FindNeighbors(epi.combined, reduction = "pca", dims = 23:53)
epi.combined <- FindClusters(epi.combined, resolution = 0.2)
DimPlot(epi.combined)
DimPlot(epi.combined,group.by = 'Dataset3')
DimPlot(epi.combined,group.by = 'Final1')
epi.combined$Tissue <- factor(epi.combined$Tissue,levels = c('Pseudo.Start','BCL5','BCEC2',
                                                             'BEF1','BEF2','BEF3',
                                                             'BEF12','BEF14','BEF15',
                                                             'BEFM1','BEFM2',
                                                             'BEFM4','BEFM5','BEFM6',
                                                             'fRat','mRat','P0-f','P0-m'))
DimPlot(epi.combined,group.by = 'Tissue')
DimPlot(epi.combined,split.by = 'Tissue',group.by = 'Final1')
epi.combined$Dataset2 <- factor(epi.combined$Dataset2,levels = c('Start','Mono','Co',
                                                                 'Tri_E','Tri_L','Quad_E',
                                                                 'Quad_L','Native'))
DimPlot(epi.combined,group.by = 'Dataset2')
DimPlot(epi.combined,split.by = 'Dataset2',group.by = 'Final1')

p1 <- DimPlot(epi.combined,group.by = 'Tissue',cols = color_pals$sample_colors)
p2 <- DimPlot(epi.combined,group.by = 'Dataset2',cols = color_pals$dataset_colors)
p3 <- DimPlot(epi.combined,group.by = 'Dataset3',cols = color_pals$starteng_colors)
p4 <- DimPlot(epi.combined,group.by = 'Final1',cols = color_pals$epi)

DefaultAssay(epi.combined) <- 'SystemToCell'
FeaturePlot(epi.combined,'Csf1—Csf1r',split.by = 'Dataset2',pt.size = 1,order=T)
FeaturePlot(epi.combined,'Wnt4—Fzd2',split.by = 'Dataset2',pt.size = 1,order=T)
FeaturePlot(epi.combined,'Tgfb1—Tgfbr3',split.by = 'Dataset2',pt.size = 1,order=T)
FeaturePlot(epi.combined,'Csf1—Csf1r',split.by = 'Dataset2',pt.size = 1,order=T)
VlnPlot(epi.combined,'Wnt4—Fzd2',group.by = 'Dataset2')


pdf(file = 'epithelium.system.to.cell.global.integration.pdf',width = 16,height = 14)
plot_grid(p1,p2,p3,p4)
plot_grid(DimPlot(epi.combined,split.by = 'Dataset2',group.by = 'Dataset2',cols = color_pals$dataset_colors),
          DimPlot(epi.combined,split.by = 'Dataset2',group.by = 'Final1',cols = color_pals$epi),
          FeaturePlot(epi.combined,'Tgfb1—Tgfbr3',split.by = 'Dataset2',pt.size = 1,order=T),
          VlnPlot(epi.combined,'Tgfb1—Tgfbr3',group.by = 'Dataset2'),nrow = 4)

dev.off()

save(epi.combined,file = 'epi.combined.system.to.cell.2023-08-15.Robj')

#####

# Mesenchyme integration by Dataset3
temp <- subset(mes,subset = nFeature_SystemToCell>100)
table(temp$Dataset3)
table(temp$Final1)
ifnb.list <- SplitObject(temp, split.by = "Dataset3")
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = ifnb.list)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
mes.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(mes.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
mes.combined <- ScaleData(mes.combined, verbose = FALSE)
mes.combined <- RunPCA(mes.combined, npcs = 100, verbose = FALSE)
ElbowPlot(mes.combined,ndims=100)
# Look carefully at PCs
pdf(file='mes.combined.PCs.pdf',width=10,height=8)
ElbowPlot(mes.combined,ndims = 100)
PCHeatmap(mes.combined,cells=200,balanced=T,dims=1:9)
PCHeatmap(mes.combined,cells=200,balanced=T,dims=10:18)
PCHeatmap(mes.combined,cells=200,balanced=T,dims=19:27)
PCHeatmap(mes.combined,cells=200,balanced=T,dims=28:36)
PCHeatmap(mes.combined,cells=200,balanced=T,dims=37:45)
PCHeatmap(mes.combined,cells=200,balanced=T,dims=46:54)
PCHeatmap(mes.combined,cells=200,balanced=T,dims=55:63)
PCHeatmap(mes.combined,cells=200,balanced=T,dims=64:72)
PCHeatmap(mes.combined,cells=200,balanced=T,dims=73:81)
PCHeatmap(mes.combined,cells=200,balanced=T,dims=82:90)
PCHeatmap(mes.combined,cells=200,balanced=T,dims=91:99)
dev.off()

mes.combined <- RunUMAP(mes.combined, reduction = "pca", dims = 16:60) # trying a new approach to get rid of non-indepednet pcs?
mes.combined <- FindNeighbors(mes.combined, reduction = "pca", dims = 16:60)
mes.combined <- FindClusters(mes.combined, resolution = 0.2)
DimPlot(mes.combined)
DimPlot(mes.combined,group.by = 'Dataset3')
DimPlot(mes.combined,group.by = 'Final1')
mes.combined$Tissue <- factor(mes.combined$Tissue,levels = c('Pseudo.Start','BCL5','BCEC2',
                                                             'BEF1','BEF2','BEF3',
                                                             'BEF12','BEF14','BEF15',
                                                             'BEFM1','BEFM2',
                                                             'BEFM4','BEFM5','BEFM6',
                                                             'fRat','mRat','P0-f','P0-m'))
DimPlot(mes.combined,group.by = 'Tissue')
DimPlot(mes.combined,split.by = 'Tissue',group.by = 'Final1')
mes.combined$Dataset2 <- factor(mes.combined$Dataset2,levels = c('Start','Mono','Co',
                                                                 'Tri_E','Tri_L','Quad_E',
                                                                 'Quad_L','Native'))
DimPlot(mes.combined,group.by = 'Dataset2')
DimPlot(mes.combined,split.by = 'Dataset2',group.by = 'Final1')

p1 <- DimPlot(mes.combined,group.by = 'Tissue',cols = color_pals$sample_colors)
p2 <- DimPlot(mes.combined,group.by = 'Dataset2',cols = color_pals$dataset_colors)
p3 <- DimPlot(mes.combined,group.by = 'Dataset3',cols = color_pals$starteng_colors)
p4 <- DimPlot(mes.combined,group.by = 'Final1',cols = color_pals$mes)

DefaultAssay(mes.combined) <- 'SystemToCell'
View(class.markers.total$Mesenchyme$Native)
FeaturePlot(mes.combined,'Dll4—Notch3',split.by = 'Dataset2',pt.size = 1,order=T,max.cutoff = 1)
VlnPlot(mes.combined,'Dll4—Notch3',group.by = 'Dataset2')

pdf(file = 'mesenchyme.system.to.cell.global.integration.pdf',width = 16,height = 14)
plot_grid(p1,p2,p3,p4)
plot_grid(DimPlot(mes.combined,split.by = 'Dataset2',group.by = 'Dataset2',cols = color_pals$dataset_colors),
          DimPlot(mes.combined,split.by = 'Dataset2',group.by = 'Final1',cols = color_pals$mes),
          FeaturePlot(mes.combined,'Dll4—Notch3',split.by = 'Dataset2',pt.size = 1,order=T,max.cutoff = 2),
          VlnPlot(mes.combined,'Dll4—Notch3',group.by = 'Dataset2'),nrow = 4)
dev.off()

save(mes.combined,file = 'mes.combined.system.to.cell.2023-08-15.Robj')

#####

# Endothelium integration by Dataset3
temp <- subset(end,subset = nFeature_SystemToCell>100)
table(temp$Dataset3)
table(temp$Final1)
ifnb.list <- SplitObject(temp, split.by = "Dataset3")
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = ifnb.list)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
end.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(end.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
end.combined <- ScaleData(end.combined, verbose = FALSE)
end.combined <- RunPCA(end.combined, npcs = 100, verbose = FALSE)
ElbowPlot(end.combined,ndims=100)
# Look carefully at PCs
pdf(file='end.combined.PCs.pdf',width=10,height=8)
ElbowPlot(end.combined,ndims = 100)
PCHeatmap(end.combined,cells=200,balanced=T,dims=1:9)
PCHeatmap(end.combined,cells=200,balanced=T,dims=10:18)
PCHeatmap(end.combined,cells=200,balanced=T,dims=19:27)
PCHeatmap(end.combined,cells=200,balanced=T,dims=28:36)
PCHeatmap(end.combined,cells=200,balanced=T,dims=37:45)
PCHeatmap(end.combined,cells=200,balanced=T,dims=46:54)
PCHeatmap(end.combined,cells=200,balanced=T,dims=55:63)
PCHeatmap(end.combined,cells=200,balanced=T,dims=64:72)
PCHeatmap(end.combined,cells=200,balanced=T,dims=73:81)
PCHeatmap(end.combined,cells=200,balanced=T,dims=82:90)
PCHeatmap(end.combined,cells=200,balanced=T,dims=91:99)
dev.off()

end.combined <- RunUMAP(end.combined, reduction = "pca", dims = 18:39) # trying a new approach to get rid of non-indepednet pcs?
end.combined <- FindNeighbors(end.combined, reduction = "pca", dims = 18:39)
end.combined <- FindClusters(end.combined, resolution = 0.2)
DimPlot(end.combined)
DimPlot(end.combined,group.by = 'Dataset3')
DimPlot(end.combined,group.by = 'Final1')
end.combined$Tissue <- factor(end.combined$Tissue,levels = c('Pseudo.Start','BCL5','BCEC2',
                                                             'BEF1','BEF2','BEF3',
                                                             'BEF12','BEF14','BEF15',
                                                             'BEFM1','BEFM2',
                                                             'BEFM4','BEFM5','BEFM6',
                                                             'fRat','mRat','P0-f','P0-m'))
DimPlot(end.combined,group.by = 'Tissue')
DimPlot(end.combined,split.by = 'Tissue',group.by = 'Final1')
end.combined$Dataset2 <- factor(end.combined$Dataset2,levels = c('Start','Mono','Co',
                                                                 'Tri_E','Tri_L','Quad_E',
                                                                 'Quad_L','Native'))
DimPlot(end.combined,group.by = 'Dataset2')
DimPlot(end.combined,split.by = 'Dataset2',group.by = 'Final1')

p1 <- DimPlot(end.combined,group.by = 'Tissue',cols = color_pals$sample_colors)
p2 <- DimPlot(end.combined,group.by = 'Dataset2',cols = color_pals$dataset_colors)
p3 <- DimPlot(end.combined,group.by = 'Dataset3',cols = color_pals$starteng_colors)
p4 <- DimPlot(end.combined,group.by = 'Final1',cols = color_pals$end)

DefaultAssay(end.combined) <- 'SystemToCell'
View(class.markers.total$Endothelium$Native)
FeaturePlot(end.combined,'Angpt1—Tek',split.by = 'Dataset2',pt.size = 1,order=T,max.cutoff = 3)
VlnPlot(end.combined,'Angpt1—Tek',group.by = 'Dataset2')

pdf(file = 'endothelium.system.to.cell.global.integration.pdf',width = 16,height = 14)
plot_grid(p1,p2,p3,p4)
plot_grid(DimPlot(end.combined,split.by = 'Dataset2',group.by = 'Dataset2',cols = color_pals$dataset_colors),
          DimPlot(end.combined,split.by = 'Dataset2',group.by = 'Final1',cols = color_pals$end),
          FeaturePlot(end.combined,'Angpt1—Tek',split.by = 'Dataset2',pt.size = 1,order=T,max.cutoff = 3),
          VlnPlot(end.combined,'Angpt1—Tek',group.by = 'Dataset2'),nrow = 4)
dev.off()

save(end.combined,file = 'end.combined.system.to.cell.2023-08-15.Robj')

#####

# Immune integration by Dataset3
temp <- subset(imm,subset = nFeature_SystemToCell>100)
table(temp$Dataset3)
table(temp$Final1)
ifnb.list <- SplitObject(temp, split.by = "Dataset3")
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = ifnb.list)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
imm.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(imm.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
imm.combined <- ScaleData(imm.combined, verbose = FALSE)
imm.combined <- RunPCA(imm.combined, npcs = 100, verbose = FALSE)
ElbowPlot(imm.combined,ndims=100)
# Look carefully at PCs
pdf(file='imm.combined.PCs.pdf',width=10,height=8)
ElbowPlot(imm.combined,ndims = 100)
PCHeatmap(imm.combined,cells=200,balanced=T,dims=1:9)
PCHeatmap(imm.combined,cells=200,balanced=T,dims=10:18)
PCHeatmap(imm.combined,cells=200,balanced=T,dims=19:27)
PCHeatmap(imm.combined,cells=200,balanced=T,dims=28:36)
PCHeatmap(imm.combined,cells=200,balanced=T,dims=37:45)
PCHeatmap(imm.combined,cells=200,balanced=T,dims=46:54)
PCHeatmap(imm.combined,cells=200,balanced=T,dims=55:63)
PCHeatmap(imm.combined,cells=200,balanced=T,dims=64:72)
PCHeatmap(imm.combined,cells=200,balanced=T,dims=73:81)
PCHeatmap(imm.combined,cells=200,balanced=T,dims=82:90)
PCHeatmap(imm.combined,cells=200,balanced=T,dims=91:99)
dev.off()

imm.combined <- RunUMAP(imm.combined, reduction = "pca", dims = 25:50) # trying a new approach to get rid of non-indepednet pcs?
imm.combined <- FindNeighbors(imm.combined, reduction = "pca", dims = 25:50)
imm.combined <- FindClusters(imm.combined, resolution = 0.2)
DimPlot(imm.combined)
DimPlot(imm.combined,group.by = 'Dataset3')
DimPlot(imm.combined,group.by = 'Final1')
imm.combined$Tissue <- factor(imm.combined$Tissue,levels = c('Pseudo.Start','BCL5','BCEC2',
                                                             'BEF1','BEF2','BEF3',
                                                             'BEF12','BEF14','BEF15',
                                                             'BEFM1','BEFM2',
                                                             'BEFM4','BEFM5','BEFM6',
                                                             'fRat','mRat','P0-f','P0-m'))
DimPlot(imm.combined,group.by = 'Tissue')
DimPlot(imm.combined,split.by = 'Tissue',group.by = 'Final1')
imm.combined$Dataset2 <- factor(imm.combined$Dataset2,levels = c('Start','Mono','Co',
                                                                 'Tri_E','Tri_L','Quad_E',
                                                                 'Quad_L','Native'))
DimPlot(imm.combined,group.by = 'Dataset2')
DimPlot(imm.combined,split.by = 'Dataset2',group.by = 'Final1')

p1 <- DimPlot(imm.combined,group.by = 'Tissue',cols = color_pals$sample_colors)
p2 <- DimPlot(imm.combined,group.by = 'Dataset2',cols = color_pals$dataset_colors)
p3 <- DimPlot(imm.combined,group.by = 'Dataset3',cols = color_pals$starteng_colors)
p4 <- DimPlot(imm.combined,group.by = 'Final1',cols = color_pals$imm)

DefaultAssay(imm.combined) <- 'SystemToCell'
View(class.markers.total$Immune$Eng_Dataset2)
FeaturePlot(imm.combined,'Tgfb1—Tgfbr3',split.by = 'Dataset2',pt.size = 1,order=T)
VlnPlot(imm.combined,'Tgfb1—Tgfbr3',group.by = 'Dataset2')
FeaturePlot(imm.combined,'Cxcl3—Cxcr1',split.by = 'Dataset2',pt.size = 1,order=T)
VlnPlot(imm.combined,'Wnt5a—Ror2',group.by = 'Dataset2')

pdf(file = 'immune.system.to.cell.global.integration.pdf',width = 16,height = 14)
plot_grid(p1,p2,p3,p4)
plot_grid(DimPlot(imm.combined,split.by = 'Dataset2',group.by = 'Dataset2',cols = color_pals$dataset_colors),
          DimPlot(imm.combined,split.by = 'Dataset2',group.by = 'Final1',cols = color_pals$imm),
          FeaturePlot(imm.combined,'Tgfb1—Tgfbr3',split.by = 'Dataset2',pt.size = 1,order=T),
          VlnPlot(imm.combined,'Tgfb1—Tgfbr3',group.by = 'Dataset2'),nrow = 4)
dev.off()

save(imm.combined,file = 'imm.combined.system.to.cell.2023-08-15.Robj')

# epi.combined$Location <- epi.combined$Final1
# Idents(epi.combined) <- epi.combined$Location
# epi.combined <- RenameIdents(epi.combined,
#                              'ATI'='Alveolar',
#                              'ATI-like'='Alveolar',
#                              'ATII'='Alveolar',
#                              'ATII-ATI'='Alveolar',
#                              'Basal-like'='Proximal',
#                              'BASC'='Alveolar',
#                              'Ciliated'='Proximal',
#                              'Secretory'='Proximal',
#                              'Tuft'='Alveolar',
#                              'Transitional'='Alveolar')
# epi.combined$Location <- Idents(epi.combined)
# DimPlot(epi.combined,split.by = 'Dataset2',group.by = 'Location')
# epi.combined$Dataset4 <- epi.combined$Dataset2
# Idents(epi.combined) <- epi.combined$Dataset4
# epi.combined <- RenameIdents(epi.combined,
#                              'Start'='Start',
#                              'Mono'='Mono',
#                              'Co'='Co',
#                              'Tri_E'='Tri',
#                              'Tri_L'='Tri',
#                              'Quad_E'='Quad',
#                              'Quad_L'='Quad',
#                              'Native'='Native')
# epi.combined$Dataset4 <- Idents(epi.combined)
# png(filename = 'epi.system.to.cell.umaps.dataset4.location.png',width = 12,height = 3,units = 'in',res=300)
# DimPlot(epi.combined,split.by = 'Dataset4',group.by = 'Location',pt.size = 0.1,shuffle = T)
# dev.off()
# DimPlot(epi.combined,split.by = 'Dataset4',group.by = 'Final1',pt.size = 0.5,shuffle = T,cols = color_pals$epi)
# 
# # legacy below
# DefaultAssay(epi.combined) <- "SystemToCell"
# VlnPlot(epi.combined,features = 'Hgf—Met',group.by = 'Tissue')
# VlnPlot(epi.combined,features = 'Wnt4—Fzd2',group.by = 'Tissue')
# VlnPlot(epi.combined,features = 'Tgfb1—Itgb6',group.by = 'Dataset4')
# VlnPlot(epi.combined,features = 'Cthrc1—Fzd6',group.by = 'Tissue')
# 
# VlnPlot(mes,features = 'Dll4—Notch3',group.by = 'Tissue')
# 
# VlnPlot(end,features = 'Vegfa—Flt1',group.by = 'Tissue')


# attempt to integrate by sample mirroring allies process ######
# # integrate by sample, exlucing native, to mirror the phenotype workflow?
# temp <- subset(epi,subset=Dataset3!='Native')
# # subset to only high information measurements
# temp <- subset(temp,nFeature_SystemToCell>100)
# 
# table(temp$Dataset3)
# table(temp$Tissue)
# ifnb.list <- SplitObject(temp, split.by = "Tissue")
# ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })
# features <- SelectIntegrationFeatures(object.list = ifnb.list)
# immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
# epi.combined <- IntegrateData(anchorset = immune.anchors)
# DefaultAssay(epi.combined) <- "integrated"
# 
# # Run the standard workflow for visualization and clustering
# epi.combined <- ScaleData(epi.combined, verbose = FALSE)
# epi.combined <- RunPCA(epi.combined, npcs = 100, verbose = FALSE)
# ElbowPlot(epi.combined,ndims=100)
# # Look carefully at PCs
# pdf(file='epi.no.native.combined.PCs.pdf',width=10,height=8)
# ElbowPlot(epi.combined,ndims = 100)
# PCHeatmap(epi.combined,cells=200,balanced=T,dims=1:9)
# PCHeatmap(epi.combined,cells=200,balanced=T,dims=10:18)
# PCHeatmap(epi.combined,cells=200,balanced=T,dims=19:27)
# PCHeatmap(epi.combined,cells=200,balanced=T,dims=28:36)
# PCHeatmap(epi.combined,cells=200,balanced=T,dims=37:45)
# PCHeatmap(epi.combined,cells=200,balanced=T,dims=46:54)
# PCHeatmap(epi.combined,cells=200,balanced=T,dims=55:63)
# PCHeatmap(epi.combined,cells=200,balanced=T,dims=64:72)
# PCHeatmap(epi.combined,cells=200,balanced=T,dims=73:81)
# PCHeatmap(epi.combined,cells=200,balanced=T,dims=82:90)
# PCHeatmap(epi.combined,cells=200,balanced=T,dims=91:99)
# dev.off()
# 
# epi.combined <- RunUMAP(epi.combined, reduction = "pca", dims = 1:50)
# epi.combined <- FindNeighbors(epi.combined, reduction = "pca", dims = 1:50)
# epi.combined <- FindClusters(epi.combined, resolution = 0.2)
# DimPlot(epi.combined)
# DimPlot(epi.combined,group.by = 'Dataset3')
# DimPlot(epi.combined,group.by = 'Final1')
# DimPlot(epi.combined,group.by = 'Tissue')
# epi.combined$Tissue <- factor(epi.combined$Tissue,levels = c('Pseudo.Start','BCL5','BCEC2',
#                                                              'BEF1','BEF2','BEF3',
#                                                              'BEF12','BEF14','BEF15',
#                                                              'BEFM1','BEFM2',
#                                                              'BEFM4','BEFM5','BEFM6',
#                                                              'fRat','mRat','P0-f','P0-m'))
# DimPlot(epi.combined,split.by = 'Tissue',group.by = 'Final1')
# 
# epi.combined$Dataset2 <- factor(epi.combined$Dataset2,levels = c('Start','Mono','Co',
#                                                                  'Tri_E','Tri_L','Quad_E',
#                                                                  'Quad_L','Native'))
# DimPlot(epi.combined,group.by = 'Dataset2')
# DimPlot(epi.combined,split.by = 'Dataset2',group.by = 'Final1')