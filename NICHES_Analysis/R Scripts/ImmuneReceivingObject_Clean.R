# This assembles and cleans previous code to make a replicable workflow generating the ImmuneReceiving object for the BEFM project.

# This script includes object definition, clustering, embedding, annotation, organization, differential expression, and saving.

# Input objects:
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Saved Objects/global.connectomics.2023-11-11.Robj")
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Allie's Object Filtration/Final Objects/all.color_palettes_2.R")

# Output objects:
# save(imm.rec,file = 'imm.rec.2024-03-15.Robj')
# save(mark.imm.rec,file = 'mark.imm.rec.2024-03-15.Robj')

# Set WD
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics")

# Require packages
require(Seurat)
require(NICHES)
require(dplyr)
require(ggthemes)
require(ggplot2)
require(cowplot)
require(knitr)

# Inspect object
global.connectomics

# Define samples of interest
samples.of.interest <- c('BEFM1','BEFM2')

# Build the type color palette
color_pals$type_colors <- c('#A40606','#9CFFFA','#B0DB43','#9C528B','#2F6690',
                            '#946846','#F1C40F','green','#0F0326','#E65F5C','#14591D','#726DA8',
                            'yellow','purple','blue','red','orange','darkgrey','magenta')
names(color_pals$type_colors) <- unique(global.connectomics$CellToCell$alra$SendingType)

# allie.colors <- c(color_pals$epi,color_pals$endo,color_pals$mes,color_pals$imm2,color_pals$imm)
# connectomic.types <- unique(global.connectomics$CellToCell$alra$SendingType)
# 
# # These types have colors assigned already
# found <- connectomic.types[connectomic.types%in% names(allie.colors)]
# # These types have subtly different names which need to be 'found'
# lost <- connectomic.types[!(connectomic.types%in% names(allie.colors))]
# lost
# lost2 <-  c("Cycling",
#            "Endo_Progenitor",
#            "Alveolar_Fibroblast","Aberrant_Fibroblast","Cycling_Mes","Remodeling_Fibroblast",  
#             "Alveolar_Macrophage_Activated","Interstitial_Macrophage","Alveolar_Macrophage_Naive","Cycling_Alveolar_Macrophage")
# colors.lost.2 <- c("#E78AC3", 
#                    "#A6D854",
#                    "#66C2A5","#8DA0CB","#FC8D62","#E78AC3",
#                    "#8DA0CB"...)
# 
# 
# 
# # Assigning the missing colors manually
# found_colors <- allie.colors[found]
# lost_colors <- allie.colors[c(4,10,11,
#                               14,12,13,
#                               )]
#   
# names(color_pals$type_colors) <- unique(global.connectomics$CellToCell$alra$SendingType)

# Isolate CellToCell edges that land on Mesenchyme within this limited dataset. 
temp <- global.connectomics$CellToCell$alra # we will use the imputed data, for now
Idents(temp) <- temp$class.Receiving
table(Idents(temp))
imm.rec <- subset(temp, idents = 'Immune')
Idents(imm.rec) <- imm.rec$orig.ident
table(Idents(imm.rec))
imm.rec <- subset(imm.rec,idents = c(samples.of.interest))
table(Idents(imm.rec))

## Clean
# Since BCL5 is not in this limited dataset, we can clean these data with reasonable thresholds to make the plots prettier
VlnPlot(imm.rec,'nFeature_CellToCell',group.by = 'orig.ident')
imm.rec <- subset(imm.rec,subset = nFeature_CellToCell>100)
VlnPlot(imm.rec,'nFeature_CellToCell',group.by = 'orig.ident')


# Excellent. Let's take a global look at this edge category before moving further:

#### Mesenchyme Receiving - Embed and Cluster
imm.rec <- ScaleData(imm.rec)
imm.rec <- FindVariableFeatures(imm.rec)
imm.rec <- RunPCA(imm.rec,npcs = 100)
ElbowPlot(imm.rec,ndims = 50)
imm.rec <- RunUMAP(imm.rec,dims = 1:10)
imm.rec <- FindNeighbors(imm.rec,dims = 1:10)
imm.rec <- FindClusters(imm.rec,resolution = 0.2)

## MesPlot (before annotation)
p1 <- DimPlot(imm.rec,group.by = 'class.Sending',shuffle = T,cols = color_pals$class_colors)+ggtitle('Sending CellClass')+NoAxes()
p2 <- DimPlot(imm.rec,group.by = 'class.Receiving',shuffle = T,cols = color_pals$class_colors)+ggtitle('Receiving CellClass')+NoAxes()
p3 <- DimPlot(imm.rec,group.by = 'Dataset2.Sending',shuffle = T,cols = color_pals$dataset_colors)+ggtitle('Dataset2')+NoAxes()
p4 <- DimPlot(imm.rec,group.by = 'SendingType',shuffle = T,cols = color_pals$type_colors)+ggtitle('Sending CellType')+NoAxes()
p5 <- DimPlot(imm.rec,group.by = 'ReceivingType',shuffle = T,cols = color_pals$type_colors)+ggtitle('Receiving CellType')+NoAxes()
p6 <- DimPlot(imm.rec,group.by = 'orig.ident',shuffle = T,cols = color_pals$sample_colors)+ggtitle('Sample')+NoAxes()

# Define plot blocks
small.mult <- cowplot::plot_grid(p1,p2,p3,p4,p5,p6,ncol = 3)
sig.arch <- DimPlot(imm.rec,group.by = 'seurat_clusters',shuffle = T,label = T)+ggtitle('Signaling Archetypes')+NoAxes()
split.arch <- DimPlot(imm.rec,split.by = 'orig.ident',group.by = 'seurat_clusters',shuffle = T,ncol = 2)+NoAxes()+ggtitle('Signaling Archetypes (Split By Sample)')

# Assemble blocks in cowplot
top <- small.mult
bottom <- cowplot::plot_grid(sig.arch,split.arch,rel_widths = c(1.2,1))
plot <- cowplot::plot_grid(top,bottom,nrow = 2)
plot

png('Immune_Receiving_2024-03-15.png',width=24,height = 22,units = 'in',res=300)
print(plot)
dev.off()

# Add annotations
Idents(imm.rec) <- imm.rec$seurat_clusters
table(Idents(imm.rec))
imm.rec <- RenameIdents(imm.rec,
                        '0'='Mes-Immune BEFM2',
                        '1'='Epi-Immune BEFM2',
                        '2'='Endo-Immune Common',
                        '3'='Immune-Auto Common',
                        '4'='Epi-Immune BEFM1',
                        '5'='Mes-Immune BEFM1')
table(Idents(imm.rec))
imm.rec$Archetype <- Idents(imm.rec)

## Imm Plot (after annotation)
p1 <- DimPlot(imm.rec,group.by = 'class.Sending',shuffle = T,cols = color_pals$class_colors)+ggtitle('Sending CellClass')+NoAxes()
p2 <- DimPlot(imm.rec,group.by = 'class.Receiving',shuffle = T,cols = color_pals$class_colors)+ggtitle('Receiving CellClass')+NoAxes()
p3 <- DimPlot(imm.rec,group.by = 'Dataset2.Sending',shuffle = T,cols = color_pals$dataset_colors)+ggtitle('Dataset2')+NoAxes()
p4 <- DimPlot(imm.rec,group.by = 'SendingType',shuffle = T,cols = color_pals$type_colors)+ggtitle('Sending CellType')+NoAxes()
p5 <- DimPlot(imm.rec,group.by = 'ReceivingType',shuffle = T,cols = color_pals$type_colors)+ggtitle('Receiving CellType')+NoAxes()
p6 <- DimPlot(imm.rec,group.by = 'orig.ident',shuffle = T,cols = color_pals$sample_colors)+ggtitle('Sample')+NoAxes()

# Define plot blocks
small.mult <- cowplot::plot_grid(p1,p2,p3,p4,p5,p6,ncol = 3)
sig.arch <- DimPlot(imm.rec,group.by = 'Archetype',shuffle = T,label = T)+NoAxes()+ggtitle('Signaling Archetypes')
split.arch <- DimPlot(imm.rec,split.by = 'orig.ident',group.by = 'seurat_clusters',shuffle = T,ncol = 2)+NoAxes()+ggtitle('Signaling Archetypes (Split By Sample)')

# Assemble blocks in cowplot
top <- small.mult
bottom <- cowplot::plot_grid(sig.arch,split.arch,rel_widths = c(1.2,1))
plot <- cowplot::plot_grid(top,bottom,nrow = 2)
plot

png('Immune_Receiving_Annotated_2024-03-15.png',width=24,height = 22,units = 'in',res=300)
print(plot)
dev.off()

# Differential expression across archetype
mark.imm.rec <- FindAllMarkers(imm.rec,min.pct = 0.1,logfc.threshold = 0.1)
mark.imm.rec$ratio <- mark.imm.rec$pct.1/mark.imm.rec$pct.2
mark.imm.rec$power <- mark.imm.rec$ratio*mark.imm.rec$avg_log2FC

# Clean up: reorder class levels for publication
imm.rec$class.Sending <- factor(imm.rec$class.Sending,levels = c('Epithelium','Endothelium','Mesenchyme','Immune'))
imm.rec$class.Receiving <- factor(imm.rec$class.Receiving,levels = c('Epithelium','Endothelium','Mesenchyme','Immune'))

## Save objects for later
save(imm.rec,file = 'imm.rec.2024-03-15.Robj')
save(mark.imm.rec,file = 'mark.imm.rec.2024-03-15.Robj')

