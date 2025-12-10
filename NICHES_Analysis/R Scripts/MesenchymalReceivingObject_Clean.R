# This assembles and cleans previous code to make a replicable workflow generating the MesenchymalReceiving object for the BEFM project.

# This script includes object definition, clustering, embedding, annotation, organization, differential expression, and saving.

# The work here is adapted from, originally, 
## "MesenchymalReceivingWorking_v2.Rmd" and
## "MesRecFeatures_v2.Rmd"

# Input objects:
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Saved Objects/global.connectomics.2023-11-11.Robj")
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Allie's Object Filtration/Final Objects/all.color_palettes_2.R")

# Output objects:
# save(mes.rec,file = 'mes.rec.2024-03-15.Robj')
# save(mark.mes.rec,file = 'mark.mes.rec.2024-03-15.Robj')

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
samples.of.interest <- c('BEF12','BEF14','BEF15','BEFM1','BEFM2')

# Add a temporary color palette allowing SendingType visualization
color_pals$type_colors <- c('#A40606','#9CFFFA','#B0DB43','#9C528B','#2F6690',
                            '#946846','#F1C40F','green','#0F0326','#E65F5C','#14591D','#726DA8',
                            'yellow','purple','blue','red','orange','darkgrey','magenta')
names(color_pals$type_colors) <- unique(global.connectomics$CellToCell$alra$SendingType)

# Isolate CellToCell edges that land on Mesenchyme within this limited dataset. 
temp <- global.connectomics$CellToCell$alra # we will use the imputed data, for now
Idents(temp) <- temp$class.Receiving
table(Idents(temp))
mes.rec <- subset(temp, idents = 'Mesenchyme')
Idents(mes.rec) <- mes.rec$orig.ident
table(Idents(mes.rec))
mes.rec <- subset(mes.rec,idents = c(samples.of.interest))
table(Idents(mes.rec))

## Clean
# Since BCL5 is not in this limited dataset, we can clean these data with reasonable thresholds to make the plots prettier
VlnPlot(mes.rec,'nFeature_CellToCell',group.by = 'orig.ident')
mes.rec <- subset(mes.rec,subset = nFeature_CellToCell>100)
VlnPlot(mes.rec,'nFeature_CellToCell',group.by = 'orig.ident')


# Excellent. Let's take a global look at this edge category before moving further:

#### Mesenchyme Receiving - Embed and Cluster
mes.rec <- ScaleData(mes.rec)
mes.rec <- FindVariableFeatures(mes.rec)
mes.rec <- RunPCA(mes.rec,npcs = 100)
ElbowPlot(mes.rec,ndims = 50)
mes.rec <- RunUMAP(mes.rec,dims = 1:20)
mes.rec <- FindNeighbors(mes.rec,dims = 1:20)
mes.rec <- FindClusters(mes.rec,resolution = 0.2)

## MesPlot (before annotation)
p1 <- DimPlot(mes.rec,group.by = 'class.Sending',shuffle = T,cols = color_pals$class_colors)+ggtitle('Sending CellClass')+NoAxes()
p2 <- DimPlot(mes.rec,group.by = 'class.Receiving',shuffle = T,cols = color_pals$class_colors)+ggtitle('Receiving CellClass')+NoAxes()
p3 <- DimPlot(mes.rec,group.by = 'Dataset2.Sending',shuffle = T,cols = color_pals$dataset_colors)+ggtitle('Dataset2')+NoAxes()
p4 <- DimPlot(mes.rec,group.by = 'SendingType',shuffle = T,cols = color_pals$type_colors)+ggtitle('Sending CellType')+NoAxes()
p5 <- DimPlot(mes.rec,group.by = 'ReceivingType',shuffle = T,cols = color_pals$type_colors)+ggtitle('Receiving CellType')+NoAxes()
p6 <- DimPlot(mes.rec,group.by = 'orig.ident',shuffle = T,cols = color_pals$sample_colors)+ggtitle('Sample')+NoAxes()

# Define plot blocks
small.mult <- cowplot::plot_grid(p1,p2,p3,p4,p5,p6,ncol = 3)
sig.arch <- DimPlot(mes.rec,group.by = 'seurat_clusters',shuffle = T,label = T)+NoAxes()+ggtitle('Signaling Archetypes')
split.arch <- DimPlot(mes.rec,split.by = 'orig.ident',group.by = 'seurat_clusters',shuffle = T,ncol = 2)+NoAxes()+ggtitle('Signaling Archetypes (Split By Sample)')

# Assemble blocks in cowplot
top <- small.mult
bottom <- cowplot::plot_grid(sig.arch,split.arch,rel_widths = c(1.2,1))
plot <- cowplot::plot_grid(top,bottom,nrow = 2)
plot

png('Mesenchymal_Receiving_2024-03-15.png',width=24,height = 22,units = 'in',res=300)
print(plot)
dev.off()

# Add annotations
Idents(mes.rec) <- mes.rec$seurat_clusters
table(Idents(mes.rec))
mes.rec <- RenameIdents(mes.rec,
                        '0'='Mes-Auto BEF12',
                        '1'='Epi-Mes Common',
                        '2'='Mes-Auto BEF14',
                        '3'='Mes-Auto BEF15',
                        '4'='Mes-Auto BEFM2',
                        '5'='Endo-Mes Common',
                        '6'='Epi-Mes BEFM1',
                        '7'='Mes-Auto BEFM1',
                        '8'='Macrophage-Mes Common')
table(Idents(mes.rec))
mes.rec$Archetype <- Idents(mes.rec)

## Epi Plot (after annotation)
p1 <- DimPlot(mes.rec,group.by = 'class.Sending',shuffle = T,cols = color_pals$class_colors)+ggtitle('Sending CellClass')+NoAxes()
p2 <- DimPlot(mes.rec,group.by = 'class.Receiving',shuffle = T,cols = color_pals$class_colors)+ggtitle('Receiving CellClass')+NoAxes()
p3 <- DimPlot(mes.rec,group.by = 'Dataset2.Sending',shuffle = T,cols = color_pals$dataset_colors)+ggtitle('Dataset2')+NoAxes()
p4 <- DimPlot(mes.rec,group.by = 'SendingType',shuffle = T,cols = color_pals$type_colors)+ggtitle('Sending CellType')+NoAxes()
p5 <- DimPlot(mes.rec,group.by = 'ReceivingType',shuffle = T,cols = color_pals$type_colors)+ggtitle('Receiving CellType')+NoAxes()
p6 <- DimPlot(mes.rec,group.by = 'orig.ident',shuffle = T,cols = color_pals$sample_colors)+ggtitle('Sample')+NoAxes()

# Define plot blocks
small.mult <- cowplot::plot_grid(p1,p2,p3,p4,p5,p6,ncol = 3)
sig.arch <- DimPlot(mes.rec,group.by = 'Archetype',shuffle = T,label = T)+NoAxes()+ggtitle('Signaling Archetypes')
split.arch <- DimPlot(mes.rec,split.by = 'orig.ident',group.by = 'seurat_clusters',shuffle = T,ncol = 2)+NoAxes()+ggtitle('Signaling Archetypes (Split By Sample)')

# Assemble blocks in cowplot
top <- small.mult
bottom <- cowplot::plot_grid(sig.arch,split.arch,rel_widths = c(1.2,1))
plot <- cowplot::plot_grid(top,bottom,nrow = 2)
plot

png('Mesenchymal_Receiving_Annotated_2024-03-15.png',width=24,height = 22,units = 'in',res=300)
print(plot)
dev.off()

# Differential expression across archetype
mark.mes.rec <- FindAllMarkers(mes.rec,min.pct = 0.1,logfc.threshold = 0.1)
mark.mes.rec$ratio <- mark.mes.rec$pct.1/mark.mes.rec$pct.2
mark.mes.rec$power <- mark.mes.rec$ratio*mark.mes.rec$avg_log2FC

# Clean up: reorder class levels for publication
mes.rec$class.Sending <- factor(mes.rec$class.Sending,levels = c('Epithelium','Endothelium','Mesenchyme','Immune'))
mes.rec$class.Receiving <- factor(mes.rec$class.Receiving,levels = c('Epithelium','Endothelium','Mesenchyme','Immune'))

## Save objects for later
save(mes.rec,file = 'mes.rec.2024-03-15.Robj')
save(mark.mes.rec,file = 'mark.mes.rec.2024-03-15.Robj')
