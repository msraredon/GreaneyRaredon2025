# This assembles and cleans previous code to make a replicable workflow generating the EndothelialReceiving object for the BEFM project.

# This script includes object definition, clustering, embedding, annotation, organization, differential expression, and saving.

# The work here is adapted from, originally, 
## "EndothelialReceivingWorking_v2.Rmd" and
## "EndoRecFeatures.Rmd"

# Input objects:
## Global connectomic object
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Saved Objects/global.connectomics.2023-11-11.Robj")
## Color palette
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Allie's Object Filtration/Final Objects/all.color_palettes_2.R")

# Output objects:
# save(end.rec,file = 'end.rec.2024-03-15.Robj')
# save(mark.end.rec,file = 'mark.end.rec.2024-03-15.Robj')

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

# Isolate CellToCell edges that land on Endothelium within this limited dataset. 
temp <- global.connectomics$CellToCell$alra # we will use the imputed data, for now
Idents(temp) <- temp$class.Receiving
table(Idents(temp))
end.rec <- subset(temp, idents = 'Endothelium')
Idents(end.rec) <- end.rec$orig.ident
table(Idents(end.rec))
end.rec <- subset(end.rec,idents = c(samples.of.interest))
table(Idents(end.rec))

## Clean
# Since BCL5 is not in this limited dataset, we can clean these data with reasonable thresholds to make the plots prettier
VlnPlot(end.rec,'nFeature_CellToCell',group.by = 'orig.ident')
end.rec <- subset(end.rec,subset = nFeature_CellToCell>100)
VlnPlot(end.rec,'nFeature_CellToCell',group.by = 'orig.ident')


# Excellent. Let's take a global look at this edge category before moving further:

#### Endothelium Receiving - Embed and Cluster
end.rec <- ScaleData(end.rec)
end.rec <- FindVariableFeatures(end.rec)
end.rec <- RunPCA(end.rec,npcs = 50)
ElbowPlot(end.rec,ndims = 50)
end.rec <- RunUMAP(end.rec,dims = 1:10)
end.rec <- FindNeighbors(end.rec,dims = 1:10)
end.rec <- FindClusters(end.rec,resolution = 0.1)


## Endo Plot (before annotation)
p1 <- DimPlot(end.rec,group.by = 'class.Sending',shuffle = T,cols = color_pals$class_colors)+ggtitle('Sending CellClass')+NoAxes()
p2 <- DimPlot(end.rec,group.by = 'class.Receiving',shuffle = T,cols = color_pals$class_colors)+ggtitle('Receiving CellClass')+NoAxes()
p3 <- DimPlot(end.rec,group.by = 'Dataset2.Sending',shuffle = T,cols = color_pals$dataset_colors)+ggtitle('Dataset2')+NoAxes()
p4 <- DimPlot(end.rec,group.by = 'SendingType',shuffle = T,cols = color_pals$type_colors)+ggtitle('Sending CellType')+NoAxes()
p5 <- DimPlot(end.rec,group.by = 'ReceivingType',shuffle = T,cols = color_pals$type_colors)+ggtitle('Receiving CellType')+NoAxes()
p6 <- DimPlot(end.rec,group.by = 'orig.ident',shuffle = T,cols = color_pals$sample_colors)+ggtitle('Sample')+NoAxes()

# Define plot blocks
small.mult <- cowplot::plot_grid(p1,p2,p3,p4,p5,p6,ncol = 3)
sig.arch <- DimPlot(end.rec,group.by = 'seurat_clusters',shuffle = T,label = T)+NoAxes()+ggtitle('Signaling Archetypes')
split.arch <- DimPlot(end.rec,split.by = 'orig.ident',group.by = 'seurat_clusters',shuffle = T,ncol = 2)+NoAxes()+ggtitle('Signaling Archetypes (Split By Sample)')

# Assemble blocks in cowplot
top <- small.mult
bottom <- cowplot::plot_grid(sig.arch,split.arch,rel_widths = c(1.2,1))
plot <- cowplot::plot_grid(top,bottom,nrow = 2)
plot

png('Endothelial_Receiving_2024-03-15.png',width=24,height = 22,units = 'in',res=300)
print(plot)
dev.off()

# Add annotations
Idents(end.rec) <- end.rec$seurat_clusters
table(Idents(end.rec))
end.rec <- RenameIdents(end.rec,
                        '0'='Baseline Mes-Endo',
                        '1'='Baseline Epi-Endo',
                        '2'='Baseline Endo-Auto',
                        '3'='Regenerative Mes-Endo',
                        '4'='Macrophage-Endo')
table(Idents(end.rec))
end.rec$Archetype <- Idents(end.rec)

## Endo Plot (after annotation)
p1 <- DimPlot(end.rec,group.by = 'class.Sending',shuffle = T,cols = color_pals$class_colors)+ggtitle('Sending CellClass')+NoAxes()
p2 <- DimPlot(end.rec,group.by = 'class.Receiving',shuffle = T,cols = color_pals$class_colors)+ggtitle('Receiving CellClass')+NoAxes()
p3 <- DimPlot(end.rec,group.by = 'Dataset2.Sending',shuffle = T,cols = color_pals$dataset_colors)+ggtitle('Dataset2')+NoAxes()
p4 <- DimPlot(end.rec,group.by = 'SendingType',shuffle = T,cols = color_pals$type_colors)+ggtitle('Sending CellType')+NoAxes()
p5 <- DimPlot(end.rec,group.by = 'ReceivingType',shuffle = T,cols = color_pals$type_colors)+ggtitle('Receiving CellType')+NoAxes()
p6 <- DimPlot(end.rec,group.by = 'orig.ident',shuffle = T,cols = color_pals$sample_colors)+ggtitle('Sample')+NoAxes()

# Define plot blocks
small.mult <- cowplot::plot_grid(p1,p2,p3,p4,p5,p6,ncol = 3)
sig.arch <- DimPlot(end.rec,group.by = 'Archetype',shuffle = T,label = T)+NoAxes()+ggtitle('Signaling Archetypes')
split.arch <- DimPlot(end.rec,split.by = 'orig.ident',group.by = 'seurat_clusters',shuffle = T,ncol = 2)+NoAxes()+ggtitle('Signaling Archetypes (Split By Sample)')

# Assemble blocks in cowplot
top <- small.mult
bottom <- cowplot::plot_grid(sig.arch,split.arch,rel_widths = c(1.2,1))
plot <- cowplot::plot_grid(top,bottom,nrow = 2)
plot

png('Endothelial_Receiving_Annotated_2024-03-15.png',width=24,height = 22,units = 'in',res=300)
print(plot)
dev.off()

# Differential expression across archetype
mark.end.rec <- FindAllMarkers(end.rec,min.pct = 0.1,logfc.threshold = 0.1)
mark.end.rec$ratio <- mark.end.rec$pct.1/mark.end.rec$pct.2
mark.end.rec$power <- mark.end.rec$ratio*mark.end.rec$avg_log2FC

# Clean up: reorder class levels for publication
end.rec$class.Sending <- factor(end.rec$class.Sending,levels = c('Epithelium','Endothelium','Mesenchyme','Immune'))
end.rec$class.Receiving <- factor(end.rec$class.Receiving,levels = c('Epithelium','Endothelium','Mesenchyme','Immune'))

## Save objects for later
save(end.rec,file = 'end.rec.2024-03-15.Robj')
save(mark.end.rec,file = 'mark.end.rec.2024-03-15.Robj')

