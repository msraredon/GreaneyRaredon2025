# This assembles and cleans previous code to make a replicable workflow generating the EpithelialReceiving object for the BEFM project.

# This script includes object definition, clustering, embedding, annotation, organization, differential expression, and saving.

# The work here is adapted from, originally, 
## "EpithelialReceivingWorking_v2.Rmd" and
## "EpiRecFeatures.Rmd"

# Input objects:
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Saved Objects/global.connectomics.2023-11-11.Robj")
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Allie's Object Filtration/Final Objects/all.color_palettes_2.R")

# Output objects:
# save(epi.rec,file = 'epi.rec.2024-03-15.Robj')
# save(mark.epi.rec,file = 'mark.epi.rec.2024-03-15.Robj')

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

# Isolate CellToCell edges that land on Epithelium within this limited dataset. 
temp <- global.connectomics$CellToCell$alra # we will use the imputed data, for now
Idents(temp) <- temp$class.Receiving
table(Idents(temp))
epi.rec <- subset(temp, idents = 'Epithelium')
Idents(epi.rec) <- epi.rec$orig.ident
table(Idents(epi.rec))
epi.rec <- subset(epi.rec,idents = c(samples.of.interest))
table(Idents(epi.rec))

## Clean
# Since BCL5 is not in this limited dataset, we can clean these data with reasonable thresholds to make the plots prettier
VlnPlot(epi.rec,'nFeature_CellToCell',group.by = 'orig.ident')
epi.rec <- subset(epi.rec,subset = nFeature_CellToCell>100)
VlnPlot(epi.rec,'nFeature_CellToCell',group.by = 'orig.ident')


# Excellent. Let's take a global look at this edge category before moving further:

#### Epithelium Receiving - Embed and Cluster

epi.rec <- ScaleData(epi.rec)
epi.rec <- FindVariableFeatures(epi.rec)
epi.rec <- RunPCA(epi.rec,npcs = 50)
ElbowPlot(epi.rec,ndims = 50)
epi.rec <- RunUMAP(epi.rec,dims = 1:20)
epi.rec <- FindNeighbors(epi.rec,dims = 1:20)
epi.rec <- FindClusters(epi.rec,resolution = 0.1)


## Epi Plot (before annotation)
p1 <- DimPlot(epi.rec,group.by = 'class.Sending',shuffle = T,cols = color_pals$class_colors)+ggtitle('Sending CellClass')+NoAxes()
p2 <- DimPlot(epi.rec,group.by = 'class.Receiving',shuffle = T,cols = color_pals$class_colors)+ggtitle('Receiving CellClass')+NoAxes()
p3 <- DimPlot(epi.rec,group.by = 'Dataset2.Sending',shuffle = T,cols = color_pals$dataset_colors)+ggtitle('Dataset2')+NoAxes()
p4 <- DimPlot(epi.rec,group.by = 'SendingType',shuffle = T,cols = color_pals$type_colors)+ggtitle('Sending CellType')+NoAxes()
p5 <- DimPlot(epi.rec,group.by = 'ReceivingType',shuffle = T,cols = color_pals$type_colors)+ggtitle('Receiving CellType')+NoAxes()
p6 <- DimPlot(epi.rec,group.by = 'orig.ident',shuffle = T,cols = color_pals$sample_colors)+ggtitle('Sample')+NoAxes()

# Define plot blocks
small.mult <- cowplot::plot_grid(p1,p2,p3,p4,p5,p6,ncol = 3)
sig.arch <- DimPlot(epi.rec,group.by = 'seurat_clusters',shuffle = T,label = T)+NoAxes()+ggtitle('Signaling Archetypes')
split.arch <- DimPlot(epi.rec,split.by = 'orig.ident',group.by = 'seurat_clusters',shuffle = T,ncol = 2)+NoAxes()+ggtitle('Signaling Archetypes (Split By Sample)')

# Assemble blocks in cowplot
top <- small.mult
bottom <- cowplot::plot_grid(sig.arch,split.arch,rel_widths = c(1.2,1))
plot <- cowplot::plot_grid(top,bottom,nrow = 2)
plot

png('Epithelial_Receiving_2024-03-15.png',width=24,height = 22,units = 'in',res=300)
print(plot)
dev.off()

# Add annotations
Idents(epi.rec) <- epi.rec$seurat_clusters
table(Idents(epi.rec))
epi.rec <- RenameIdents(epi.rec,
                        '0'='Mes-Epi Common',
                        '1'='Epi-Auto Common',
                        '2'='Endo-Epi Common',
                        '3'='Epi-Auto Regenerative',
                        '4'='Mes-Epi Regenerative',
                        '5'='Macrophage-Epi')
table(Idents(epi.rec))
epi.rec$Archetype <- Idents(epi.rec)

## Epi Plot (after annotation)
p1 <- DimPlot(epi.rec,group.by = 'class.Sending',shuffle = T,cols = color_pals$class_colors)+ggtitle('Sending CellClass')+NoAxes()
p2 <- DimPlot(epi.rec,group.by = 'class.Receiving',shuffle = T,cols = color_pals$class_colors)+ggtitle('Receiving CellClass')+NoAxes()
p3 <- DimPlot(epi.rec,group.by = 'Dataset2.Sending',shuffle = T,cols = color_pals$dataset_colors)+ggtitle('Dataset2')+NoAxes()
p4 <- DimPlot(epi.rec,group.by = 'SendingType',shuffle = T,cols = color_pals$type_colors)+ggtitle('Sending CellType')+NoAxes()
p5 <- DimPlot(epi.rec,group.by = 'ReceivingType',shuffle = T,cols = color_pals$type_colors)+ggtitle('Receiving CellType')+NoAxes()
p6 <- DimPlot(epi.rec,group.by = 'orig.ident',shuffle = T,cols = color_pals$sample_colors)+ggtitle('Sample')+NoAxes()

# Define plot blocks
small.mult <- cowplot::plot_grid(p1,p2,p3,p4,p5,p6,ncol = 3)
sig.arch <- DimPlot(epi.rec,group.by = 'Archetype',shuffle = T,label = T)+NoAxes()+ggtitle('Signaling Archetypes')
split.arch <- DimPlot(epi.rec,split.by = 'orig.ident',group.by = 'seurat_clusters',shuffle = T,ncol = 2)+NoAxes()+ggtitle('Signaling Archetypes (Split By Sample)')

# Assemble blocks in cowplot
top <- small.mult
bottom <- cowplot::plot_grid(sig.arch,split.arch,rel_widths = c(1.2,1))
plot <- cowplot::plot_grid(top,bottom,nrow = 2)
plot

png('Epithelial_Receiving_Annotated_2024-03-15.png',width=24,height = 22,units = 'in',res=300)
print(plot)
dev.off()

# Differential expression across archetype
mark.epi.rec <- FindAllMarkers(epi.rec,min.pct = 0.1,logfc.threshold = 0.1)
mark.epi.rec$ratio <- mark.epi.rec$pct.1/mark.epi.rec$pct.2
mark.epi.rec$power <- mark.epi.rec$ratio*mark.epi.rec$avg_log2FC

# Clean up: reorder class levels for publication
epi.rec$class.Sending <- factor(epi.rec$class.Sending,levels = c('Epithelium','Endothelium','Mesenchyme','Immune'))
epi.rec$class.Receiving <- factor(epi.rec$class.Receiving,levels = c('Epithelium','Endothelium','Mesenchyme','Immune'))

## Save objects for later
save(epi.rec,file = 'epi.rec.2024-03-15.Robj')
save(mark.epi.rec,file = 'mark.epi.rec.2024-03-15.Robj')
