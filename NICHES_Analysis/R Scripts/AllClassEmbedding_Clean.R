# 2024-03-15
# MSBR

# This script is a clean, replicable workflow that makes a custom embedding of 5 specific samples of interest in the BEFM project.

# Compiled from Markdown document originally written 2024-02-28.
# 'AllClass_Embedding_EngOnly_TriQuad_Dev.Rmd'

# It takes the global phenotype object as input.

# It saves to disc a small subset phenotype object, with a custom embedding, and three contaminating macrophages removed.

# Set WD:
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics")

# Input object(s):
load('Saved Objects/downsampled.imputed.Robj')

# Output object(s):
# save(sub.fine,file = 'transcriptome.TriQuad.customEmbed.fineFilter.proofRead.2024-03-15.Robj')
# save(ctc.sub.fine,file = 'connectome.TriQuad.finefilter.proofRead.2024-03-15.Robj')

# Require packages
require(ggplot2)
require(Seurat)
require(RColorBrewer)
require(cowplot)
require(dplyr)
library(scales)
library(viridis)

# Load color palette
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Allie's Object Filtration/Final Objects/all.color_palettes_2.R")
color_pals

# Organize object metadata
downsampled$class <- factor(downsampled$class,levels = c('Epithelium','Endothelium','Mesenchyme','Immune'))

# Limit just to samples of interest (create object named "sub")
Idents(downsampled) <- downsampled$Dataset2
table(Idents(downsampled))
sub <- subset(downsampled,idents=c('Tri_L','Quad_E'),invert=F)
table(Idents(sub))
table(sub$orig.ident)

# Excellent. Next, we identify those genes that segregate classes, trying to get the calculation to not take very seriously those features that separate batches:
table(sub$class)
Idents(sub) <- sub$class
DefaultAssay(sub) <- 'RNA'

# Old approach (not great here...)
# mark <- FindAllMarkers(sub,min.pct = 0.75,logfc.threshold = 0.5,min.cells.feature = 1000,return.thresh = 0.0000000001) # These high threshold parameters are designed to (mostly) force it to only consider features common to multiple samples
# mark$ratio <- mark$pct.1/mark$pct.2
# mark$power <- mark$ratio*mark$avg_log2FC
# goi <- unique(mark$gene)
# goi

# New experiment for allie
Idents(sub) <- sub$orig.ident
bef12 <- subset(sub,idents = 'BEF12')
bef14 <- subset(sub,idents = 'BEF14')
bef15 <- subset(sub,idents = 'BEF15')
befm1 <- subset(sub,idents = 'BEFM1')
befm2 <- subset(sub,idents = 'BEFM2')
Idents(bef12) <- bef12$class
Idents(bef14) <- bef14$class
Idents(bef15) <- bef15$class
Idents(befm1) <- befm1$class
Idents(befm2) <- befm2$class
table(Idents(befm1))
table(Idents(befm2))
mark.bef12 <- FindAllMarkers(bef12,min.pct = 0.75,logfc.threshold = 1) # experimental
mark.bef14 <- FindAllMarkers(bef14,min.pct = 0.75,logfc.threshold = 1) # experimental
mark.bef15 <- FindAllMarkers(bef15,min.pct = 0.75,logfc.threshold = 1) # experimental
mark.befm1 <- FindAllMarkers(befm1,min.pct = 0.75,logfc.threshold = 1) # experimental
mark.befm2 <- FindAllMarkers(befm2,min.pct = 0.75,logfc.threshold = 1) # experimental

# Concatenate
temp <- c(unique(mark.bef12$gene),unique(mark.bef14$gene),unique(mark.bef15$gene),unique(mark.befm1$gene),unique(mark.befm2$gene))
# Filter
goi <- names(table(temp)>1) # must be found in at least two sample to be used for embedding

# Inspect goi
goi
length(goi)

# Seems reasonable as a starting point. Let's see what a default UMAP looks like using these goi:
sub <- ScaleData(sub,features = goi)
sub <- RunUMAP(sub,features = goi,seed.use =104657) # playing with seed.use here will modify how the clusters are laid out on the page relative to one another
p1 <- DimPlot(sub,group.by = 'class',cols = color_pals$class_colors)
p2 <- DimPlot(sub,group.by = 'orig.ident',cols = color_pals$sample_colors,shuffle = T)
plot_grid(p1,p2,nrow=1)

# Not bad, but serious mesenchymal-specific separation (real signal). We don't care for this though, here. Let's exclude the variance that separates the mesenchyme samples from one another:
Idents(sub) <- sub$class
mes.sub <- subset(sub,idents = 'Mesenchyme')
Idents(mes.sub) <- mes.sub$orig.ident
mark.to.remove <- FindAllMarkers(mes.sub,min.pct = 0.5,logfc.threshold = 0.5)
goi.2 <- goi[!(goi %in% mark.to.remove$gene)]

# Add some features that I want in there (help endothelium to separate from mesenchyme)
goi.2 <- c(goi.2,'Cdh5') # Must include this major lineage marker, I am demanding it
mark.end <- FindMarkers(sub,ident.1 = 'Endothelium',min.pct = 0.5,logfc.threshold = 0.5) # Find some others in the data
mark.end$gene <- rownames(mark.end)
to.add <- mark.end$gene[!(mark.end$gene %in% mark.to.remove$gene)]
goi.3 <- unique(c(goi.2,to.add))

# Inspect goi.3
goi.3
length(goi)

# Look at embedding with goi.3
sub <- ScaleData(sub,features = goi.3)
sub <- RunUMAP(sub,features = goi.3,seed.use = 42,min.dist = 0.1) # playing with seed.use here will modify how the clusters are laid out on the page relative to one another
p1 <- DimPlot(sub,group.by = 'class',cols = color_pals$class_colors)
p2 <- DimPlot(sub,group.by = 'orig.ident',cols = color_pals$sample_colors,shuffle = T)
plot_grid(p1,p2,nrow=1)

# Looks great!

# Clean-up: # remove the 2(!) contaminating macrophages in tri-culture samples (make "sub.fine" object)
table(Idents(bef12))
barcodes.remove <- rownames(bef12@meta.data[which(bef12@meta.data$class == 'Immune'),])
sub.fine <- subset(sub,cells = barcodes.remove,invert=T)
table(Idents(sub.fine),sub.fine$orig.ident)

# Clean up: # set up metadata slot 'Condition' and order for publication
sub.fine$Condition <- sub.fine$Dataset2
sub.fine$Condition <- factor(sub.fine$Condition,levels = c('Tri_L','Quad_E'))

# plot for posterity
plot_grid(p1,p2,nrow=1)

## Save to disc 
save(sub.fine,file = 'transcriptome.TriQuad.customEmbed.fineFilter.proofRead.2024-03-15.Robj')


#### Make connectivity sub object that matches the above phenotype object

# Load global connectomics from disc
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Saved Objects/global.connectomics.2023-11-11.Robj")

# Pull the NICHES data made using imputed slot
ctc <- global.connectomics$CellToCell$alra

# Subset to just the 5 samples of interest
Idents(ctc) <- ctc$Dataset2.Sending
table(Idents(ctc))
ctc.sub <- subset(ctc,idents = c('Tri_L','Quad_E'),invert=F)
table(Idents(ctc.sub))
ctc.sub <- ScaleData(ctc.sub)
table(ctc.sub$orig.ident)

## Fine tuning (removing the two macrophages in BEF12!)
Idents(ctc.sub) <- ctc.sub$orig.ident
bef12.ctc <- subset(ctc.sub,idents = 'BEF12')
barcodes.remove.ctc <- rownames(bef12.ctc@meta.data[which(bef12.ctc@meta.data$class.Sending == 'Immune'),])
barcodes.remove.ctc <- c(barcodes.remove.ctc,
                         rownames(bef12.ctc@meta.data[which(bef12.ctc@meta.data$class.Receiving == 'Immune'),]))
ctc.sub.fine <- subset(ctc.sub,cells = barcodes.remove.ctc,invert=T)

# Inspect ctc.sub.fine
table(Idents(ctc.sub.fine),ctc.sub.fine$class.Sending)
table(Idents(ctc.sub.fine),ctc.sub.fine$class.Receiving)

# Looks good!

## Save to disc
save(ctc.sub.fine,file = 'connectome.TriQuad.finefilter.proofRead.2024-03-15.Robj')
