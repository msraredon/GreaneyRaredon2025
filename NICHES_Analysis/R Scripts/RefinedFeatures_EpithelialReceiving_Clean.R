# 2024-03-21
# MSBR

# This script assembles and cleans the work to date on Epithelial Receiving Analysis for BEFM.

# Copies and clean exploratory scripting originally from:
## 'EpiRecFeatures.Rmd' (heatmaps and circuit plots) and 
## 'AllClassEmbedding_Clean.Rmd' (network plots)

# Set WD:
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics")

# Global parameters
update_geom_defaults("violin", aes(linewidth = 0)) # Removes outlines on violin plots
min.connectivity = 0.00001

# Local functions
source('/Users/msbr/GitHub/NICHESMethods/ToolBuilding/CustomHeatmapOnly3.R')
source('/Users/msbr/GitHub/NICHESMethods/ToolBuilding/NetworkPlot.R')

# Input object(s):
## 1.
## Custom embedding of phenotype data for only 5 samples of interest, cleaned to remove contaminating macrophages in tri-culture conditions, and with organized metadata
## Created by 'R Scripts/AllClassEmbedding_Clean.R'
## the object is called "sub.fine" when loaded
load('Saved Objects/transcriptome.TriQuad.customEmbed.fineFilter.proofRead.2024-03-21.Robj') 

## 2.
## Connectivity data that exactly matches the above phenotype object, for plotting
## Created by 'R Scripts/AllClassEmbedding_Clean.R'
## the object is called "ctc.sub.fine" when loaded
load('Saved Objects/connectome.TriQuad.finefilter.proofRead.2024-03-21.Robj')

## 3.
## Epithelial Receiving Connectivity Object, clustered, annotated, organized for plotting
## Created by 'R Scripts/EpithelialReceivingObject_Clean.R
load("Saved Objects/epi.rec.2024-03-21.Robj")

## 4.
## Color palette from Allie
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Allie's Object Filtration/Final Objects/all.color_palettes_2.R")

## 5. Global connectomics to make full color palette (inefficient but works)
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Saved Objects/global.connectomics.2023-11-11.Robj")

## 6. Ultra marker list for epithelium
## Created by GlobalStatisticalAnalysis.R
load("mark.epi.ultra.Robj")
load("mark.end.ultra.Robj")
load("mark.mes.ultra.Robj")
load("mark.imm.ultra.Robj")

# Make exclusive
mark.epi.final <- mark.epi.ultra[!(mark.epi.ultra$gene %in% mark.end.ultra$gene)&
                                   !(mark.epi.ultra$gene %in% mark.mes.ultra$gene)&
                                   !(mark.epi.ultra$gene %in% mark.imm.ultra$gene),]

# Require packages
require(Seurat)
require(NICHES)
require(dplyr)
require(ggthemes)
require(ggplot2)
require(cowplot)
require(knitr)

# Add a temporary color palette allowing SendingType visualization
color_pals$type_colors <- c('#A40606','#9CFFFA','#B0DB43','#9C528B','#2F6690',
                            '#946846','#F1C40F','green','#0F0326','#E65F5C','#14591D','#726DA8',
                            'yellow','purple','blue','red','orange','darkgrey','magenta')
names(color_pals$type_colors) <- unique(global.connectomics$CellToCell$alra$SendingType)

# Modify color palette to have a slot that aligns exactly with metadata handles we will use below
color_pals$class <- color_pals$class_colors

## Update to v5
epi.rec <- UpdateSeuratObject(epi.rec)
epi.rec

# Inspect
table(epi.rec$Archetype,epi.rec$seurat_clusters)

# Clean up: re-order SendingType and ReceivingType for plotting
names(table(epi.rec$SendingType))
epi.rec$SendingType <- factor(epi.rec$SendingType,
                              levels=c('Basal-like','Transitional','BASC','ATI-like','Cycling_Epi',
                                       'Endo_Progenitor','Microvascular','Cycling_Microvascular','Lymphatic','Cycling_Lymphatic',
                                       'Aberrant_Fibroblast','Remodeling_Fibroblast','Alveolar_Fibroblast','Pericytes','Cycling_Mes',
                                       'Alveolar_Macrophage_Naive','Alveolar_Macrophage_Activated','Cycling_Alveolar_Macrophage','Interstitial_Macrophage'))
table(epi.rec$SendingType)
table(epi.rec$ReceivingType)
epi.rec$ReceivingType <- factor(epi.rec$ReceivingType,
                                levels=c('Basal-like','Transitional','BASC','ATI-like','Cycling_Epi',
                                         'Endo_Progenitor','Microvascular','Cycling_Microvascular','Lymphatic','Cycling_Lymphatic',
                                         'Aberrant_Fibroblast','Remodeling_Fibroblast','Alveolar_Fibroblast','Pericytes','Cycling_Mes',
                                         'Alveolar_Macrophage_Naive','Alveolar_Macrophage_Activated','Cycling_Alveolar_Macrophage','Interstitial_Macrophage'))
table(epi.rec$ReceivingType)

# Clean up: check work
sum(is.na(epi.rec$SendingType))
sum(is.na(epi.rec$ReceivingType))

## All features as FeaturePlots [approved]
archetypes <- as.character(unique(mark.epi.final$cluster))
for(j in 1:length(archetypes)){
  print(j)
  sub.markers <- mark.epi.final[mark.epi.final$cluster==archetypes[j],]
 for(i in 1:nrow(sub.markers)){
   print(i)
  png(file = paste(archetypes[j],sub.markers$gene[i],'FeaturePlot_EpitheliumReceiving_2024-03-21.png'),width = 8,height = 8,units = 'in',res=300)
  print(FeaturePlot(epi.rec,
              features = sub.markers$gene[i],
              order = T,
              pt.size=1,
              label=T,
              label.size = 4,
              slot = 'data')+
    scale_colour_gradientn(colours = c('lightgray','#348AA7','#FAA916','#EF233C'),name='Connectivity')+
    NoAxes())
  dev.off()
 }
}

## Feature Violin Combo Plots for Top Sixteen markers in each archetype
for(j in 1:length(archetypes)){
  print(j)
  sub.markers <- mark.epi.final[mark.epi.final$cluster==archetypes[j],]
  sub.markers <- sub.markers[sub.markers$avg_log2FC>1,] # At least one fold change difference
  sub.markers <- sub.markers %>% top_n(16,ratio) # top 16 by ratio

  # Define goi
  goi <- data.frame(sub.markers)
  goi <- goi$gene
  
  # Feature Plots
  plot.list <- list()
  for(i in 1:length(goi)){
    plot.list[[i]] = FeaturePlot(epi.rec,
                                 features = goi[i],
                                 order = T,
                                 pt.size=1,
                                 label=T,
                                 label.size = 4,
                                 slot = 'data')+
      scale_colour_gradientn(colours = c('lightgray','#348AA7','#FAA916','#EF233C'),name='Connectivity')+
      NoAxes()
  }
  total.plot.features <- cowplot::plot_grid(plotlist = plot.list)
  
  # Export (features)
  png(filename = paste(archetypes[j],'Top Sixteen Features 2024-03-17.png'),width = 28,height = 22,units = 'in',res=600)
  print(total.plot.features)
  dev.off()
  
  # Violin Plots (by sample)
  plot.list.vln <- list()
  for(i in 1:length(goi)){
    plot.list.vln[[i]] = VlnPlot(epi.rec,
                                 features = goi[i],
                                 group.by = 'Sample.Sending',
                                 cols = color_pals$sample_colors,
                                 pt.size = 0,
                                 adjust=2)+ylab('Connectivity')+NoLegend()+ylim(min.connectivity,NA)
    
  }
  total.plot.vlns <- cowplot::plot_grid(plotlist = plot.list.vln,ncol = 1)
  
  # Violin Plots (by sending celltype)
  plot.list.vln.sending <- list()
  for(i in 1:length(goi)){
    plot.list.vln.sending[[i]] = VlnPlot(epi.rec,
                                         features = goi[i],
                                         group.by = 'SendingType',
                                         cols = color_pals$type_colors,
                                         pt.size = 0,
                                         adjust=2)+ylab('Connectivity')+NoLegend()+ylim(min.connectivity,NA)
    
  }
  total.plot.vlns.sending <- cowplot::plot_grid(plotlist = plot.list.vln.sending,ncol = 1)
  
  # Violin Plots (by receiving celltype)
  plot.list.vln.receiving <- list()
  for(i in 1:length(goi)){
    plot.list.vln.receiving[[i]] = VlnPlot(epi.rec,
                                           features = goi[i],
                                           group.by = 'ReceivingType',
                                           cols = color_pals$type_colors,
                                           pt.size = 0,
                                           adjust=2)+ylab('Connectivity')+NoLegend()+ylim(min.connectivity,NA)
    
  }
  total.plot.vlns.receiving <- cowplot::plot_grid(plotlist = plot.list.vln.receiving,ncol = 1)
  
  
  # Remake feature plots vertical
  total.plot.features.vertical <- cowplot::plot_grid(plotlist = plot.list,ncol=1)
  
  # Export (features & violins)
  total.plot <- cowplot::plot_grid(total.plot.features.vertical,
                                   total.plot.vlns,
                                   total.plot.vlns.sending,
                                   total.plot.vlns.receiving,
                                   ncol = 4,
                                   rel_widths = c(1,1,2,1))
  png(filename = paste(archetypes[j],'Top Sixteen Features with Violins 2024-03-17.png'),width = 20,height = 72,units = 'in',res=600)
  print(total.plot)
  dev.off()

# j Loop end
}
  
 

## Heatmap of thresholded features (yuck)
goi <- c()
row.labels.epi <- c()
for(j in 1:length(archetypes)){
  print(j)
  sub.markers <- mark.epi.final[mark.epi.final$cluster==archetypes[j],]
  #sub.markers <- sub.markers[sub.markers$avg_log2FC>1,] # At least one fold change difference
  sub.markers <- sub.markers %>% top_n(30,ratio) # top 16 by ratio
  
  # define labels
  to.label <- sub.markers %>% top_n(5,ratio)
  to.label <- data.frame(to.label)
  row.labels.epi <- c(row.labels.epi,to.label$gene)
  
  # Define goi
  goi.temp <- data.frame(sub.markers)
  goi.temp <- goi.temp$gene
  # concatenate
  goi <- c(goi,goi.temp)
}
# Make unique
goi <- unique(goi)

png('EpiRec.Heatmap.2024-03-21.png',width=8,height=12,units='in',res=300)
CustomHeatmapOnly3(epi.rec,
                   data.type = 'CellToCell',
                   primary = 'Archetype' ,
                   secondary = 'class.Sending' ,
                   tertiary = 'orig.ident' ,
                   #quarternary = 'orig.ident' ,
                   primary.cols = NULL ,
                   secondary.cols = color_pals$class_colors,
                   tertiary.cols = color_pals$sample_colors,#color_pals$dataset_colors[5:6],
                   #quarternary.cols = color_pals$sample_colors,
                   features = goi,
                   labels = c('Signaling Archetype','Sending Cell Class','Sample'),
                   selected.row.anotations=row.labels.epi,
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.2)
dev.off()

## Handcrafted Heatmap


# End of script


##### move to other scripts
## Network Plots for all features

# Define features to plot
all.features <- goi.epi

# PreProcess: add Condition handle to ctc.sub.fine
ctc.sub.fine$Condition <- ctc.sub.fine$Dataset2.Sending
ctc.sub.fine$Condition <- factor(ctc.sub.fine$Condition,
                                 levels = c('Tri_L','Quad_E'))

# PreProcess: for each feature, find the maximum edge value
max.values <- qlcMatrix::rowMax(as.matrix(ctc.sub.fine@assays$CellToCell@counts[all.features,]))
max.values <- data.frame(feature = all.features,
                         max.value = as.matrix(max.values))
max.values$log.transformed <- log1p(max.values$max.value)

# Set connectivity threshold for edge pruning
min.connectivity.thresh <- 0.25

# Epithelium NetworkPlots for all features
for(i in 1:length(all.features)){
  feature = all.features[i]
  png(file = paste(feature,'NetworkPlot',min.connectivity.thresh,'2024-03-21.png',sep='_'),width = 6,height = 5,units = 'in',res=600)
  print(NetworkPlot(transcriptome.object = sub.fine,
                    connectome.object = ctc.sub.fine,
                    mechanism.of.interest = feature,
                    legends.to.plot = 'class',
                    legend.palettes = color_pals,
                    split.by = 'Condition',
                    min.connectivity.thresh = min.connectivity.thresh,
                    connectivity.color.min = 0,
                    connectivity.color.max = max.values[max.values$feature==feature,]$log.transformed,
                    line.thickness = 0.1,
                    black.points = T))
  dev.off()
}

# Circuit Plots for all features

