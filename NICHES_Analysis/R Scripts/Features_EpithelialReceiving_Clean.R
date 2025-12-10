# 2024-03-15
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

# Local functions3
source('/Users/msbr/GitHub/NICHESMethods/ToolBuilding/CustomHeatmapOnly3.R')
source('/Users/msbr/GitHub/NICHESMethods/ToolBuilding/NetworkPlot.R')

# Input object(s):
## 1.
## Custom embedding of phenotype data for only 5 samples of interest, cleaned to remove contaminating macrophages in tri-culture conditions, and with organized metadata
## Created by 'R Scripts/AllClassEmbedding_Clean.R'
## the object is called "sub.fine" when loaded
load('Saved Objects/transcriptome.TriQuad.customEmbed.fineFilter.proofRead.2024-03-15.Robj') 

## 2.
## Connectivity data that exactly matches the above phenotype object, for plotting
## Created by 'R Scripts/AllClassEmbedding_Clean.R'
## the object is called "ctc.sub.fine" when loaded
load('Saved Objects/connectome.TriQuad.finefilter.proofRead.2024-03-15.Robj')

## 3.
## Epithelial Receiving Connectivity Object, clustered, annotated, organized for plotting
## Marker list associated with the above
## Created by 'R Scripts/EpithelialReceivingObject_Clean.R
load("Saved Objects/epi.rec.2024-03-15.Robj")
load("Saved Objects/mark.epi.rec.2024-03-15.Robj")

## 4.
## Color palette from Allie
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Allie's Object Filtration/Final Objects/all.color_palettes_2.R")
## Global connectomics to make full color palette (inefficient but works)
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Saved Objects/global.connectomics.2023-11-11.Robj")

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

## Cluster 0 (Mes-Epi Common)

goi.0 <-  c('Dlk1—Notch1','Bmp5—Bmpr1b','Rspo3—Lrp6',
            'Angpt1—Itgb1','Bmp5—Acvr2a','Thbs2—Itgb1',
            'Bmp5—Bmpr1a','Thbs2—Cd47','Ptn—Sdc1','Mst1—Mst1r',
            'Fgf7—Fgfr1','Igf2—Igf2r','Bmp5—Bmpr1b','Ptn—Sdc1',
            'Fgf10—Fgfr2','Dkk2—Lrp6','Bmp5—Bmpr2')
# Feature Plots
plot.list <- list()
for(i in 1:length(goi.0)){
  plot.list[[i]] = FeaturePlot(epi.rec,
                               features = goi.0[i],
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
png(filename = paste("Mes-Epi Common",'Selected Features 2024-03-17.png'),width = 28,height = 22,units = 'in',res=600)
total.plot.features
dev.off()

# Violin Plots (by sample)
plot.list.vln <- list()
for(i in 1:length(goi.0)){
  plot.list.vln[[i]] = VlnPlot(epi.rec,
                               features = goi.0[i],
                               group.by = 'Sample.Sending',
                               cols = color_pals$sample_colors,
                               pt.size = 0,
                               adjust=2)+ylab('Connectivity')+NoLegend()+ylim(min.connectivity,NA)

}
total.plot.vlns <- cowplot::plot_grid(plotlist = plot.list.vln,ncol = 1)

# Violin Plots (by sending celltype)
plot.list.vln.sending <- list()
for(i in 1:length(goi.0)){
  plot.list.vln.sending[[i]] = VlnPlot(epi.rec,
                               features = goi.0[i],
                               group.by = 'SendingType',
                               cols = color_pals$type_colors,
                               pt.size = 0,
                               adjust=2)+ylab('Connectivity')+NoLegend()+ylim(min.connectivity,NA)
  
}
total.plot.vlns.sending <- cowplot::plot_grid(plotlist = plot.list.vln.sending,ncol = 1)

# Violin Plots (by receiving celltype)
plot.list.vln.receiving <- list()
for(i in 1:length(goi.0)){
  plot.list.vln.receiving[[i]] = VlnPlot(epi.rec,
                                       features = goi.0[i],
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
png(filename = paste("Mes-Epi Common",'Selected Feature with Violins 2024-03-17.png'),width = 20,height = 72,units = 'in',res=600)
total.plot
dev.off()

## Cluster 1 (Epi-Autocrine Common)

goi.1 <-  c('Efna3—Epha1','Efna3—Epha2','Efna3—Epha4','Ntf4—Ngfr','Il10—Il10rb',
            'Fgf9—Fgfr1','Dlk2—Notch1','Wnt3a—Fzd6','Wnt3a—Lrp6','Wnt3a—Ryk',
            'Fgf9—Fgfr2','Bmp7—Eng','Bmp7—Bmpr1b','Wnt7b—Fzd1','Bmp7—Bmpr2',
            'Areg—Erbb3','Ereg—Erbb3','Areg—Egfr')
# Feature Plots
plot.list <- list()
for(i in 1:length(goi.1)){
  plot.list[[i]] = FeaturePlot(epi.rec,
                               features = goi.1[i],
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
png(filename = paste("Epi-Auto Common",'Selected Features 2024-03-17.png'),width = 28,height = 22,units = 'in',res=600)
total.plot.features
dev.off()

# Violin Plots (by sample)
plot.list.vln <- list()
for(i in 1:length(goi.1)){
  plot.list.vln[[i]] = VlnPlot(epi.rec,
                               features = goi.1[i],
                               group.by = 'Sample.Sending',
                               cols = color_pals$sample_colors,
                               pt.size = 0,
                               adjust=2)+ylab('Connectivity')+NoLegend()+ylim(min.connectivity,NA)
  
}
total.plot.vlns <- cowplot::plot_grid(plotlist = plot.list.vln,ncol = 1)

# Violin Plots (by sending celltype)
plot.list.vln.sending <- list()
for(i in 1:length(goi.1)){
  plot.list.vln.sending[[i]] = VlnPlot(epi.rec,
                                       features = goi.1[i],
                                       group.by = 'SendingType',
                                       cols = color_pals$type_colors,
                                       pt.size = 0,
                                       adjust=2)+ylab('Connectivity')+NoLegend()+ylim(min.connectivity,NA)
  
}
total.plot.vlns.sending <- cowplot::plot_grid(plotlist = plot.list.vln.sending,ncol = 1)

# Violin Plots (by receiving celltype)
plot.list.vln.receiving <- list()
for(i in 1:length(goi.1)){
  plot.list.vln.receiving[[i]] = VlnPlot(epi.rec,
                                         features = goi.1[i],
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
png(filename = paste("Epi-Auto Common",'Selected Feature with Violins 2024-03-17.png'),width = 20,height = 72,units = 'in',res=600)
total.plot
dev.off()


## Cluster 2 (End-Epi Common)

goi.2 <-  c('Dll1—Notch2','Dll1—Notch1','Jag1—Notch4','Dll1—Notch3','Bmp3—Bmpr1b',
            'Dll1—Notch4','Gas6—Axl','Tnfsf12—Tnfrsf25','Jag1—Notch3','Bmp6—Bmpr2',
            'Bmp3—Bmpr2','Vegfc—Nrp2','Vip—Sctr','Tnfsf4—Tnfrsf4','Dll4—Notch1',
            'Dll4—Notch2','Pdgfb—Pdgfrb','Angpt2—Tie1','Dll4—Notch3','Dll4—Notch4',
            'Dhh—Ptch1','Sema6d—Plxna1')
# Feature Plots
plot.list <- list()
for(i in 1:length(goi.2)){
  plot.list[[i]] = FeaturePlot(epi.rec,
                               features = goi.2[i],
                               order = T,
                               pt.size=1,
                               label=T,
                               label.size = 4,
                               slot='data')+
    scale_colour_gradientn(colours = c('lightgray','#348AA7','#FAA916','#EF233C'),name='Connectivity')+
    NoAxes()
}
total.plot.features <- cowplot::plot_grid(plotlist = plot.list)

# Export (features)
png(filename = paste("Endo-Epi Common",'Selected Features 2024-03-17.png'),width = 28,height = 22,units = 'in',res=600)
total.plot.features
dev.off()

# Violin Plots (by sample)
plot.list.vln <- list()
for(i in 1:length(goi.2)){
  plot.list.vln[[i]] = VlnPlot(epi.rec,
                               features = goi.2[i],
                               group.by = 'Sample.Sending',
                               cols = color_pals$sample_colors,
                               pt.size = 0,
                               adjust=2)+ylab('Connectivity')+NoLegend()+ylim(min.connectivity,NA)
  
}
total.plot.vlns <- cowplot::plot_grid(plotlist = plot.list.vln,ncol = 1)

# Violin Plots (by sending celltype)
plot.list.vln.sending <- list()
for(i in 1:length(goi.2)){
  plot.list.vln.sending[[i]] = VlnPlot(epi.rec,
                                       features = goi.2[i],
                                       group.by = 'SendingType',
                                       cols = color_pals$type_colors,
                                       pt.size = 0,
                                       adjust=2)+ylab('Connectivity')+NoLegend()+ylim(min.connectivity,NA)
  
}
total.plot.vlns.sending <- cowplot::plot_grid(plotlist = plot.list.vln.sending,ncol = 1)

# Violin Plots (by receiving celltype)
plot.list.vln.receiving <- list()
for(i in 1:length(goi.2)){
  plot.list.vln.receiving[[i]] = VlnPlot(epi.rec,
                                         features = goi.2[i],
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
png(filename = paste("Endo-Epi Common",'Selected Feature with Violins 2024-03-17.png'),width = 20,height = 72,units = 'in',res=600)
total.plot
dev.off()


## Cluster 3 Epi-Auto Regenerative

goi.3 <- c('Tnf—Tnfrsf1b','Sct—Sctr','Wnt7a—Lrp6','Npnt—Itgb1','Ereg—Egfr',
           'Wnt7b—Lrp5','Hbegf—Egfr','Lamb3—Itga3','Npnt—Itgb1','Tnf—Tnfrsf21',
           'Tnf—Tnfrsf21','Tnf—Ltbr','Lamb3—Itga2','Mmp7—Cd44','Lamb3—Itga2','Btc—Egfr')
# Feature Plots
plot.list <- list()
for(i in 1:length(goi.3)){
  plot.list[[i]] = FeaturePlot(epi.rec,
                               features = goi.3[i],
                               order = T,
                               pt.size=1,
                               label=T,
                               label.size = 4,
                               slot='data')+
    scale_colour_gradientn(colours = c('lightgray','#348AA7','#FAA916','#EF233C'),name='Connectivity')+
    NoAxes()
}
total.plot.features <- cowplot::plot_grid(plotlist = plot.list)

# Export (features)
png(filename = paste("Epi-Auto Regenerative",'Selected Features 2024-03-17.png'),width = 28,height = 22,units = 'in',res=600)
total.plot.features
dev.off()

# Violin Plots (by sample)
plot.list.vln <- list()
for(i in 1:length(goi.3)){
  plot.list.vln[[i]] = VlnPlot(epi.rec,
                               features = goi.3[i],
                               group.by = 'Sample.Sending',
                               cols = color_pals$sample_colors,
                               pt.size = 0,
                               adjust=2)+ylab('Connectivity')+NoLegend()+ylim(min.connectivity,NA)
  
}
total.plot.vlns <- cowplot::plot_grid(plotlist = plot.list.vln,ncol = 1)

# Violin Plots (by sending celltype)
plot.list.vln.sending <- list()
for(i in 1:length(goi.3)){
  plot.list.vln.sending[[i]] = VlnPlot(epi.rec,
                                       features = goi.3[i],
                                       group.by = 'SendingType',
                                       cols = color_pals$type_colors,
                                       pt.size = 0,
                                       adjust=2)+ylab('Connectivity')+NoLegend()+ylim(min.connectivity,NA)
  
}
total.plot.vlns.sending <- cowplot::plot_grid(plotlist = plot.list.vln.sending,ncol = 1)

# Violin Plots (by receiving celltype)
plot.list.vln.receiving <- list()
for(i in 1:length(goi.3)){
  plot.list.vln.receiving[[i]] = VlnPlot(epi.rec,
                                         features = goi.3[i],
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
png(filename = paste("Epi-Auto Regenerative",'Selected Feature with Violins 2024-03-17.png'),width = 20,height = 72,units = 'in',res=600)
total.plot
dev.off()


## Cluster 4 Mes-Epi Regenerative

goi.4 <- c('Fgf1—Fgfr3','Fgf1—Egfr','Fgf1—Fgfr2','Fgf2—Fgfr3','Fgf12—Fgfr3','Lama1—Itga2',
           'Csf2—Csf2ra','Ngf—Maged1','Wnt16—Fzd6','Fgf7—Fgfr3','Fn1—Itgb6','Vegfa—Ephb2',
           'Wnt5a—Lrp5','Col4a1—Itga2','Tnc—Egfr','Tnc—Egfr')
# Feature Plots
plot.list <- list()
for(i in 1:length(goi.4)){
  plot.list[[i]] = FeaturePlot(epi.rec,
                               features = goi.4[i],
                               order = T,
                               pt.size=1,
                               label=T,
                               label.size = 4,
                               slot='data')+
    scale_colour_gradientn(colours = c('lightgray','#348AA7','#FAA916','#EF233C'),name='Connectivity')+
    NoAxes()
}
total.plot.features <- cowplot::plot_grid(plotlist = plot.list)

# Export (features)
png(filename = paste('Mes-Epi Regenerative','Selected Features 2024-03-17.png'),width = 28,height = 22,units = 'in',res=600)
total.plot.features
dev.off()

# Violin Plots (by sample)
plot.list.vln <- list()
for(i in 1:length(goi.4)){
  plot.list.vln[[i]] = VlnPlot(epi.rec,
                               features = goi.4[i],
                               group.by = 'Sample.Sending',
                               cols = color_pals$sample_colors,
                               pt.size = 0,
                               adjust=2)+ylab('Connectivity')+NoLegend()+ylim(min.connectivity,NA)
  
}
total.plot.vlns <- cowplot::plot_grid(plotlist = plot.list.vln,ncol = 1)

# Violin Plots (by sending celltype)
plot.list.vln.sending <- list()
for(i in 1:length(goi.4)){
  plot.list.vln.sending[[i]] = VlnPlot(epi.rec,
                                       features = goi.4[i],
                                       group.by = 'SendingType',
                                       cols = color_pals$type_colors,
                                       pt.size = 0,
                                       adjust=2)+ylab('Connectivity')+NoLegend()+ylim(min.connectivity,NA)
  
}
total.plot.vlns.sending <- cowplot::plot_grid(plotlist = plot.list.vln.sending,ncol = 1)

# Violin Plots (by receiving celltype)
plot.list.vln.receiving <- list()
for(i in 1:length(goi.4)){
  plot.list.vln.receiving[[i]] = VlnPlot(epi.rec,
                                         features = goi.4[i],
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
png(filename = paste("Mes-Epi Regenerative",'Selected Feature with Violins 2024-03-17.png'),width = 20,height = 72,units = 'in',res=600)
total.plot
dev.off()



## Cluster 5 Macrophage-Epi

goi.5 <- c('Osm—Osmr','Il1b—Il1r1','Il1b—Adrb2','Tnfsf15—Tnfrsf25','Ebi3—Il6st',
           'Il18—Il18r1','Camp—Igf1r','Tnfsf13—Tnfrsf1a','Nrg4—Egfr','Nrg4—Erbb2',
           'Efna2—Epha1','Wnt5a—Ryk','Amh—Egfr','Hgf—Met','Igf1—Igf1r','Cxcl2—Cxcr2')
# Feature Plots
plot.list <- list()
for(i in 1:length(goi.5)){
  plot.list[[i]] = FeaturePlot(epi.rec,
                               features = goi.5[i],
                               order = T,
                               pt.size=1,
                               label=T,
                               label.size = 4,
                               slot='data')+
    scale_colour_gradientn(colours = c('lightgray','#348AA7','#FAA916','#EF233C'),name='Connectivity')+
    NoAxes()
}
total.plot.features <- cowplot::plot_grid(plotlist = plot.list)

# Export (features)
png(filename = paste('Macrophage-Epi','Selected Features 2024-03-17.png'),width = 28,height = 22,units = 'in',res=600)
total.plot.features
dev.off()

# Violin Plots (by sample)
plot.list.vln <- list()
for(i in 1:length(goi.5)){
  plot.list.vln[[i]] = VlnPlot(epi.rec,
                               features = goi.5[i],
                               group.by = 'Sample.Sending',
                               cols = color_pals$sample_colors,
                               pt.size = 0,
                               adjust=2)+ylab('Connectivity')+NoLegend()+ylim(min.connectivity,NA)
  
}
total.plot.vlns <- cowplot::plot_grid(plotlist = plot.list.vln,ncol = 1)

# Violin Plots (by sending celltype)
plot.list.vln.sending <- list()
for(i in 1:length(goi.5)){
  plot.list.vln.sending[[i]] = VlnPlot(epi.rec,
                                       features = goi.5[i],
                                       group.by = 'SendingType',
                                       cols = color_pals$type_colors,
                                       pt.size = 0,
                                       adjust=2)+ylab('Connectivity')+NoLegend()+ylim(min.connectivity,NA)
  
}
total.plot.vlns.sending <- cowplot::plot_grid(plotlist = plot.list.vln.sending,ncol = 1)

# Violin Plots (by receiving celltype)
plot.list.vln.receiving <- list()
for(i in 1:length(goi.5)){
  plot.list.vln.receiving[[i]] = VlnPlot(epi.rec,
                                         features = goi.5[i],
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
png(filename = paste("Macrophage-Epi",'Selected Feature with Violins 2024-03-17.png'),width = 20,height = 72,units = 'in',res=600)
total.plot
dev.off()

## Heatmap for Janelia poster, updated for publication

goi.0 <-  c('Angpt1—Itgb1','Bmp5—Acvr2a','Thbs2—Itgb1',
            'Fgf10—Fgfr2','Bmp5—Bmpr1a','Thbs2—Cd47','Dlk1—Notch1','Bmp5—Bmpr1b','Rspo3—Lrp6',
            'Angpt1—Itgb1','Bmp5—Acvr2a','Thbs2—Itgb1',
            'Bmp5—Bmpr1a','Thbs2—Cd47','Ptn—Sdc1','Mst1—Mst1r','Fgf7—Fgfr1','Igf2—Igf2r',
            'Bmp5—Bmpr1b','Ptn—Sdc1','Dkk2—Lrp6','Bmp5—Bmpr2')

goi.1 <-  c('Efna3—Epha1','Efna3—Epha2','Efna3—Epha4','Ntf4—Ngfr','Il10—Il10rb',
            'Fgf9—Fgfr1','Dlk2—Notch1','Wnt3a—Fzd6','Wnt3a—Lrp6','Wnt3a—Ryk','Wnt7b—Fzd1',
            'Fgf9—Fgfr2','Bmp7—Eng','Bmp7—Bmpr1b','Bmp7—Bmpr2','Areg—Erbb3','Ereg—Erbb3','Areg—Egfr')

goi.2 <-  c('Jag1—Notch4','Dll1—Notch3','Bmp3—Bmpr1b','Dll1—Notch4','Gas6—Axl',
            'Tnfsf12—Tnfrsf25','Jag1—Notch3','Bmp6—Bmpr2','Bmp3—Bmpr2','Vegfc—Nrp2',
            'Vip—Sctr','Dll1—Notch2','Dll1—Notch1','Dhh—Ptch1','Tnfsf4—Tnfrsf4',
            'Dll4—Notch1','Dll4—Notch2','Pdgfb—Pdgfrb','Angpt2—Tie1','Dll4—Notch3',
            'Dll4—Notch4','Sema6d—Plxna1')

goi.3 <- c('Tnf—Tnfrsf1b','Sct—Sctr','Wnt7a—Lrp6',
           'Npnt—Itgb1','Ereg—Egfr','Wnt7b—Lrp5',
           'Hbegf—Egfr','Lamb3—Itga3','Npnt—Itgb1','Tnf—Tnfrsf21','Tnf—Tnfrsf21',
           'Tnf—Ltbr','Lamb3—Itga2','Mmp7—Cd44','Lamb3—Itga2','Btc—Egfr')

goi.4 <- c('Fgf1—Fgfr3','Fgf1—Egfr','Fgf1—Fgfr2','Fgf2—Fgfr3','Fgf12—Fgfr3','Lama1—Itga2',
           'Csf2—Csf2ra','Ngf—Maged1','Wnt16—Fzd6',
           'Fgf7—Fgfr3','Fn1—Itgb6','Vegfa—Ephb2',
           'Wnt5a—Lrp5','Wnt5a—Ryk','Col4a1—Itga2','Tnc—Egfr')

goi.5 <- c('Amh—Egfr','Hgf—Met','Igf1—Igf1r','Cxcl2—Cxcr2','Ebi3—Il6st','Il18—Il18r1',
           'Tnfsf13—Tnfrsf1a','Nrg4—Egfr','Nrg4—Erbb2','Efna2—Epha1')

goi.epi <- unique(c(goi.0,goi.1,goi.2,goi.3,goi.4,goi.5))

row.labels.epi <- c('Bmp5—Bmpr1b','Rspo3—Lrp6','Fgf10—Fgfr2',
                    'Il10—Il10rb','Wnt3a—Fzd6','Wnt7b—Fzd1',
                    'Dll1—Notch2','Vip—Sctr','Dhh—Ptch1',
                    'Sct—Sctr','Wnt7a—Lrp6','Wnt7b—Lrp5',"Lamb3—Itga3",'Tnf—Ltbr',
                    'Fgf7—Fgfr3','Fn1—Itgb6','Col4a1—Itga2','Wnt5a—Ryk' ,'Il18—Il18r1','Nrg4—Egfr')

png('EpiRec.Heatmap.2024-03-16.png',width=7,height=9,units='in',res=300)
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
                   features = goi.epi,
                   labels = c('Signaling Archetype','Sending Cell Class','Sample'),
                   selected.row.anotations=row.labels.epi,
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.2)
dev.off()


#End


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
  png(file = paste(feature,'NetworkPlot',min.connectivity.thresh,'2024-03-16.png',sep='_'),width = 6,height = 5,units = 'in',res=600)
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

