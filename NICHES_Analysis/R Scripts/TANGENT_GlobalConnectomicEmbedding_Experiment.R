# Quick Experiment: Gloabl connectivity embeddings

# Load data
## Epithelial Receiving Connectivity Object, clustered, annotated, organized for plotting
## Marker list associated with the above
## Created by 'R Scripts/EpithelialReceivingObject_Clean.R
load("Saved Objects/epi.rec.2024-03-15.Robj")
load("Saved Objects/mark.epi.rec.2024-03-15.Robj")
## Endothelial Receiving Connectivity Object, clustered, annotated, organized for plotting
## Marker list associated with the above
## Created by 'R Scripts/EndothelialReceivingObject_Clean.R
load("Saved Objects/end.rec.2024-03-15.Robj")
load("Saved Objects/mark.end.rec.2024-03-15.Robj")
## Mesenchyme Receiving Connectivity Object, clustered, annotated, organized for plotting
## Marker list associated with the above
## Created by 'R Scripts/MesenchymeReceivingObject_Clean.R
load("Saved Objects/mes.rec.2024-03-15.Robj")
load("Saved Objects/mark.mes.rec.2024-03-15.Robj")
## Immune Receiving Connectivity Object, clustered, annotated, organized for plotting
## Marker list associated with the above
## Created by 'R Scripts/ImmuneReceivingObject_Clean.R
load("Saved Objects/imm.rec.2024-03-15.Robj")
load("Saved Objects/mark.imm.rec.2024-03-15.Robj")

# Merge together to get global object with Archetype labels
merge <- merge(epi.rec,list(end.rec,mes.rec,imm.rec))

#Scale
merge <- ScaleData(merge)

# Find features that differentiate receiving class and nothing else
Idents(merge) <- merge$class.Receiving
foi <- FindAllMarkers(merge,min.cells.feature = 1000,only.pos=F)

# Exclude features that allow intra-class segregation
Idents(epi.rec) <- epi.rec$Sample.Receiving
epi.remove <- FindAllMarkers(epi.rec,min.cells.feature = 250,min.pct = 0.25,only.pos=F)
Idents(mes.rec) <- mes.rec$Sample.Receiving
mes.remove <- FindAllMarkers(mes.rec,min.cells.feature = 250,min.pct = 0.25,only.pos=F)
Idents(end.rec) <- end.rec$Sample.Receiving
end.remove <- FindAllMarkers(end.rec,min.cells.feature = 250,min.pct = 0.25,only.pos=F)
Idents(imm.rec) <- imm.rec$Sample.Receiving
imm.remove <- FindAllMarkers(imm.rec,min.cells.feature = 250,min.pct = 0.25,only.pos=F)
to.remove <- unique(c(epi.remove$gene,mes.remove$gene,end.remove$gene,imm.remove$gene))
foi.2 <- foi[!(foi$gene %in% to.remove),]

# Embed
merge <- RunUMAP(merge,features = unique(foi$gene))

# Plot
DimPlot(merge,group.by='class.Receiving',cols = color_pals$class)+NoAxes()
DimPlot(merge,group.by='Sample.Receiving',cols = color_pals$sample_colors)+NoAxes()
DimPlot(merge,group.by='Dataset2.Sending',cols = color_pals$dataset_colors)+NoAxes()

DimPlot(merge,group.by='Archetype',label=T)+NoAxes()+NoLegend() # temp colors for viz
DimPlot(merge,group.by='Archetype',label=T,split.by = 'Dataset2.Sending')+NoAxes()+NoLegend() # temp colors for viz

names(table(merge$Archetype))

regenerative.archetypes <- c("Epi-Auto Regenerative","Epi-Mes BEFM1","Mes-Epi Regenerative","Regenerative Mes-Endo",
                             "Mes-Auto BEFM1" ,"Mes-Auto BEFM2", "Mes-Immune BEFM1","Mes-Immune BEFM2",
                             "Epi-Immune BEFM1","Epi-Immune BEFM2","Macrophage-Endo","Macrophage-Epi")
common.archetypes <- c("Baseline Endo-Auto","Baseline Epi-Endo","Baseline Mes-Endo","Endo-Epi Common","Endo-Immune Common",
                       "Endo-Mes Common","Epi-Auto Common","Epi-Mes Common","Immune-Auto Common", "Macrophage-Mes Common","Mes-Epi Common")
aberrant.archetypes <- c( "Mes-Auto BEF12","Mes-Auto BEF14","Mes-Auto BEF15")

Idents(merge) <- merge$Archetype
edge.regen <- WhichCells(merge, idents = regenerative.archetypes)
edge.common <- WhichCells(merge, idents = common.archetypes)
edge.aberrant <- WhichCells(merge, idents = aberrant.archetypes)

merge <- SetIdent(merge, cells = edge.regen, value = 'Regenerative Archetypes')
merge <- SetIdent(merge, cells = edge.common, value = 'Common Archetypes')
merge <- SetIdent(merge, cells = edge.aberrant, value = 'Aberrant Archetypes')

merge$Archetype.Category <- Idents(merge)
Idents(merge) <- merge$class.Receiving

p1 <- DimPlot(merge,group.by='Archetype',label=T,split.by = 'Archetype.Category')+NoAxes()+NoLegend() # temp colors for viz
p2 <- DimPlot(merge,group.by='class.Receiving',cols = color_pals$class,label=T,split.by = 'Archetype.Category')+NoAxes()+NoLegend() # temp colors for viz
p3 <- DimPlot(merge,group.by='Sample.Receiving',cols = color_pals$sample_colors,label=T,split.by = 'Archetype.Category')+NoAxes()+NoLegend() # temp colors for viz
cowplot::plot_grid(p1,p2,p3,nrow = 3)
png(filename = paste("Global Connectomic Embedding Experiment 2024-03-17.png"),width = 28,height = 22,units = 'in',res=600)
print(cowplot::plot_grid(p1,p2,p3,nrow = 3))
dev.off()

merge.standard <- merge
merge.standard <- FindVariableFeatures(merge.standard)
merge.standard <- RunPCA(merge.standard)
ElbowPlot(merge.standard,ndims = 50)
merge.standard <- RunUMAP(merge.standard,dims = 1:50)

p1 <- DimPlot(merge.standard,group.by='Archetype',label=T,split.by = 'Archetype.Category')+NoAxes()+NoLegend() # temp colors for viz
p2 <- DimPlot(merge.standard,group.by='class.Receiving',cols = color_pals$class,label=T,split.by = 'Archetype.Category')+NoAxes()+NoLegend() # temp colors for viz
p3 <- DimPlot(merge.standard,group.by='Sample.Receiving',cols = color_pals$sample_colors,label=T,split.by = 'Archetype.Category')+NoAxes()+NoLegend() # temp colors for viz
cowplot::plot_grid(p1,p2,p3,nrow = 3)
png(filename = paste("Global Connectomic Embedding Experiment Standard 2024-03-17.png"),width = 28,height = 22,units = 'in',res=600)
print(cowplot::plot_grid(p1,p2,p3,nrow = 3))
dev.off()

DimPlot(merge.standard,group.by='class.Sending',label=T)+NoAxes()+NoLegend() # temp colors for viz
DimPlot(merge.standard,group.by='class.Receiving',label=T)+NoAxes()+NoLegend() # temp colors for viz

# What if we just used this global to narrow in on really specific findings?
