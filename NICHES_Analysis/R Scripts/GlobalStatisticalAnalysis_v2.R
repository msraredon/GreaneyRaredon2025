# Global Statistical Analysis - Global Connectomics BEFM

# This script pulls all of the BEFM connectomic data together in one place. 
# We then compute a bunch of marker lists that will be useful later on, because they describe specific variance of interest both for study and exclusion.

# Set WD
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics")

# Packages
require(Seurat)
require(ggplot2)

### Load data 

## Epithelial Receiving Connectivity Object, clustered, annotated, organized for plotting
## Created by 'R Scripts/EpithelialReceivingObject_Clean.R
load("Saved Objects/epi.rec.2024-03-15.Robj")

## Endothelial Receiving Connectivity Object, clustered, annotated, organized for plotting
## Created by 'R Scripts/EndothelialReceivingObject_Clean.R
load("Saved Objects/end.rec.2024-03-15.Robj")

## Mesenchyme Receiving Connectivity Object, clustered, annotated, organized for plotting
## Created by 'R Scripts/MesenchymeReceivingObject_Clean.R
load("Saved Objects/mes.rec.2024-03-15.Robj")

## Immune Receiving Connectivity Object, clustered, annotated, organized for plotting
## Created by 'R Scripts/ImmuneReceivingObject_Clean.R
load("Saved Objects/imm.rec.2024-03-15.Robj")

## PreProcessing

# Merge together to get global object that includes class-specific Archetype labels
global <- merge(epi.rec,list(end.rec,mes.rec,imm.rec))

# Remove edges from the 2 contaminating macrophages

Idents(global) <- global$orig.ident
bef12.ctc <- subset(global,idents = 'BEF12')
barcodes.remove.ctc <- rownames(bef12.ctc@meta.data[which(bef12.ctc@meta.data$class.Sending == 'Immune'),])
barcodes.remove.ctc <- c(barcodes.remove.ctc,
                         rownames(bef12.ctc@meta.data[which(bef12.ctc@meta.data$class.Receiving == 'Immune'),]))
global <- subset(global,cells = barcodes.remove.ctc,invert=T)


#Scale
global <- ScaleData(global)

# Define metadata handles
# a. Condition (Dataset 2)
global$Condition <- global$Dataset2.Receiving
global$Condition <- factor(global$Condition,levels=c('Tri_L','Quad_E'))

# b. Vector (class.Joint)
global$Vector <- global$class.Joint

# c. Niche (class.Receiving)
global$Niche <- global$class.Receiving

# d. Influence (class.Sending)
global$Influence <- global$class.Sending

# e. Archetype (Archetype)
global$Archetype <- global$Archetype

#### Initialize a mega marker list of lists
marks <- list()
marks$global <- list()
marks$condition <- list()
marks$niche <- list()
marks$influence <- list()
marks$vector <- list()

#### Compute Marker Lists

# 1. Markers within the Global Object (marks$global)

# a. Archetypes
Idents(global) <- global$Archetype
archetype.markers <- FindAllMarkers(global,min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
archetype.markers$ratio <- archetype.markers$pct.1/archetype.markers$pct.2
archetype.markers$power <- archetype.markers$ratio*archetype.markers$avg_log2FC
marks$global$archetype <- archetype.markers

# b. Condition
Idents(global) <- global$Condition
condition.markers <- FindAllMarkers(global,min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
condition.markers$ratio <- condition.markers$pct.1/condition.markers$pct.2
condition.markers$power <- condition.markers$ratio*condition.markers$avg_log2FC
marks$global$condition <- condition.markers

# c. Vector
Idents(global) <- global$Vector
temp.markers <- FindAllMarkers(global,min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
temp.markers$ratio <- temp.markers$pct.1/temp.markers$pct.2
temp.markers$power <- temp.markers$ratio*temp.markers$avg_log2FC
marks$global$vector  <- temp.markers

# d. Niche
Idents(global) <- global$Niche
temp.markers <- FindAllMarkers(global,min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
temp.markers$ratio <- temp.markers$pct.1/temp.markers$pct.2
temp.markers$power <- temp.markers$ratio*temp.markers$avg_log2FC
marks$global$niche  <- temp.markers

# e. Influence
Idents(global) <- global$Influence
temp.markers <- FindAllMarkers(global,min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
temp.markers$ratio <- temp.markers$pct.1/temp.markers$pct.2
temp.markers$power <- temp.markers$ratio*temp.markers$avg_log2FC
marks$global$influence  <- temp.markers

# 2. Markers within each Condition (marks$condition)
Idents(global) <- global$Condition
temp.split <- SplitObject(global)
names(temp.split)

# a. Vector
temp.mark.list <- list()
for(i in 1:length(temp.split)){
  Idents(temp.split[[i]]) <- temp.split[[i]]$Vector
  temp.markers <- FindAllMarkers(temp.split[[i]],min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
  temp.markers$ratio <- temp.markers$pct.1/temp.markers$pct.2
  temp.markers$power <- temp.markers$ratio*temp.markers$avg_log2FC
  temp.mark.list[[i]]  <- temp.markers
}
names(temp.mark.list) <- names(temp.split)
marks$condition$vector <- temp.mark.list

# b. Niche
temp.mark.list <- list()
for(i in 1:length(temp.split)){
  Idents(temp.split[[i]]) <- temp.split[[i]]$Niche
  temp.markers <- FindAllMarkers(temp.split[[i]],min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
  temp.markers$ratio <- temp.markers$pct.1/temp.markers$pct.2
  temp.markers$power <- temp.markers$ratio*temp.markers$avg_log2FC
  temp.mark.list[[i]]  <- temp.markers
}
names(temp.mark.list) <- names(temp.split)
marks$condition$niche <- temp.mark.list

# c. Influence
temp.mark.list <- list()
for(i in 1:length(temp.split)){
  Idents(temp.split[[i]]) <- temp.split[[i]]$Influence
  temp.markers <- FindAllMarkers(temp.split[[i]],min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
  temp.markers$ratio <- temp.markers$pct.1/temp.markers$pct.2
  temp.markers$power <- temp.markers$ratio*temp.markers$avg_log2FC
  temp.mark.list[[i]]  <- temp.markers
}
names(temp.mark.list) <- names(temp.split)
marks$condition$influence <- temp.mark.list

# d. Archetype
temp.mark.list <- list()
for(i in 1:length(temp.split)){
  Idents(temp.split[[i]]) <- temp.split[[i]]$Archetype
  temp.markers <- FindAllMarkers(temp.split[[i]],min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
  temp.markers$ratio <- temp.markers$pct.1/temp.markers$pct.2
  temp.markers$power <- temp.markers$ratio*temp.markers$avg_log2FC
  temp.mark.list[[i]]  <- temp.markers
}
names(temp.mark.list) <- names(temp.split)
marks$condition$archetype <- temp.mark.list

# 3. Markers within each Niche (marks$niche)
Idents(global) <- global$Niche
temp.split <- SplitObject(global)
names(temp.split)

# a. Archetype
temp.mark.list <- list()
for(i in 1:length(temp.split)){
  Idents(temp.split[[i]]) <- temp.split[[i]]$Archetype
  temp.markers <- FindAllMarkers(temp.split[[i]],min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
  temp.markers$ratio <- temp.markers$pct.1/temp.markers$pct.2
  temp.markers$power <- temp.markers$ratio*temp.markers$avg_log2FC
  temp.mark.list[[i]]  <- temp.markers
}
names(temp.mark.list) <- names(temp.split)
marks$niche$archetype <- temp.mark.list

# b. Vector
temp.mark.list <- list()
for(i in 1:length(temp.split)){
  Idents(temp.split[[i]]) <- temp.split[[i]]$Vector
  temp.markers <- FindAllMarkers(temp.split[[i]],min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
  temp.markers$ratio <- temp.markers$pct.1/temp.markers$pct.2
  temp.markers$power <- temp.markers$ratio*temp.markers$avg_log2FC
  temp.mark.list[[i]]  <- temp.markers
}
names(temp.mark.list) <- names(temp.split)
marks$niche$vector <- temp.mark.list

# c. Influence
temp.mark.list <- list()
for(i in 1:length(temp.split)){
  Idents(temp.split[[i]]) <- temp.split[[i]]$Influence
  temp.markers <- FindAllMarkers(temp.split[[i]],min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
  temp.markers$ratio <- temp.markers$pct.1/temp.markers$pct.2
  temp.markers$power <- temp.markers$ratio*temp.markers$avg_log2FC
  temp.mark.list[[i]]  <- temp.markers
}
names(temp.mark.list) <- names(temp.split)
marks$niche$influence <- temp.mark.list

# d. Condition
temp.mark.list <- list()
for(i in 1:length(temp.split)){
  Idents(temp.split[[i]]) <- temp.split[[i]]$Condition
  temp.markers <- FindAllMarkers(temp.split[[i]],min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
  temp.markers$ratio <- temp.markers$pct.1/temp.markers$pct.2
  temp.markers$power <- temp.markers$ratio*temp.markers$avg_log2FC
  temp.mark.list[[i]]  <- temp.markers
}
names(temp.mark.list) <- names(temp.split)
marks$niche$condition <- temp.mark.list

# 4. Markers within each Influence (marks$influence)
Idents(global) <- global$Influence
temp.split <- SplitObject(global)
names(temp.split)

# a. Archetype
temp.mark.list <- list()
for(i in 1:length(temp.split)){
  Idents(temp.split[[i]]) <- temp.split[[i]]$Archetype
  temp.markers <- FindAllMarkers(temp.split[[i]],min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
  temp.markers$ratio <- temp.markers$pct.1/temp.markers$pct.2
  temp.markers$power <- temp.markers$ratio*temp.markers$avg_log2FC
  temp.mark.list[[i]]  <- temp.markers
}
names(temp.mark.list) <- names(temp.split)
marks$influence$archetype <- temp.mark.list

# b. Vector
temp.mark.list <- list()
for(i in 1:length(temp.split)){
  Idents(temp.split[[i]]) <- temp.split[[i]]$Vector
  temp.markers <- FindAllMarkers(temp.split[[i]],min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
  temp.markers$ratio <- temp.markers$pct.1/temp.markers$pct.2
  temp.markers$power <- temp.markers$ratio*temp.markers$avg_log2FC
  temp.mark.list[[i]]  <- temp.markers
}
names(temp.mark.list) <- names(temp.split)
marks$influence$vector <- temp.mark.list

# c. Niche
temp.mark.list <- list()
for(i in 1:length(temp.split)){
  Idents(temp.split[[i]]) <- temp.split[[i]]$Niche
  temp.markers <- FindAllMarkers(temp.split[[i]],min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
  temp.markers$ratio <- temp.markers$pct.1/temp.markers$pct.2
  temp.markers$power <- temp.markers$ratio*temp.markers$avg_log2FC
  temp.mark.list[[i]]  <- temp.markers
}
names(temp.mark.list) <- names(temp.split)
marks$influence$niche <- temp.mark.list

# d. Condition
temp.mark.list <- list()
for(i in 1:length(temp.split)){
  Idents(temp.split[[i]]) <- temp.split[[i]]$Condition
  temp.markers <- FindAllMarkers(temp.split[[i]],min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
  temp.markers$ratio <- temp.markers$pct.1/temp.markers$pct.2
  temp.markers$power <- temp.markers$ratio*temp.markers$avg_log2FC
  temp.mark.list[[i]]  <- temp.markers
}
names(temp.mark.list) <- names(temp.split)
marks$influence$condition <- temp.mark.list

# 5. Markers within each Vector (marks$vector)
Idents(global) <- global$Vector
temp.split <- SplitObject(global)
names(temp.split)

# a. Archetype
temp.mark.list <- list()
for(i in 1:length(temp.split)){
  Idents(temp.split[[i]]) <- temp.split[[i]]$Archetype
  temp.markers <- FindAllMarkers(temp.split[[i]],min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
  temp.markers$ratio <- temp.markers$pct.1/temp.markers$pct.2
  temp.markers$power <- temp.markers$ratio*temp.markers$avg_log2FC
  temp.mark.list[[i]]  <- temp.markers
}
names(temp.mark.list) <- names(temp.split)
marks$vector$archetype <- temp.mark.list

# b. Influence
temp.mark.list <- list()
for(i in 1:length(temp.split)){
  Idents(temp.split[[i]]) <- temp.split[[i]]$Influence
  temp.markers <- FindAllMarkers(temp.split[[i]],min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
  temp.markers$ratio <- temp.markers$pct.1/temp.markers$pct.2
  temp.markers$power <- temp.markers$ratio*temp.markers$avg_log2FC
  temp.mark.list[[i]]  <- temp.markers
}
names(temp.mark.list) <- names(temp.split)
marks$vector$influence <- temp.mark.list

# c. Niche
temp.mark.list <- list()
for(i in 1:length(temp.split)){
  Idents(temp.split[[i]]) <- temp.split[[i]]$Niche
  temp.markers <- FindAllMarkers(temp.split[[i]],min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
  temp.markers$ratio <- temp.markers$pct.1/temp.markers$pct.2
  temp.markers$power <- temp.markers$ratio*temp.markers$avg_log2FC
  temp.mark.list[[i]]  <- temp.markers
}
names(temp.mark.list) <- names(temp.split)
marks$vector$niche <- temp.mark.list

# d. Condition
temp.mark.list <- list()
for(i in 1:length(temp.split)){
  Idents(temp.split[[i]]) <- temp.split[[i]]$Condition
  temp.markers <- FindAllMarkers(temp.split[[i]],min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
  temp.markers$ratio <- temp.markers$pct.1/temp.markers$pct.2
  temp.markers$power <- temp.markers$ratio*temp.markers$avg_log2FC
  temp.mark.list[[i]]  <- temp.markers
}
names(temp.mark.list) <- names(temp.split)
marks$vector$condition <- temp.mark.list

# Save for later
save(marks,file = 'marks.2024-03-24.Robj')
# Save for later
save(global,file = 'global.connectomics.TriQuad.2024-03-29.Robj')

