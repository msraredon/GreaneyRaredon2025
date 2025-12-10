# Global connectomics Part 1: Turning data into NICHES data
# This script takes data from starting cell populations, engineered populations, and native populations, and calculates connectomic data for downstream analysis.
# We only consider macrophage populations within the immune compartment
# We exclude transplant experiment samples from this workflow
# Currently choosing to perform NICHES on raw RNA data without imputation, given the very high level of multi-factorial cross-sample variance in this project

# Require packages
require(Seurat)
require(NICHES)

# Set WD
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics")

# Load color palette (does not include native colors yet)
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Allie's Object Filtration/Final Objects/all.color_palettes_2.R")
color_pals

# Load engineered lung and starting cell data
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Allie's Object Filtration/Final Objects/epi.seurat.objs.int.2023-07-27.Robj")
epi <- object
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Allie's Object Filtration/Final Objects/endo.seurat.objs.int.2023-07-27.Robj")
end <- object
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Allie's Object Filtration/Final Objects/mes.seurat.objs.int.2023-07-27.Robj")
mes <- object
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Allie's Object Filtration/Final Objects/mac.seurat.objs.int.2023-07-27.Robj")
imm <- object
rm(object)

# Investigate metadata
table(mes$Final1)
table(mes$Final2)
table(mes$Dataset2)
table(mes$Dataset3)
table(mes$orig.ident)
table(epi$Final1)
table(epi$Final2)
table(epi$Dataset2)
table(epi$Dataset3)
table(epi$orig.ident)
table(imm$Final1)
table(imm$Final2)
table(imm$Dataset2)
table(imm$Dataset3)
table(imm$orig.ident)
table(end$Final1)
table(end$Final2)
table(end$Dataset2)
table(end$Dataset3)
table(end$orig.ident)

# Set Idents to "Final2"
Idents(epi) <- epi$Final2
Idents(end) <- end$Final2
Idents(mes) <- mes$Final2
Idents(imm) <- imm$Final2
length(unique(Idents(epi))) # 5
length(unique(Idents(end))) # 5
length(unique(Idents(mes))) # 5
length(unique(Idents(imm))) # 4

# Merge into one big object and make idents unique
merge <- merge(epi,list(end,mes,imm))
length(unique(merge$Final2)) # 19 good
Idents(merge) <- merge$Final2

# Split into samples
sample.list <- SplitObject(merge,split.by = 'orig.ident')
names(sample.list)

# Clean up workspace
rm(epi,end,mes,imm)
gc()

# Create starting "tissue" by sampling even numbers of cells from each starting class object, cleaned to only include relevant class
# Epi
table(sample.list$BC1P3$class)
table(sample.list$BC1P6$class)
epi.1 <- sample.list$BC1P3
epi.2 <- sample.list$BC1P6
# End
table(sample.list$RLMVEC$class)
end <- sample.list$RLMVEC
# Mes
table(sample.list$FB13$class)
table(sample.list$FB14$class)
mes.1 <- subset(sample.list$FB13,subset = class =='Mesenchyme')
mes.2 <- subset(sample.list$FB14,subset = class =='Mesenchyme')
# Imm
table(sample.list$MacAlv$class)
imm <- sample.list$MacAlv
# Set Idents to class to allow even downsampling
Idents(epi.1) <- epi.1$class
Idents(epi.2) <- epi.2$class
Idents(end) <- end$class
Idents(mes.1) <- mes.1$class
Idents(mes.2) <- mes.2$class
Idents(imm) <- imm$class

# Downsample
epi.1 <- subset(epi.1,cells = WhichCells(epi.1,downsample = 1500))
epi.2 <- subset(epi.2,cells = WhichCells(epi.2,downsample = 1500))
end <- subset(end,cells = WhichCells(end,downsample = 3000))
mes.1 <- subset(mes.1,cells = WhichCells(mes.1,downsample = 1500))
mes.2 <- subset(mes.2,cells = WhichCells(mes.2,downsample = 1500))
imm <- subset(imm,cells = WhichCells(imm,downsample = 3000))

# Assemble 'pseudo.start' tissue
pseudo.start <- merge(epi.1,list(epi.2,end,mes.1,mes.2,imm))
table(pseudo.start$class)
table(pseudo.start$Final2)

# Swap out starting populations for the pseudo.start tissue
names(sample.list)
to.remove <- c('BC1P3','BC1P6','RLMVEC','MacAlv','FB13','FB14')
tissue.list <- sample.list[!(names(sample.list) %in% to.remove)]
names(tissue.list)
final.names <- c(names(tissue.list),'Pseudo.Start')
tissue.list <- c(tissue.list,pseudo.start)
names(tissue.list) <- final.names
names(tissue.list)

# Add native reference data here
# Load Native Blueprint
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Sam's Object Filtration/Native_Blueprint/lung.combined.clean.classed.annotated.final.2022-07-24.Robj")
# Revised Annotations
load("~/T5/BEFM Swing Space Collected/Native_Blueprint_Connectomics_Revised/meta.data.2022-08-22.Robj")
lung.combined@meta.data <- meta.data
table(lung.combined$CellType_Final)

# format native metadat to work with other metadata
lung.combined$Final1 <- lung.combined$CellType_Final
lung.combined$Final2 <- paste(lung.combined$CellType_Final,'Native',sep='_')
lung.combined$Dataset2 <- 'Native'
lung.combined$Dataset3 <- 'Native'

# remove non-macrophage immune cells
table(lung.combined$CellType_Final)
Idents(lung.combined) <- lung.combined$CellType_Final
immune.cells.to.remove <- c('B','Mo_Act','Neutro','Mo_C','Mo_NC','T_Killer','T','DC','NK','Eos','ILC','Mast','Plasma','Baso')
lung.sub <- subset(lung.combined,idents = immune.cells.to.remove,invert=T)
table(Idents(lung.sub))

# split by tissue
table(lung.sub$orig.ident)
split <- SplitObject(lung.sub,split.by = 'orig.ident')
names(split)

# concatenate native in with the other tissues
tissue.list <- c(tissue.list,split)
names(tissue.list)

# Add Tissue meta.data
for(i in 1:length(tissue.list)){
  tissue.list[[i]]@meta.data$Tissue <- names(tissue.list)[i]
}

# Make consistent size? Might not be necessary....not sure how it affects the average values....lossy...avoid for now?

# Run NICHES by tissue
con.list <- list()
for (i in 1:length(tissue.list)){
  print(i)
  Idents(tissue.list[[i]]) <- tissue.list[[i]]$Final2
  con.list[[i]] <- RunNICHES(tissue.list[[i]],
                        species = 'rat',
                        assay = 'RNA', # Note that we are choosing here to not use imputed data
                        LR.database = 'fantom5',
                        meta.data.to.map = names(tissue.list[[i]]@meta.data),
                        CellToCell = T,
                        CellToSystem = T,
                        SystemToCell = T,
                        cell_types = 'Final2',
                        blend = 'mean.adj')
}

# Extract System to Cell Measurements, Format, and Save
# extract
temp.list <- list()
for(i in 1:length(con.list)){
  temp.list[[i]] <- con.list[[i]]$SystemToCell # Isolate SystemToCell Signaling
  gc()
}
system.to.cell <- merge(temp.list[[1]],temp.list[2:length(temp.list)])
# assess
system.to.cell
table(system.to.cell$Tissue)
table(system.to.cell$Final1)
table(system.to.cell$Final2)
table(system.to.cell$Dataset2)
table(system.to.cell$Dataset3)

# order metadata
system.to.cell$Dataset2 <- factor(system.to.cell$Dataset2,
                                  levels = c('Start','Mono','Co', 'Tri_E','Tri_L','Quad_E','Quad_L','Native'))
system.to.cell$Dataset3 <- factor(system.to.cell$Dataset3,
                                  levels = c('Start','Eng','Native'))
table(system.to.cell$Dataset2)
table(system.to.cell$Dataset3)

system.to.cell$Tissue <- factor(system.to.cell$Tissue,
                                  levels = c('Pseudo.Start','BCL5','BCEC2',
                                             'BEF1','BEF2','BEF3',
                                             'BEF12','BEF14','BEF15',
                                             'BEFM1','BEFM2',
                                             'BEFM4','BEFM5','BEFM6',
                                             'fRat','mRat','P0-f','P0-m'))
table(system.to.cell$Tissue)
# scale
system.to.cell <- ScaleData(system.to.cell,features = rownames(system.to.cell))

# save
save(system.to.cell,file='system.to.cell.BEFM.project.2023-08-03.Robj')
