# Require packages
require(Seurat)
require(NICHES)
require(dplyr)

# Set WD
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics")

# Load color palette (does not include native colors yet)
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Allie's Object Filtration/Final Objects/all.color_palettes_2.R")
color_pals

# load system to cell data
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/system.to.cell.BEFM.project.2023-08-03.Robj")
system.to.cell$Tissue <- factor(system.to.cell$Tissue,
                                levels = c('Pseudo.Start','BCL5','BCEC2',
                                           'BEF1','BEF2','BEF3',
                                           'BEF12','BEF14','BEF15',
                                           'BEFM1','BEFM2',
                                           'BEFM4','BEFM5','BEFM6',
                                           'fRat','mRat','P0-f','P0-m'))

# Set idents to class and subset out each class on its own
Idents(system.to.cell) <- system.to.cell$class
epi <- subset(system.to.cell,subset=class=='Epithelium')
end <- subset(system.to.cell,subset=class=='Endothelium')
mes <- subset(system.to.cell,subset=class=='Mesenchyme')
imm <- subset(system.to.cell,subset=class=='Immune')
epi <- ScaleData(epi)
end <- ScaleData(end)
mes <- ScaleData(mes)
imm <- ScaleData(imm)

# for each class, find three distinct useful marker lists and store
class.list <- unique(system.to.cell$class)
class.markers.total <- list()
for(i in 1:length(class.list)){
  # pull one class
  temp <- subset(system.to.cell,subset=class==class.list[i])
  
  # subset to only high information measurements
  temp <- subset(temp,nFeature_SystemToCell>50)
  
  # determine native niche markers for each cell type compared to each other
  native <- subset(temp,subset = Dataset3=='Native')
  Idents(native) <- native$Final2
  mark.native <- FindAllMarkers(native,min.pct = 0.01,logfc.threshold = 0.01,only.pos = T)
  mark.native$ratio <- mark.native$pct.1/mark.native$pct.2
  mark.native$power <- mark.native$ratio*mark.native$avg_log2FC
  mark.native.top <- mark.native %>% group_by(cluster) %>% top_n(20,power)
  
  # determine starting vs. eng vs. native markers
  Idents(temp) <- temp$Dataset3
  mark.dataset3 <- FindAllMarkers(temp,min.pct = 0.01,logfc.threshold = 0.01,only.pos = F)
  mark.dataset3$ratio <- mark.dataset3$pct.1/mark.dataset3$pct.2
  mark.dataset3$power <- mark.dataset3$ratio*mark.dataset3$avg_log2FC
  mark.dataset3.top <- mark.dataset3 %>% group_by(cluster) %>% top_n(20,power)
  
  # determine engineered mono / co / tri / quad markers within the engineered chunk
  eng <- subset(temp,subset = Dataset3=='Eng')
  Idents(eng) <- eng$Dataset2
  mark.eng.Dataset2 <- FindAllMarkers(eng,min.pct = 0.01,logfc.threshold = 0.01,only.pos = T)
  mark.eng.Dataset2$ratio <- mark.eng.Dataset2$pct.1/mark.eng.Dataset2$pct.2
  mark.eng.Dataset2$power <- mark.eng.Dataset2$ratio*mark.eng.Dataset2$avg_log2FC
  mark.eng.Dataset2.top <- mark.eng.Dataset2 %>% group_by(cluster) %>% top_n(20,power)
  
  # determine celltype markers within the engineered chunk
  eng <- subset(temp,subset = Dataset3=='Eng')
  Idents(eng) <- eng$Final1
  mark.eng.Final1 <- FindAllMarkers(eng,min.pct = 0.01,logfc.threshold = 0.01,only.pos = T)
  mark.eng.Final1$ratio <- mark.eng.Final1$pct.1/mark.eng.Final1$pct.2
  mark.eng.Final1$power <- mark.eng.Final1$ratio*mark.eng.Final1$avg_log2FC
  mark.eng.Final1.top <- mark.eng.Final1 %>% group_by(cluster) %>% top_n(20,power)
  
  # determine celltype markers within the whole dataset
  total <- temp
  Idents(total) <- total$Final1
  mark.total.Final1 <- FindAllMarkers(total,min.pct = 0.01,logfc.threshold = 0.01,only.pos = T)
  mark.total.Final1$ratio <- mark.total.Final1$pct.1/mark.total.Final1$pct.2
  mark.total.Final1$power <- mark.total.Final1$ratio*mark.total.Final1$avg_log2FC
  mark.total.Final1.top <- mark.total.Final1 %>% group_by(cluster) %>% top_n(20,power)
  
  # output
  class.marker.set <- list(mark.native,mark.native.top,
                           mark.dataset3,mark.dataset3.top,
                           mark.eng.Dataset2,mark.eng.Dataset2.top,
                           mark.eng.Final1,mark.eng.Final1.top)
  names(class.marker.set) <- c('Native','Native_Top',
                               'Dataset3','Dataset3_Top',
                               'Eng_Dataset2','Eng_Dataset2_Top',
                               'Eng_Final1','Eng_Final1_Top')
  class.markers.total[[i]] <- class.marker.set
}

names(class.markers.total) <- class.list
View(class.markers.total$Epithelium$Dataset3_Top)
View(class.markers.total$Immune$Dataset3)

save(class.markers.total,file = 'class.markers.total.2023-08-15.Robj')

