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

### PreProcessing
# Compute intra-class archetype markers
Idents(mes.rec) <- mes.rec$Archetype
Idents(epi.rec) <- epi.rec$Archetype
Idents(end.rec) <- end.rec$Archetype
Idents(imm.rec) <- imm.rec$Archetype
mark.mes.rec <- FindAllMarkers(mes.rec,min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
mark.epi.rec <- FindAllMarkers(epi.rec,min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
mark.end.rec <- FindAllMarkers(end.rec,min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
mark.imm.rec <- FindAllMarkers(imm.rec,min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
mark.mes.rec$ratio <- mark.mes.rec$pct.1/mark.mes.rec$pct.2
mark.epi.rec$ratio <- mark.epi.rec$pct.1/mark.epi.rec$pct.2
mark.end.rec$ratio <- mark.end.rec$pct.1/mark.end.rec$pct.2
mark.imm.rec$ratio <- mark.imm.rec$pct.1/mark.imm.rec$pct.2
mark.mes.rec$power <- mark.mes.rec$ratio*mark.mes.rec$avg_log2FC
mark.epi.rec$power <- mark.epi.rec$ratio*mark.epi.rec$avg_log2FC
mark.end.rec$power <- mark.end.rec$ratio*mark.end.rec$avg_log2FC
mark.imm.rec$power <- mark.imm.rec$ratio*mark.imm.rec$avg_log2FC

### Assembly
# Merge together to get global object that includes class-specific Archetype labels
merge <- merge(epi.rec,list(end.rec,mes.rec,imm.rec))

#Scale
merge <- ScaleData(merge)

# Clean up: Setup metdata for plotting later
merge$Condition <- merge$Dataset2.Receiving
merge$Condition <- factor(merge$Condition,
                          levels=c('Tri_L','Quad_E'))
table(merge$Condition)


### Analysis

# 1. Global markers of archetypes
Idents(merge) <- merge$Archetype
archetype.markers <- FindAllMarkers(merge,min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
archetype.markers$VarianceType <- 'GlobalArchetype'

#2. Global markers of receiving class
Idents(merge) <- merge$class.Receiving
class.receiving.markers <- FindAllMarkers(merge,min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
class.receiving.markers$VarianceType <- 'ClassReceiving'

#3. Global markers of sending class
Idents(merge) <- merge$class.Sending
class.sending.markers <- FindAllMarkers(merge,min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
class.sending.markers$VarianceType <- 'ClassSending'

#4. Global markers of sample (Sample)
Idents(merge) <- merge$Sample.Receiving
sample.markers <- FindAllMarkers(merge,min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
sample.markers$VarianceType <- 'Sample'

#5. Global markers of dataset (Condition)
Idents(merge) <- merge$Condition
condition.markers <- FindAllMarkers(merge,min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
condition.markers$VarianceType <- 'Condition'

# 6. Define the "bounds" of each class, globally - those features which DO mark the niches of the 4 classes
epi.bounds <- class.receiving.markers[class.receiving.markers$cluster=='Epithelium',]
end.bounds <- class.receiving.markers[class.receiving.markers$cluster=='Endothelium',]
mes.bounds <- class.receiving.markers[class.receiving.markers$cluster=='Mesenchyme',]
imm.bounds <- class.receiving.markers[class.receiving.markers$cluster=='Immune',]

# 7. Identify a limited list of features which describe intra-class variance AND also mark the same class globally
mark.epi.limited <- mark.epi.rec[mark.epi.rec$gene %in% epi.bounds$gene,]
mark.end.limited <- mark.end.rec[mark.end.rec$gene %in% end.bounds$gene,]
mark.mes.limited <- mark.mes.rec[mark.mes.rec$gene %in% mes.bounds$gene,]
mark.imm.limited <- mark.imm.rec[mark.imm.rec$gene %in% imm.bounds$gene,]

# 8.  Identify a limited list of features which describe intra-class variance AND also mark the same class globally AND
# does not mark any other class globally

mark.epi.super <- mark.epi.limited
epi.curve <- data.frame(Criteria = NA,NumFeatures = NA)
epi.curve[1,]$NumFeatures = nrow(mark.epi.super)
mark.epi.super <- mark.epi.super[!(mark.epi.super$gene %in% end.bounds$gene),]
epi.curve[2,]$NumFeatures = nrow(mark.epi.super)
mark.epi.super <- mark.epi.super[!(mark.epi.super$gene %in% mes.bounds$gene),]
epi.curve[3,]$NumFeatures = nrow(mark.epi.super)
mark.epi.super <- mark.epi.super[!(mark.epi.super$gene %in% imm.bounds$gene),]
epi.curve[4,]$NumFeatures = nrow(mark.epi.super)
epi.curve$Criteria <- c('Limited IntraClass Archetype Markers',
                        'that do not mark End',
                        'and do not mark Mes',
                        'and do not mark Imm')
rownames(epi.curve) <- c('Limited IntraClass Archetype Markers',
                         'that do not mark End',
                         'and do not mark Mes',
                         'and do not mark Imm')
epi.curve$Criteria <- factor(epi.curve$Criteria,
                             levels = c('Limited IntraClass Archetype Markers',
                                        'that do not mark End',
                                        'and do not mark Mes',
                                        'and do not mark Imm') )
png(filename = 'Epithelial Niche Archtype Markers Filtration.png',width = 7,height = 7,units = 'in',res=300)
ggplot(data = epi.curve,
       aes(x = Criteria,y=NumFeatures))+
  geom_point()+
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylim(0,600)+
  ylab('Number of Differential Features')+
  ggtitle('Epithelial Niche Archetype Markers')
dev.off()

# Endothelium 
mark.end.super <- mark.end.limited
end.curve <- data.frame(Criteria = NA,NumFeatures = NA)
end.curve[1,]$NumFeatures = nrow(mark.end.super)
mark.end.super <- mark.end.super[!(mark.end.super$gene %in% epi.bounds$gene),]
end.curve[2,]$NumFeatures = nrow(mark.end.super)
mark.end.super <- mark.end.super[!(mark.end.super$gene %in% mes.bounds$gene),]
end.curve[3,]$NumFeatures = nrow(mark.end.super)
mark.end.super <- mark.end.super[!(mark.end.super$gene %in% imm.bounds$gene),]
end.curve[4,]$NumFeatures = nrow(mark.end.super)
end.curve$Criteria <- c('Limited IntraClass Archetype Markers',
                        'that do not mark Epi',
                        'and do not mark Mes',
                        'and do not mark Imm')
rownames(end.curve) <- c('Limited IntraClass Archetype Markers',
                         'that do not mark Epi',
                         'and do not mark Mes',
                         'and do not mark Imm')
end.curve$Criteria <- factor(end.curve$Criteria,
                             levels = c('Limited IntraClass Archetype Markers',
                                        'that do not mark Epi',
                                        'and do not mark Mes',
                                        'and do not mark Imm') )
png(filename = 'Endothelial Niche Archtype Markers Filtration.png',width = 7,height = 7,units = 'in',res=300)
ggplot(data = end.curve,
       aes(x = Criteria,y=NumFeatures))+
  geom_point()+
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylim(0,600)+
  ylab('Number of Differential Features')+
  ggtitle('Endothelial Niche Archetype Markers')
dev.off()

# Mesenchyme
mark.mes.super <- mark.mes.limited
mes.curve <- data.frame(Criteria = NA,NumFeatures = NA)
mes.curve[1,]$NumFeatures = nrow(mark.mes.super)
mark.mes.super <- mark.mes.super[!(mark.mes.super$gene %in% end.bounds$gene),]
mes.curve[2,]$NumFeatures = nrow(mark.mes.super)
mark.mes.super <- mark.mes.super[!(mark.mes.super$gene %in% epi.bounds$gene),]
mes.curve[3,]$NumFeatures = nrow(mark.mes.super)
mark.mes.super <- mark.mes.super[!(mark.mes.super$gene %in% imm.bounds$gene),]
mes.curve[4,]$NumFeatures = nrow(mark.mes.super)
mes.curve$Criteria <- c('Limited IntraClass Archetype Markers',
                        'that do not mark End',
                        'and do not mark Epi',
                        'and do not mark Imm')
rownames(mes.curve) <- c('Limited IntraClass Archetype Markers',
                         'that do not mark End',
                         'and do not mark Epi',
                         'and do not mark Imm')
mes.curve$Criteria <- factor(mes.curve$Criteria,
                             levels = c('Limited IntraClass Archetype Markers',
                                        'that do not mark End',
                                        'and do not mark Epi',
                                        'and do not mark Imm') )
png(filename = 'Mesenchyme Niche Archtype Markers Filtration.png',width = 7,height = 7,units = 'in',res=300)
ggplot(data = mes.curve,
       aes(x = Criteria,y=NumFeatures))+
  geom_point()+
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylim(0,1500)+
  ylab('Number of Differential Features')+
  ggtitle('Mesenchyme Niche Archetype Markers')
dev.off()

# Immune
mark.imm.super <- mark.imm.limited
imm.curve <- data.frame(Criteria = NA,NumFeatures = NA)
imm.curve[1,]$NumFeatures = nrow(mark.imm.super)
mark.imm.super <- mark.imm.super[!(mark.imm.super$gene %in% epi.bounds$gene),]
imm.curve[2,]$NumFeatures = nrow(mark.imm.super)
mark.imm.super <- mark.imm.super[!(mark.imm.super$gene %in% mes.bounds$gene),]
imm.curve[3,]$NumFeatures = nrow(mark.imm.super)
mark.imm.super <- mark.imm.super[!(mark.imm.super$gene %in% end.bounds$gene),] #0
imm.curve[4,]$NumFeatures = nrow(mark.imm.super)
imm.curve$Criteria <- c('Limited IntraClass Archetype Markers',
                        'that do not mark Epi',
                        'and do not mark Mes',
                        'and do not mark End')
rownames(imm.curve) <- c('Limited IntraClass Archetype Markers',
                         'that do not mark Epi',
                         'and do not mark Mes',
                         'and do not mark End')
imm.curve$Criteria <- factor(imm.curve$Criteria,
                             levels = c('Limited IntraClass Archetype Markers',
                                       'that do not mark Epi',
                                       'and do not mark Mes',
                                       'and do not mark End') )
png(filename = 'Immune Niche Archtype Markers Filtration.png',width = 7,height = 7,units = 'in',res=300)
ggplot(data = imm.curve,
       aes(x = Criteria,y=NumFeatures))+
  geom_point()+
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylim(0,600)+
  ylab('Number of Differential Features')+
  ggtitle('Immune Niche Archetype Markers')
dev.off()

# 9. Experiment: Class + IntraClass Condition + Archetype

# Compute intra-class condition markers
Idents(epi.rec) <- epi.rec$Dataset2.Receiving
Idents(mes.rec) <- mes.rec$Dataset2.Receiving
Idents(end.rec) <- end.rec$Dataset2.Receiving # Immune not required, no conditional
intra.epi.mark <- FindAllMarkers(epi.rec,min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
intra.mes.mark <- FindAllMarkers(mes.rec,min.pct=0.1,logfc.threshold = 0.1,only.pos = T)
intra.end.mark <- FindAllMarkers(end.rec,min.pct=0.1,logfc.threshold = 0.1,only.pos = T) # Immune not required, no conditional

# Intersection between class markers and intraclass condition markers
epi.selected <- intersect(epi.bounds$gene,intra.epi.mark$gene)
mes.selected <- intersect(mes.bounds$gene,intra.mes.mark$gene)
end.selected <- intersect(end.bounds$gene,intra.end.mark$gene)
imm.selected <- imm.bounds$gene # Does not apply to immune

# Internal archetype markers that are also internal condition markers that are also external class markers
epi.selected <- mark.epi.rec[mark.epi.rec$gene %in% epi.selected,] # Internal archetype markers that are also markers of the class as a whole and are markers of one condition over the other
mes.selected <- mark.mes.rec[mark.mes.rec$gene %in% mes.selected,] # Internal archetype markers that are also markers of the class as a whole and are markers of one condition over the other
end.selected <- mark.end.rec[mark.end.rec$gene %in% end.selected,] # Internal archetype markers that are also markers of the class as a whole and are markers of one condition over the other
imm.selected <- mark.imm.rec[mark.imm.rec$gene %in% imm.selected,] # Condition is not applicable for immune niche

# Limit to just the home class
mark.epi.ultra <- epi.selected
mark.epi.ultra <- mark.epi.ultra[!(mark.epi.ultra$gene %in% mes.bounds$gene),]
mark.epi.ultra <- mark.epi.ultra[!(mark.epi.ultra$gene %in% end.bounds$gene),]
mark.epi.ultra <- mark.epi.ultra[!(mark.epi.ultra$gene %in% imm.bounds$gene),]

mark.mes.ultra <- mes.selected
mark.mes.ultra <- mark.mes.ultra[!(mark.mes.ultra$gene %in% epi.bounds$gene),]
mark.mes.ultra <- mark.mes.ultra[!(mark.mes.ultra$gene %in% imm.bounds$gene),]
mark.mes.ultra <- mark.mes.ultra[!(mark.mes.ultra$gene %in% end.bounds$gene),]

mark.end.ultra <- end.selected
mark.end.ultra <- mark.end.ultra[!(mark.end.ultra$gene %in% epi.bounds$gene),]
mark.end.ultra <- mark.end.ultra[!(mark.end.ultra$gene %in% imm.bounds$gene),]
mark.end.ultra <- mark.end.ultra[!(mark.end.ultra$gene %in% mes.bounds$gene),]

mark.imm.ultra <- imm.selected
mark.imm.ultra <- mark.imm.ultra[!(mark.imm.ultra$gene %in% epi.bounds$gene),]
mark.imm.ultra <- mark.imm.ultra[!(mark.imm.ultra$gene %in% mes.bounds$gene),]
mark.imm.ultra <- mark.imm.ultra[!(mark.imm.ultra$gene %in% end.bounds$gene),]

# Save for later
save(mark.epi.ultra,file = 'mark.epi.ultra.Robj')
save(mark.end.ultra,file = 'mark.end.ultra.Robj')
save(mark.imm.ultra,file = 'mark.imm.ultra.Robj')
save(mark.mes.ultra,file = 'mark.mes.ultra.Robj')

