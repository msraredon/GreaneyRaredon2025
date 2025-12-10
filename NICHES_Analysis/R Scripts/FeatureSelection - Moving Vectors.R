#### LOADS

# Load annotated marker lists
## MADE BY "AnnotateGlobalMarkerLists.R"
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/marks.annotated.2024-03-29.Robj")

# Load global connectivity object for this project (5 samples only)
## MADE BY "
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Saved Objects/global.connectomics.TriQuad.2024-03-29.Robj")

# Allie's comments - this has solved all our problems and creates a highly useful marker list
# Likely sufficient for our paper

# Packages
require(Seurat)
require(ggplot2)

# Set WD
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics")

# Pull out the two marker lists we would like to study closely - split by condition, and then in each condition, markers of VECTOR
temp1 <- marks$condition$vector$Tri_L
temp2 <- marks$condition$vector$Quad_E

# Identify all of the mechanisms that we will be investigating
to.assess <- unique(c(temp1$gene,temp2$gene))

# Assess conservation status
stash.list1 <- list()
stash.list2 <- list()
for(i in 1:length(to.assess)){
  # Assess conservation
  sub1 <- temp1[temp1$gene==to.assess[i],]
  sub2 <- temp2[temp2$gene==to.assess[i],]
  sub1 <- sub1[order(sub1$avg_log2FC, decreasing = TRUE),] # This ordering parameter matters, p_val is perhaps not ideal. Log FC?
  sub2 <- sub2[order(sub2$avg_log2FC, decreasing = TRUE),] # This ordering parameter matters, p_val is perhaps not ideal. Log FC?
  stash.list1[[i]] <- as.character(sub1$cluster)
  stash.list2[[i]] <- as.character(sub2$cluster)
}
names(stash.list1) <- to.assess
names(stash.list2) <- to.assess

# Ask if the ranking is the same between Tri and Quad
stash.list.compared <- list()
for(i in 1:length(to.assess)){
  temp <- stash.list1[[to.assess[i]]] == stash.list2[[to.assess[i]]]
  stash.list.compared[[i]] <- temp
}
names(stash.list.compared) <- to.assess

# Clean up and reorganize in a way that I like
befm.moi <- data.frame(MOI = names(stash.list.compared))
befm.moi$topology.conservation <- NA
for(i in 1:nrow(befm.moi)){
  befm.moi[i,]$topology.conservation <- paste(as.character(stash.list.compared[[i]]),collapse = "-")
}

# Add feature annotations
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/fantom5.rat.annotations.Robj")
mechanisms.to.annotate <- befm.moi$MOI
annotations.to.add <- fantom5.rat.annotations[mechanisms.to.annotate,]
befm.moi.annotated <- cbind(befm.moi,annotations.to.add)

# Figure out where the signal came from
temp1 <- marks$condition$vector$Tri_L
temp2 <- marks$condition$vector$Quad_E

# Does it mark any vector in Tri? Any vector in Quad
befm.moi.annotated$in.tri <- befm.moi.annotated$MOI %in% temp1$gene
befm.moi.annotated$in.quad <- befm.moi.annotated$MOI %in% temp2$gene

# Top rank vector in Tri? In Quad?
befm.moi.annotated$top.in.tri <- NA
befm.moi.annotated$top.in.quad <- NA
for(i in 1:nrow(befm.moi.annotated)){
  befm.moi.annotated[i,]$top.in.tri <- stash.list1[[befm.moi.annotated[i,]$MOI]][1]
  befm.moi.annotated[i,]$top.in.quad <- stash.list2[[befm.moi.annotated[i,]$MOI]][1]
}

# Add condition annotations
befm.moi.annotated$top.conserved <- befm.moi.annotated$top.in.tri == befm.moi.annotated$top.in.quad
# Add more
befm.moi.annotated$is.gained <- befm.moi.annotated$in.quad == TRUE & befm.moi.annotated$in.tri == FALSE
befm.moi.annotated$is.lost <- befm.moi.annotated$in.quad == FALSE & befm.moi.annotated$in.tri == TRUE
befm.moi.annotated$always.present <- befm.moi.annotated$in.quad == TRUE & befm.moi.annotated$in.tri == TRUE

# Fix the NAs
befm.moi.annotated[is.na(befm.moi.annotated$top.conserved),]$top.conserved  <- 'NA'
befm.moi.annotated[is.na(befm.moi.annotated$top.in.tri),]$top.in.tri  <- 'NA'
befm.moi.annotated[is.na(befm.moi.annotated$top.in.quad),]$top.in.quad  <- 'NA'

# Fix the empties
befm.moi.annotated[befm.moi.annotated$topology.conservation=="",]$topology.conservation <- FALSE

# Add a global delta metric
Idents(global) <- global$Condition
temp <- FindAllMarkers(global,features = befm.moi.annotated$MOI,min.pct = 0,logfc.threshold = 0,min.cells.feature = 0,min.cells.group = 0,return.thresh = 1,only.pos = T)
temp$ratio <- temp$pct.1/temp$pct.2
temp$power <- temp$ratio*temp$avg_log2FC
befm.moi.annotated$global.delta <- NA
for(i in 1:nrow(befm.moi.annotated)){
  befm.moi.annotated[i,]$global.delta <- temp[temp$gene==befm.moi.annotated[i,]$MOI,]$avg_log2FC
}

# Set idents to class vector
Idents(global) <- global$Vector

# Reorder for save
names(befm.moi.annotated)
new.column.order <- c( "MOI","in.tri","in.quad","top.in.tri","top.in.quad","global.delta",
                       "LIGAND.CATEGORY","RECEPTOR.CATEGORY","LIGAND.FAMILY","RECEPTOR.FAMILY",
                       "is.gained","is.lost","always.present","top.conserved",
                       "MECHANISM","LIGAND","RECEPTOR","QUERY.SET","DEFINITION.LEVEL","topology.conservation")
befm.moi.annotated <- befm.moi.annotated[,new.column.order]

#### 2024-04-15 MSBR Updating the Feature Annotations ####

# Load local annotation function, in development
source('/Users/msbr/GitHub/NICHESMethods/ToolBuilding/RuleSetFunction.R')

# Load human FANTOM5 database, which is what I'm using to do annotations, because gene symbols are better known
fantom5 <- NICHES::ncomms8866$Pair.Name

# Annotate
fantom5.annotations <- RuleSetFunction(query.set = fantom5,
                                       query.set.name = 'FANTOM5',
                                       special.separating.character = "_")
# Load rat FANTOM5
fantom5.rat <- NICHES::ncomms8866_rat

# Because they are still index aligned, we can simply transfer the labels
fantom5.rat.annotations <- fantom5.annotations
fantom5.rat.annotations$LIGAND <- fantom5.rat$Ligand.ApprovedSymbol
fantom5.rat.annotations$RECEPTOR <- fantom5.rat$Receptor.ApprovedSymbol
fantom5.rat.annotations$MECHANISM <- paste(fantom5.rat.annotations$LIGAND,fantom5.rat.annotations$RECEPTOR,sep = '—') # NICHES separator, will align with marker lists this way
rownames(fantom5.rat.annotations) <- fantom5.rat.annotations$MECHANISM

# Update befm.moi.annotations
befm.moi.annotated$LIGAND.FAMILY <- fantom5.rat.annotations[befm.moi.annotated$MOI,]$LIGAND.FAMILY
befm.moi.annotated$RECEPTOR.FAMILY <- fantom5.rat.annotations[befm.moi.annotated$MOI,]$RECEPTOR.FAMILY
befm.moi.annotated$LIGAND.CATEGORY <- fantom5.rat.annotations[befm.moi.annotated$MOI,]$LIGAND.CATEGORY
befm.moi.annotated$RECEPTOR.CATEGORY <- fantom5.rat.annotations[befm.moi.annotated$MOI,]$RECEPTOR.CATEGORY

# View to fill
View(befm.moi.annotated[is.na(befm.moi.annotated$LIGAND.CATEGORY),])

#########
# Organizing for Allie 2024-04-03
# Break into topologic patterns
markers.gained <- befm.moi.annotated[befm.moi.annotated$is.gained==TRUE,]
markers.lost <- befm.moi.annotated[befm.moi.annotated$is.lost==TRUE,]
markers.moved <- befm.moi.annotated[befm.moi.annotated$always.present==TRUE &
                                          befm.moi.annotated$top.conserved==FALSE,]
markers.conserved <- befm.moi.annotated[befm.moi.annotated$always.present==TRUE &
                                          befm.moi.annotated$top.conserved==TRUE,]
# Break into categorical pieces
# Gained
growth.factors.gained <- markers.gained[markers.gained$LIGAND.CATEGORY%in%'Growth Factor',]
cytokines.gained <- markers.gained[markers.gained$LIGAND.CATEGORY%in%'Cytokine',]
matrix.gained <- markers.gained[markers.gained$LIGAND.CATEGORY%in%'Matrix',]
spatial.gained <- markers.gained[markers.gained$LIGAND.CATEGORY%in%'Spatial Guidance',]
misc.gained <- markers.gained[!(markers.gained$LIGAND.CATEGORY%in%'Matrix') &
                                !(markers.gained$LIGAND.CATEGORY%in%'Cytokine') &
                                !(markers.gained$LIGAND.CATEGORY%in%'Growth Factor') &
                                !(markers.gained$LIGAND.CATEGORY%in%'Spatial Guidance'),]
# Lost
growth.factors.lost <- markers.lost[markers.lost$LIGAND.CATEGORY%in%'Growth Factor',]
cytokines.lost <- markers.lost[markers.lost$LIGAND.CATEGORY%in%'Cytokine',]
matrix.lost <- markers.lost[markers.lost$LIGAND.CATEGORY%in%'Matrix',]
spatial.lost <- markers.lost[markers.lost$LIGAND.CATEGORY%in%'Spatial Guidance',]
misc.lost <- markers.lost[!(markers.lost$LIGAND.CATEGORY%in%'Matrix') &
                                !(markers.lost$LIGAND.CATEGORY%in%'Cytokine') &
                                !(markers.lost$LIGAND.CATEGORY%in%'Growth Factor') &
                                !(markers.lost$LIGAND.CATEGORY%in%'Spatial Guidance'),]

# Moved
growth.factors.moved <- markers.moved[markers.moved$LIGAND.CATEGORY%in%'Growth Factor',]
cytokines.moved <- markers.moved[markers.moved$LIGAND.CATEGORY%in%'Cytokine',]
matrix.moved <- markers.moved[markers.moved$LIGAND.CATEGORY%in%'Matrix',]
spatial.moved <- markers.moved[markers.moved$LIGAND.CATEGORY%in%'Spatial Guidance',]
misc.moved <- markers.moved[!(markers.moved$LIGAND.CATEGORY%in%'Matrix') &
                            !(markers.moved$LIGAND.CATEGORY%in%'Cytokine') &
                            !(markers.moved$LIGAND.CATEGORY%in%'Growth Factor') &
                            !(markers.moved$LIGAND.CATEGORY%in%'Spatial Guidance'),]

# Conserved
growth.factors.conserved <- markers.conserved[markers.conserved$LIGAND.CATEGORY%in%'Growth Factor',]
cytokines.conserved <- markers.conserved[markers.conserved$LIGAND.CATEGORY%in%'Cytokine',]
matrix.conserved <- markers.conserved[markers.conserved$LIGAND.CATEGORY%in%'Matrix',]
spatial.conserved <- markers.conserved[markers.conserved$LIGAND.CATEGORY%in%'Spatial Guidance',]
misc.conserved <- markers.conserved[!(markers.conserved$LIGAND.CATEGORY%in%'Matrix') &
                              !(markers.conserved$LIGAND.CATEGORY%in%'Cytokine') &
                              !(markers.conserved$LIGAND.CATEGORY%in%'Growth Factor') &
                              !(markers.conserved$LIGAND.CATEGORY%in%'Spatial Guidance'),]

# Quick test
VlnPlot(global,'Tgfb1—Tgfbr1',group.by = 'Vector',split.by = 'Condition',pt.size = 0.1) # GOOD!
VlnPlot(global,'Tgfb1—Itgb6',group.by = 'Vector',split.by = 'Condition',pt.size = 0.1) # GOOD!
VlnPlot(global,'Cthrc1—Fzd5',group.by = 'Vector',split.by = 'Condition',pt.size = 0.1) # GOOD!
VlnPlot(global,'Csf2—Csf1r',group.by = 'Vector',split.by = 'Condition',pt.size = 0.1) # GOOD!
VlnPlot(global,'Efna3—Epha1',group.by = 'Vector',split.by = 'Condition',pt.size = 0.1) # GOOD!
VlnPlot(global,'Vegfa—Kdr',group.by = 'Vector',split.by = 'Condition',pt.size = 0.1) # GOOD!
VlnPlot(global,'Wnt7a—Lrp6',group.by = 'Vector',split.by = 'Condition',pt.size = 0) # GOOD!
VlnPlot(global,'Igf1—Insr',group.by = 'Vector',split.by = 'Condition',pt.size = 0) # GOOD!

# Save for later
save(befm.moi.annotated,file='befm.moi.annotated.Robj')

save(markers.gained,file='markers.gained.Robj')
save(markers.lost,file='markers.lost.Robj')
save(markers.conserved,file='markers.conserved.Robj')
save(markers.moved,file='markers.moved.Robj')

save(growth.factors.gained,file='growth.factors.gained.Robj')
save(growth.factors.lost,file='growth.factors.lost.Robj')
save(growth.factors.conserved,file='growth.factors.conserved.Robj')
save(growth.factors.moved,file='growth.factors.moved.Robj')

save(matrix.gained,file='matrix.gained.Robj')
save(matrix.lost,file='matrix.lost.Robj')
save(matrix.conserved,file='matrix.conserved.Robj')
save(matrix.moved,file='matrix.moved.Robj')

save(cytokines.gained,file='cytokines.gained.Robj')
save(cytokines.lost,file='cytokines.lost.Robj')
save(cytokines.conserved,file='cytokines.conserved.Robj')
save(cytokines.moved,file='cytokines.moved.Robj')

save(spatial.gained,file='spatial.gained.Robj')
save(spatial.lost,file='spatial.lost.Robj')
save(spatial.conserved,file='spatial.conserved.Robj')
save(spatial.moved,file='spatial.moved.Robj')

save(misc.gained,file='misc.gained.Robj')
save(misc.lost,file='misc.lost.Robj')
save(misc.conserved,file='misc.conserved.Robj')
save(misc.moved,file='misc.moved.Robj')

# convert to excel
# convert special character
befm.moi.annotated.formatted <- apply(befm.moi.annotated, 2, function(x) gsub("—", " - ", x))
write.csv(befm.moi.annotated.formatted,file='befm.moi.annotated.formatted.csv')
write.table(befm.moi.annotated.formatted,file='befm.moi.annotated.formatted.txt',sep='\t',col.names = NA)

# 2024-04-15

# Count the number of features in each category

num.matrix.conserved <- nrow(matrix.conserved)
num.matrix.gained <- nrow(matrix.gained)
num.matrix.lost <- nrow(matrix.lost)
num.matrix.moved <- nrow(matrix.moved)

num.cytokines.conserved <- nrow(cytokines.conserved)
num.cytokines.gained <- nrow(cytokines.gained)
num.cytokines.lost <- nrow(cytokines.lost)
num.cytokines.moved <- nrow(cytokines.moved)

num.growth.factors.conserved <- nrow(growth.factors.conserved)
num.growth.factors.gained <- nrow(growth.factors.gained)
num.growth.factors.lost <- nrow(growth.factors.lost)
num.growth.factors.moved <- nrow(growth.factors.moved)

num.spatial.conserved <- nrow(spatial.conserved)
num.spatial.gained <- nrow(spatial.gained)
num.spatial.lost <- nrow(spatial.lost)
num.spatial.moved <- nrow(spatial.moved)

num.misc.conserved <- nrow(misc.conserved)
num.misc.gained <- nrow(misc.gained)
num.misc.lost <- nrow(misc.lost)
num.misc.moved <- nrow(misc.moved)

# Make dataframe for ggplotting
data <- as.data.frame(t(data.frame(num.matrix.conserved=num.matrix.conserved,
           num.matrix.gained=num.matrix.gained,
           num.matrix.lost=num.matrix.lost,
           num.matrix.moved=num.matrix.moved,
           
           num.misc.conserved=num.misc.conserved,
           num.misc.gained=num.misc.gained,
           num.misc.lost=num.misc.lost,
           num.misc.moved=num.misc.moved,
           
           num.cytokines.conserved=num.cytokines.conserved,
           num.cytokines.gained=num.cytokines.gained,
           num.cytokines.lost=num.cytokines.lost,
           num.cytokines.moved=num.cytokines.moved,
           
           num.spatial.conserved=num.spatial.conserved,
           num.spatial.gained=num.spatial.gained,
           num.spatial.lost=num.spatial.lost,
           num.spatial.moved=num.spatial.moved,
           
           num.growth.factors.conserved=num.growth.factors.conserved,
           num.growth.factors.gained=num.growth.factors.gained,
           num.growth.factors.lost=num.growth.factors.lost,
           num.growth.factors.moved=num.growth.factors.moved)))
# Add category variable
data$ID <- rownames(data)

# Add class
data$Class <- NA
data[data$ID=='num.cytokines.conserved',]$Class <- 'Cytokines'
data[data$ID=='num.cytokines.lost',]$Class <- 'Cytokines'
data[data$ID=='num.cytokines.gained',]$Class <- 'Cytokines'
data[data$ID=='num.cytokines.moved',]$Class <- 'Cytokines'

data[data$ID=='num.matrix.conserved',]$Class <- 'Matrix'
data[data$ID=='num.matrix.lost',]$Class <- 'Matrix'
data[data$ID=='num.matrix.gained',]$Class <- 'Matrix'
data[data$ID=='num.matrix.moved',]$Class <- 'Matrix'

data[data$ID=='num.misc.conserved',]$Class <- 'Everything Else'
data[data$ID=='num.misc.lost',]$Class <- 'Everything Else'
data[data$ID=='num.misc.gained',]$Class <- 'Everything Else'
data[data$ID=='num.misc.moved',]$Class <- 'Everything Else'

data[data$ID=='num.growth.factors.conserved',]$Class <- 'Growth Factors'
data[data$ID=='num.growth.factors.lost',]$Class <- 'Growth Factors'
data[data$ID=='num.growth.factors.gained',]$Class <- 'Growth Factors'
data[data$ID=='num.growth.factors.moved',]$Class <- 'Growth Factors'

data[data$ID=='num.spatial.conserved',]$Class <- 'Spatial Guidance'
data[data$ID=='num.spatial.lost',]$Class <- 'Spatial Guidance'
data[data$ID=='num.spatial.gained',]$Class <- 'Spatial Guidance'
data[data$ID=='num.spatial.moved',]$Class <- 'Spatial Guidance'

# Add Category
data$Category <- NA
data[data$ID=='num.cytokines.conserved',]$Category <- 'Conserved'
data[data$ID=='num.cytokines.lost',]$Category <- 'Lost'
data[data$ID=='num.cytokines.gained',]$Category <- 'Gained'
data[data$ID=='num.cytokines.moved',]$Category <- 'Moved'

data[data$ID=='num.matrix.conserved',]$Category <- 'Conserved'
data[data$ID=='num.matrix.lost',]$Category <- 'Lost'
data[data$ID=='num.matrix.gained',]$Category <- 'Gained'
data[data$ID=='num.matrix.moved',]$Category <- 'Moved'

data[data$ID=='num.misc.conserved',]$Category <- 'Conserved'
data[data$ID=='num.misc.lost',]$Category <- 'Lost'
data[data$ID=='num.misc.gained',]$Category <- 'Gained'
data[data$ID=='num.misc.moved',]$Category <- 'Moved'

data[data$ID=='num.growth.factors.conserved',]$Category <- 'Conserved'
data[data$ID=='num.growth.factors.lost',]$Category <- 'Lost'
data[data$ID=='num.growth.factors.gained',]$Category <- 'Gained'
data[data$ID=='num.growth.factors.moved',]$Category <- 'Moved'

data[data$ID=='num.spatial.conserved',]$Category <- 'Conserved'
data[data$ID=='num.spatial.lost',]$Category <- 'Lost'
data[data$ID=='num.spatial.gained',]$Category <- 'Gained'
data[data$ID=='num.spatial.moved',]$Category <- 'Moved'

# See https://stackoverflow.com/questions/72411477/how-to-reorder-bars-by-multiple-variables
to.plot <- data %>% arrange(desc(across(.cols=c("Class", "V1")))) %>%
  rowid_to_column()

totals <- to.plot %>%
  group_by(Class) %>%
  summarize(total = sum(V1))

  ggplot(data = to.plot,
         aes(x =  reorder(ID, rowid),  
             y = V1, 
             fill = Category))+
  geom_bar(stat = "identity",  width = 0.8)+
  theme_minimal()+
  ylab('Count')+
  xlab(NULL)+
  #scale_x_discrete(labels = reorder(to.plot$Class,to.plot$rowid))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ggtitle('NICHES Signficant Findings')


####### Make heatmaps for publication #######
  
# Make some epic heatmaps for publication
source("~/GitHub/NICHESMethods/ToolBuilding/CustomHeatmapOnly3.R")
source("~/GitHub/NICHESMethods/ToolBuilding/CustomHeatmap.R")
  source("~/GitHub/NICHESMethods/ToolBuilding/CustomHeatmapDendrogram.R")

  # Identify metadat slots of interest for organization  
names(global@meta.data)


######### COlOR WORK #########
## Color palette from Allie
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Allie's Object Filtration/Final Objects/all.color_palettes_2.R")
## 5. Global connectomics to make full color palette (inefficient but works)
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Saved Objects/global.connectomics.2023-11-11.Robj")
# Add a temporary color palette allowing SendingType visualization
color_pals$type_colors <- c('#A40606','#9CFFFA','#B0DB43','#9C528B','#2F6690',
                            '#946846','#F1C40F','green','#0F0326','#E65F5C','#14591D','#726DA8',
                            'yellow','purple','blue','red','orange','darkgrey','magenta')
names(color_pals$type_colors) <- unique(global.connectomics$CellToCell$alra$SendingType)

# Modify color palette to have a slot that aligns exactly with metadata handles we will use below
color_pals$class <- color_pals$class_colors
color_pals$Condition <- c("#40A6D1","#08306B")
names(color_pals$Condition) <- c('Tri_L','Quad_E')

########

DefaultAssay(global)
global <- ScaleData(global)

global$class.Sending <- factor(global$class.Sending,
                                levels = c('Epithelium','Endothelium','Mesenchyme','Immune'))

global$class.Receiving <- factor(global$class.Receiving,
                                levels = c('Epithelium','Endothelium','Mesenchyme','Immune'))

Idents(global) <- global$class.Receiving

# Pull each niche
epi.rec <- subset(global,idents = 'Epithelium')
end.rec <- subset(global,idents = 'Endothelium')
mes.rec <- subset(global,idents = 'Mesenchyme')
imm.rec <- subset(global,idents = 'Immune')

# Re scale each
epi.rec <- ScaleData(epi.rec)
end.rec <- ScaleData(end.rec)
mes.rec <- ScaleData(mes.rec)
imm.rec <- ScaleData(imm.rec)

 ## Epi
# Extract epi specific lists
epi.gained <- markers.gained[grep('- Epithelium',markers.gained$top.in.quad),]
epi.conserved <- markers.conserved[grep('- Epithelium',markers.conserved$top.in.quad),]
epi.lost <- markers.lost[grep('- Epithelium',markers.lost$top.in.tri),]
epi.moved <- markers.moved[grep('- Epithelium',markers.moved$top.in.quad),]

# Set a threshold for plotting?
epi.moved <- epi.moved[epi.moved$global.delta>2,]
epi.gained <- epi.gained[epi.gained$global.delta>2,]
epi.conserved <- epi.conserved[epi.conserved$global.delta>2,]
epi.lost <- epi.lost[epi.lost$global.delta>2,]

# Define rows
moi <- c(epi.lost$MOI,
         epi.conserved$MOI,
         epi.moved$MOI,
         epi.gained$MOI)

png(file = 'epi.test.png',width = 10,height = 10,units = 'in',res=300)
CustomHeatmapDendrogram(object = epi.rec,
                   data.type = 'CellToCell',
                   primary = 'Condition' ,
                   secondary = 'class.Sending',
                   tertiary = 'class.Receiving' ,
                   #quarternary = 'orig.ident' ,
                   primary.cols = color_pals$Condition,
                   secondary.cols = color_pals$class_colors, # Need to be a named list of colors
                   tertiary.cols = color_pals$class_colors,
                   #quarternary.cols = NULL,
                   features = moi,
                   labels = c('Condition','Sending Class','Receiving Class'),
                   selected.row.anotations=NULL,#row.annotations,
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.7)
dev.off()

## Endo
# Extract endo specific lists
end.gained <- markers.gained[grep('- Endothelium',markers.gained$top.in.quad),]
end.conserved <- markers.conserved[grep('- Endothelium',markers.conserved$top.in.quad),]
end.lost <- markers.lost[grep('- Endothelium',markers.lost$top.in.tri),]
end.moved <- markers.moved[grep('- Endothelium',markers.moved$top.in.quad),]

# Set a threshold for plotting?
end.moved <- end.moved[end.moved$global.delta>2,]
end.gained <- end.gained[end.gained$global.delta>2,]
end.conserved <- end.conserved[end.conserved$global.delta>2,]
end.lost <- end.lost[end.lost$global.delta>2,]

# Define rows
moi <- c(end.lost$MOI,
         end.conserved$MOI,
         end.moved$MOI,
         end.gained$MOI)

png(file = 'end.test.png',width = 10,height = 10,units = 'in',res=300)
CustomHeatmapDendrogram(object = end.rec,
                        data.type = 'CellToCell',
                        primary = 'Condition' ,
                        secondary = 'class.Sending',
                        tertiary = 'class.Receiving' ,
                        #quarternary = 'orig.ident' ,
                        primary.cols = color_pals$Condition,
                        secondary.cols = color_pals$class_colors, # Need to be a named list of colors
                        tertiary.cols = color_pals$class_colors,
                        #quarternary.cols = NULL,
                        features = moi,
                        labels = c('Condition','Sending Class','Receiving Class'),
                        selected.row.anotations=NULL,#row.annotations,
                        selected.label.size = 10,
                        use.scale.data = T,
                        range.frac = 0.7)
dev.off()

## Mes
# Extract mes specific lists
mes.gained <- markers.gained[grep('- Mesenchyme',markers.gained$top.in.quad),]
mes.conserved <- markers.conserved[grep('- Mesenchyme',markers.conserved$top.in.quad),]
mes.lost <- markers.lost[grep('- Mesenchyme',markers.lost$top.in.tri),]
mes.moved <- markers.moved[grep('- Mesenchyme',markers.moved$top.in.quad),]

# Set a threshold for plotting?
mes.moved <- mes.moved[mes.moved$global.delta>2,]
mes.gained <- mes.gained[mes.gained$global.delta>2,]
mes.conserved <- mes.conserved[mes.conserved$global.delta>2,]
mes.lost <- mes.lost[mes.lost$global.delta>2,]

# Define rows
moi <- c(mes.lost$MOI,
         mes.conserved$MOI,
         mes.moved$MOI,
         mes.gained$MOI)

png(file = 'mes.test.png',width = 10,height = 10,units = 'in',res=300)
CustomHeatmapDendrogram(object = mes.rec,
                        data.type = 'CellToCell',
                        primary = 'Condition' ,
                        secondary = 'class.Sending',
                        tertiary = 'class.Receiving' ,
                        #quarternary = 'orig.ident' ,
                        primary.cols = color_pals$Condition,
                        secondary.cols = color_pals$class_colors, # Need to be a named list of colors
                        tertiary.cols = color_pals$class_colors,
                        #quarternary.cols = NULL,
                        features = moi,
                        labels = c('Condition','Sending Class','Receiving Class'),
                        selected.row.anotations=NULL,#row.annotations,
                        selected.label.size = 10,
                        use.scale.data = T,
                        range.frac = 0.7)
dev.off()

## Imm
# Extract imm specific lists
imm.gained <- markers.gained[grep('- Immune',markers.gained$top.in.quad),]
imm.conserved <- markers.conserved[grep('- Immune',markers.conserved$top.in.quad),]
imm.lost <- markers.lost[grep('- Immune',markers.lost$top.in.tri),]
imm.moved <- markers.moved[grep('- Immune',markers.moved$top.in.quad),]

# Set a threshold for plotting?
imm.moved <- imm.moved[imm.moved$global.delta>2,]
imm.gained <- imm.gained[imm.gained$global.delta>2,]
imm.conserved <- imm.conserved[imm.conserved$global.delta>2,]
imm.lost <- imm.lost[imm.lost$global.delta>2,]

# Define rows
moi <- c(imm.lost$MOI,
         imm.conserved$MOI,
         imm.moved$MOI,
         imm.gained$MOI)

png(file = 'imm.test.png',width = 10,height = 10,units = 'in',res=300)
CustomHeatmapDendrogram(object = imm.rec,
                        data.type = 'CellToCell',
                        primary = 'Condition' ,
                        secondary = 'class.Sending',
                        tertiary = 'class.Receiving' ,
                        #quarternary = 'orig.ident' ,
                        primary.cols = color_pals$Condition,
                        secondary.cols = color_pals$class_colors, # Need to be a named list of colors
                        tertiary.cols = color_pals$class_colors,
                        #quarternary.cols = NULL,
                        features = moi,
                        labels = c('Condition','Sending Class','Receiving Class'),
                        selected.row.anotations=NULL,#row.annotations,
                        selected.label.size = 10,
                        use.scale.data = T,
                        range.frac = 0.7)
dev.off()

## Remove integrins?
no.int <- befm.moi.annotated[!(befm.moi.annotated$RECEPTOR.CATEGORY%in%'Integrin'),]

## Apply threshold for plotting (ugh)
no.int.thresh <- no.int[no.int$global.delta>2,]
no.int.thresh.for.epi <- no.int[no.int$global.delta>1.5,] # preserves more things of interest
### Epi
## Identify things with top vector marker landing on epithelium, in either tri or quad
epi.indices <- unique(union(grep('- Epithelium',no.int.thresh.for.epi$top.in.tri),
                            grep('- Epithelium',no.int.thresh.for.epi$top.in.quad)))
epi.mark <- no.int.thresh.for.epi[epi.indices,]

## Epi heatmap gen 2
epi.moi <- epi.mark$MOI

png(file = 'epi.test.2.png',width = 10,height = 10,units = 'in',res=300)
CustomHeatmapDendrogram(object = epi.rec,
                        data.type = 'CellToCell',
                        primary = 'Condition' ,
                        secondary = 'class.Sending',
                        tertiary = 'class.Receiving' ,
                        #quarternary = 'orig.ident' ,
                        primary.cols = color_pals$Condition,
                        secondary.cols = color_pals$class_colors, # Need to be a named list of colors
                        tertiary.cols = color_pals$class_colors,
                        #quarternary.cols = NULL,
                        features = epi.moi,
                        labels = c('Condition','Sending Class','Receiving Class'),
                        selected.row.anotations=NULL,#row.annotations,
                        selected.label.size = 10,
                        use.scale.data = T,
                        range.frac = 0.5)
dev.off()

### End
## Identify things with top vector marker landing on endothelium, in either tri or quad
end.indices <- unique(union(grep('- Endothelium',no.int.thresh$top.in.tri),
                            grep('- Endothelium',no.int.thresh$top.in.quad)))
end.mark <- no.int.thresh[end.indices,]

## End heatmap gen 2
end.moi <- end.mark$MOI
png(file = 'end.test.2.png',width = 12,height = 12,units = 'in',res=300)
CustomHeatmapDendrogram(object = end.rec,
                        data.type = 'CellToCell',
                        primary = 'Condition' ,
                        secondary = 'class.Sending',
                        tertiary = 'class.Receiving' ,
                        #quarternary = 'orig.ident' ,
                        primary.cols = color_pals$Condition,
                        secondary.cols = color_pals$class_colors, # Need to be a named list of colors
                        tertiary.cols = color_pals$class_colors,
                        #quarternary.cols = NULL,
                        features = end.moi,
                        labels = c('Condition','Sending Class','Receiving Class'),
                        selected.row.anotations=NULL,#row.annotations,
                        selected.label.size = 10,
                        use.scale.data = T,
                        range.frac = 0.5)
dev.off()

### Mes
## Identify things with top vector marker landing on mesenchyme, in either tri or quad
mes.indices <- unique(union(grep('- Mesenchyme',no.int.thresh$top.in.tri),
                            grep('- Mesenchyme',no.int.thresh$top.in.quad)))
mes.mark <- no.int.thresh[mes.indices,]

## Mes heatmap gen 2
mes.moi <- mes.mark$MOI

png(file = 'mes.test.2.png',width = 12,height = 12,units = 'in',res=300)
CustomHeatmapDendrogram(object = mes.rec,
                        data.type = 'CellToCell',
                        primary = 'Condition' ,
                        secondary = 'class.Sending',
                        tertiary = 'class.Receiving' ,
                        #quarternary = 'orig.ident' ,
                        primary.cols = color_pals$Condition,
                        secondary.cols = color_pals$class_colors, # Need to be a named list of colors
                        tertiary.cols = color_pals$class_colors,
                        #quarternary.cols = NULL,
                        features = mes.moi,
                        labels = c('Condition','Sending Class','Receiving Class'),
                        selected.row.anotations=NULL,#row.annotations,
                        selected.label.size = 10,
                        use.scale.data = T,
                        range.frac = 0.5)
dev.off()

### Imm
## Identify things with top vector marker landing on immune, in either tri or quad
imm.indices <- unique(union(grep('- Immune',no.int.thresh$top.in.tri),
                            grep('- Immune',no.int.thresh$top.in.quad)))
imm.mark <- no.int.thresh[imm.indices,]

## Imm heatmap gen 2
imm.moi <- imm.mark$MOI

png(file = 'imm.test.2.png',width = 12,height = 12,units = 'in',res=300)
CustomHeatmapDendrogram(object = imm.rec,
                        data.type = 'CellToCell',
                        primary = 'Condition' ,
                        secondary = 'class.Sending',
                        tertiary = 'class.Receiving' ,
                        #quarternary = 'orig.ident' ,
                        primary.cols = color_pals$Condition,
                        secondary.cols = color_pals$class_colors, # Need to be a named list of colors
                        tertiary.cols = color_pals$class_colors,
                        #quarternary.cols = NULL,
                        features = imm.moi,
                        labels = c('Condition','Sending Class','Receiving Class'),
                        selected.row.anotations=NULL,#row.annotations,
                        selected.label.size = 10,
                        use.scale.data = T,
                        range.frac = 0.5)
dev.off()

 ### GEn 3 sending first, then receving, then condition

## Epi heatmap gen 3
png(file = 'epi.test.3.png',width = 10,height = 10,units = 'in',res=300)
CustomHeatmapDendrogram(object = epi.rec,
                        data.type = 'CellToCell',
                        primary = 'class.Sending' ,
                        secondary = 'Condition',
                        tertiary = 'class.Receiving' ,
                        #quarternary = 'orig.ident' ,
                        primary.cols = color_pals$class_colors,
                        secondary.cols = color_pals$Condition, # Need to be a named list of colors
                        tertiary.cols = color_pals$class_colors,
                        #quarternary.cols = NULL,
                        features = epi.moi,
                        labels = c('Sending Class','Condition','Receiving Class'),
                        selected.row.anotations=NULL,#row.annotations,
                        selected.label.size = 10,
                        use.scale.data = T,
                        range.frac = 0.5)
dev.off()

## End heatmap gen 3
png(file = 'end.test.3.png',width = 10,height = 10,units = 'in',res=300)
CustomHeatmapDendrogram(object = end.rec,
                        data.type = 'CellToCell',
                        primary = 'class.Sending' ,
                        secondary = 'Condition',
                        tertiary = 'class.Receiving' ,
                        #quarternary = 'orig.ident' ,
                        primary.cols = color_pals$class_colors,
                        secondary.cols = color_pals$Condition, # Need to be a named list of colors
                        tertiary.cols = color_pals$class_colors,
                        #quarternary.cols = NULL,
                        features = end.moi,
                        labels = c('Sending Class','Condition','Receiving Class'),
                        selected.row.anotations=NULL,#row.annotations,
                        selected.label.size = 10,
                        use.scale.data = T,
                        range.frac = 0.5)
dev.off()


## Mes heatmap gen 3
png(file = 'mes.test.3.png',width = 10,height = 10,units = 'in',res=300)
CustomHeatmapDendrogram(object = mes.rec,
                        data.type = 'CellToCell',
                        primary = 'class.Sending' ,
                        secondary = 'Condition',
                        tertiary = 'class.Receiving' ,
                        #quarternary = 'orig.ident' ,
                        primary.cols = color_pals$class_colors,
                        secondary.cols = color_pals$Condition, # Need to be a named list of colors
                        tertiary.cols = color_pals$class_colors,
                        #quarternary.cols = NULL,
                        features = mes.moi,
                        labels = c('Sending Class','Condition','Receiving Class'),
                        selected.row.anotations=NULL,#row.annotations,
                        selected.label.size = 10,
                        use.scale.data = T,
                        range.frac = 0.5)
dev.off()

## Imm heatmap gen 3
png(file = 'imm.test.3.png',width = 10,height = 10,units = 'in',res=300)
CustomHeatmapDendrogram(object = imm.rec,
                        data.type = 'CellToCell',
                        primary = 'class.Sending' ,
                        secondary = 'Condition',
                        tertiary = 'class.Receiving' ,
                        #quarternary = 'orig.ident' ,
                        primary.cols = color_pals$class_colors,
                        secondary.cols = color_pals$Condition, # Need to be a named list of colors
                        tertiary.cols = color_pals$class_colors,
                        #quarternary.cols = NULL,
                        features = imm.moi,
                        labels = c('Sending Class','Condition','Receiving Class'),
                        selected.row.anotations=NULL,#row.annotations,
                        selected.label.size = 10,
                        use.scale.data = T,
                        range.frac = 0.5)
dev.off()

VlnPlot(global,'Hbegf—Egfr',split.by = 'Condition',group.by = 'Vector')


### Gen 4: split columns by sending class, tight range.frac, custom row order

#### Epi heatmap gen 4 ####
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics")

# subset to flattened number per vector
Idents(epi.rec) <- epi.rec$Vector
table(Idents(epi.rec))

epi.rec.sub <- subset(epi.rec,downsample = 4160)
table(Idents(epi.rec.sub))

# make preliminary, to get order
png(file = 'epi.test.4.png',width = 10,height = 10,units = 'in',res=300)
CustomHeatmapDendrogram(object = epi.rec.sub,
                        data.type = 'CellToCell',
                        primary = 'class.Sending' ,
                        secondary = 'class.Receiving',
                        tertiary = 'Condition' ,
                        #quarternary = 'orig.ident' ,
                        primary.cols = color_pals$class_colors,
                        secondary.cols = color_pals$class_colors, # Need to be a named list of colors
                        tertiary.cols = color_pals$Condition,
                        #quarternary.cols = NULL,
                        features = epi.moi,
                        labels = c('Sending Class','Receiving Class','Condition'),
                        selected.row.anotations=NULL,#row.annotations,
                        selected.label.size = 10,
                        use.scale.data = T,
                        range.frac = 0.3)
dev.off()
epi.moi.ordered <- epi.mark[row.order.output,]$MOI # this pulls the row order from the above plot, for modification below:

epi.moi.ordered <- c( "Efna3—Epha1","Il24—Il20ra","Ntf4—Ngfr", 
                      "Efna3—Epha7","Tgfa—Egfr",  
 "Shbg—Cd177", "Plau—St14",  "Tnf—Tnfrsf21" , 
  "Cxcl11—Ackr3"  , "Sct—Sctr",   "Hbegf—Cd9",  "Adipoq—Adipor2" ,"Adipoq—Adipor1",   
"Hbegf—Egfr", "Ereg—Egfr",  "Cdh1—Igf1r", "Cdh1—Egfr", 
  "Areg—Egfr",  "C3—Cd81","S100a8—Tlr4","S100a9—Tlr4",
  
  "Sema4c—Plxnb2" , "Tnfsf4—Traf2"  , "Igfbp4—Lrp6" ,   "Efna1—Epha1"  ,  "Tnfsf4—Tnfrsf4" ,"Mfng—Notch1"   ,
                     "Dll4—Notch1","Dll1—Notch1","Efna1—Ephb6","Efna1—Epha7","Vwf—Sirpa",  "Ccl5—Sdc1", 
 "Ccl5—Sdc4","Bmp3—Bmpr1b","Bmp4—Bmpr1b","Bmp6—Bmpr1b","Rgmb—Bmpr1b",
 
  
  "Dlk1—Notch1","Slit1—Sdc1", "Mst1—Mst1r", "Bdnf—Ngfr","Ptn—Sdc1",   "Ptn—Plxnb2",
"Fn1—Il17rc",  "Col1a1—Cd44","Rtn4—Gjb2", 
 "Mdk—Tspan1", "C3—C5ar2",  
  "Vegfa—Egfr", "Tg—Asgr1","Efnb2—Ephb6","Efnb1—Ephb6", 
  "Efnb2—Ephb2","Fgf2—Fgfr3", "Mmp9—Ephb2", "Il6—F3", "Vegfa—Sirpa","Fgf1—Fgfr2", "Fgf1—Fgfr3",
  "Fgf1—Egfr",  


  "Amh—Egfr",
"Nrg4—Egfr",  "Sema4d—Plxnb2" , "Sema4d—Erbb2"  , "Sema4d—Plxnb1"  ,"Il18—Il18r1","Sema4d—Met",
"Nrg4—Erbb2", "Tnfsf9—Traf2"  , "Il1b—Adrb2","Il1b—Il1r2", "Pf4—Fgfr2",  "Psen1—Notch3",   "Psen1—Notch1" , 
 "Efna2—Epha7","Rtn4—Rtn4r", "Psap—Gpr37l1"  , "Psap—Celsr1","Gnai2—Egfr", "Anxa1—Egfr")

# Remove specific mechanisms that do not look like real signal by this lens
edit.out <- c("Efemp2—Lingo1","Hbegf—Cd9",  "Efna5—Epha7","Efna5—Ephb6","Btc—Egfr",  "Efna1—Epha1","Tg—Asgr1","Rtn4—Gjb2","B2m—Cd3g","Calcb—Ramp1",   "Vegfb—Tyro3",
              "Dll3—Notch1","Jag2—Notch3","Hdc—Hrh4","Hsp90b1—Asgr1",  "Ngf—Ngfr", "Gpi—Amfr",    "Sema3f—Plxna1" ,"Hsp90b1—Tlr4","Hsp90aa1—Cftr", "Hsp90b1—Tlr2","Cxcl3—Cxcr2","Cxcl1—Cxcr2")
epi.moi.ordered <- epi.moi.ordered[!(epi.moi.ordered%in% edit.out)]

png(file = 'epi.test.4.png',width = 10,height = 10,units = 'in',res=300)
CustomHeatmapDendrogram(object = epi.rec.sub,
                        data.type = 'CellToCell',
                        primary = 'class.Sending' ,
                        secondary = 'class.Receiving',
                        tertiary = 'Condition' ,
                        #quarternary = 'orig.ident' ,
                        primary.cols = color_pals$class_colors,
                        secondary.cols = color_pals$class_colors, # Need to be a named list of colors
                        tertiary.cols = color_pals$Condition,
                        #quarternary.cols = NULL,
                        features = epi.moi.ordered,
                        labels = c('Sending Class','Receiving Class','Condition'),
                        selected.row.anotations=NULL,#row.annotations,
                        selected.label.size = 10,
                        use.scale.data = T,
                        range.frac = 0.3,
                        row.dendrogram = F)
dev.off()



#### End heatmap gen 4 ####
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics")

# subset to flattened number per vector
Idents(end.rec) <- end.rec$Vector
table(Idents(end.rec))
end.rec.sub <- subset(end.rec,downsample = 4200)

# make preliminary, to get order
png(file = 'end.test.4.png',width = 10,height = 10,units = 'in',res=300)
CustomHeatmapDendrogram(object = end.rec.sub,
                        data.type = 'CellToCell',
                        primary = 'class.Sending' ,
                        secondary = 'Condition',
                        tertiary = 'class.Receiving' ,
                        #quarternary = 'orig.ident' ,
                        primary.cols = color_pals$class_colors,
                        secondary.cols = color_pals$Condition, # Need to be a named list of colors
                        tertiary.cols = color_pals$class_colors,
                        #quarternary.cols = NULL,
                        features = end.moi,
                        labels = c('Sending Class','Condition','Receiving Class'),
                        selected.row.anotations=NULL,#row.annotations,
                        selected.label.size = 10,
                        use.scale.data = T,
                        range.frac = 0.3)
dev.off()
end.moi.ordered <- end.mark[row.order.output,]$MOI # this pulls the row order from the above plot, for modification below:
end.moi.ordered <- c(
  
  "Lamc2—Cd151", "Lamb3—Cd151", "Sema3e—Plxnd1",   
  "Sema3e—Nrp1", "Calcb—Calcrl","Wnt7a—Fzd9",
  "Adm2—Calcrl", "Sema4a—Plxnd1"  ,
  "Vegfa—Nrp2",  "Vegfa—Kdr",    "Lin7c—Abca1", "Dusp18—Cd151","Psen1—Notch4" , 
  "Anxa1—Dysf",  "Cd24—Selp",   "Cgn—Tgfbr2",  "Efna3—Ephb1",

  
 
  "Ecm1—Cachd1", "Sema3f—Nrp1", "Efna1—Ephb1", 
  "Igfbp4—Lrp6", "Sema3g—Nrp2", "Sema6a—Plxna2", 
  "Cubn—Lrp2",   "Tcn2—Lrp2",   "Cd34—Selp",  
  "Angpt2—Tie1", "Fat4—Dchs1",  "Dll4—Notch4",  "Tgfb1—Eng",     
  "Tgfb1—Tgfbr2",  "Tgfb1—Cxcr4",   "Rims2—Abca1","Psap—Sort1",
  "Tnfsf12—Tnfrsf25", "Tg—Lrp2",
  
  "Lrpap1—Lrp2","Apln—Aplnr","Efna5—Ephb1","Efnb1—Ephb1",
                   "Serpine1—Lrp2", "Lpl—Gpihbp1","Lpl—Lrp2","Pltp—Abca1", "Col1a1—Flt4",  
                  "Wnt16—Fzd6", "Col5a1—Sdc3","Col5a3—Sdc3","Fgf2—Sdc3", "Ccl2—Ccr10", "Ccl7—Ccr10",
                   "Eda—Eda2r",  "C4b—Cd46", "Efnb3—Ephb1", 
                  "Mdk—Gpc2","Mdk—Lrp2",    "Mdk—Ptprb",   "Mdk—Sdc3",    "Apoe—Lrp2",   "Apoe—Vldlr", 
                  "Angptl4—Tie1","Dll3—Notch4",
                  
                  "Il7—Il2rg",   "Selplg—Esam", "Selplg—Selp", "Selplg—Sele",  "Mmp12—Plaur",
                   "Agt—Lrp2",    "Tnfsf13—Tnfrsf1a", "Osm—Il6st",   "Ebi3—Il6st",  "Il1b—Il1rap", "Ebi3—Il27ra",
                  "Pf4—Ldlr",    "Pf4—Thbd",    "Pf4—Procr", "Tnf—Tnfrsf1a","Tnf—Traf2") 

# Remove specific mechanisms that do not look like real signal by this lens
end.edit.out <- c("Agrp—Sdc3", "Plau—Lrp2","Sema6a—Plxna4", "Col18a1—Kdr","Ccl27—Ccr10", "Wnt5a—Fzd5","Scgb1a1—Lrp2","Cxcl13—Ccr10","Gpha2—Tshr","Nampt—Insr" ,"C3—Cd46")
end.moi.ordered <- end.moi.ordered[!(end.moi.ordered%in% end.edit.out)]

png(file = 'end.test.4.png',width = 10,height = 10,units = 'in',res=300)
CustomHeatmapDendrogram(object = end.rec.sub,
                        data.type = 'CellToCell',
                        primary = 'class.Sending' ,
                        secondary = 'class.Receiving',
                        tertiary = 'Condition' ,
                        #quarternary = 'orig.ident' ,
                        primary.cols = color_pals$class_colors,
                        secondary.cols = color_pals$class_colors, # Need to be a named list of colors
                        tertiary.cols = color_pals$Condition,
                        #quarternary.cols = NULL,
                        features = end.moi.ordered,
                        labels = c('Sending Class','Receiving Class','Condition'),
                        selected.row.anotations=NULL,#row.annotations,
                        selected.label.size = 10,
                        use.scale.data = T,
                        range.frac = 0.3,
                        row.dendrogram = F)
dev.off()



#### Mes heatmap gen 4 ####
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics")

# subset to flattened number per vector
Idents(mes.rec) <- mes.rec$Vector
table(Idents(mes.rec))
mes.rec.sub <- subset(mes.rec,downsample = 5000)
# make preliminary, to get order
png(file = 'mes.test.4.png',width = 10,height = 10,units = 'in',res=300)
CustomHeatmapDendrogram(object = mes.rec.sub,
                        data.type = 'CellToCell',
                        primary = 'class.Sending' ,
                        secondary = 'Condition',
                        tertiary = 'class.Receiving' ,
                        #quarternary = 'orig.ident' ,
                        primary.cols = color_pals$class_colors,
                        secondary.cols = color_pals$Condition, # Need to be a named list of colors
                        tertiary.cols = color_pals$class_colors,
                        #quarternary.cols = NULL,
                        features = mes.moi,
                        labels = c('Sending Class','Condition','Receiving Class'),
                        selected.row.anotations=NULL,#row.annotations,
                        selected.label.size = 10,
                        use.scale.data = T,
                        range.frac = 0.5)
dev.off()
mes.moi.ordered <- mes.mark[row.order.output,]$MOI # this pulls the row order from the above plot, for modification below:
mes.moi.ordered <- c(
  
  "Hbegf—Cd44",    "Tnf—Ltbr",      "Wnt7a—Lrp6",    "Mmp7—Cd44",       "Sema3e—Nrp1",  
  
 "Rims1—Slc18a3","Psen1—Notch3",  "Dll1—Notch3",   
  "Psap—Lrp1",     "Bmp6—Acvr2a",   "Rgmb—Neo1",     "Nxph1—Nrxn3",   "Il16—Kcna3",   
  "Sema3g—Nrp2",   "Ace—Bdkrb2",    "Vwf—Tnfrsf11b", "Dhh—Cdon",      "Vip—Npr3",     
  "Efna1—Epha3",   "Tac4—Tacr3",   "Dll4—Notch3",  "Jag1—Notch3",  
  
    "Apoe—Lrp1",      "Cxcl12—Sdc4",   "C3—Ifitm1",    
    "Col14a1—Cd44",  "Thbs2—Cd47",   
  "A2m—Lrp1",      "Ccl19—Ackr4",   "Sfrp1—Fzd2",    "Slit1—Gpc1",    "Slit1—Robo1",  
 "Flt3lg—Flt3",   "Il15—Il15ra",   "Fn1—Nt5e",      "Gdf10—Acvr1b",    
 "Fn1—Col13a1",   "Nid1—Col13a1",   "Ccl2—Ackr4",     "Ccl7—Ackr4",   "Il15—Il2ra",    
 
  "Vegfa—Nrp2",    "Vegfb—Nrp1",   
  "Pgf—Nrp1",      "Vegfa—Nrp1",    "P4hb—Gpr162",   "Tf—Gpr162",     
  "Adam10—Epha3",  "Sema3f—Nrp1",
 "Serpine2—Lrp1", "Plat—Lrp1",     "Col1a1—Cd44",   "Ntn1—Unc5b",   
  "Ntn1—Adora2b",  "Ntn1—Neo1",     "Mmp9—Lrp1",   
  "Lama1—Sdc4",    "Lama1—Sdc2",    "Lama1—Gpc1",    "Efna2—Epha3",   "Tnfsf13—Sdc2", 
  "Tnfsf13—Fas",   "Tnfsf13—Tnfrsf11b"  ,  "Nrg4—Egfr",     "Osm—Osmr",      "Il1b—Il1r1",   
 "Osm—Lifr",      "Pf4—Sdc2",      "C1qb—Lrp1",     "C1qa—Cspg4",    "Agt—Agtr2"    
  
)

# Remove specific mechanisms that do not look like real signal by this lens
mes.edit.out <- c( "Csf2—Il3ra",   
                   "Csf2—Sdc2",  
                   "Gnai2—Adcy7", "Hsp90b1—Lrp1","Ncam1—AABR07007068.1", "Gnai2—Adora1", "Gnai2—Ptpru","Artn—AABR07007068.1",  "Cxcl13—Ackr4",  "Nppc—Npr3",    "Ccl11—Ackr4" , "Timp1—Cd63", "Ccl5—Gpr75")
mes.moi.ordered <- mes.moi.ordered[!(mes.moi.ordered%in% mes.edit.out)]

png(file = 'mes.test.4.png',width = 10,height = 10,units = 'in',res=300)
CustomHeatmapDendrogram(object = mes.rec.sub,
                        data.type = 'CellToCell',
                        primary = 'class.Sending' ,
                        secondary = 'class.Receiving',
                        tertiary = 'Condition' ,
                        #quarternary = 'orig.ident' ,
                        primary.cols = color_pals$class_colors,
                        secondary.cols = color_pals$class_colors, # Need to be a named list of colors
                        tertiary.cols = color_pals$Condition,
                        #quarternary.cols = NULL,
                        features = mes.moi.ordered,
                        labels = c('Sending Class','Receiving Class','Condition'),
                        selected.row.anotations=NULL,#row.annotations,
                        selected.label.size = 10,
                        use.scale.data = T,
                        range.frac = 0.3,
                        row.dendrogram = F)
dev.off()

#### Imm heatmap gen 4 ####
# subset to flattened number per vector
Idents(imm.rec) <- imm.rec$Vector
table(Idents(imm.rec))
imm.rec.sub <- subset(imm.rec,downsample = 640)

# make preliminary, to get order
png(file = 'imm.test.4.png',width = 10,height = 10,units = 'in',res=300)
CustomHeatmapDendrogram(object = imm.rec.sub,
                        data.type = 'CellToCell',
                        primary = 'class.Sending' ,
                        secondary = 'Condition',
                        tertiary = 'class.Receiving' ,
                        #quarternary = 'orig.ident' ,
                        primary.cols = color_pals$class_colors,
                        secondary.cols = color_pals$Condition, # Need to be a named list of colors
                        tertiary.cols = color_pals$class_colors,
                        #quarternary.cols = NULL,
                        features = imm.moi,
                        labels = c('Sending Class','Condition','Receiving Class'),
                        selected.row.anotations=NULL,#row.annotations,
                        selected.label.size = 10,
                        use.scale.data = T,
                        range.frac = 0.5)
dev.off()

imm.moi.ordered <- imm.mark[row.order.output,]$MOI # this pulls the row order from the above plot, for modification below:
imm.moi.ordered <- c("Inhbb—Acvr1c", "Il10—Il10ra","Tnf—Tnfrsf1b",  "Cxcl6—Cxcr1",  
  "Cxcl3—Cxcr1",   "Ppbp—Cxcr1",    "Agrn—Atp1a3",   "Adipoq—Adipor1","C3—C5ar2",     
  "C3—C3ar1",      "Csf3—Csf1r",    "Csf3—Csf3r",    "Cxcl6—Cxcr2",   "Cxcl3—Cxcr2",  
  "Slpi—Cd4",     
  
  "Sema6d—Trem2",  "Egfl8—Ncr3",    "Podxl—Sell",    "Ccl5—Ccr1",
  "Cd34—Sell",     "F8—Asgr2",
  "Tslp—Il7r",  "Vwf—Sirpa",      "Sema6d—Tyrobp", "Psap—Sort1",    "Rims2—Abca1",  
 

"Fgf2—Sdc3",    
 "Mmp9—Cd44",     "Csf1—Csf1r",    "Fn1—Nt5e",      "Fn1—C5ar1",     "Timp1—Cd63",   
"Il15—Il15ra",   "Il34—Csf1r",    "Mdk—Sdc3",      "Col5a1—Sdc3",   "Col5a3—Sdc3",  
"Agrp—Sdc3",     "Cxcl1—Cxcr1", "Vcan—Sell",     "Podxl2—Sell",   "Ptgs2—Alox5","Vegfa—Sirpa", "Hsp90b1—Tlr7",  

"Il16—Cd4", "Ly86—Cd180",        "Cfh—Sell",      "Cxcl2—Cxcr2",   "Pltp—Abca1",   
"Il7—Il7r",      "Ccl4—Ccr1",     "Ccl3—Ccr1",     "Il18—Cd48",     "Sema4d—Cd72",  
  "Alox5ap—Alox5", "Camp—P2rx7",   "Cxcl2—Cxcr1" )

# Remove specific mechanisms that do not look like real signal by this lens
imm.edit.out <- c("Lin7c—Htr2c",  "Gnai2—Ccr5",        "Sema7a—Plxnc1",   "Gnai2—Adcy7",     "Il16—Ccr5", "Il23a—Il23r",  "Gnai2—Cxcr1",   "Gnai2—P2ry12",  "Gnai2—C5ar1",  
                  "App—Cd74",      "Ccl2—Ccr5",     "Ccl7—Ccr5",       "Scgb3a2—Marco" ,
                  "Cxcl12—Cd4",    "AABR07064312.1—C5ar1" ,"Csf2—Csf2ra",   "Csf2—Csf2rb",   
                  "Csf2—Csf1r",    "Csf2—Csf3r","Ccl3—Ccr5",     "Ccl4—Ccr5",    
                  "Orm1—Ccr5", "Ccl12—Ccr1","Tnfsf13b—Tnfrsf13c","Selplg—Sell")
imm.moi.ordered <- imm.moi.ordered[!(imm.moi.ordered%in% imm.edit.out)]
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics")

png(file = 'imm.test.4.png',width = 10,height = 10,units = 'in',res=300)
CustomHeatmapDendrogram(object = imm.rec.sub,
                        data.type = 'CellToCell',
                        primary = 'class.Sending' ,
                        secondary = 'class.Receiving',
                        tertiary = 'Condition' ,
                        #quarternary = 'orig.ident' ,
                        primary.cols = color_pals$class_colors,
                        secondary.cols = color_pals$class_colors, # Need to be a named list of colors
                        tertiary.cols = color_pals$Condition,
                        #quarternary.cols = NULL,
                        features = imm.moi.ordered,
                        labels = c('Sending Class','Receiving Class','Condition'),
                        selected.row.anotations=NULL,#row.annotations,
                        selected.label.size = 10,
                        use.scale.data = T,
                        range.frac = 0.3,
                        row.dendrogram = F)
dev.off()

## 2024-04-24

# Make circuit plots aligning with the above heatmap rows

# 0. Load transcriptomic data
load('Saved Objects/transcriptome.TriQuad.customEmbed.fineFilter.proofRead.2024-03-15.Robj') # sub.fine
 
# 1. Load CircuitPlot and dependts
source("~/GitHub/NICHESMethods/ToolBuilding/CircuitPlot.R")
source("~/GitHub/NICHESMethods/ToolBuilding/DefineNodeObject.R")
source("~/GitHub/NICHESMethods/ToolBuilding/DefineEdgeObject.R")
source("~/GitHub/NICHESMethods/ToolBuilding/AggregateNodeData.R")
source("~/GitHub/NICHESMethods/ToolBuilding/AggregateEdgeData.R")
source("~/GitHub/NICHESMethods/ToolBuilding/ggCircuit.R")

# 2. Define the TRI and QUAD datasets (allows splitting, below)
Idents(global) <- global$Condition
table(Idents(global))
tri.connectivity <- subset(global,idents = "Tri_L")
quad.connectivity <- subset(global,idents = "Quad_E")

Idents(sub.fine) <- sub.fine$Condition
table(Idents(sub.fine))
tri.phenotype <- subset(sub.fine,idents = 'Tri_L')
quad.phenotype <- subset(sub.fine, idents = 'Quad_E')


# 2. Define temp plotting function
TempFunc <- function(feature,
                     max.edge.value = 0.1){
  title <- title_theme <- cowplot::ggdraw() +
    cowplot::draw_label(feature, 
               #fontfamily = theme_georgia()$text$family, 
               #fontface = theme_georgia()$plot.title$face, 
               x = 0.05, hjust = 0)
  p4 <- CircuitPlot(transcr.obj = quad.phenotype,
                    connect.obj = quad.connectivity,
                    feature = feature,
                    group.by = 'class',
                    edge.fixed.size = T,min.edge.value = 0,max.edge.value = max.edge.value,cols.use = cols.use)
  p3 <- CircuitPlot(transcr.obj = tri.phenotype,
                    connect.obj = tri.connectivity,
                    feature = feature,
                    group.by = 'class',
                    edge.fixed.size = T,min.edge.value = 0,max.edge.value = max.edge.value,cols.use = cols.use)
  # p2 <- CircuitPlot(transcr.obj = tran.co,
  #             connect.obj = con.co,
  #             feature = feature,
  #             group.by = 'class',
  #             edge.fixed.size = T,min.edge.value = 0,max.edge.value = max.edge.value,cols.use = cols.use)
  # p1 <- CircuitPlot(transcr.obj = tran.mono,
  #             connect.obj = con.mono,
  #             feature = feature,
  #             group.by = 'class',
  #             edge.fixed.size = T,min.edge.value = 0,max.edge.value = max.edge.value,cols.use = cols.use)
  plot <- cowplot::plot_grid(p3,p4,nrow = 1)
  print(cowplot::plot_grid(title,plot,ncol=1,rel_heights = c(0.15, 1)))
}

# Define new metadata slot for finding maximum edge value
merge <- merge(tri.connectivity,quad.connectivity)
merge$tempSlot <- paste(merge$Dataset2.Sending,merge$class.Joint)
temp <- AverageExpression(merge,group.by = 'tempSlot',slot = 'counts')
temp <- data.frame(temp$CellToCell)
temp$max <- apply(temp, 1, max, na.rm=TRUE)

# Set colors
cols.use <- color_pals$class_colors

# height and width of circuit plots
circuit.height <- 2
circuit.width <- 3.6

# Epithelium
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Epi Circuit Plots")
for(i in 1:length(epi.moi.ordered)){
  feature = epi.moi.ordered[i]
  png(filename = paste(feature,'_CircuitPlot.png',sep=''),width = circuit.width,height = circuit.height,units = 'in',res=600)
  print(TempFunc(feature = feature,
           max.edge.value = temp[epi.moi.ordered[i],]$max))
  dev.off()
}

# Endothelium
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Endo Circuit Plots")
for(i in 1:length(end.moi.ordered)){
  feature = end.moi.ordered[i]
  png(filename = paste(feature,'_CircuitPlot.png',sep=''),width = circuit.width,height = circuit.height,units = 'in',res=600)
  print(TempFunc(feature = feature,
                 max.edge.value = temp[end.moi.ordered[i],]$max))
  dev.off()
}

# Mesenchyme
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Mes Circuit Plots")
for(i in 1:length(mes.moi.ordered)){
  feature = mes.moi.ordered[i]
  png(filename = paste(feature,'_CircuitPlot.png',sep=''),width = circuit.width,height = circuit.height,units = 'in',res=600)
  print(TempFunc(feature = feature,
                 max.edge.value = temp[mes.moi.ordered[i],]$max))
  dev.off()
}

# Immune
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Imm Circuit Plots")
for(i in 1:length(imm.moi.ordered)){
  feature = imm.moi.ordered[i]
  png(filename = paste(feature,'_CircuitPlot.png',sep=''),width = circuit.width,height = circuit.height,units = 'in',res=600)
  print(TempFunc(feature = feature,
                 max.edge.value = temp[imm.moi.ordered[i],]$max))
  dev.off()
}
# Allie Requested Additions 2024-06-12
befm_moi_annotated_formatted_AMG.immune_ligands <- read_excel("befm.moi.annotated.formatted_AMG.xlsm", sheet = "immune ligand")
befm_moi_annotated_formatted_AMG.immune_receptors <- read_excel("befm.moi.annotated.formatted_AMG.xlsm", sheet = "immune receptor")
# Requested from "immune ligand" sheet
allie.request.immune.ligand <- befm_moi_annotated_formatted_AMG.immune_ligands$MOI
allie.request.immune.ligand <- gsub(' - ','—',allie.request.immune.ligand)
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Allie_Selected_Immune_Ligand")
for(i in 1:length(allie.request.immune.ligand)){
  feature = allie.request.immune.ligand[i]
  png(filename = paste(feature,'_CircuitPlot.png',sep=''),width = circuit.width,height = circuit.height,units = 'in',res=600)
  print(TempFunc(feature = feature,
                 max.edge.value = temp[allie.request.immune.ligand[i],]$max))
  dev.off()
}
# Requested from "immune receptor" sheet
allie.request.immune.receptor <- befm_moi_annotated_formatted_AMG.immune_receptors$MOI
allie.request.immune.receptor <- gsub(' - ','—',allie.request.immune.receptor)

setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Allie_Selected_Immune_Receptor")
for(i in 1:length(allie.request.immune.receptor)){
  feature = allie.request.immune.receptor[i]
  png(filename = paste(feature,'_CircuitPlot.png',sep=''),width = circuit.width,height = circuit.height,units = 'in',res=600)
  print(TempFunc(feature = feature,
                 max.edge.value = temp[allie.request.immune.receptor[i],]$max))
  dev.off()
}

## Violin plots of phenotype RNA Expression
update_geom_defaults("violin", aes(linewidth = 0)) # Removes outlines on violin plots

TempVlnFunc <- function(phenotype.object,
                        mechanism = 'Il10—Il10ra',
                        group.by = 'class',
                        split.by='Condition',
                        cols = color_pals$Condition,
                        pt.size = 0.01,
                        adjust=0,
                        min.cut = 0){
      
      # split mechanism into ligand and receptor 
    mech.split <- strsplit(mechanism,split = '—')
    # make first plot (ligand)
    p1 <- VlnPlot(phenotype.object,
            features = mech.split[[1]][1],
            group.by = group.by,
            split.by=split.by,
            cols = cols,
            pt.size = pt.size,
            adjust=adjust)+ylim(min.cut,NA)+ylab('Ligand Expression')+NoLegend()
    # make second plot (receptor)
    p2 <- VlnPlot(phenotype.object,
                  features = mech.split[[1]][2],
                  group.by = group.by,
                  split.by=split.by,
                  cols = cols,
                  pt.size = pt.size,
                  adjust=adjust)+ylim(min.cut,NA)+ylab('Receptor Expression')+NoLegend()
    output.plot <- cowplot::plot_grid(p1,p2)
    print(output.plot)
    return(output.plot)
    }

# Epithelium Vln Plots
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Epi Vln Plots")
for(i in 1:length(epi.moi.ordered)){
  feature = epi.moi.ordered[i]
  png(filename = paste(feature,'_VlnPhenotypePlot.png',sep=''),width = 10,height = 3,units = 'in',res=600)
  print(TempVlnFunc(phenotype.object = sub.fine,
                    mechanism = feature))
  dev.off()
}

# Endothelium Vln Plots
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Endo Vln Plots")
for(i in 1:length(end.moi.ordered)){
  feature = end.moi.ordered[i]
  png(filename = paste(feature,'_VlnPhenotypePlot.png',sep=''),width = 10,height = 3,units = 'in',res=600)
  print(TempVlnFunc(phenotype.object = sub.fine,
                    mechanism = feature))
  dev.off()
}

# Mesenchyme Vln Plots
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Mes Vln Plots")
for(i in 1:length(mes.moi.ordered)){
  feature = mes.moi.ordered[i]
  png(filename = paste(feature,'_VlnPhenotypePlot.png',sep=''),width = 10,height = 3,units = 'in',res=600)
  print(TempVlnFunc(phenotype.object = sub.fine,
                    mechanism = feature))
  dev.off()
}

# Immune Vln Plots
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Imm Vln Plots")
for(i in 1:length(imm.moi.ordered)){
  feature = imm.moi.ordered[i]
  png(filename = paste(feature,'_VlnPhenotypePlot.png',sep=''),width = 10,height = 3,units = 'in',res=600)
  print(TempVlnFunc(phenotype.object = sub.fine,
                    mechanism = feature))
  dev.off()
}

## Network plots
source("~/GitHub/NICHESMethods/ToolBuilding/NetworkPlot.R")

## Epithelium

# Define features to plot
all.features <- epi.moi.ordered

# PreProcess: for each feature, find the maximum edge value and log transform (for plotting)
max.values <- qlcMatrix::rowMax(as.matrix(global@assays$CellToCell@counts[all.features,]))
max.values <- data.frame(feature = all.features,
                         max.value = as.matrix(max.values))
max.values$log.transformed <- log1p(max.values$max.value)

# Set connectivity threshold for edge pruning
min.connectivity.thresh <- 0.25

# Epithelium NetworkPlots for all features
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Epi Network Plots")
for(i in 1:length(all.features)){
  feature = all.features[i]
  png(file = paste(feature,'NetworkPlot',min.connectivity.thresh,'2024-04-26.png',sep='_'),width = 6,height = 5,units = 'in',res=600)
  print(NetworkPlot(transcriptome.object = sub.fine,
                    connectome.object = global,
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

## Endothelium

# Define features to plot
all.features <- end.moi.ordered

# PreProcess: for each feature, find the maximum edge value and log transform (for plotting)
max.values <- qlcMatrix::rowMax(as.matrix(global@assays$CellToCell@counts[all.features,]))
max.values <- data.frame(feature = all.features,
                         max.value = as.matrix(max.values))
max.values$log.transformed <- log1p(max.values$max.value)

# Set connectivity threshold for edge pruning
min.connectivity.thresh <- 0.25

# Endothelial NetworkPlots for all features
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Endo Network Plots")
for(i in 1:length(all.features)){
  feature = all.features[i]
  png(file = paste(feature,'NetworkPlot',min.connectivity.thresh,'2024-04-26.png',sep='_'),width = 6,height = 5,units = 'in',res=600)
  print(NetworkPlot(transcriptome.object = sub.fine,
                    connectome.object = global,
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


## Mesenchyme

# Define features to plot
all.features <- mes.moi.ordered

# PreProcess: for each feature, find the maximum edge value and log transform (for plotting)
max.values <- qlcMatrix::rowMax(as.matrix(global@assays$CellToCell@counts[all.features,]))
max.values <- data.frame(feature = all.features,
                         max.value = as.matrix(max.values))
max.values$log.transformed <- log1p(max.values$max.value)

# Set connectivity threshold for edge pruning
min.connectivity.thresh <- 0.25

# Mesenchymal NetworkPlots for all features
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Mes Network Plots")
for(i in 1:length(all.features)){
  feature = all.features[i]
  png(file = paste(feature,'NetworkPlot',min.connectivity.thresh,'2024-04-26.png',sep='_'),width = 6,height = 5,units = 'in',res=600)
  print(NetworkPlot(transcriptome.object = sub.fine,
                    connectome.object = global,
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

## Immune

# Define features to plot
all.features <- imm.moi.ordered

# PreProcess: for each feature, find the maximum edge value and log transform (for plotting)
max.values <- qlcMatrix::rowMax(as.matrix(global@assays$CellToCell@counts[all.features,]))
max.values <- data.frame(feature = all.features,
                         max.value = as.matrix(max.values))
max.values$log.transformed <- log1p(max.values$max.value)

# Set connectivity threshold for edge pruning
min.connectivity.thresh <- 0.25

# Immune NetworkPlots for all features
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Imm Network Plots")
for(i in 1:length(all.features)){
  feature = all.features[i]
  png(file = paste(feature,'NetworkPlot',min.connectivity.thresh,'2024-04-26.png',sep='_'),width = 6,height = 5,units = 'in',res=600)
  print(NetworkPlot(transcriptome.object = sub.fine,
                    connectome.object = global,
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

## Requested by Allie 2024-06-12
source("~/GitHub/NICHESMethods/ToolBuilding/NetworkPlot.R")

# "Immune Ligand" request
all.features <- allie.request.immune.ligand

# PreProcess: for each feature, find the maximum edge value and log transform (for plotting)
max.values <- qlcMatrix::rowMax(as.matrix(global@assays$CellToCell@counts[all.features,]))
max.values <- data.frame(feature = all.features,
                         max.value = as.matrix(max.values))
max.values$log.transformed <- log1p(max.values$max.value)

# Set connectivity threshold for edge pruning
min.connectivity.thresh <- 0.25

# Requested Allie "Immune Ligand" NetworkPlots for all features
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Allie_ImmuneLigand_NetworkPlots")
for(i in 1:length(all.features)){
  feature = all.features[i]
  png(file = paste(feature,'NetworkPlot',min.connectivity.thresh,'2024-06-12.png',sep='_'),width = 6,height = 5,units = 'in',res=600)
  print(NetworkPlot(transcriptome.object = sub.fine,
                    connectome.object = global,
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

# "Immune Receptor" request
all.features <- allie.request.immune.receptor

# PreProcess: for each feature, find the maximum edge value and log transform (for plotting)
max.values <- qlcMatrix::rowMax(as.matrix(global@assays$CellToCell@counts[all.features,]))
max.values <- data.frame(feature = all.features,
                         max.value = as.matrix(max.values))
max.values$log.transformed <- log1p(max.values$max.value)

# Set connectivity threshold for edge pruning
min.connectivity.thresh <- 0.25

# Requested Allie "Immune Receptor" NetworkPlots for all features
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Allie_ImmuneReceptor_NetworkPlots")
for(i in 1:length(all.features)){
  feature = all.features[i]
  png(file = paste(feature,'NetworkPlot',min.connectivity.thresh,'2024-06-12.png',sep='_'),width = 6,height = 5,units = 'in',res=600)
  print(NetworkPlot(transcriptome.object = sub.fine,
                    connectome.object = global,
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






