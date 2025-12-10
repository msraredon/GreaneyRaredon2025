# Require packages
require(Seurat)
require(NICHES)
require(dplyr)
require(cowplot)
require(ggplot2)
require(ggrepel)

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


sc.ent <- function(X){
  dist.x <- dist(X,method = 'euclidean')
  signet.x <- data.frame(cmdscale(dist.x)) # signaling network (can be visualized)
  dist.x.mat <- as.matrix(dist.x) # distance
  ss.x <- sum(dist.x.mat^2) # sum of squares
  nfeat.x <- sum(X>0)
  ssff.x <- ss.x/nfeat.x^2
  
  output <- list(dist.x,signet.x,dist.x.mat,
                 ss.x,nfeat.x,ssff.x)
  
  names(output) <- c('dist.x','signet.x','dist.x.mat',
                     'ss.x','nfeat.x','ssff.x')
  return(output)
}

library(future.apply)
plan("multisession", workers = 8)  #adjust based on how many cores you have. I have 10 total so i am reserving 2 for other tasks
plan()
options(future.globals.maxSize= 10000*1024^2) # 10 Gb futures export limit, bug fix for below

# Apply entropic single cell function
Idents(epi) <- epi$Tissue
table(Idents(epi))
little.demo <- subset(epi,cells=WhichCells(epi,downsample = 100))
test <- future_apply(little.demo@assays$SystemToCell@data,
                     MARGIN = 2,
                     FUN = sc.ent) # holy shit it worked

save(test,file = 'test.little.demo.epi.Robj')

# ok so now we want to...
#test$`28x14`$ssff.x
ssff <- c()
for(i in 1:length(test)){
  ssff <- c(ssff,test[[i]][['ssff.x']])
}
little.demo$ssff <- log(ssff)

VlnPlot(little.demo[,!is.na(little.demo$ssff)],features = 'ssff')

# information?
little.demo$info <- -log(ssff)
VlnPlot(little.demo[,!is.na(little.demo$ssff)],'info')

# does it correlate with feature depth? should it?
FeatureScatter(little.demo[,!is.na(little.demo$ssff)],'nFeature_SystemToCell','ssff')

# BCL5 Signet Viz
names(test)
bcl5 <- test[1:100]
bcl5.signets <- list()
for(i in 1:length(bcl5)){
  if(ncol(bcl5[[i]]$signet.x)>0){
  bcl5.signets[[i]] <- bcl5[[i]]$signet.x
  
  }else{
    bcl5.signets[[i]] <- NA
  }
}
# get a mean matrix for visualization
bcl5.signets <- bcl5.signets[!is.na(bcl5.signets)]
bcl5.signet.mean <- Reduce("+", bcl5.signets) / length(bcl5.signets)
ggplot(bcl5.signet.mean,
       aes(x=X1,y=X2,label=rownames(bcl5.signet.mean)))+
  geom_point()+
  geom_text_repel()+theme_classic()

 # BEFM1 Signet Viz
names(test)
befm1 <- test[801:900]
befm1.signets <- list()
for(i in 1:length(befm1)){
  if(ncol(befm1[[i]]$signet.x)>0){
    befm1.signets[[i]] <- befm1[[i]]$signet.x
  }else{
    befm1.signets[[i]] <- NA
  }
}
# get a mean matrix for visualization
befm1.signets <- befm1.signets[!is.na(befm1.signets)]
befm1.signet.mean <- Reduce("+", befm1.signets) / length(befm1.signets)
ggplot(befm1.signet.mean,
       aes(x=X1,y=X2,label=rownames(befm1.signet.mean)))+
  geom_point()+
  geom_text_repel()+theme_classic()


