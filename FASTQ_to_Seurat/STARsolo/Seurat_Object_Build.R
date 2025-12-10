# Set WD
setwd('/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/DGEs_30000_Rank_Cutoff')

# Packages
require(Seurat)
require(ggplot2)

# Get Filenames
file.names <- list.files(file.path("BEFM_Project"))
file.names

# Get Sample Names
sample.names <- strsplit(file.names, split = "[.]")
for (i in 1:length(sample.names)){
  sample.names[[i]] <- sample.names[[i]][1]
}

# Define File Paths
file.paths <- paste("BEFM_Project",file.names,sep='/')
file.paths

# Load Data
data <- list()
for (i in 1:length(file.paths)){
  rm(output)
  load(file = file.paths[i])
  data[[i]] <- output
  gc()
}

# Convert to Seurat object
seurat.data <- list()
for (i in 1:length(data)){
  seurat.data[[i]] <- CreateSeuratObject(counts = data[[i]][['Gene']])
  seurat.data[[i]]$Sample <- sample.names[[i]][1] # Tag with sample metadata
  seurat.data[[i]][['GeneFull']] <- CreateAssayObject(data[[i]][['GeneFull']])
  seurat.data[[i]][['spliced']] <- CreateAssayObject(data[[i]][['spliced']])
  seurat.data[[i]][['unspliced']] <- CreateAssayObject(data[[i]][['unspliced']])
  gc()
}
names(seurat.data) <- sample.names

# Clean
rm(data)
gc()

# Merge
merge <- merge(seurat.data[[1]],seurat.data[-1])
gc()

# Define key metrics
merge[["percent.mt"]] <- PercentageFeatureSet(merge, pattern = "^Mt-")
perc.spliced <- 100*(colSums(merge@assays$spliced)/(colSums(merge@assays$spliced)+colSums(merge@assays$unspliced)))
merge$perc.spliced <- perc.spliced


# Plot QC metrics
pdf(file = 'QC_First_Extraction.pdf',width=24,height=12)
VlnPlot(merge, features = c("nCount_RNA"), ncol = 1,pt.size=0.1)+scale_y_log10()+geom_hline(yintercept=800, linetype='dashed', color='red', size=1)
VlnPlot(merge, features = c("nCount_RNA"), ncol = 1,pt.size=0.1)+scale_y_log10()+geom_hline(yintercept=800, linetype='dashed', color='red', size=1)
VlnPlot(merge, features = c("nFeature_RNA"), ncol = 1,group.by = 'Sample',pt.size=0.1)+scale_y_log10()+geom_hline(yintercept=500, linetype='dashed', color='red', size=1)
VlnPlot(merge, features = c("nFeature_RNA"), ncol = 1,pt.size=0.1)+scale_y_log10()+geom_hline(yintercept=500, linetype='dashed', color='red', size=1)
VlnPlot(merge, features = c("percent.mt"), ncol = 1,group.by = 'Sample',pt.size=0.1)+geom_hline(yintercept=7.5, linetype='dashed', color='red', size=1)+geom_hline(yintercept=0.1, linetype='dashed', color='red', size=1)
VlnPlot(merge, features = c("percent.mt"), ncol = 1,pt.size=0.1)+geom_hline(yintercept=7.5, linetype='dashed', color='red', size=1)+geom_hline(yintercept=0.1, linetype='dashed', color='red', size=1)
dev.off()

# Secondary Filtering
merge2 <- subset(merge,percent.mt > 0.1 & percent.mt < 7.5 & nFeature_RNA > 500 & nCount_RNA > 800)

# Plot QC metrics after filtering
pdf(file = 'QC_After_Filtration.pdf',width=24,height=12)
VlnPlot(merge2, features = c("nCount_RNA"), ncol = 1,pt.size=0.1)+scale_y_log10()+geom_hline(yintercept=800, linetype='dashed', color='red', size=1)
VlnPlot(merge2, features = c("nCount_RNA"), ncol = 1,pt.size=0.1)+scale_y_log10()+geom_hline(yintercept=800, linetype='dashed', color='red', size=1)
VlnPlot(merge2, features = c("nFeature_RNA"), ncol = 1,group.by = 'Sample',pt.size=0.1)+scale_y_log10()+geom_hline(yintercept=500, linetype='dashed', color='red', size=1)
VlnPlot(merge2, features = c("nFeature_RNA"), ncol = 1,pt.size=0.1)+scale_y_log10()+geom_hline(yintercept=500, linetype='dashed', color='red', size=1)
VlnPlot(merge2, features = c("percent.mt"), ncol = 1,group.by = 'Sample',pt.size=0.1)+geom_hline(yintercept=7.5, linetype='dashed', color='red', size=1)+geom_hline(yintercept=0.1, linetype='dashed', color='red', size=1)
VlnPlot(merge2, features = c("percent.mt"), ncol = 1,pt.size=0.1)+geom_hline(yintercept=7.5, linetype='dashed', color='red', size=1)+geom_hline(yintercept=0.1, linetype='dashed', color='red', size=1)
dev.off()
