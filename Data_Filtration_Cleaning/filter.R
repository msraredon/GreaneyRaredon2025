#### FILTER
#### Code for publication unifying all 'filtration' scripts
#### Replicable method for taking STARSolo DGEs to filtered objects
####    including data filtration supplemental figures




# make_initial_objects.R
#### Loading and assembling list of Seurat objects and metadata object
## Total objects to make (n=19):
## BC1P3, BC1P6, RLMVEC, FB13, FB14, BAL, BCL5, BCEC2, BEF1, BEF2, BEF3,
## BEF12, BEF14, BEF15, BEFM1, BEFM2, BEFM4, BEFM5, BEFM6

# Get File names
file.names <- list.files(path = ".", pattern = "min.nUMI")

# Get Sample names
sample.names <- strsplit(file.names, split = "[.]")
for (i in 1:length(sample.names)){
  sample.names[[i]] <- sample.names[[i]][1]
}
sample.names <- unlist(sample.names)
sample.names <- gsub("-", "_", sample.names)

# Load Seurat Objects Function - Creates Seurat objects w all assays/metadata
LoadSeuratData <- function(file.names,sample.names){
    data <- list()
    for (i in 1:length(file.names)){
      message(sample.names[i])
      load(file = file.names[i])
      data[[i]] <- CreateSeuratObject(counts = output[['Gene']])
      data[[i]][['GeneFull']] <- CreateAssayObject(output[['GeneFull']])
      data[[i]][['spliced']] <- CreateAssayObject(output[['spliced']])
      data[[i]][['unspliced']] <- CreateAssayObject(output[['unspliced']])
      data[[i]][['percent.mt']] <- PercentageFeatureSet(data[[i]], pattern = '^Mt-')
      data[[i]]$perc.spliced <- 100*(colSums(data[[i]]@assays$spliced)/(colSums(data[[i]]@assays$spliced)+colSums(data[[i]]@assays$unspliced)))
      data[[i]]$Sample <- sample.names[i]
      gc()
    }
    return(data)
}

seurat.data <- LoadSeuratData(file.names = file.names,sample.names = sample.names)
names(seurat.data) <- sample.names

# Load MetaData Function - Creates dataframe for some plots (less data)
LoadMetaData <- function(seurat.data){
  data <- data.frame()
  for (i in 1:length(seurat.data)){
    data <- rbind(data,seurat.data[[i]]@meta.data)
    gc()
  }
  return(data)
}

bef.metadata <- LoadMetaData(seurat.data = seurat.data)
bef.metadata$Sample <- factor(bef.metadata$Sample,levels=c('BC1P3','BC1P6','RLMVEC','FB13','FB14','BAL',
    'BCL5','BCEC2','BEF1','BEF2','BEF3','BEF12','BEF14','BEF15','BEFM1','BEFM2','BEFM4','BEFM5','BEFM6'))




# prefiltration_plots.R
#### Initial QC plots post-30k rank cutoff, pre-filtration and pre-clean
## Plot full dataset in single QC violin plot
png(filename = 'BEF_nCount_RNA_Unfiltered.png',width = 20,height = 10,res = 200,units = 'in')
ggplot(bef.metadata, aes(x=Sample, y=nCount_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nCount_RNA')
dev.off()

png(filename = 'BEF_nFeature_RNA_Unfiltered.png',width = 20,height = 10,res = 200,units = 'in')
ggplot(bef.metadata, aes(x=Sample, y=nFeature_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nFeature_RNA')
dev.off()

png(filename = 'BEF_percent.mt_Unfiltered.png',width = 20,height = 10,res = 200,units = 'in')
ggplot(bef.metadata, aes(x=Sample, y=percent.mt, fill=Sample)) +
  geom_violin(trim=FALSE)+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('percent.mt')
dev.off()


## Plot individual object QC plots
# Standard violin plots of nCount_RNA, nFeature_RNA, percent.mt
plots <- vector('list',length(seurat.data))
for (i in 1:length(seurat.data)) {
    message(names(seurat.data[i]))
    object <- seurat.data[[i]]
    p <- VlnPlot(object, features = c('nCount_RNA'),pt.size=0.1,log=T) +
        theme(legend.position = 'none',axis.title.x = element_blank())
    q <- VlnPlot(object, features = c('nFeature_RNA'),pt.size=0.1,log=T) +
        theme(legend.position = 'none',axis.title.x = element_blank())
    r <- VlnPlot(object, features = c('percent.mt'),pt.size=0.1,y.max=50) +
        theme(legend.position = 'none',axis.title.x = element_blank())
    plots[[i]] <- plot_grid(p,q,r,ncol=3)
    png(paste0(names(seurat.data[i]),'_qc_1.png'),width = 1500,height=500)
    print(plots[[i]])
    dev.off()
}
pdf(file = paste('BEF_QC_First_Extraction.pdf',sep='_'),width=24,height=12)
print(plots)
dev.off()

# QC plots based on Luecken 2019, Figure 2D - nFeature_RNA vs nCount_RNA, colored by percent.mt
samples <- levels(bef.metadata$Sample)
plots <- vector('list',length(samples))
for (i in 1:length(samples)) {
    message(samples[i])
    object <- subset(bef.metadata,Sample %in% samples[i])
    p <- ggplot(object, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point(aes(color=percent.mt)) +
        scale_color_viridis(limits = c(0, 50), oob = scales::squish) +
        labs(title = samples[i]) +
        theme(plot.title = element_text(size = 24, hjust = 0.5),
            axis.title.x = element_text(size = 20,color='black'),
            axis.title.y = element_text(size = 20,color='black'),
            axis.text.x = element_text(size=12,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.line = element_line(size = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())
    plots[[i]] <- p
    png(paste0(samples[i],'_qc_2.png'),width = 600,height=500)
    print(p)
    dev.off()
}
pdf(file = paste('BEF_QC_First_Extraction_2.pdf',sep='_'),width=12,height=10)
print(plots)
dev.off()

# Same plot with all cells together, colored by percent.mt
png('BEF_qc_2.png',width = 700,height=500)
ggplot(bef.metadata, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point(aes(color=percent.mt)) +
    scale_color_viridis(limits = c(0, 50), oob = scales::squish) +
    labs(title = 'All BEF Samples') +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size = 12,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

# Luecken 2019, Figure 2D - nFeature_RNA vs nCount_RNA, colored by perc.spliced
samples <- levels(bef.metadata$Sample)
plots <- vector('list',length(samples))
for (i in 1:length(samples)) {
    message(samples[i])
    object <- subset(bef.metadata,Sample %in% samples[i])
    p <- ggplot(object, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point(aes(color=perc.spliced)) +
        scale_color_viridis(option="magma",direction=-1,limits = c(50, 100), oob = scales::squish) +
        labs(title = samples[i]) +
        theme(plot.title = element_text(size = 24, hjust = 0.5),
            axis.title.x = element_text(size = 20,color='black'),
            axis.title.y = element_text(size = 20,color='black'),
            axis.text.x = element_text(size=12,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.line = element_line(size = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())
    plots[[i]] <- p
    png(paste0(samples[i],'_qc_3.png'),width = 600,height=500)
    print(p)
    dev.off()
}
pdf(file = paste('BEF_QC_First_Extraction_3.pdf',sep='_'),width=12,height=10)
print(plots)
dev.off()

# Same plot with all cells together, colored by perc.spliced
png('BEF_qc_3.png',width = 700,height=500)
ggplot(bef.metadata, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point(aes(color=perc.spliced)) +
    scale_color_viridis(option="magma",direction=-1,limits = c(50, 100), oob = scales::squish) +
    labs(title = 'All BEF Samples') +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size = 12,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()




# filtration_1.R
#### Choosing filtration cutoffs for each object 
## Functions for iterating on nCount_RNA, nFeature_RNA & percent.mt filtration parameters

# Three violin plots of QC features w cutoff line
QCViolinPlotter <- function(sample.name,min.counts,min.features,max.mt){
    x <- which(sample.name == names(seurat.data))
    object <- seurat.data[[x]]
    filtered <- subset(object,nFeature_RNA > min.features & nCount_RNA > min.counts & percent.mt < max.mt)
    tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
    p <- VlnPlot(object, features = c("nCount_RNA"),pt.size=0.1,log=T) +
        geom_hline(yintercept=min.counts, linetype='dashed', color='red', size=1) +
        labs(subtitle = paste("nUMI >",min.counts)) + labs(caption = paste(" ")) +
        theme(legend.position = 'none',axis.title.x = element_blank(),axis.text.x = element_blank(),plot.subtitle = element_text(color='red'))
    q <- VlnPlot(object, features = c("nFeature_RNA"),pt.size=0.1,log=T) +
        geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
        labs(subtitle = paste("nGene >",min.features)) + labs(caption = paste(" ")) +
        theme(legend.position = 'none',axis.title.x = element_blank(),axis.text.x = element_blank(),plot.subtitle = element_text(color='red'))
    r <- VlnPlot(object, features = c("percent.mt"),pt.size=0.1,y.max=50) +
        geom_hline(yintercept=max.mt, linetype='dashed', color='red', size=1) +
        labs(subtitle = paste("pt.mt <",max.mt)) +
        labs(caption = paste("Barcodes retained with filtration / starting total:",ncol(filtered),"/",ncol(object))) +
        theme(legend.position = 'none',axis.title.x = element_blank(),axis.text.x = element_blank(),plot.subtitle = element_text(color='red'),plot.caption = element_text(color='red'))
    png(paste0(sample.name,tag,'_filterselect1.png'),width = 1500,height=500)
    print(plot_grid(p,q,r,ncol=3))
    dev.off()
    p <<- plot_grid(p,q,r,ncol=3)
}

# Scatter plot of nFeature_RNA vs nCount_RNA, colored by percent.mt
TriplePercentMtPlotter <- function(sample.name,min.counts,min.features,max.mt){
    object <- subset(bef.metadata,Sample %in% sample.name)
    lowpass <- subset(object,percent.mt < max.mt)
    highmt <- subset(object,percent.mt > max.mt)
    filtered <- subset(object,nFeature_RNA > min.features & nCount_RNA > min.counts & percent.mt < max.mt)
    tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
    p <- ggplot() +
        geom_point(data=highmt, aes(x=nCount_RNA, y=nFeature_RNA), color = 'gray') +
        geom_point(data=lowpass, aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
        scale_color_viridis(limits = c(0, 50), oob = scales::squish) +
        geom_vline(xintercept=min.counts, linetype='dashed', color='red', size=1) +
        geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
        labs(title = sample.name) +
        labs(subtitle = paste("nCount_RNA >",min.counts,"& nFeature_RNA >",min.features,"& percent.mt <",max.mt)) +
        labs(caption = paste("Barcodes retained with filtration / starting total:",nrow(filtered),"/",nrow(object))) +
        theme(plot.title = element_text(size = 24, hjust = 0.5),
            plot.subtitle = element_text(color='red'),
            plot.caption = element_text(color='red'),
            axis.title.x = element_text(size = 20,color='black'),
            axis.title.y = element_text(size = 20,color='black'),
            axis.text.x = element_text(size=12,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.line = element_line(size = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())
    png(paste0(sample.name,tag,'_filterselect2.png'),width = 600,height=500)
    print(p)
    dev.off()
    p <<- p
}

# Scatter plot of nFeature_RNA vs nCount_RNA, colored by perc.spliced
TriplePercentSplPlotter <- function(sample.name,min.counts,min.features,max.mt){
    object <- subset(bef.metadata,Sample %in% sample.name)
    lowpass <- subset(object,percent.mt < max.mt)
    highmt <- subset(object,percent.mt > max.mt)
    filtered <- subset(object,nFeature_RNA > min.features & nCount_RNA > min.counts & percent.mt < max.mt)
    tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
    p <- ggplot() +
        geom_point(data=highmt, aes(x=nCount_RNA, y=nFeature_RNA), color = 'gray') +
        geom_point(data=lowpass, aes(x=nCount_RNA, y=nFeature_RNA, color=perc.spliced)) +
        scale_color_viridis(option="magma",direction=-1,limits = c(50, 100), oob = scales::squish) +
        geom_vline(xintercept=min.counts, linetype='dashed', color='red', size=1) +
        geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
        labs(title = sample.name) +
        labs(subtitle = paste("nCount_RNA >",min.counts,"& nFeature_RNA >",min.features,"& percent.mt <",max.mt)) +
        labs(caption = paste("Barcodes retained with filtration / starting total:",nrow(filtered),"/",nrow(object))) +
        theme(plot.title = element_text(size = 24, hjust = 0.5),
            plot.subtitle = element_text(color='red'),
            plot.caption = element_text(color='red'),
            axis.title.x = element_text(size = 20,color='black'),
            axis.title.y = element_text(size = 20,color='black'),
            axis.text.x = element_text(size=12,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.line = element_line(size = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())
    png(paste0(sample.name,tag,'_filterselect3.png'),width = 600,height=500)
    print(p)
    dev.off()
    p <<- p
}

# Iterating filtration on individual samples
sample.name <- 'BEFM1'
min.counts <- 200
min.features <- 200
max.mt <- 10

# Plot QC metrics
QCViolinPlotter(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt)
TriplePercentMtPlotter(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt)
TriplePercentSplPlotter(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt)




# filtration_2.R
#### Choosing filtration cutoffs for each object
## Iterate filtration parameters and visualize effects on data quality and clustering

# Titrating QC parameters on existing UMAP
titrationUMAP <- function(object,min.counts,min.features,max.mt){
  object@meta.data$Passing <- 'NoPass'
  cells.Pass <- WhichCells(object,expression = nCount_RNA > min.counts & nFeature_RNA > min.features & percent.mt < max.mt)
  meta.data.pass <- data.frame(Barcode = cells.Pass,Passing = 'Pass')
  rownames(meta.data.pass) <- cells.Pass
  if(length(cells.Pass) < ncol(object)){
    cells.NoPass <- WhichCells(object,cells = cells.Pass,invert = T)
    meta.data.nopass <- data.frame(Barcode = cells.NoPass,Passing = 'NoPass')
    rownames(meta.data.nopass) <- cells.NoPass
  }else{meta.data.nopass <- data.frame()}
  meta.data <- rbind(meta.data.pass,meta.data.nopass)
  object <- AddMetaData(object,metadata = meta.data)
  tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
  p <- DimPlot(object,group.by = 'Passing',cols = c('Pass'='Grey','NoPass'='Red'),shuffle=T) +
        labs(title = sample.name,subtitle = paste('nCount_RNA >',min.counts,'& nFeature_RNA >',min.features,'& percent.mt <',max.mt)) +
        theme(plot.subtitle = element_text(color='red'))
  png(paste(sample.name,tag,'_titration.png',sep=""), width = 500, height = 500)
  print(p)
  dev.off()
  p <<- p
}


## Iterate filtration parameters, plot and evaluate
# Select filters for individual samples
sample.name <- 'BEFM1'
setwd(paste("/nobackup1/greaneya/DGEs_30K_Rank_Cutoff/",sample.name,sep=''))
min.counts <- 1000
min.features <- 400
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')

# Plot QC metrics
QCViolinPlotter(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt)
p1 <- p
TriplePercentMtPlotter(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt)
p2 <- p
TriplePercentSplPlotter(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt)
p3 <- p

# Pull object, apply cutoffs & dimensional reduction
x <- which(sample.name == names(seurat.data))
object <- seurat.data[[x]]
object <- subset(object, nCount_RNA > min.counts & nFeature_RNA > min.features & percent.mt < max.mt)
object <- NormalizeData(object)
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, npcs = 50, verbose = F)
png(paste0(sample.name,tag,'_elbow.png'), width = 800, height = 500)
ElbowPlot(object,ndims = 50)+labs(title = sample.name)
dev.off()
pdf(paste0(sample.name,tag,'_pcheatmap1.pdf'), width = 12, height = 16)
DimHeatmap(object, dims = 1:15, cells = 200, balanced = T)
dev.off()
pdf(paste0(sample.name,tag,'_pcheatmap2.pdf'), width = 12, height = 16)
DimHeatmap(object, dims = 16:30, cells = 200, balanced = T)
dev.off()
pdf(paste0(sample.name,tag,'_pcheatmap3.pdf'), width = 12, height = 16)
DimHeatmap(object, dims = 31:45, cells = 200, balanced = T)
dev.off()

# Cluster, plot & check QC
pcs <- 40
res <- 0.5
object <- FindNeighbors(object, dims = 1:pcs)
object <- FindClusters(object, resolution = res)
object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
tag2 <- paste(tag,'pcs',pcs,'res',res,sep='.')
p4 <- UMAPPlot(object,label = T) + labs(title = sample.name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoLegend() + NoAxes()
p5 <- FeaturePlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'), label = T,repel=T)
p6 <- VlnPlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'),ncol=2,pt.size=0.1)
p7 <- FeaturePlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','perc.spliced','Mki67','Top2a'),ncol=2,label = T,repel=T)
p8 <- VlnPlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','perc.spliced','Mki67','Top2a'),ncol=2,pt.size=0.1)

# Backtrack to original filtration plots
nclusters <- length(levels(object$seurat_clusters))

# ViolinPlots by cluster
p <- VlnPlot(object, features = c("nCount_RNA"),pt.size=0.1,log=T) +
    geom_hline(yintercept=min.counts, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("nCount_RNA >",min.counts)) + theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'))
q <- VlnPlot(object, features = c("nFeature_RNA"),pt.size=0.1,log=T) +
    geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("nFeature_RNA >",min.features)) + theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'))
r <- VlnPlot(object, features = c("percent.mt"),pt.size=0.1,y.max=(max.mt+0.5)) +
    geom_hline(yintercept=max.mt, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("percent.mt <",max.mt)) +
    theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'),plot.caption = element_text(color='red'))
p9 <- plot_grid(p,q,r,ncol=3)

# ViolinPlots whole sample, jitter colored by cluster
p <- VlnPlot(object, features = c("nCount_RNA"),pt.size=0,log=T,group.by='orig.ident',cols=1) +
    geom_jitter(aes(color = factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    geom_hline(yintercept=min.counts, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("nCount_RNA >",min.counts)) + theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'))
q <- VlnPlot(object, features = c("nFeature_RNA"),pt.size=0,log=T,group.by='orig.ident',cols=1) +
    geom_jitter(aes(color = factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("nFeature_RNA >",min.features)) + theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'))
r <- VlnPlot(object, features = c("percent.mt"),pt.size=0,group.by='orig.ident',cols=1,y.max=(max.mt+0.5)) +
    geom_jitter(aes(color = factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    geom_hline(yintercept=max.mt, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("percent.mt <",max.mt)) +
    theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'))
p10 <- plot_grid(p,q,r,ncol=3)

# Scatter plot colored by cluster
list <- subset(bef.metadata,Sample %in% sample.name)
filtered <- subset(list,nCount_RNA > min.counts & nFeature_RNA > min.features & percent.mt < max.mt)
p11 <- ggplot() +
    geom_point(data=filtered, aes(x=nCount_RNA, y=nFeature_RNA, color=factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    geom_vline(xintercept=min.counts, linetype='dashed', color='red', size=1) +
    geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
    labs(title = sample.name) +
    labs(subtitle = paste("nCount_RNA >",min.counts,"& nFeature_RNA >",min.features,"& percent.mt <",max.mt)) +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
        plot.subtitle = element_text(color='red'),
        plot.caption = element_text(color='red'),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size = 12,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Class Characterization
# Epi: 'Krt5','Sftpc','Hopx','Aqp5','Sftpb','Defb4','Lyz2','Scgb1a1'
# Endo: 'Plvap','Lyve1','Prx','Car4','Apln','Aplnr','Prox1'
# Mes: 'Col13a1','Itga8','Col14a1','Lgr5','Dcn','Myh11','Acta2','Wnt11','Gucy1a1','Gucy2a1','Lgr6'
# Imm: 'Mrc1','Cd68','Naaa','C1qc'
p12 <- FeaturePlot(object, c('Krt5','Sftpc','Hopx','Aqp5'), label = T, repel = T)
p13 <- FeaturePlot(object, c('Sftpc','Sftpb','Defb4','Lyz2'), label = T, repel = T)

# Assemble plots into a cowplot
plots <- list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13)
ht <- ceiling(length(plots)/3)*600
png(paste0(sample.name,tag2,'_allplots.png'), width = 2000,height=ht)
plot_grid(plotlist=plots,ncol=3)
dev.off()


## Titration of QC parameters from existing UMAP
# Set equal number of tested parameters for each for most useful plotting
countplots <- list()
count.test <- c(150,200,250,300)  # SET VALUES
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
  countplots[[i]] <- p
}
featureplots <- list()
feature.test <- c(150,200,250,300)  # SET VALUES
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
  featureplots[[i]] <- p
}
mtplots <- list()
mt.test <- c(20,18,15,10)  # SET VALUES
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
  mtplots[[i]] <- p
}
# Assemble plots into a cowplot
wd <- max(length(countplots),length(featureplots),length(mtplots))*500
png(paste0(sample.name,tag2,'_alltitrations.png'),width = wd,height = 1500)
plot_grid(plotlist=c(countplots,featureplots,mtplots),nrow=3)
dev.off()


## Additional Optional Characterization:
# Ribosomal Characterization
object[['percent.ribo']] <- PercentageFeatureSet(object, pattern = "^Rp")
png(paste0(sample.name,tag2,'_fp_ribo.png'), width = 500, height = 500)
FeaturePlot(object, c('percent.ribo'), label = T)
dev.off()
png(paste0(sample.name,tag2,'_vln_ribo.png'), width = 500, height = 500)
VlnPlot(object, c('percent.ribo'),pt.size=0.1)
dev.off()





# filtration_3.R
#### Choosing filtration cutoffs for each object
## Generate supplemental figure to visualize initial filtration

# Kneeplot function
Barcode_rank_plot <- function(nBarcodes, xmin, xmax){
  # Load UMIperCellSorted
  UMIperCellSorted <- read.table(list.files(path = ".", pattern = "UMIperCellSorted"), quote="\"", comment.char="")
  # Make a knee plot
  counts <- UMIperCellSorted$V1
  counts <- sort(counts, decreasing = TRUE)
  ranks <- seq(length(counts))
  counts[duplicated(counts)] <- NA
  ranks[duplicated(counts)] <- NA
  datf <- data.frame(Rank = ranks, Count = counts)
  datf <- datf[!is.na(datf[, 2]), , drop = FALSE]
  # Define cutoff point in terms of nUMI
  cutoff <- UMIperCellSorted[nBarcodes,]
  # Plot
  p <- ggplot(datf, aes_string(x = "Rank", y = "Count")) +
    geom_point() +
    scale_x_continuous(trans = "log10",
                       breaks=c(10, 100, 500, 1000, 5000, 10000, 30000, 50000, 100000, 500000),
                       limits=c(xmin, xmax)) +
    scale_y_continuous(name = "nCount_RNA", trans = "log10",
                       breaks=c(1, 10, 100, 150, 200, 500, 1000, 5000, 10000, 25000, 50000),
                       limits=c(1, 500000)) +
    geom_vline(xintercept=nBarcodes, linetype="dashed", color="red")  +
    geom_hline(yintercept=cutoff, linetype="dashed", color="blue")  +
    theme_minimal(base_size = 18) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(sample.name)+
    annotate("text", label = paste("Rank Threshold =",nBarcodes), x = xmin, y = 4,hjust=0)+
    annotate("text", label = paste("nUMI Threshold =",cutoff), x = xmin, y = 2,hjust=0)
  return(p)
}

# Input selected cutoffs for each object
sample.name <- 'BEFM1'
setwd(paste("/nobackup1/greaneya/DGEs_30K_Rank_Cutoff/",sample.name,sep=''))
min.counts <- 1000
min.features <- 400
max.mt <- 25
pcs <- 40
res <- 0.5
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')

# Kneeplot
nBarcodes <- 30000 ; xmin <- 1 ; xmax <- 100000
p1 <- Barcode_rank_plot(nBarcodes, xmin, xmax)
p1 <- plot_grid(NULL,p1,NULL,nrow=3,rel_heights = c(1,3,1),labels = c('A','',''),label_size = 40)

# QC Plots
x <- which(sample.name == names(seurat.data))
object <- seurat.data[[x]]
p <- VlnPlot(object, features = c("nCount_RNA"),pt.size=0.1,log=T) +
    theme(legend.position = 'none',axis.title.x = element_blank(),axis.text.x = element_blank())
q <- VlnPlot(object, features = c("nFeature_RNA"),pt.size=0.1,log=T) +
    theme(legend.position = 'none',axis.title.x = element_blank(),axis.text.x = element_blank())
r <- VlnPlot(object, features = c("percent.mt"),pt.size=0.1,y.max=50) +
    theme(legend.position = 'none',axis.title.x = element_blank(),axis.text.x = element_blank())
p2 <- plot_grid(p,q,r,ncol=3)

object <- subset(bef.metadata,Sample %in% sample.name)
p3 <- ggplot(object, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point(aes(color=percent.mt)) +
    scale_color_viridis(limits = c(0, 50), oob = scales::squish) +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size = 12,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

p4 <- ggplot(object, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point(aes(color=perc.spliced)) +
    scale_color_viridis(option="magma",direction=-1,limits = c(50, 100), oob = scales::squish) +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size = 12,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

p234 <- plot_grid(p2,p3,p4,ncol=3,labels = c('B','C','D'),label_size = 40)

# Filtration Selects
QCViolinPlotter(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt)
p5 <- p
TriplePercentMtPlotter(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt)
p6 <- p
TriplePercentSplPlotter(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt)
p7 <- p
p567 <- plot_grid(p5,p6,p7,ncol=3,labels = c('E','F','G'),label_size = 40)
p234567 <- plot_grid(p234,p567,nrow=2)
toprow <- plot_grid(p1,p234567,ncol=2,rel_widths = c(1, 2))

# UMAP
object <- seurat.data[[x]]
object <- subset(object, nCount_RNA > min.counts & nFeature_RNA > min.features & percent.mt < max.mt)
object <- NormalizeData(object)
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, npcs = 50, verbose = F)
object <- FindNeighbors(object, dims = 1:pcs)
object <- FindClusters(object, resolution = res)
object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
tag2 <- paste(tag,'pcs',pcs,'res',res,sep='.')
p8 <- UMAPPlot(object,label = T) + labs(title = sample.name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoLegend() + NoAxes() +
        theme(plot.title = element_text(size = 24),plot.subtitle = element_text(size = 20))
p8 <- plot_grid(NULL,p8,NULL,nrow=3,rel_heights = c(1,4,1),labels = c('H','',''),label_size = 40)

# QC fp & vln
p9 <- FeaturePlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','perc.spliced','Mki67','Top2a'),ncol=2,label = T,repel=T)
p10 <- VlnPlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','perc.spliced','Mki67','Top2a'),ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())
p910 <- plot_grid(p9,p10,ncol=2,labels = c('I','J'),label_size = 40)

# Lineage fp & vln
p11 <- FeaturePlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'),label=T,repel=T)
p12 <- VlnPlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'),ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())
p1112 <- plot_grid(p11,p12,ncol=2,labels = c('K','L'),label_size = 40)
p9101112 <- plot_grid(p910,p1112,nrow=2)
midrow <- plot_grid(p8,p9101112,ncol=2,rel_widths = c(1, 2))

# Filterback
nclusters <- length(levels(object$seurat_clusters))
# ViolinPlots not broken apart, jitter points colored by cluster
p <- VlnPlot(object, features = c("nCount_RNA"),pt.size=0,log=T,group.by='orig.ident',cols=1) +
    geom_jitter(aes(color = factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    geom_hline(yintercept=min.counts, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("nCount_RNA >",min.counts)) + theme(legend.position = 'none',axis.title.x = element_blank(),axis.text.x = element_blank(),plot.subtitle = element_text(color='red'))
q <- VlnPlot(object, features = c("nFeature_RNA"),pt.size=0,log=T,group.by='orig.ident',cols=1) +
    geom_jitter(aes(color = factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("nFeature_RNA >",min.features)) + theme(legend.position = 'none',axis.title.x = element_blank(),axis.text.x = element_blank(),plot.subtitle = element_text(color='red'))
r <- VlnPlot(object, features = c("percent.mt"),pt.size=0,group.by='orig.ident',cols=1,y.max=(max.mt+0.5)) +
    geom_jitter(aes(color = factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    geom_hline(yintercept=max.mt, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("percent.mt <",max.mt)) +
    theme(legend.position = 'none',axis.title.x = element_blank(),axis.text.x = element_blank(),plot.subtitle = element_text(color='red'))
p13 <- plot_grid(p,q,r,ncol=3)

# Scatter Plot colored by cluster
list <- subset(bef.metadata,Sample %in% sample.name)
filtered <- subset(list,nCount_RNA > min.counts & nFeature_RNA > min.features & percent.mt < max.mt)
p14 <- ggplot() +
    geom_point(data=filtered, aes(x=nCount_RNA, y=nFeature_RNA, color=factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    geom_vline(xintercept=min.counts, linetype='dashed', color='red', size=1) +
    geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("nCount_RNA >",min.counts,"& nFeature_RNA >",min.features,"& percent.mt <",max.mt)) +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
        plot.subtitle = element_text(color='red'),
        plot.caption = element_text(color='red'),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size = 12,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

# Class characterization
# Epi: ['Krt5','Hopx'] | 'Krt5','Sftpc','Hopx','Aqp5' ; 'Sftpc','Sftpb','Defb4','Lyz2'
# Endo: ['Lyve1','Prx'] | 'Lyve1','Prx','Plvap',('Ca4') ; 'Apln',('Aplnr')
# Mes: ['Dcn','Acta2'] | 'Dcn',('Col13a1','Itga8','Col14a1') ; 'Acta2',('Wnt11','Myh11') ; ('Gucy1a1','Gucy2a1' ; 'Lgr6')
# Imm: ['Mrc1','Naaa'] | 'Cd3d','S100a8'
# BEF: 'Krt5','Hopx','Lyve1','Prx','Dcn','Acta2' | BEFM: 'Krt5','Hopx','Lyve1','Dcn','Mrc1','Naaa'
p15 <- FeaturePlot(object, c('Lyve1','Prx','Plvap','Apln'), label = T, repel = T)
lastrow <- plot_grid(p13,p14,p15,ncol=3,labels = c('M','N','O'),label_size = 40)

# Assemble into a cowplot figure
plots <- list(toprow,midrow,lastrow)
png(paste0(sample.name,tag2,'_suppfig.png'), width = 2000,height=2750)
plot_grid(plotlist=plots,ncol=1,rel_heights = c(2,2,1))
dev.off()




# filtration_4.R
#### Choosing filtration cutoffs for each object
## Save new object with selected filtration parameters

## Selected cutoffs:
#   BC1P3: 1000/600/5  20  DONE
#   BC1P6: 1000/1000/10  15  DONE
#   RLMVEC: 2000/1500/10  15  DONE
#   FB13: 1000/800/7.5  20  DONE
#   FB14: 1000/1000/5  20  DONE
#   BAL: 500/400/15  20  DONE
#   BCL5: 140/90/7.5  20  DONE
#   BCEC2: 1600/500/25  30  DONE
#   BEF1: 300/200/10  30  DONE
#   BEF2: 400/250/10  30  DONE
#   BEF3: 300/250/10  30  DONE
#   BEF12: 1000/500/20  30  DONE
#   BEF14: 1500/700/20  30  DONE
#   BEF15: 1400/600/20  30  DONE
#   BEFM1: 1000/400/25  40  DONE
#   BEFM2: 1100/600/20  40  DONE
#   BEFM4: 600/400/15  40  DONE
#   BEFM5: 400/250/15  40  DONE
#   BEFM6: 800/500/18  40  DONE

# Data frame of selected filtration parameters
objects <- c('BC1P3','BC1P6','RLMVEC','FB13','FB14','BAL','BCL5','BCEC2','BEF1','BEF2','BEF3','BEF12','BEF14','BEF15',
        'BEFM1','BEFM2','BEFM4','BEFM5','BEFM6')
min.counts <- c(1000,1000,2000,1000,1000,500,140,1600,300,400,300,1000,1500,1400,1000,1100,600,400,800)
min.features <- c(600,1000,1500,800,1000,400,90,500,200,250,250,500,700,600,400,600,400,250,500)
max.mt <- c(5,10,10,7.5,5,15,7.5,25,10,10,10,20,20,20,25,20,15,15,18)
pcs <- c(20,15,15,20,20,20,20,30,30,30,30,30,30,30,40,40,40,40,40)

filter1 <- data.frame(min.counts,min.features,max.mt,pcs)
rownames(filter1) <- objects

# Function to apply filtration to seurat.data
UpdateSeuratData <- function(data,filter) {
    new.data <- list() ; final.names <- c()
    myplots <- vector('list',length(data))
    for (i in 1:length(data)) {
        object <- data[[i]]
        sample.name <- names(data[i])
        message(sample.name)
        params <- filter[sample.name,]
        min.counts <- params$min.counts
        min.features <- params$min.features
        max.mt <- params$max.mt
        pcs <- params$pcs
        # make new object
        object <- subset(object,nFeature_RNA > min.features & nCount_RNA > min.counts & percent.mt < max.mt)
        object <- NormalizeData(object)
        object <- FindVariableFeatures(object)
        object <- ScaleData(object)
        object <- RunPCA(object, npcs = 50, verbose = F)
        object <- FindNeighbors(object, dims = 1:pcs)
        object <- FindClusters(object, resolution = 0.5)
        object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
        p <- UMAPPlot(object,label = T) + labs(title = sample.name,subtitle = paste("PC's =",pcs,"\nRes = 0.5",'\nCells =',ncol(object))) + NoLegend() + NoAxes()
        myplots[[i]] <- p
        new.data[[i]] <- object
        final.names[[i]] <- sample.name
    }
    names(new.data) <- final.names
    pdf(file = paste('UpdateSeuratData_umapchecks_1.pdf',sep=''),width=5,height=5)
    print(myplots)
    dev.off()
    return(new.data)
}

new.seurat.data <- UpdateSeuratData(data=seurat.data,filter=filter1)

save(new.seurat.data,file = paste("bef.seurat.objs.filter1.",Sys.Date(),".Robj",sep=""))

# Create new bef.metadata from new seurat.data (same as above)
bef.metadata <- LoadMetaData(seurat.data = new.seurat.data)
bef.metadata$Sample <- factor(bef.metadata$Sample,levels=c('BC1P3','BC1P6','RLMVEC','FB13','FB14','BAL','BCL5','BCEC2','BEF1','BEF2','BEF3',
    'BEF12','BEF14','BEF15','BEFM1','BEFM2','BEFM4','BEFM5','BEFM6'))

save(bef.metadata,file = paste("bef.metadata.filter1.",Sys.Date(),".Robj",sep=""))







