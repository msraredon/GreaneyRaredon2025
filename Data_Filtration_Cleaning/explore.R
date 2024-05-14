#### EXPLORE
#### Code for publication unifying all 'explore' scripts
#### Replicable method for integrating cleaned objects by class




## Note: We tried generating unified class objects by merge, mapping, or CCA integration, in that order.
##  Sufficient sample alignment was only achieved by integration. Integration was iterated across multiple
##  sample combinations, including with samples that are not included in this manuscript. All scripts have
##  been edited to run only with published samples. Only code that contributed to final processed objects is shown.




#### Assembling relevant datasets and lists of objects
## Generate list of all individual cleaned objects
allobjects <- c("BC1P3","BC1P6","RLMVEC","FB13","FB14","BAL","BCL5","BCEC2","BEF1","BEF2","BEF3",
    "BEF12","BEF14","BEF15","BEFM1","BEFM2","BEFM4","BEFM5","BEFM6")
seurat.data <- list()
for (i in 1:length(allobjects)){
    sample.name <- allobjects[i]
    message(sample.name)
    setwd(paste0("../",sample.name))
    load(list.files(path = ".", pattern = "clean.2022-07-23"))
    seurat.data[[i]] <- object
}
names(seurat.data) <- allobjects

save(seurat.data,file = paste0("bef.seurat.objs.",Sys.Date(),".Robj"))


## Object lists by cell class - Epithelium, Endothelium, Mesenchyme, Immune
epi <- list() ; epi.names <- c()
endo <- list() ; endo.names <- c()
mes <- list() ; mes.names <- c()
imm <- list() ; imm.names <- c()
for (i in 1:length(seurat.data)){
    object <- seurat.data[[i]]
    sample.name <- names(seurat.data[i])
    message(sample.name)
    Idents(object) <- object$class
    if ("Epithelium" %in% object$class){
        epi[i] <- subset(object,idents = c('Epithelium'))
        epi.names[i] <- sample.name
    }
    if ("Endothelium" %in% object$class){
        endo[i] <- subset(object,idents = c('Endothelium'))
        endo.names[i] <- sample.name
    }
    if ("Mesenchyme" %in% object$class){
        mes[i] <- subset(object,idents = c('Mesenchyme'))
        mes.names[i] <- sample.name
    }
    if ("Immune" %in% object$class){
        imm[i] <- subset(object,idents = c('Immune'))
        imm.names[i] <- sample.name
    }
}
names(epi) <- epi.names
names(endo) <- endo.names
names(mes) <- mes.names
names(imm) <- imm.names

# Remove empty slots from list
epi <- epi[lengths(epi) != 0]
endo <- endo[lengths(endo) != 0]
mes <- mes[lengths(mes) != 0]
imm <- imm[lengths(imm) != 0]

save(epi,file = paste0("epi.seurat.objs.",Sys.Date(),".Robj"))
save(endo,file = paste0("endo.seurat.objs.",Sys.Date(),".Robj"))
save(mes,file = paste0("mes.seurat.objs.",Sys.Date(),".Robj"))
save(imm,file = paste0("imm.seurat.objs.",Sys.Date(),".Robj"))


## Table of metadata info for each object
df <- data.frame()
Name <- c() ; Chemistry <- c() ; Cells <- c() ; Epithelium <- c() ; Endothelium <- c() ; Mesenchyme <- c() ; Immune <- c()
Mean_Counts <- c() ; Mean_Features <- c()
for (i in 1:length(seurat.data)){
    object <- seurat.data[[i]]
    sample.name <- names(seurat.data[i])
    message(sample.name)

    Name[i] <- sample.name
    Chemistry[i] <- unique(object$Chemistry)
    Mean_Counts[i] <- round(mean(object$nCount_RNA))
    Mean_Features[i] <- round(mean(object$nFeature_RNA))
    Cells[i] <- ncol(object)
    Epithelium[i] <- length(which(object$class == 'Epithelium'))
    Endothelium[i] <- length(which(object$class == 'Endothelium'))
    Mesenchyme[i] <- length(which(object$class == 'Mesenchyme'))
    Immune[i] <- length(which(object$class == 'Immune'))

}

df <- cbind(Name,Chemistry,Mean_Counts,Mean_Features,Cells,Epithelium,Endothelium,Mesenchyme,Immune)

p <- tableGrob(df,theme = ttheme_minimal(base_size = 20))

png(paste("bef.seurat.objs_clean_table.png"), width = 1200,height=800)
plot(p)
dev.off()




#### Generating Epithelial object

## Edit this elsewhere first...



#### EPITHELIUM ####




# explore_epi_2.R
#### Initial exploration of all Epithelium
## CCA Integration & Merge
# Keep? or obsolete integration?

## Load all epithelium
setwd("/nobackup1/greaneya/DGEs_30K_Rank_Cutoff/Epi")
load("./epi.seurat.objs.2022-07-27.Robj")


# explore_epi_3.R
#### Initial exploration of all Epithelium
## Integration by Reference Mapping
# Omit bc didn't work/use \\


# explore_epi_4.R
#### Exploration of all Epithelium
## Integration by Reference Mapping
# Omit bc didn't work/use \\


# explore_epi_5.R
#### Exploration of Engineered Epithelium
## Re-do CCA Integration
# Keep? or obsolete integration?


# explore_epi_6.R
#### Exploration of Engineered Epithelium
## Embedding-agnostic marker comparisons
# Keep? or non-essential exploration?


# focus_epi_1.R
#### Focused Analysis of Engineered Epithelium
## BEFM1 Gold-standard Evaulation
# Non-essential exploration \\


# focus_epi_2.R
#### Focused Analysis of Engineered Epithelium
## Archetype Analysis of Integrated Cells

## Try grouping clusters
Idents(object) <- object$seurat_clusters
object <- RenameIdents(object,
        '0'='Krt5',
        '1'='Trans',
        '2'='Hopx',
        '3'='Hopx',
        '4'='Krt5',
        '5'='Trans',
        '6'='Krt5',
        '7'='Trans',
        '8'='Cycling',
        '9'='Trans',
        '10'='Cycling',
        '11'='Trans',
        '12'='Other',
        '13'='Hopx',
        '14'='Secretory',
        '15'='Trans')
object$grouping <- Idents(object)

Idents(cl12) <- cl12$seurat_clusters
cl12 <- RenameIdents(cl12,
        '0'='EMT',
        '1'='BASC',
        '2'='Krt5',
        '3'='Endothelium',
        '4'='Endothelium',
        '5'='Immune',
        '6'='BASC',
        '7'='Hopx')
cl12$grouping <- Idents(cl12)


# focus_epi_3.R
#### Focused Analysis of Engineered Epithelium
## Trying to Redo Archetype Analysis of Integrated Cells

## Remove Endo/Imm doublets, recluster (in focus_epi_2.R)
Idents(object) <- object$grouping2
object <- subset(object, idents = c('Endothelium','Immune'), invert = T)
object$grouping2 <- Idents(object)
object$putative_labels <- object$grouping2


# focus_epi_4.R
#### Focused Analysis of Engineered Epithelium
## Archetype Analysis of Integrated Cells, Continued

## Function to query genes of a GO term
GOgenes2 <- function(goid,totalgenelist,totalmappedidlist){
    goid_locale <- sapply(totalgenelist, function(x) grep(paste('"id": "GO:',goid,sep=''), totalgenelist[ ,1]))
    label_locale <- goid_locale + 1
    label <- gsub("                    \"label\": \"", "", totalgenelist[label_locale,])
    label <- gsub(" ","_",label)
    label <- gsub("\"","",label)
    label <- gsub("label:_","",label)
    label <- gsub("-","_",label)
    mappedidlist_locale <- totalmappedidlist[totalmappedidlist<goid_locale, ]
    mappedidlist_locale <- tail(mappedidlist_locale, n = 1)

    mappedgenes <- c()
    for (i in (mappedidlist_locale+1):nrow(totalgenelist)) {
        if (totalgenelist[i, ] != paste("]},")) {
            newgene <- gsub("[^a-zA-Z0-9]", "", totalgenelist[i, ])
            mappedgenes <- c(mappedgenes,newgene)
        } else {
            break
        }
    }
    mappedgenes <- gsub('^Mt', 'Mt-', mappedgenes)
    mappedgenes <<- mappedgenes
    label <<- label
}


## Categorical Investigation - by Archetype
load('./epi.seurat.objs.int.marks.2022-12-15.Robj')
object <- ScaleData(object,features = rownames(object))
cluster <- c('Krt5')
# Single GO Category (per cluster)
go_cat <- c('0031099')

totalgenelist <- as.data.frame(read_excel("epi.seurat.objs.int.marks.2022-12-15.xlsm", sheet = paste(cluster,'(2)')))
totalmappedidlist <- sapply(totalgenelist, function(x) grep(c('mapped_id_list'), totalgenelist[ ,1]))

mygenes <- list() ; totalgenes <- c() ; names <- c()
for (i in 1:length(go_cat)) {
    goid <- go_cat[i]
    GOgenes2(goid=goid,totalgenelist=totalgenelist,totalmappedidlist=totalmappedidlist)
    mygenes[[i]] <- mappedgenes
    totalgenes <- unique(c(totalgenes,mappedgenes))
    names <- c(names,label)
}
names(mygenes) <- names

genes.use <- intersect(totalgenes,rownames(object))
object <- AddModuleScore(object,features=list(genes.use),name=paste(label,'score',sep="."))
Idents(object) <- object$putative_labels

# Prepare plots
go <- as.data.frame(read_excel("epi.seurat.objs.int.marks.2022-12-15.xlsm", sheet = cluster))
go <- go[grep(go_cat, go[,1]), ]
go <- go[,c(1,6,7,8)] 
colnames(go)[2:4] = c('fold Enrichment','raw P-value','FDR')
p1 <- tableGrob(go,theme = ttheme_minimal(base_size = 16))

gomarks <- markers[which(markers$cluster == cluster),]
gomarks <- gomarks[gomarks$gene %in% totalgenes, ]
gomarks <- gomarks[1:25, ]
topmarks <- round(gomarks[, 1:5], digits = 4)
topmarks$cluster <- gomarks$cluster
topmarks$gene <- gomarks$gene
topmarks$ratio <- round(gomarks[1:25, ]$ratio, digits = 4)
topmarks <- na.omit(as.data.frame(topmarks))
p2 <- tableGrob(topmarks,theme = ttheme_minimal(base_size = 16))

toprow <- plot_grid(p1,p2,ncol=2,
    labels=c(paste('Select GO Category of Differentially Upregulated Genes in',cluster),
        paste('Upregulated Genes in',cluster)),label_size = 25,hjust = 0, label_x = 0.01)

p3 <- FeaturePlot(object,paste(label,'score1',sep="."),label=T,repel=T) + labs(title = paste('Engineered Epi -',label,'score')) +
    theme(plot.title = element_text(size = 24))
p4 <- VlnPlot(object,paste(label,'score1',sep="."),pt.size = 0,group.by='Dataset2') + 
    theme(axis.title.x = element_blank()) + NoLegend()
p <- VlnPlot(object,paste(label,'score1',sep="."),pt.size = 0,group.by='putative_labels') + 
    theme(axis.title.x = element_blank()) + NoLegend()
q <- VlnPlot(object,paste(label,'score1',sep="."),pt.size = 0,group.by='Sample') + 
    theme(axis.title.x = element_blank()) + NoLegend()
p5 <- plot_grid(p,q,ncol=1)
midrow <- plot_grid(p3,p4,p5,ncol=3)

# Violin plots of all genes in a Category of GO's
gogenes <- intersect(totalgenes,rownames(object))
cl_marks <- markers[which(markers$cluster == cluster),]
gomarks <- cl_marks[cl_marks$gene %in% gogenes, ]
gomarks <- gomarks[order(-gomarks$avg_log2FC),]
topmarks <- na.omit(gomarks[1:36, ])
genelist <- topmarks$gene

myplots <- vector('list',length(genelist))
for (i in 1:length(genelist)) {
    gene <- genelist[i]
    message(gene)
    myplots[[i]] <- local({
        i <- i
        p1 <- VlnPlot(object,gene,pt.size = 0,group.by='Sample') + xlab('') + ylab('') + NoLegend() +
            theme(axis.text.x=element_blank(),
            axis.text.y = element_text(size = 28, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line(color = "black"),
            plot.title = element_text(size = 38, hjust = 0.5))
    })
}
lastrow <- plot_grid(plotlist=myplots,ncol=6,nrow=6)

# Cowplot all panels
plots <- list(toprow,midrow,lastrow)
png(paste(dir_name,cluster,label,'allplots.png',sep="_"), width = 2000,height=2750)
plot_grid(plotlist=plots,ncol=1,rel_heights=c(1,1,2))
dev.off()


# focus_epi_5.R
#### Focused Analysis of Engineered Epithelium
## Native Benchmark Analyses, to align with other cell classes
# Omit bc didn't work/use \\


# focus_epi_6.R
#### Focused Analysis of Engineered Epithelium
## Pseudotime via Slingshot
# Non-essential exploration \\


# focus_epi_7.R
#### Focused Analysis of Engineered Epithelium
## Remove TXP3_L Sample (output: epi.seurat.objs.int.2023-01-26.Robj (wo TXP_L, *new putative_labels))

## Load integrated object with putative labels
setwd("/nobackup1/greaneya/DGEs_30K_Rank_Cutoff/Epi")
dir_name <- c('epi')
load('./epi.seurat.objs.int.2022-12-15.Robj')  # grouped, w starting cells, regressed

# Remove 'TXP3_L' sample
Idents(object) <- object$Sample
object <- subset(object,idents = 'TXP3_L',invert = T)

object$Sample <- droplevels(object$Sample)
object$Sample <- factor(object$Sample,levels=c('BC1P3','BC1P6','BCL5','BCEC2','BEF1',
    'BEF2','BEF3','BEF12','BEF14','BEF15','BEFM1','BEFM2','BEFM4','BEFM5','BEFM6'))
object$Dataset2 <- droplevels(object$Dataset2)
object$Dataset2 <- factor(object$Dataset2,levels=c('Start','Mono','Co','Tri_E','Tri_L',
    'Quad_E','Quad_L'))


## Perform an integrated analysis on remaining cells
DefaultAssay(object) <- "integrated"
object <- ScaleData(object)
object <- RunPCA(object, npcs = 100, verbose = F)

pcs <- 37
res <- 0.6
tag <- paste('pcs',pcs,'res',res,sep='.')
DefaultAssay(object) <- "integrated"
object <- FindNeighbors(object, dims = 1:pcs)
object <- FindClusters(object, resolution = res)
object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
p1 <- UMAPPlot(object,group.by = 'seurat_clusters',label = T) + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res),caption = paste("nCells =",ncol(object))) +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20),plot.caption = element_text(size = 18,color='red')) + NoAxes()
p2 <- UMAPPlot(object,group.by = 'Sample') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
p3 <- UMAPPlot(object,group.by = 'putative_labels') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
Idents(object) <- object$Sample
toprow <- plot_grid(p1,p2,p3,ncol=3)

data <- prop.table(table(object$seurat_clusters,object$Sample),1)*100
data <- as.data.frame(data)
p5 <- ggplot(data, aes(x = Var1, y = Freq, fill = Var2)) + NoLegend() +
  geom_col(position = "dodge") + scale_fill_manual(values = hue_pal()(length(unique(object$Sample)))) +
  labs(title = c('Cluster Composition by Sample')) + xlab("Cluster") + ylab("% cluster per sample") +
  theme(plot.title = element_text(size = 24),
        legend.title = element_blank(),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
midrow1 <- plot_grid(p5)

DefaultAssay(object) <- "RNA"
Idents(object) <- object$seurat_clusters
p6 <- FeaturePlot(object, c('Krt5','Hopx','Aqp5','Sftpc'),ncol=2,label = F)  # Epi
p7 <- FeaturePlot(object, c('Nkx2-1','Sftpb','Defb4','Lyz2'),ncol=2,label = F)  # ATII
p8 <- FeaturePlot(object, c('Ager','Aqp5','Hopx','Pdpn'),ncol=2,label = F)  # ATI
p9 <- FeaturePlot(object, c('Krt5','Krt8','Igfbp2','Sox2'),ncol=2,label = F)  # Basal
p10 <- FeaturePlot(object, c('Scgb1a1','Rarres1','Muc20','Clca1'),ncol=2,label = F)  # Secretory
p11 <- FeaturePlot(object, c('Dclk1','Espn','Sox9','Pou2f3'),ncol=2,label = F)  # Tuft1
p12 <- FeaturePlot(object, c('Top2a','Trpm5','Sox9','Lgr5'),ncol=2,label = F)  # Tuft2
p13 <- FeaturePlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','Mki67'),ncol=2,label = F)  # QC
p14 <- FeaturePlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'),ncol=2,label = F)  # Lineage
midrow2 <- plot_grid(p6,p7,p8,ncol=3)
midrow3 <- plot_grid(p9,p10,p11,ncol=3)
lastrow <- plot_grid(p12,p13,p14,ncol=3)

# Cowplot all panels
plots <- list(toprow,midrow1,midrow2,midrow3,lastrow)
png(paste(dir_name,tag,'_integrated7_allplots.png',sep=""), width = 2000,height=2750)
plot_grid(plotlist=plots,ncol=1,rel_heights = c(2,1,2,2,2))
dev.off()


## Evaluate new cluster alignment with old putative_labels
table(object$seurat_clusters,object$putative_labels)
# new clusters and old putative_labels align somewhat but not perfectly
# redefine Krt5, Hopx, Trans, but keep Cycling, EMT, BASC, Secretory
# REMAKE PUTATIVE LABELS
Idents(object) <- object$putative_labels
cycling <- subset(object,idents='Cycling')
emt <- subset(object,idents='EMT')
basc <- subset(object,idents='BASC')
# Secretory aligns perfectly with cluster 13
Idents(object) <- object$seurat_clusters
# create new column in metadata
object$sub_cluster <- as.character(object$seurat_clusters)
object$sub_cluster[Cells(cycling)] <- paste(Idents(cycling))
object$sub_cluster[Cells(emt)] <- paste(Idents(emt))
object$sub_cluster[Cells(basc)] <- paste(Idents(basc))

table(object$sub_cluster,object$putative_labels)

p1 <- UMAPPlot(object,group.by = 'sub_cluster') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20)) + NoAxes()
p2 <- UMAPPlot(object,group.by = 'putative_labels') + labs(title = dir_name,subtitle = paste("PC's =",pcs,"\nRes =",res)) + NoAxes() +
        theme(plot.title = element_text(size = 24,hjust = 0),plot.subtitle = element_text(size = 20))
toprow <- plot_grid(p1,p2,ncol=2)
p3 <- FeaturePlot(object, c('Krt5','Hopx'),ncol=2,order=T)
p4 <- VlnPlot(object,c('Krt5','Hopx'),group.by='sub_cluster',ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())

png(paste(dir_name,tag,'_integrated7_subcluster.png',sep=""), width = 1000,height=1500)
plot_grid(plotlist=list(toprow,p3,p4),ncol=1)
dev.off()

Idents(object) <- object$sub_cluster
object <- RenameIdents(object,
                '0'='Hopx',
                '1'='Krt5',
                '2'='Trans',
                '3'='Hopx',
                '4'='Trans',
                '5'='Krt5',
                '6'='Trans',
                '7'='Trans',
                '8'='Krt5',
                '9'='Trans',
                '10'='Krt5',
                '11'='Krt5',
                '12'='Krt5',
                '13'='Secretory')
object$putative_labels <- Idents(object)

table(object$putative_labels)  # confirm
object$putative_labels <- factor(object$putative_labels,levels=c('Krt5','Hopx',
    'Trans','Cycling','EMT','BASC','Secretory'))


# focus_epi_8.R
#### Focused Analysis of Engineered Epithelium
## Re-do Slingshot & Run Velocyto w Pagoda2

### Pseudotime using Slingshot (v2.2.1)
DefaultAssay(object) <- "integrated"  # recommended
slobject <- as.SingleCellExperiment(object,assay = c('integrated','RNA'))

slobject <- slingshot(slobject,reducedDim = reducedDim(slobject,'PCA')[,pcs],
    approx_points = 40,start.clus = 'Cycling',clusterLabels = 'putative_labels')

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(slobject$slingPseudotime_1, breaks=100)]

png(paste(dir_name,'sling_int3_umap_1.png',sep='_'), width = 800, height = 800)
plot(reducedDims(slobject)$UMAP, col = plotcol, pch=16, asp = 1)
dev.off()

save(slobject,file = paste("epi.seurat.objs.sling",Sys.Date(),"Robj",sep="."))


## Pull slingshot pseudotime metadata
sling.metadata <- slingPseudotime(slobject)

sling.metadata2 <- slingCurveWeights(slobject)

object <- AddMetaData(object,metadata = sling.metadata,col.name=colnames(sling.metadata))

png(paste(dir_name,'sling_int3_fp_1.png',sep='_'), width = 800, height = 400)
FeaturePlot(object,c('Lineage1','Lineage2'),ncol=2,label = T,repel=T,order=T)
dev.off()

# Save metadata
sling.metadata <- slingPseudotime(slobject)
colnames(sling.metadata) <- paste0(colnames(sling.metadata),'_ps')

sling.metadata2 <- slingCurveWeights(slobject)
colnames(sling.metadata2) <- paste0(colnames(sling.metadata2),'_cw')

sling.metadata <- data.frame(sling.metadata,sling.metadata2)

save(sling.metadata,file = paste0("sling.metadata.",Sys.Date(),".Robj"))



### Velocyto on Engineered Epithelium
library(velocyto.R)
library('pagoda2')

## Setup with Pagoda2
DefaultAssay(object) <- "integrated"
Idents(object) <- object$Final1
object <- subset(object, downsample = 2000)
object <- FindVariableFeatures(object)
length(object@assays$integrated@var.features)
object <- ScaleData(object)
nrow(GetAssayData(object,slot="scale.data"))

r <- Pagoda2$new(as(GetAssayData(object,slot="scale.data"),"dgCMatrix"),modelType="raw")
r$reductions[['PCA']] <- object@reductions$pca
r$reductions[['umap']] <- object@reductions$umap
r$embeddings[['PCA']][['umap']] <- Embeddings(object[['umap']])

cluster.label <- object$Final1
cell.colors <- sccore:::fac2col(cluster.label,
        level.colors=type_colors[levels(object$Final1)])

png(paste0(dir_name,'_pagumap_3.png'), width = 800, height = 800)
r$plotEmbedding(type='PCA',embeddingType='umap',colors=cell.colors,mark.clusters=T,
        min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,
        main='cell clusters')
dev.off()


## Calculating Velocity
emat <- object$spliced; nmat <- object$unspliced
emat <- emat[,rownames(r$counts)]; nmat <- nmat[,rownames(r$counts)]

cluster.label <- object$Final1
cell.colors <- sccore:::fac2col(cluster.label,
        level.colors=type_colors[levels(object$Final1)])

cell.dist <- as.dist(1-armaCor(t(r$embeddings[['PCA']][['umap']])))

emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.08)
length(intersect(rownames(emat),rownames(nmat)))
# [1] 912

fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,
        cell.dist=cell.dist,fit.quantile=fit.quantile)
emb <- r$embeddings$PCA$umap

x <- show.velocity.on.embedding.cor(emb,rvel.cd,n=200,n.cores=1,scale='sqrt',
        cell.colors=ac(cell.colors,alpha=0.5),cex=0.8,arrow.scale=8,show.grid.flow=T,
        min.grid.cell.mass=5,grid.n=20,arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)

res <- 400
png(paste(dir_name,'_pagvel_3.png',sep=''), width=650*res/72, height=650*res/72,res=res)
show.velocity.on.embedding.cor(emb,rvel.cd,n=800,n.cores=1,scale='sqrt',cc=x$cc,
        cell.colors=ac(cell.colors),cex=1,arrow.scale=15,show.grid.flow=T,
        min.grid.cell.mass=30,grid.n=20,arrow.lwd=4,do.par=F,cell.border.alpha = 0.1,
        xlab = c('UMAP_1'),ylab= c('UMAP_2'))
dev.off()



### Velocyto on Engineered Epithelium - BASC subset focus
# Subset based on embedding
cut <- Embeddings(object[['umap']])[which(Embeddings(object[['umap']])[,'UMAP_1'] > -3),]
cut <- cut[which(cut[,'UMAP_1'] < 0.5),]
cut <- cut[which(cut[,'UMAP_2'] > -4),]
cut <- cut[which(cut[,'UMAP_2'] < 0),]

sub <- subset(object,cells=rownames(cut))  # 628 cells

## Setup with Pagoda2
DefaultAssay(sub) <- "integrated"
Idents(sub) <- sub$Final1
sub <- FindVariableFeatures(sub)
length(sub@assays$integrated@var.features)
sub <- ScaleData(sub)
nrow(GetAssayData(sub,slot="scale.data"))

r <- Pagoda2$new(as(GetAssayData(sub,slot="scale.data"),"dgCMatrix"),modelType="raw")
r$reductions[['PCA']] <- sub@reductions$pca
r$reductions[['umap']] <- sub@reductions$umap
r$embeddings[['PCA']][['umap']] <- Embeddings(sub[['umap']])

cluster.label <- sub$Final1
cell.colors <- sccore:::fac2col(cluster.label,
        level.colors=type_colors[levels(sub$Final1)])

png(paste(dir_name,'_pagumap_sub_1.png',sep=''), width = 800, height = 800)
r$plotEmbedding(type='PCA',embeddingType='umap',colors=cell.colors,mark.clusters=T,
        min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,
        main='cell clusters')
dev.off()


## Calculating Velocity
emat <- sub$spliced; nmat <- sub$unspliced
emat <- emat[,rownames(r$counts)]; nmat <- nmat[,rownames(r$counts)]

cluster.label <- sub$Final1
cell.colors <- sccore:::fac2col(cluster.label,
        level.colors=type_colors[levels(sub$Final1)])

cell.dist <- as.dist(1-armaCor(t(r$embeddings[['PCA']][['umap']])))

emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.08)
length(intersect(rownames(emat),rownames(nmat)))
# [1] 667

fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,
        cell.dist=cell.dist,fit.quantile=fit.quantile)
emb <- r$embeddings$PCA$umap

x <- show.velocity.on.embedding.cor(emb,rvel.cd,n=200,n.cores=1,scale='sqrt',
        cell.colors=ac(cell.colors,alpha=0.5),cex=0.8,arrow.scale=8,show.grid.flow=T,
        min.grid.cell.mass=5,grid.n=20,arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)

res <- 400
png(paste(dir_name,'_pagvel_sub_1.png',sep=''), width=400*res/72, height=400*res/72,res=res)
show.velocity.on.embedding.cor(emb,rvel.cd,n=25,n.cores=1,scale='sqrt',cc=x$cc,
        cell.colors=ac(cell.colors),cex=1,arrow.scale=5,show.grid.flow=T,
        min.grid.cell.mass=5,grid.n=18,arrow.lwd=4,do.par=F,cell.border.alpha = 0.1,
        min.arrow.size = 0.25,xlab = c('UMAP_1'),ylab= c('UMAP_2'))
dev.off()



### Plot velocity on slingshot pseudotime coloring
object <- AddMetaData(object,metadata = sling.metadata,col.name=names(sling.metadata))
lins <- data.frame(object$Lineage1_ps,object$Lineage2_ps)
object$Lineage_avg <- rowMeans(lins,na.rm=TRUE)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(slobject$slingPseudotime_1, breaks=100)]

png(paste(dir_name,'_pagvelpseudo_1.png',sep=''), width = 800, height = 800)
show.velocity.on.embedding.cor(emb,rvel.cd,n=200,n.cores=1,scale='sqrt',
        cell.colors=ac(cell.colors,alpha=0.5),cex=0.8,arrow.scale=8,show.grid.flow=TRUE,
        min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
dev.off()

# Plot velocyto colored by slingshot
cluster.label <- sort(object$Lineage1_ps)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(ncol(object))
cell.colors <- sccore:::fac2col(cluster.label,shuffle=FALSE,
        level.colors=colors)

png(paste(dir_name,'_pagvel_4.png',sep=''), width = 800, height = 800)
show.velocity.on.embedding.cor(emb,rvel.cd,n=1000,n.cores=1,scale='sqrt',cc=x$cc,
        cell.colors=ac(cell.colors,alpha=0.5),cex=1,arrow.scale=12,show.grid.flow=T,
        min.grid.cell.mass=50,grid.n=15,arrow.lwd=4,do.par=F,cell.border.alpha = 0.1)
dev.off()





#### ENDOTHELIUM ####

# explore_endo_1.R
#### Initial exploration of all Endothelium
## Re-clustering individual objects & Gene Expression Heatmaps


# explore_endo_2.R
#### Initial exploration of all Endothelium
## CCA Integration & Merge


# explore_endo_3.R
#### Initial exploration of all Endothelium
## Integration by Reference Mapping


# explore_endo_4.R
#### Exploration of all Endothelium
## Integration by Reference Mapping wo TXP_R's


# explore_endo_5.R
#### Exploration of Engineered Endothelium
## Re-do CCA Integration







#### MESENCHYME ####

# explore_mes_1.R
#### Initial exploration of all Mesenchyme
## Re-clustering individual objects & Gene Expression Heatmaps


# explore_mes_2.R
#### Initial exploration of all Mesenchyme
## CCA Integration & Merge


# explore_mes_3.R
#### Initial exploration of all Mesenchyme
## Integration by Reference Mapping


# explore_mes_4.R
#### Exploration of all Mesenchyme
## Integration by Reference Mapping


# explore_mes_5.R
#### Exploration of Engineered Mesenchyme
## Re-do CCA Integration




#### IMMUNE ####

# explore_imm_1.R
#### Initial exploration of all Immune
## Re-clustering individual objects & Gene Expression Heatmaps


# explore_imm_2.R
#### Initial exploration of all Immune
## Update Cell Labels & CCA Integration


# explore_imm_3.R
#### Initial exploration of all Immune
## Integration by Reference Mapping


# explore_imm_4.R
#### Exploration of all Immune
## Integration by Reference Mapping wo TXP_R's


# explore_imm_5.R
#### Exploration of Engineered Immune
## Re-do CCA Integration








#### Extra investigation, DO NOT INCLUDE

#### Epithelium
# Investigation to trace integration
# Figure object = epi.seurat.objs.int.2023-07-27.Robj
load('./Epi/epi.seurat.objs.int.2023-01-26.Robj')  # made in focus_epi_7.R
opt1 <- object ; rm(object)
load('./Epi/epi.refined.2023-03-06.Robj')  # Sam reclustered and re-defined Final1+2 metadata?
opt2 <- epi ; rm(epi)
load('./Epi/epi.seurat.objs.int.2023-07-27.Robj')  # I brought all together in update_eng_objects.R
opt3 <- object ; rm(object)
# all appear to be same size (52878 cells), so I'll use 2023-01-26...
# except it doesn't have integration, check 2022-12-15 (made in focus_epi_3.R)
load('./Epi/epi.seurat.objs.int.2022-12-15.Robj')  # 52936 cells - contains 58 TXP3_L cells
opt4 <- object ; rm(object)
# still doesn't have integration, check 2022-09-19 (made in explore_epi_5.R)
load('./Epi/epi.seurat.objs.int.2022-09-19.Robj')  # 53176 cells - contains 58 TXP3_L cells and...?
opt5 <- object ; rm(object)
# This suggests the last time an integration was performed on epithelium was explore_epi_5.R,
#   which included TXP3_L and some contaminating endo and imm cells, which were later removed
#   Do we care? Do we just bury it in the code?

