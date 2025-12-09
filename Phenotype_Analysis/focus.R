#### FOCUS
#### Code for publication unifying all 'focus' scripts
#### Replicable methods for focused analyses of classed objects,
####    including archetype labels, slingshot pseudotime, velocyto RNA velocity, ...




## Note: Engineered cell archetypes and cell types were identified through iterative clustering of integrated class
##  objects, combined with targeted transcriptomic characterizations. Iterations of archetype definition are not
##  presented here, but rather the final steps in assigning archetype labels, which should yield metadata that is
##  consistent with the final processed objects uploaded to GEO.




#### EPITHELIUM ####


# focus_epi_7.R
#### Focused Analysis of Engineered Epithelium
## Remove TXP3_L Sample

## OR Just present final archetype definition in processed objects...

### Archetype and cell type definition of integrated object clustering
## Re-integrate epithelial object including starting cells and engineered populations
# Load object containing putative labels from previous analyses
dir_name <- c('epi')
load('./epi.seurat.objs.int.2022-12-15.Robj')

object$Sample <- factor(object$Sample,levels=c('BC1P3','BC1P6','BCL5','BCEC2','BEF1',
    'BEF2','BEF3','BEF12','BEF14','BEF15','BEFM1','BEFM2','BEFM4','BEFM5','BEFM6'))
object$Dataset2 <- factor(object$Dataset2,levels=c('Start','Mono','Co','Tri_E','Tri_L',
    'Quad_E','Quad_L'))

# Re-embed and -cluster on existing integration
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


## Remake putative labels
table(object$seurat_clusters,object$putative_labels)
Idents(object) <- object$seurat_clusters
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

table(object$Sample)  # confirm
table(object$Dataset2)  # confirm


save(object,file = paste("epi.seurat.objs.int.",Sys.Date(),".Robj",sep=""))
load('./epi.seurat.objs.int.2023-01-26.Robj')  # wo TXP_L, *new putative_labels




# focus_epi_8.R
#### Focused Analysis of Engineered Epithelium
## Re-do Slingshot & Run Velocyto w Pagoda2




#### ENDOTHELIUM ####

# focus_endo_1.R
#### Focused Analysis of Engineered Endothelium
## Archetype Analysis of Integrated Cells


# focus_endo_2.R
#### Focused Analysis of Engineered Endothelium
## Archetype Analysis of Integrated Cells


# focus_endo_3.R
#### Focused Analysis of Engineered Endothelium
## Analysis of Endothelium over 


# focus_endo_4.R
#### Focused Analysis of Engineered Endothelium
## Remove Lowinfo cluster from Archetype Analysis


# focus_endo_5.R
#### Focused Analysis of Engineered Endothelium
## Pseudotime via Slingshot


# focus_endo_6.R
#### Focused Analysis of Engineered Endothelium
## Run Velocyto w Pagoda2




#### MESENCHYME ####

# focus_mes_1.R
#### Focused Analysis of Engineered Mesenchyme
## Archetype Analysis of Integrated Cells


# focus_mes_2.R
#### Focused Analysis of Engineered Mesenchyme
## Archetype Analysis of Integrated Cells


# focus_mes_3.R
#### Focused Analysis of Engineered Mesenchyme
## Pseudotime on Archetype Analysis


# focus_mes_4.R
#### Focused Analysis of Engineered Mesenchyme
## Velocity on Archetype Analysis


# focus_mes_5.R
#### Focused Analysis of Engineered Mesenchyme
## Fibrotic Evaluation


# focus_mes_6.R
#### Focused Analysis of Engineered Mesenchyme
## Remove TXP3_L Sample


# focus_mes_7.R
#### Focused Analysis of Engineered Mesenchyme
## Pseudotime via Slingshot


# focus_mes_8.R
#### Focused Analysis of Engineered Mesenchyme
## Run Velocyto w Pagoda2




#### IMMUNE ####

# focus_imm_1.R
#### Focused Analysis of Engineered Immune
## Defining Engineered Cell Types


# focus_imm_2.R
#### Focused Analysis of Engineered Immune
## Immune Phenotype Archetype Analysis to align with other cell classes (attempted w native)


# focus_imm_3.R
#### Focused Analysis of Engineered Immune
## Re-Integrating & Clustering Engineered Immune Cells without any Native


# focus_imm_4.R
#### Focused Analysis of Engineered Immune
## Immune Phenotype Archetype plots, to align with other cell classes


# focus_imm_5.R
#### Focused Analysis of Engineered Immune
## Gold Standard and Native Benchmark Analyses, to align with other cell classes


# focus_imm_6.R
#### Focused Analysis of Engineered Immune
## Remove TXP3_L Sample


# focus_imm_7.R
#### Focused Analysis of Engineered Immune
## Parse Macrophage Heterogeneity


# focus_imm_8.R
#### Focused Analysis of Engineered Immune
## Pseudotime by Slingshot & Velocyto


# focus_imm_9.R
#### Focused Analysis of Engineered Immune
## Re-do Velocyto




#### ADDITIONAL ####

# focus_rus_1.R
#### Focused Analysis guided by Ruslan
## May 2024, Regional ECM expression in engineered epithelium & ECM gene list


# focus_rus_2.R
#### Focused Analysis guided by Ruslan
## May 2024, Generating native v2 object to evaluate regional ECM expression in native


# focus_rus_3.R
#### Focused Analysis guided by Ruslan
## May 2024, Regional ECM expression in native & engineered epithelium


# focus_rus_4.R
#### Focused Analysis guided by Ruslan
## June 2024, Regional ECM expression in native & engineered epithelium WO starting peBC


# focus_rus_5.R
#### Focused Analysis guided by Ruslan
## June 2024, ECM Figure


# focus_rus_6.R
#### Focused Analysis guided by Ruslan
## June 2024, Functional Delegation
