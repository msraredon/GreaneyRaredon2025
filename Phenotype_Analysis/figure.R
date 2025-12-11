#### FIGURE
#### Code for published figure plots
#### Replicable scripts for generating scRNAseq plots in published figures
####    including Figures 2E, 3, 4B, 5, 6D-E, S1-19, & S21-31




#### FIGURE 2 ####
### Panel E - Mmp7 expression across engingeered lungs and in Co
# in engineered lungs, by dataset
p <- VlnPlot(object,'Mmp7',group.by='Dataset4',assay='RNA',slot='data',cols=dataset_colors[levels(object$Dataset4)],pt.size=0) + 
            labs(title = 'Mmp7') + NoLegend() + 
            theme(plot.title = element_text(size = 30,hjust = 0,face = "italic"),
            axis.text.x=element_text(size=14,color='black',angle=0,hjust = 0.5),
            axis.text.y = element_text(size = 14,color='black'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.title.x = element_blank())
res <- 400
png(paste0('all_vln_Mmp7.png'),width = 230*res/72, height = 230*res/72,res=res)
print(p)
dev.off()

# just in Co
Idents(object) <- object$Dataset2
object.sub <- subset(object,idents=c('Co'))
object.sub$Dataset2 <- droplevels(object.sub$Dataset2)
object.sub$class <- droplevels(object.sub$class)

p <- VlnPlot(object.sub,'Mmp7',group.by='class',assay='RNA',slot='data',cols=class_colors[levels(object$class)],pt.size=0) + 
            labs(title = 'Mmp7') + NoLegend() + 
            theme(plot.title = element_text(size = 30,hjust = 0,face = "italic"),
            axis.text.x=element_text(size=14,color='black',angle=0,hjust = 0.5),
            axis.text.y = element_text(size = 14,color='black'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.title.x = element_blank())
res <- 400
png(paste0('all_vln_Co_Mmp7.png'),width = 135*res/72, height = 230*res/72,res=res)
print(p)
dev.off()




#### FIGURE 3 ####
### Phenotype overview by class
## Load samples, colors, set-up
# Load up-to-date epithelium object
load('./Epi/epi.seurat.objs.int.2023-07-27.Robj')
epi <- object ; rm(object)
epi$Dataset4 <- factor(epi$Dataset4,levels=c('Start','Mono','Co','Tri','Quad'))
load('./all.color_palettes_2.R')
sample_colors <- color_pals$sample_colors
dataset_colors <- color_pals$dataset_colors
class_colors <- color_pals$class_colors
epi_colors <- color_pals[['epi']]

# Load up-to-date endothelium object
load('./Endo/endo.seurat.objs.int.2023-07-27.Robj')
endo <- object ; rm(object)
endo$Dataset4 <- factor(endo$Dataset4,levels=c('Start','Co','Tri','Quad'))
Idents(endo) <- endo$Final1
endo$Final1 <- factor(endo$Final1,levels=c('Progenitor','Lymphatic',
    'Cycling_Lymphatic','Cycling_Microvascular','Microvascular'))
endo <- RenameIdents(endo,
        "Cycling_Lymphatic"='Cyc_Lymph',
        "Cycling_Microvascular"='Cyc_MV')
endo$Final3 <- Idents(endo)
endo$Final3 <- factor(endo$Final3,levels=c('Progenitor','Lymphatic',
    'Cyc_Lymph','Cyc_MV','Microvascular'))
endo <- endo
endo_colors <- color_pals[['endo']]
endo_colors2 <- endo_colors
names(endo_colors2) <- c("Microvascular","Cyc_MV","Lymphatic","Cyc_Lymph","Progenitor")

# Load up-to-date mesenchyme object
load('./Mes/mes.seurat.objs.int.2023-07-27.Robj')
mes <- object ; rm(object)
mes$Dataset4 <- factor(mes$Dataset4,levels=c('Start','Tri','Quad'))
mes_colors <- color_pals[['mes']]

# Load up-to-date macrophage object
load('./Imm/mac.seurat.objs.int.2023-07-27.Robj')
imm <- object ; rm(object)
imm$Dataset4 <- factor(imm$Dataset4,levels=c('Start','Tri','Quad'))
starteng_colors <- color_pals$starteng_colors
imm_colors <- color_pals[['imm']]


## Epithelium plots
# UMAP plot of integrated archetypes
if (ncol(epi) < 10000){pt.size=1.5} else if (ncol(epi) >= 10000){pt.size=NULL}
p1 <- UMAPPlot(epi,group.by = 'Final1',cols = epi_colors,pt.size=pt.size,label = T,repel=F,label.size = 6) + NoLegend() +
    labs(title = 'Engineered Epithelium',caption = paste("nCells =",ncol(epi))) + 
    theme(plot.title = element_text(size = 30,hjust = 0,face = "plain"),
        plot.caption = element_text(size = 18,color='red'),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'))

# DotPlot of engineered marker genes
genes <- c('Mki67','Krt5','Krt14','Krt8','Sftpc','Scgb1a1','Hopx','Aqp5')
Idents(epi) <- epi$Final1
epi$Final1 <- factor(epi$Final1,levels=c('Cycling','Basal-like','Transitional','BASC','ATI-like'))
p2 <- DotPlot(epi,features=genes,group.by='Final1',cols=c("lightgrey",class_colors[[1]]),
    dot.scale=10,scale.min=10,scale.max=90,assay='RNA') +
    scale_y_discrete(label=rev(levels(epi$Final1)),limits=rev) +
    guides(color = guide_colorbar(title = 'Avg. Exp.',order=1)) +
    guides(size = guide_legend(order = 2,title = 'Percent\nExp.')) +
    theme(axis.text.x=element_text(size = 18,angle = 45,vjust=1,hjust=1,face = 'italic',color = 'black'),
        axis.text.y = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(color = "black"))

# Featureplot of Slingshot pseudotime
metadata <- epi@meta.data
Pseudotime <- metadata %>%
    group_by(Final1) %>%
    summarize(average_pseudotime = mean(Lineage_avg, na.rm = TRUE))
Pseudotime <- setNames(Pseudotime$average_pseudotime, Pseudotime$Final1)[levels(epi$Final1)]
p3 <- FeaturePlot(epi,'Lineage_avg',min.cutoff=min(Pseudotime),max.cutoff=max(Pseudotime),cols = colorRampPalette(brewer.pal(11,'Spectral')[-6])(60),pt.size=pt.size,order=F) +
    labs(title = 'Pseudotime') +
    theme(plot.title = element_text(size = 30,hjust = 0.2,face = "plain"),
        legend.text = element_text(size = 12)) + NoAxes()

# Featureplot of key archetype genes
p4_1 <- FeaturePlot(epi, c('Krt5','Hopx'),ncol=2,label = F) & NoAxes() & 
    theme(plot.title = element_text(size = 20,face='bold.italic'),
            axis.title.x = element_blank(),
            legend.text = element_text(size = 16))
bascs <- subset(epi,subset = UMAP_1 > -3 & UMAP_1 < 1)
bascs <- subset(bascs,subset = UMAP_2 > -4 & UMAP_2 < 0)
p4_2 <- FeaturePlot(bascs, c('Sftpc','Scgb1a1'),ncol=2,label = F) & 
    theme(plot.title = element_text(size = 20,face='bold.italic'),
            axis.title.x = element_blank(),
            legend.text = element_text(size = 16))
p4 <- plot_grid(p4_1,p4_2,ncol=1)

topleft <- plot_grid(p1,p2,p3,p4,nrow=2)


## Endothelium plots
# UMAP plot of integrated archetypes
if (ncol(endo) < 10000){pt.size=1.5} else if (ncol(endo) >= 10000){pt.size=NULL}
p5 <- UMAPPlot(endo,group.by = 'Final3',cols = endo_colors2,pt.size=pt.size,label = T,repel=F,label.size = 6) + NoLegend() +
    labs(title = 'Engineered Endothelium',caption = paste("nCells =",ncol(endo))) + 
    theme(plot.title = element_text(size = 30,hjust = 0,face = "plain"),
        plot.caption = element_text(size = 18,color='red'),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'))

# DotPlot of engineered marker genes
genes <- c('Vegfa','Atf3','Mmrn1','Mki67','Kdr','Kit','Dll4')
Idents(endo) <- endo$Final3
endo$Final3 <- factor(endo$Final3,levels=c('Progenitor','Lymphatic',
    'Cyc_Lymph','Cyc_MV','Microvascular'))
p6 <- DotPlot(endo,features=genes,group.by='Final3',cols=c("lightgrey",class_colors[[2]]),
    dot.scale=10,scale.min=10,scale.max=90,assay='RNA') +  
    scale_y_discrete(label=rev(levels(endo$Final3)),limits=rev) +
    guides(color = guide_colorbar(title = 'Avg. Exp.',order=1)) +
    guides(size = guide_legend(order = 2,title = 'Percent\nExp.')) +
    theme(axis.text.x=element_text(size = 18,angle = 45,vjust=1,hjust=1,face = 'italic',color = 'black'),
        axis.text.y = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(color = "black"))

# Featureplot of Slingshot pseudotime
metadata <- endo@meta.data
Pseudotime <- metadata %>%
    group_by(Final3) %>%
    summarize(average_pseudotime = mean(Lineage1_ps, na.rm = TRUE))
Pseudotime <- setNames(Pseudotime$average_pseudotime, Pseudotime$Final3)[levels(endo$Final3)]
if (ncol(endo) < 10000){pt.size=1.5} else if (ncol(endo) >= 10000){pt.size=NULL}
p7 <- FeaturePlot(endo,'Lineage1_ps',min.cutoff=min(Pseudotime),max.cutoff=max(Pseudotime),cols = colorRampPalette(brewer.pal(11,'Spectral')[-6])(60),pt.size=pt.size,order=F) +
    labs(title = 'Pseudotime') +
    theme(plot.title = element_text(size = 30,hjust = 0.2,face = "plain"),
        legend.text = element_text(size = 12)) + NoAxes()

# Featureplot of key archetype genes
p8 <- FeaturePlot(endo, c('Kdr','Mmrn1','Ptprb','Atf3'),ncol=2,label = F) & NoAxes() & 
    theme(plot.title = element_text(size = 20,face='bold.italic'),
            axis.title.x = element_blank(),
            legend.text = element_text(size = 16))

topright <- plot_grid(p5,p6,p7,p8,nrow=2)


## Mesenchyme plots
# UMAP plot of integrated archetypes
if (ncol(mes) < 10000){pt.size=1.5} else if (ncol(mes) >= 10000){pt.size=NULL}
p9 <- UMAPPlot(mes,group.by = 'Final1',cols = mes_colors,pt.size=pt.size,label = T,repel=F,label.size = 6) + NoLegend() +
    labs(title = 'Engineered Mesenchyme',caption = paste("nCells =",ncol(mes))) + 
    theme(plot.title = element_text(size = 30,hjust = 0,face = "plain"),
        plot.caption = element_text(size = 18,color='red'),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'))

# DotPlot of engineered marker genes
genes <- c('Mki67','Angpt1','Meox2','Pdgfra','Dcn','Mmp11','Gucy1b1','Notch3','Hsp90b1')
Idents(mes) <- mes$Final1
mes$Final1 <- factor(mes$Final1,levels=c('Cycling','Alveolar','Remodeling',
        'Pericytes','Aberrant'))
p10 <- DotPlot(mes,features=genes,group.by='Final1',cols=c("lightgrey",class_colors[[3]]),
    dot.scale=10,scale.min=10,scale.max=90,assay='RNA') +
    scale_y_discrete(label=rev(levels(mes$Final1)),limits=rev) +
    guides(color = guide_colorbar(title = 'Avg. Exp.',order=1)) +
    guides(size = guide_legend(order = 2,title = 'Percent\nExp.')) +
    theme(axis.text.x=element_text(size = 18,angle = 45,vjust=1,hjust=1,face = 'italic',color = 'black'),
        axis.text.y = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(color = "black"))

# Featureplot of Slingshot pseudotime
metadata <- mes@meta.data
Pseudotime <- metadata %>%
    group_by(Final1) %>%
    summarize(average_pseudotime = mean(Lineage_avg, na.rm = TRUE))
Pseudotime <- setNames(Pseudotime$average_pseudotime, Pseudotime$Final1)[levels(mes$Final1)]
if (ncol(mes) < 10000){pt.size=1.5} else if (ncol(mes) >= 10000){pt.size=NULL}
p11 <- FeaturePlot(mes,'Lineage_avg',min.cutoff=min(Pseudotime),max.cutoff=max(Pseudotime),cols = colorRampPalette(brewer.pal(11,'Spectral')[-6])(60),pt.size=pt.size,order=F) +
    labs(title = 'Pseudotime') +
    theme(plot.title = element_text(size = 30,hjust = 0.2,face = "plain"),
        legend.text = element_text(size = 12)) + NoAxes()

# Featureplot of key archetype genes
p12_1 <- FeaturePlot(mes, c('Pdgfra','Axin2'),ncol=2,label = F,order=T) & NoAxes() & 
    theme(plot.title = element_text(size = 20,face='bold.italic'),
            axis.title.x = element_blank(),
            legend.text = element_text(size = 16))
pericyte <- subset(mes,subset = UMAP_1 > 1 & UMAP_1 < 5)
pericyte <- subset(pericyte,subset = UMAP_2 > 1 & UMAP_2 < 5)
p12_2 <- FeaturePlot(pericyte, c('Gucy1b1'),label = F) +
    theme(plot.title = element_text(size = 20,face='bold.italic'),
            axis.title.x = element_blank(),
            legend.text = element_text(size = 16))
p12_3 <- FeaturePlot(mes, c('Thy1'),label = F,order=T) + NoAxes() + 
    theme(plot.title = element_text(size = 20,face='bold.italic'),
            axis.title.x = element_blank(),
            legend.text = element_text(size = 16))
p12_23 <- plot_grid(p12_2,p12_3,ncol=2)
p12 <- plot_grid(p12_1,p12_23,ncol=1)

botleft <- plot_grid(p9,p10,p11,p12,nrow=2)


## Immune plots
# UMAP plot of integrated archetypes
if (ncol(imm) < 10000){pt.size=1.5} else if (ncol(imm) >= 10000){pt.size=NULL}
p13 <- UMAPPlot(imm,group.by = 'Final35',cols = imm_colors,pt.size=pt.size,shuffle=T,label = F) + 
    labs(title = 'Engineered Macrophages',caption = paste("nCells =",ncol(imm))) + 
    theme(plot.title = element_text(size = 30,hjust = 0,face = "plain"),
        plot.caption = element_text(size = 18,color='red'),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'))

# DotPlot of engineered marker genes
genes <- c('Mcemp1','Naaa','Lipa','Mki67','C1qb','C1qa','C1qc')
Idents(imm) <- imm$Final35
imm$Final35 <- factor(imm$Final35,levels=c('Alveolar','Cycling','Interstitial'))
p14 <- DotPlot(imm,features=genes,group.by='Final35',cols=c("lightgrey",class_colors[[4]]),
    dot.scale=10,scale.min=10,scale.max=90,assay='RNA') + 
    scale_y_discrete(label=rev(levels(imm$Final35)),limits=rev) +
    guides(color = guide_colorbar(title = 'Avg. Exp.',order=1)) +
    guides(size = guide_legend(order = 2,title = 'Percent\nExp.')) +
    theme(axis.text.x=element_text(size = 18,angle = 45,vjust=1,hjust=1,face = 'italic',color = 'black'),
        axis.text.y = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(color = "black"))

# Featureplot of Slingshot pseudotime
if (ncol(imm) < 10000){pt.size=1.5} else if (ncol(imm) >= 10000){pt.size=NULL}
p15 <- FeaturePlot(imm,'Lineage_avg',min.cutoff=min(Pseudotime),max.cutoff=max(Pseudotime),cols = colorRampPalette(brewer.pal(11,'Spectral')[-6])(60),pt.size=pt.size,order=F) + NoAxes() +
    labs(title = 'Pseudotime') +
    theme(plot.title = element_text(size = 30,hjust = 0.2,face = "plain"),
        legend.text = element_text(size = 12)) + NoAxes()

# Featureplot of key archetype genes
p16 <- FeaturePlot(imm, c('Naaa','C1qb','Mt2A','Lcn2'),ncol=2,order=T,label = F) & NoAxes() & 
    theme(plot.title = element_text(size = 20,face='bold.italic'),
            axis.title.x = element_blank(),
            legend.text = element_text(size = 16))

botright <- plot_grid(p13,p14,p15,p16,nrow=2)

# Heatmap of macrophage activation markers (horizontal) - Cells not Average
Idents(imm) <- imm$Final35
imm2 <- subset(imm, downsample = 100)
imm2$Final35 <- factor(imm2$Final35,levels=c('Alveolar','Cycling','Interstitial'))
imm2$Sample <- factor(imm2$Sample,levels=c('FB13','FB14','BAL','BEF12','BEFM1',
        'BEFM2','BEFM4','BEFM5','BEFM6'))
imm2$Dataset3 <- factor(imm2$Dataset3,levels=c('Start','Eng'))
imm2 <- ScaleData(imm2,features = rownames(imm2))
Idents(imm2) <- imm2$Final35
genes2 <- c('Bst2','Cd9','Xrcc5','Cdc42ep3','Pla2g2d','Slc39a2','Lpl','Rgcc','Ly6al',
        'Crip1','S100a4','Scgb1a1',  # Start_Alveolar
        'Stmn1','Tuba1b','Lgals1','Tmem176b','Tmem176a','Limd2','Ccl2','Cd14','Ccl7',
        'Ccl12','Spp1','Ms4a7',  # Start_Interstitial
        'Slc11a1','Fabp4','Mt1','Mt2A','Cxcl2','Cxcl1','Fos','Clec10a','Vsig4','Retn',
        'Clca4l','S100a9','Lcn2','Mgp',  # Eng_Alveolar
        'Ctsl','Ifitm1','Hopx','Pla2g7','Tf','Slpi','Gsta1','F13a1','Gpx3','Sparc',
        'Col1a1')  # Eng_Interstitial 
gene_split2 <- c(rep('Start Alveolar',12),rep('Start Interstitial',12),
    rep('Engineered Alveolar',14),rep('Engineered Interstitial',11))
genes.use2 <- intersect(rownames(GetAssayData(imm2,slot = 'scale.data')),genes2)
imm.cells2 <- GetAssayData(imm2, slot = "scale.data")[genes.use2, ]
imm.cells2 <- t(scale(t(imm.cells2)))
unityNormalize <- function(x){(x-min(x))/(max(x)-min(x))}
imm.cells2 <- t(apply(as.matrix(imm.cells2), 1, unityNormalize))
cellorder2 <- data.frame(Dataset=imm2$Dataset3,Final35=imm2$Final35,Sample=imm2$Sample)
cellorder2 <- cellorder2[order(cellorder2[,1],cellorder2[,2],cellorder2[,3]),]
imm.cells2 <- imm.cells2[,rownames(cellorder2)]
imm.cells2 <- imm.cells2[genes2,]
dim(imm.cells2)

Dataset <- as.matrix(imm2$Dataset3)
Dataset <- as.matrix(Dataset[rownames(cellorder2),])
dataset_split <- c(rep('Start',301),rep('Eng',99))
dataset_split2 <- c(rep('Start_Alveolar',150),rep('Start_Interstitial',76),
    rep('Eng_Alveolar',50),rep('Eng_Interstitial',24))
CellType <- as.matrix(imm2$Final35)
CellType <- as.matrix(CellType[rownames(cellorder2),])
Sample <- as.matrix(imm2$Sample)
Sample <- as.matrix(Sample[rownames(cellorder2),])
colors.inferno2 <- colorRamp2(breaks = c(seq(min(imm.cells2),max(imm.cells2),length.out=60)), inferno(n=60), space = "RGB")
ra2 = rowAnnotation(df = list(data.frame(Dataset),
                            data.frame(CellType),
                            data.frame(Sample)),
    col = list(Dataset=starteng_colors,
                CellType=imm_colors,
                Sample=sample_colors),
    show_annotation_name = c(T),show_legend = c(F),
    annotation_name_gp = gpar(fontsize = 14),
    annotation_name_side = 'top',
    annotation_name_rot = 60)
lgd1 = Legend(at = levels(imm2$Dataset3),
    legend_gp = gpar(fill=starteng_colors),
    title = "Dataset",
    title_gp = gpar(fontface = "bold",fontsize = 16),
    labels_gp = gpar(fontsize = 14),
    grid_width = unit(0.75, "cm"),
    grid_height = unit(0.75, "cm"))
lgd2 = Legend(at = levels(imm2$Final35),
    legend_gp = gpar(fill=imm_colors[levels(imm2$Final35)]),
    title = "Cell Type",
    title_gp = gpar(fontface = "bold",fontsize = 16),
    labels_gp = gpar(fontsize = 14),
    grid_width = unit(0.75, "cm"),
    grid_height = unit(0.75, "cm"))
lgd3 = Legend(at = levels(imm2$Sample),
    legend_gp = gpar(fill=sample_colors[levels(imm2$Sample)]),
    title = "Sample",
    title_gp = gpar(fontface = "bold",fontsize = 16),
    labels_gp = gpar(fontsize = 14),
    grid_width = unit(0.75, "cm"),
    grid_height = unit(0.75, "cm"))
lgd4 = Legend(col_fun = colors.inferno2,
    title = "Gene\nExp.", 
    title_gp = gpar(fontsize = 16),
    labels = c("Low","High"), at = c(0,1),
    labels_gp = gpar(fontsize = 14),
    grid_width = unit(0.75, "cm"),
    legend_height = unit(3, "cm"))
pd2 = packLegend(lgd1,lgd2,lgd3,column_gap = unit(1, "cm"),max_height = unit(20, "cm"))

imm.heatmap2 <- Heatmap(as.matrix(t(imm.cells2)),
        col = colors.inferno2,
        left_annotation = ra2,
        row_split = factor(dataset_split2,levels=levels(imm2$Final4)),
        column_split = factor(gene_split2,levels=unique(gene_split2)),
        row_title_gp = gpar(fontsize = 0),
        column_title_gp = gpar(fontsize = 16),
        column_names_gp = gpar(fontsize = 12,fontface='italic'),
        column_names_rot = 45,
        cluster_rows=F,
        cluster_columns=F,
        cluster_column_slices=F,
        show_row_names=F,
        show_column_names = T,
        show_heatmap_legend = F)
p14 = grid.grabExpr(draw(imm.heatmap2,annotation_legend_list = pd2,annotation_legend_side="left",
    heatmap_legend_list = lgd4, heatmap_legend_side="right", padding = unit(c(2, 2, 10, 2), "mm")))

# High resolution png
res <- 400
png(paste('imm_figure_act_1.png'),width=1100*res/72,height=500*res/72,res=res)
plot_grid(p14)
dev.off()


## Cowplot all panels - with space for images
toprow <- plot_grid(topleft,NULL,topright,ncol=3,rel_widths=c(9,1,9))
botrow <- plot_grid(botleft,NULL,botright,ncol=3,rel_widths=c(9,1,9))
plots <- plot_grid(toprow,NULL,NULL,botrow,nrow=4,rel_heights=c(16,8,1,16))

res <- 400
png(paste('all_figure_pheno_12.png'),width=1800*res/72,height=1875*res/72,res=res)
plot_grid(plots)
dev.off()




#### FIGURE 5 ####
### Functional Delegation plots
## Load engineered object of Tri and Quad engineered lungs and moved gene list
load('transcriptome.TriQuad.customEmbed.fineFilter.proofRead.2024-03-15.Robj')
object <- sub.fine ; rm(sub.fine)
object$Dataset2 <- factor(object$Dataset2,levels=c('Tri_L','Quad_E'))
load('befm.goi.annotated.Robj')
befm.goi <- befm.goi.annotated ; rm(befm.goi.annotated)
load('./all.color_palettes_2.R')
sample_colors <- color_pals$sample_colors
dataset_colors <- color_pals$dataset_colors
class_colors <- color_pals$class_colors


## Volcano plots of moved genes
# Total genes moved to Immune (left Tri, right Quad)
moved <- befm.goi[which(befm.goi$moved == 'TRUE'),]
moved <- moved[which(moved$top.in.quad == 'Immune'),]
nrow(moved)  # 1918

# Add unified avg_log2FC
moved.sub <- moved[c('GOI','top.in.tri','top.in.quad','avg_log2FC_quad','avg_log2FC_tri','p.val.top.node_quad',
    'p.val.top.node_tri','global.delta')]
moved.quad <- moved.sub
moved.quad$avg_log2FC <- moved.quad$avg_log2FC_quad
moved.quad$p.val.top.node <- moved.quad$p.val.top.node_quad
moved.quad$class <- moved.quad$top.in.quad
moved.tri <- moved.sub
moved.tri$avg_log2FC <- -moved.tri$avg_log2FC_tri
moved.tri$p.val.top.node <- moved.tri$p.val.top.node_tri
moved.tri$class <- moved.tri$top.in.tri

moved.volc <- rbind(moved.quad,moved.tri)
moved.volc$class <- factor(moved.volc$class,levels=names(class_colors))
moved.volc <- moved.volc[order(moved.volc$class),]
nrow(moved.volc)  # 3836

# Pull GOI by GO category
GOgenes2 <- function(goid,totalgenelist,totalmappedidlist){
    goid_locale <- sapply(totalgenelist, function(x) grep(paste('"id": "GO:',goid,sep=''), totalgenelist[ ,1]))
    label_locale <- goid_locale + 1
    label <- gsub("                    \"label\": \"", "", totalgenelist[label_locale,])
    label <- gsub(" ","_",label)
    label <- gsub("\"","",label)
    label <- gsub("label:_","",label)  # added
    label <- gsub("-","_",label)  # added
    mappedidlist_locale <- totalmappedidlist[totalmappedidlist<goid_locale, ]
    mappedidlist_locale <- tail(mappedidlist_locale, n = 1)

    mappedgenes <- c()
    for (i in (mappedidlist_locale+1):nrow(totalgenelist)) {
        if (totalgenelist[i, ] != paste("]},")) {   # shortened
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


## Lipid Metabolism
go_cat <- c('0006629','0044255','0008610','0016042','0010876','0034440','0006869','0030258','0010883',
   '0055088','0009395','1905952','0045834','0045834','0097384','0006638','0032365','0046460','0006664','1903509',
   '0006692','0046890','0097006','0008611','0010884','0019915','0015914')

totalgenelist <- as.data.frame(read_excel("befm.goi.annotated.formatted_AMG.xlsm", sheet = c('to_imm (2)')))
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
length(totalgenes)

goi <- intersect(moved.volc$GOI,totalgenes)
moved.volc.goi <- moved.volc[moved.volc$GOI %in% goi,]
nrow(moved.volc.goi)/2 == length(goi)

# significant quad labels
moved.volc.goi <- moved.volc.goi[order(moved.volc.goi$p.val.top.node_quad),]
labels <- unique(head(moved.volc.goi$GOI,20))

random.goi <- moved.volc.goi[sample(nrow(moved.volc.goi)),]
res <- 400
png(paste('befm.goi_volcano_22.png'), width = 400*res/72,height=400*res/72,res=res)
EnhancedVolcano(random.goi,
    title = 'Lipid Metabolism',
    subtitle = ' ',
    caption = paste('nGenes =',nrow(random.goi)/2),
    lab=random.goi$GOI,
    selectLab=labels,
    labSize = 4.0,
    labFace = 'italic',
    x='avg_log2FC',
    y='p.val.top.node',
    pCutoff = 10e-50,
    FCcutoff = 2,
    cutoffLineType = 'blank',
    cutoffLineCol = 'black',
    colCustom=class_colors[random.goi$class],
    vline = 0,
    drawConnectors = T,
    legendPosition = 'none',
    gridlines.major = F,
    gridlines.minor = F) +
    ggplot2::theme(plot.title = element_text(size = 24,hjust = 0,face = "plain"),
        plot.caption = element_text(size = 14,color='red'),
        axis.text.x=element_text(size=14,color='black'),
        axis.text.y = element_text(size = 14,color='black'),
        axis.title.x = element_text(size = 16,color='black'),
        axis.title.y = element_text(size = 16,color='black'))
dev.off()


## Endocytosis
go_cat <- c('0006897','0098657','0032456','0007032','0034058','0006898','0006909','0030100','0001845','0010508',
   '0045807','0090385','0043277','0050764')

totalgenelist <- as.data.frame(read_excel("befm.goi.annotated.formatted_AMG.xlsm", sheet = c('to_imm (2)')))
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
length(totalgenes)

goi <- intersect(moved.volc$GOI,totalgenes)
moved.volc.goi <- moved.volc[moved.volc$GOI %in% goi,]
nrow(moved.volc.goi)/2 == length(goi)

# significant quad labels
moved.volc.goi <- moved.volc.goi[order(moved.volc.goi$p.val.top.node_quad),]
labels <- unique(head(moved.volc.goi$GOI,20))

random.goi <- moved.volc.goi[sample(nrow(moved.volc.goi)),]
res <- 400
png(paste('befm.goi_volcano_31.png'), width = 400*res/72,height=400*res/72,res=res)
EnhancedVolcano(random.goi,
    title = 'Endocytosis',
    subtitle = ' ',
    caption = paste('nGenes =',nrow(random.goi)/2),
    lab=random.goi$GOI,
    selectLab=labels,
    labSize = 4.0,
    labFace = 'italic',
    x='avg_log2FC',
    y='p.val.top.node',
    pCutoff = 10e-50,
    FCcutoff = 2,
    cutoffLineType = 'blank',
    cutoffLineCol = 'black',
    colCustom=class_colors[random.goi$class],
    vline = 0,
    drawConnectors = T,
    max.overlaps = 100,
    legendPosition = 'none',
    gridlines.major = F,
    gridlines.minor = F) +
    ggplot2::theme(plot.title = element_text(size = 24,hjust = 0,face = "plain"),
        plot.caption = element_text(size = 14,color='red'),
        axis.text.x=element_text(size=14,color='black'),
        axis.text.y = element_text(size = 14,color='black'),
        axis.title.x = element_text(size = 16,color='black'),
        axis.title.y = element_text(size = 16,color='black'))
dev.off()


## Violin plots of moved genes
Idents(object) <- object$Condition
object <- RenameIdents(object,'Tri_L'='Tri','Quad_E'='Quad')
object$Condition <- Idents(object)

moved.volc.goi <- moved.volc.goi[order(moved.volc.goi$p.val.top.node_quad),]
moved.volc.goi <- moved.volc.goi[order(-moved.volc.goi$global.delta),]
labels <- unique(head(moved.volc.goi$GOI,20))

# updated environment to seurat 4.3.2 for plotting vlnplots
labels <- c('Alox5ap','Coro1a','Il18','Msr1','Myo1f','Naaa','Syk','Lpl','Igf1','Hgf','Il6r','Il10rb','Il20rb')
res <- 400
for (i in 1:length(labels)){
    gene <- labels[i]
    message(gene)

    p <- VlnPlot(object,gene,group.by='Condition',split.by='class',log=F,cols=class_colors[levels(object$class)],pt.size=0.1,add.noise=F) + 
        labs(title = gene) + NoLegend() + 
        theme(plot.title = element_text(size = 30,hjust = 0,face = "italic"),
        axis.text.x=element_text(size=14,color='black',angle=0,hjust = 0.5),
        axis.text.y = element_text(size = 14,color='black'),
        axis.title.y = element_text(size = 16,color='black'),
        axis.title.x = element_blank()) +
        geom_vline(xintercept=1.5,linetype="longdash",color = "black")

    png(paste0('befm.goi_movedvln_',gene,'_3.png'), width = 300*res/72, height = 300*res/72,res=res)
    print(p)
    dev.off()
}

png(paste0('befm.goi_movedvln_legend_2.png'), width = 500*res/72, height = 250*res/72,res=res)
VlnPlot(object,gene,group.by='Condition',split.by='class',log=T,cols=class_colors[levels(object$class)],pt.size=0.1,add.noise=F)
dev.off()




#### FIGURE 6 ####
### Panels D-E - Endothelial and Epithelial barrier gene expression
## Load up-to-date Epi & Endo objects
load('./all.color_palettes_2.R')
load('./Epi/epi.seurat.objs.int.2023-07-27.Robj')
epi <- object ; rm(object)
load('./Endo/endo.seurat.objs.int.2023-07-27.Robj')
endo <- object ; rm(object)
Idents(endo) <- endo$Final1
endo <- RenameIdents(endo,
        "Cycling_Lymphatic"='Cyc_Lymph',
        "Cycling_Microvascular"='Cyc_MV')
endo$Final3 <- Idents(endo)
endo$Final3 <- factor(endo$Final3,levels=c('Progenitor','Lymphatic',
    'Cyc_Lymph','Cyc_MV','Microvascular'))
dataset_colors <- color_pals$dataset_colors


## Violin Plot Epi
Idents(epi) <- epi$Dataset2
epi2 <- subset(epi,idents=c('Tri_E','Quad_L'),invert=T)
epi2$Dataset2 <- droplevels(epi2$Dataset2)
Idents(epi2) <- epi2$Dataset2
epi2 <- RenameIdents(epi2,'Tri_L'='Tri','Quad_E'='Quad')
epi2$Dataset4 <- Idents(epi2)
epi2$Dataset4 <- factor(epi2$Dataset4,levels=c('Start','Mono','Co',
    'Tri','Quad'))
names(dataset_colors) <- c('Start','Mono','Co','Tri_E','Tri','Quad','Quad_L','TXP_L')

epi2 <- ScaleData(epi2,features = rownames(epi2))
epi.barrier <- c('Cldn3','Cldn4','Cldn7','Cldn23','Ocln','Tjp1','Tjp2','Tjp3','Afdn',
    'Rhoc','Actn4','F11r','Anxa1','Anxa2','Anxa3','Anxa9','Rac1','Cd24','Lyn','Erbb3',
    'Ceacam1','Epcam','Fn1','Iqgap1','Col1a1','Pdpn','Myadm','Flna','Nedd9','Nedd4',
    'Gja1','Cdh1','Pkp1','Tln1','Cast','C3','Col4a1','Ramp2','Ndrg1','Ccl2','Itgb1','Cd63')
genes.use <- intersect(rownames(GetAssayData(epi2,slot = 'scale.data')),epi.barrier)
epi2 <- AddModuleScore(epi2, features = list(genes.use),name = 'barrier.score')
Idents(epi2) <- epi2$Dataset4

res <- 400
png('epi_barriervln_1.png', width = 500*res/72, height = 400*res/72,res=res)
VlnPlot(epi2,'barrier.score1',split.by='Dataset4',log=F,cols=dataset_colors[levels(epi2$Dataset4)],pt.size=0,add.noise=F) + 
        labs(title = c('Epithelial Barrier Gene Score')) + ylab('Barrier Gene Score') + NoLegend() + 
        theme(plot.title = element_text(size = 28,hjust = 0,face = 'plain'),
        axis.text.x=element_text(size=16,color='black',angle=0,hjust = 0.5),
        axis.text.y = element_text(size = 16,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.title.x = element_blank()) +
        geom_boxplot(aes(color=epi2$class),width=0.4,outlier.color='black') +
        scale_color_manual(values = c('white','white','white','white','white'))
dev.off()


## Violin Plot Endo
Idents(endo) <- endo$Dataset2
endo2 <- subset(endo,idents=c('Tri_E','Quad_L'),invert=T)
endo2$Dataset2 <- droplevels(endo2$Dataset2)
Idents(endo2) <- endo2$Dataset2
endo2 <- RenameIdents(endo2,'Tri_L'='Tri','Quad_E'='Quad')
endo2$Dataset4 <- Idents(endo2)
endo2$Dataset4 <- factor(endo2$Dataset4,levels=c('Start','Co',
    'Tri','Quad'))
names(dataset_colors) <- c('Start','Mono','Co','Tri_E','Tri','Quad','Quad_L','TXP_L')

endo2 <- ScaleData(endo2,features = rownames(endo2))
endo.barrier <- c('Myl12a','Myl12b','Afdn','Tjp1','Ptprm','Ets1','Tie1','Ramp2','Ramp3',
    'S1pr1','Myct1','Gpr4','Lyl1','Rasip1','Cdh5','Thsd1','Gja5','Esam','Fzd4','Foxf1',
    'Erg','Podxl','Itgb1','Pxn','Gja4','Dll4','Notch4','Cd34','Fgd5','Angpt2','Clec14a',
    'Cldn15','Adam15','Cd93','Egfl7','Ephb4','Jam2','Lama5','Lyn','Msn','Notch1',
    'Tek')  # Curated from main investigation & GO (all support vascular barrier)
genes.use <- intersect(rownames(GetAssayData(endo2,slot = 'scale.data')),endo.barrier)
endo2 <- AddModuleScore(endo2, features = list(genes.use),name = 'barrier.score')
Idents(endo2) <- endo2$Dataset4

res <- 400
png('endo_barriervln_1.png', width = 500*res/72, height = 400*res/72,res=res)
VlnPlot(endo2,'barrier.score1',split.by='Dataset4',log=F,cols=dataset_colors[levels(endo2$Dataset4)],pt.size=0,add.noise=T) + 
        labs(title = c('Endothelial Barrier Gene Score')) + ylab('Barrier Gene Score') + NoLegend() + 
        theme(plot.title = element_text(size = 28,hjust = 0,face='plain'),
        axis.text.x=element_text(size=16,color='black',angle=0,hjust = 0.5),
        axis.text.y = element_text(size = 16,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.title.x = element_blank()) +
        geom_boxplot(aes(color=endo2$class),width=0.4,outlier.color='black') +
        scale_color_manual(values = c('white','white','white','white'))
dev.off()




#### FIGURES S1-19 ####
### Data Cleaning Supplement for each object
# Load up list of original, uncleaned objects & final integrated object containing all cell classes
load('./bef.seurat.objs.2022-06-01.Robj')
load('./bef.metadata2022-06-01.Robj')
load('./all.color_palettes_2.R')
sample_colors <- color_pals$sample_colors
dataset_colors <- color_pals$dataset_colors
class_colors <- color_pals$class_colors
load('./bef.seurat.objs.int.2023-04-06.Robj')
Idents(object) <- object$Final2
object <- RenameIdents(object,
    'Cycling_Alveolar_Macrophage'='Cycling_Mac_Alv',
    'Alveolar_Macrophage'='Mac_Alv',
    'Interstitial_Macrophage'='Mac_Int',
    'Cycling_Lymphatic'='Cycling_Endo',
    'Cycling_Microvascular'='Cycling_Endo')
object$Final3 <- Idents(object)
object$Final3 <- factor(object$Final3,levels=c(
    'Basal-like','Transitional','ATI-like','BASC','Cycling_Epi',
    'Microvascular','Endo_Progenitor','Lymphatic','Cycling_Endo',
    'Alveolar_Fibroblast','Aberrant_Fibroblast','Remodeling_Fibroblast','Pericytes','Cycling_Mes',
    'Mac_Alv','Mac_Int','Cycling_Mac_Alv'))
object$class <- factor(object$class,levels=c('Epithelium','Endothelium','Mesenchyme','Immune'))


## Custom plotting functions
Barcode_rank_plot <- function(sample.name,file.name){
    xmin <- 1 ; xmax <- 100000
    # Load UMIperCellSorted
    UMIperCellSorted <- read.table(paste0(file.name,"_Gene_UMIperCellSorted.txt"), quote="\"", comment.char="")
    # Make a knee plot
    counts <- UMIperCellSorted$V1
    counts <- sort(counts, decreasing = TRUE)
    ranks <- seq(length(counts))
    counts[duplicated(counts)] <- NA
    ranks[duplicated(counts)] <- NA
    datf <- data.frame(Rank = ranks, Count = counts)
    datf <- datf[!is.na(datf[, 2]), , drop = FALSE]
    # Define cutoff point in terms of nUMI
    cutoff <- UMIperCellSorted[30000,]
    # Plot
    p <- ggplot(datf, aes_string(x = "Rank", y = "Count")) +
        geom_point() +
        scale_x_continuous(trans = "log10",
                            breaks=c(10, 100, 500, 1000, 5000, 10000, 30000, 50000, 100000, 500000),
                            limits=c(xmin, xmax)) +
        scale_y_continuous(name = "nCount_RNA", trans = "log10",
                            breaks=c(1, 10, 100, 150, 200, 500, 1000, 5000, 10000, 25000, 50000),
                            limits=c(1, 500000)) +
        geom_vline(xintercept=30000, linetype="dashed", color="red")  +
        geom_hline(yintercept=cutoff, linetype="dashed", color="blue")  +
        theme_minimal(base_size = 18) +
        theme(plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(sample.name)+
        annotate("text", label = paste("Rank Threshold =",30000), x = xmin, y = 4,hjust=0)+
        annotate("text", label = paste("nUMI Threshold =",cutoff), x = xmin, y = 2,hjust=0)
    p1 <<- p
}

QCViolinPlotter <- function(sample.name,file.name,min.counts,min.features,max.mt){
    meta.object <- subset(bef.metadata,Sample %in% file.name)
    x <- which(file.name == names(seurat.data))
    object <- seurat.data[[x]]
    filtered <- subset(object,nFeature_RNA > min.features & nCount_RNA > min.counts & percent.mt < max.mt)
    p1 <- ggplot(meta.object, aes(x=Sample, y=nCount_RNA, fill=Sample)) + geom_violin(trim=F) + scale_fill_manual(values=sample_colors[[sample.name]]) +
        ggtitle('nCount_RNA') + geom_jitter(size=0.1) + scale_y_log10() + theme_classic() + NoLegend() +
        geom_hline(yintercept=min.counts, linetype='dashed', color='red', size=1) +
        labs(subtitle = paste("nUMI >",min.counts),caption = paste(" ")) +
        theme(plot.title = element_text(size = 24),
            plot.subtitle = element_text(size=20,color='red'),
            plot.caption = element_text(size=20,color='red'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.title.x = element_blank(),
            axis.text.x = element_blank())
    p2 <- ggplot(meta.object, aes(x=Sample, y=nFeature_RNA, fill=Sample)) + geom_violin(trim=F) + scale_fill_manual(values=sample_colors[[sample.name]]) +
        ggtitle('nFeature_RNA') + geom_jitter(size=0.1) + scale_y_log10() + theme_classic() + NoLegend() +
        geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
        labs(subtitle = paste("nGene >",min.features),caption = paste(" ")) +
        theme(plot.title = element_text(size = 24),
            plot.subtitle = element_text(size=20,color='red'),
            plot.caption = element_text(size=20,color='red'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.title.x = element_blank(),
            axis.text.x = element_blank())
    p3 <- ggplot(meta.object, aes(x=Sample, y=percent.mt, fill=Sample)) + geom_violin(trim=F) + scale_fill_manual(values=sample_colors[[sample.name]]) +
        ggtitle('percent.mt') + geom_jitter(size=0.1) + ylim(0,30) + theme_classic() + NoLegend() +
        geom_hline(yintercept=max.mt, linetype='dashed', color='red', size=1) +
        labs(subtitle = paste("pt.mt <",max.mt)) +
        labs(caption = paste("Barcodes retained with filtration / starting total:",ncol(filtered),"/",ncol(object))) +
        theme(plot.title = element_text(size = 24),
            plot.subtitle = element_text(size=20,color='red'),
            plot.caption = element_text(size=20,color='red'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.title.x = element_blank(),
            axis.text.x = element_blank())
    p2 <<- p1 ; p3 <<- p2 ; p4 <<- p3
}

TriplePercentMtPlot <- function(sample.name){
    meta.object <- subset(bef.metadata,Sample %in% file.name)
    p <- ggplot() +
        geom_point(data=meta.object, aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
        scale_color_viridis(limits = c(0, 30), oob = scales::squish) + theme_classic() +
        labs(title = 'percent.mt \u2014 pre-filter') +
        labs(subtitle = paste('')) +
        theme(plot.title = element_text(size = 24),
            plot.subtitle = element_text(size=20,color='red'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.title.x = element_text(size = 16,color='black'),
            axis.text.x = element_text(size=12,color='black'))
    p5 <<- p
}

TriplePercentMtPlotFilter <- function(sample.name,min.counts,min.features,max.mt){
    meta.object <- subset(bef.metadata,Sample %in% file.name)
    lowpass <- subset(meta.object,percent.mt < max.mt)
    highmt <- subset(meta.object,percent.mt > max.mt)
    filtered <- subset(meta.object,nFeature_RNA > min.features & nCount_RNA > min.counts & percent.mt < max.mt)
    p <- ggplot() +
        geom_point(data=highmt, aes(x=nCount_RNA, y=nFeature_RNA), color = 'gray') +
        geom_point(data=lowpass, aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
        scale_color_viridis(limits = c(0, 30), oob = scales::squish) + theme_classic() +
        geom_vline(xintercept=min.counts, linetype='dashed', color='red', linewidth=1) +
        geom_hline(yintercept=min.features, linetype='dashed', color='red', linewidth=1) +
        labs(title = 'percent.mt \u2014 filtered') +
        labs(subtitle = paste("nUMI_RNA >",min.counts,"& nGene_RNA >",min.features,"& pt.mt <",max.mt)) +
        theme(plot.title = element_text(size = 24),
            plot.subtitle = element_text(size=18,color='red'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.title.x = element_text(size = 16,color='black'),
            axis.text.x = element_text(size=12,color='black'))
    p6 <<- p
}

TriplePercentSplPlot <- function(sample.name){
    meta.object <- subset(bef.metadata,Sample %in% file.name)
    p <- ggplot() +
        geom_point(data=meta.object, aes(x=nCount_RNA, y=nFeature_RNA, color=perc.spliced)) +
        scale_color_viridis(option="magma",direction=-1,limits = c(50, 100), oob = scales::squish) + theme_classic() +
        labs(title = 'percent.spliced \u2014 pre-filter') +
        labs(subtitle = paste('')) +
        theme(plot.title = element_text(size = 24),
            plot.subtitle = element_text(size=20,color='red'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.title.x = element_text(size = 16,color='black'),
            axis.text.x = element_text(size=12,color='black'))
    p7 <<- p
}

TriplePercentSplPlotFilter <- function(sample.name,min.counts,min.features,max.mt){
    meta.object <- subset(bef.metadata,Sample %in% file.name)
    lowpass <- subset(meta.object,percent.mt < max.mt)
    highmt <- subset(meta.object,percent.mt > max.mt)
    filtered <- subset(meta.object,nFeature_RNA > min.features & nCount_RNA > min.counts & percent.mt < max.mt)
    p <- ggplot() +
        geom_point(data=highmt, aes(x=nCount_RNA, y=nFeature_RNA), color = 'gray') +
        geom_point(data=lowpass, aes(x=nCount_RNA, y=nFeature_RNA, color=perc.spliced)) +
        scale_color_viridis(option="magma",direction=-1,limits = c(50, 100), oob = scales::squish) + theme_classic() +
        geom_vline(xintercept=min.counts, linetype='dashed', color='red', size=1) +
        geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
        labs(title = 'percent.spliced \u2014 filtered') +
        labs(subtitle = paste("nUMI_RNA >",min.counts,"& nGene_RNA >",min.features,"& pt.mt <",max.mt)) +
        theme(plot.title = element_text(size = 24),
            plot.subtitle = element_text(size=18,color='red'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.title.x = element_text(size = 16,color='black'),
            axis.text.x = element_text(size=12,color='black'))
    p8 <<- p
}

SeuratCluster <- function(sample.name,object,pcs){
    DefaultAssay(object) <- 'RNA'
    object <- NormalizeData(object)
    object <- FindVariableFeatures(object)
    object <- ScaleData(object)
    object <- RunPCA(object, npcs = 50, verbose = F)
    object <- FindNeighbors(object, dims = 1:pcs)
    object <- FindClusters(object, resolution = 0.5)
    object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
    p <- UMAPPlot(object,group.by='class',cols=class_colors,label = T) + NoLegend() + 
        labs(title = paste(sample.name,'final \u2014 class'),subtitle = paste("PC's =",pcs),caption = paste("nCells =",ncol(object))) + 
        theme(plot.title = element_text(size = 24,hjust=0),
            plot.subtitle = element_text(size=20,color='red'),
            plot.caption = element_text(size=20,color='red'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.title.x = element_text(size = 16,color='black'),
            axis.text.x = element_text(size=12,color='black'))
    cluster.object <<- object
    p9 <<- p
}

MultiPlotLineage <- function(cluster.object){
    Idents(cluster.object) <- cluster.object$class
    a <- FeaturePlot(cluster.object,c('Epcam','Cdh5','Col1a1','Ptprc'),ncol=2,label = F)
    b <- VlnPlot(cluster.object,c('Epcam','Cdh5','Col1a1','Ptprc'),cols=class_colors,ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())
    p <- plot_grid(a,b,ncol=2)
    p10 <<- p
}

BackUMAP <- function(sample.name,cluster.object){
    p <- UMAPPlot(cluster.object,group.by='Final3',label = F) + 
        labs(title = paste(sample.name,'final \u2014 archetypes'),caption = paste("nCells =",ncol(cluster.object))) + 
        theme(plot.title = element_text(size = 24,hjust=0),
            plot.caption = element_text(size=20,color='red'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.title.x = element_text(size = 16,color='black'),
            axis.text.x = element_text(size=12,color='black'))
    p11 <<- p
}

MultiPlotQC <- function(cluster.object){
    Idents(cluster.object) <- cluster.object$Final3
    a <- FeaturePlot(cluster.object,c('percent.mt','nCount_RNA','nFeature_RNA','Mki67'),ncol=2,label = F)
    b <- VlnPlot(cluster.object,c('percent.mt','nCount_RNA','nFeature_RNA','Mki67'),ncol=2,pt.size=0.1) & theme(axis.title.x = element_blank())
    p <- plot_grid(a,b,ncol=2)
    p12 <<- p
}

FilterbackViolin <- function(sample.name,object,min.counts,min.features,max.mt){
    meta.object <- object@meta.data
    nclusters <- length(levels(meta.object$Final3))
    a <- ggplot(meta.object, aes(x=Sample, y=nCount_RNA, fill=Sample)) + 
        geom_violin(fill='black',color='black',trim=F) +
        geom_jitter(aes(color = factor(Final3)), size = 0.1) +
        scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
        ggtitle('nCount_RNA') + scale_y_log10() + theme_classic() + NoLegend() +
        geom_hline(yintercept=min.counts, linetype='dashed', color='red', size=1) +
        labs(subtitle = paste("nUMI >",min.counts)) +
        theme(plot.title = element_text(size = 24),
            plot.subtitle = element_text(size=20,color='red'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.title.x = element_blank(),
            axis.text.x = element_blank())
    b <- ggplot(meta.object, aes(x=Sample, y=nFeature_RNA, fill=Sample)) + 
        geom_violin(fill='black',color='black',trim=F) +
        geom_jitter(aes(color = factor(Final3)), size = 0.1) +
        scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
        ggtitle('nFeature_RNA') + scale_y_log10() + theme_classic() + NoLegend() +
        geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
        labs(subtitle = paste("nGene >",min.features)) +
        theme(plot.title = element_text(size = 24),
            plot.subtitle = element_text(size=20,color='red'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.title.x = element_blank(),
            axis.text.x = element_blank())
    c <- ggplot(meta.object, aes(x=Sample, y=percent.mt, fill=Sample)) + 
        geom_violin(fill='black',color='black',trim=F) +
        geom_jitter(aes(color = factor(Final3)), size = 0.1) +
        scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
        ggtitle('percent.mt') + ylim(0,30) + theme_classic() + NoLegend() +
        geom_hline(yintercept=max.mt, linetype='dashed', color='red', size=1) +
        labs(subtitle = paste("pt.mt <",max.mt)) +
        theme(plot.title = element_text(size = 24),
            plot.subtitle = element_text(size=20,color='red'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.title.x = element_blank(),
            axis.text.x = element_blank())
    p <- plot_grid(a,b,c,ncol=3)
    p13 <<- p
}

FilterbackScatter <- function(sample.name,object,min.counts,min.features,max.mt){
    meta.object <- object@meta.data
    nclusters <- length(levels(meta.object$Final3))
    p <- ggplot() + geom_point(data=meta.object, aes(x=nCount_RNA, y=nFeature_RNA, color=factor(Final3)), size = 1) +
        scale_colour_manual(name="color", values=hue_pal()(nclusters)) + theme_classic() + NoLegend() +
        geom_vline(xintercept=min.counts, linetype='dashed', color='red', linewidth=1) +
        geom_hline(yintercept=min.features, linetype='dashed', color='red', linewidth=1) +
        labs(title = 'archetype \u2014 final') +
        labs(subtitle = paste("nUMI_RNA >",min.counts,"& nGene_RNA >",min.features,"& pt.mt <",max.mt)) +
        theme(plot.title = element_text(size = 24),
            plot.subtitle = element_text(size=20,color='red'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.title.x = element_text(size = 16,color='black'),
            axis.text.x = element_text(size=12,color='black'))
    p14 <<- p
}

##   BC1P3: 1000/600/5  20  v2  DONE
##   BC1P6: 1000/1000/10  15  v2  DONE
##   BCL5: 140/90/7.5  20  v2  DONE
##   BCEC2: 1600/500/25  30  v3  DONE
##   BEF1: 300/200/10  30  v2  DONE
##   BEF2: 400/250/10  30  v2  DONE
##   BEF3: 300/250/10  30  v3  DONE
##   BEF12: 1000/500/20  30  v3  DONE
##   BEF14: 1500/700/20  30  v3  DONE
##   BEF15: 1400/600/20  30  v3  DONE
##   BEFM1: 1000/400/25  40  v3  DONE
##   BEFM2: 1100/600/20  40  v3  DONE
##   BEFM4: 600/400/15  40  v3.1  DONE
##   BEFM5: 400/250/15  40  v3.1  DONE
##   BEFM6: 800/500/18  40  v3.1  DONE
##   FB13: 1000/800/7.5*  20  v2  DONE
##   FB14: 1000/1000/5*  20  v2  DONE
##   MacAlv: 500/400/15*  20  v3.1  DONE
##   RLMVEC: 2000/1500/10*  15  v3.1  DONE

# Setup and select cutoffs for each object
sample.name <- 'BAL'
file.name <- 'MacAlv'
#file.name <- sample.name
min.counts <- 500
min.features <- 400
max.mt <- 15
pcs <- 20
chemistry <- 'v3.1'
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
setwd(paste0("./",file.name))
paste0(sample.name,tag,'_suppfig_2.png')

Idents(object) <- object$Sample
sub <- subset(object,idents = sample.name)
sub$Final3 <- droplevels(sub$Final3)

title <- ggdraw() + draw_label(paste0('Sample: ',sample.name,' (',chemistry,' chemistry)'),
    size=30,fontface = 'bold',x = 0,hjust = 0)
title.1 <- ggdraw() + draw_label(paste(sample.name,'initital cleaning'),
    size=24,fontface = 'bold',x = 0,hjust = -0.1)
Barcode_rank_plot(sample.name=sample.name,file.name=file.name)
QCViolinPlotter(sample.name=sample.name,file.name=file.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt)
TriplePercentMtPlot(sample.name=sample.name)
TriplePercentMtPlotFilter(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt)
TriplePercentSplPlot(sample.name=sample.name)
TriplePercentSplPlotFilter(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt)
title.2 <- ggdraw() + draw_label(paste(sample.name,'final object'),
    size=24,fontface = 'bold',x = 0,hjust = -0.1)
SeuratCluster(sample.name=sample.name,object=sub,pcs=pcs)
MultiPlotLineage(cluster.object=cluster.object)
BackUMAP(sample.name=sample.name,cluster.object=cluster.object)
MultiPlotQC(cluster.object=cluster.object)
FilterbackViolin(sample.name=sample.name,object=sub,min.counts=min.counts,min.features=min.features,max.mt=max.mt)
FilterbackScatter(sample.name=sample.name,object=sub,min.counts=min.counts,min.features=min.features,max.mt=max.mt)

top_row <- plot_grid(p1,p2,p3,p4,ncol=4,rel_widths=c(2,1,1,1),labels=c('A','B','',''),label_size=40)
mid_row1 <- plot_grid(p5,p6,p7,p8,ncol=4,labels=c('C','','D',''),label_size=40)
mid_row2 <- plot_grid(p9,p10,ncol=2,rel_widths=c(1,2),labels=c('E','F'),label_size=40)
mid_row3 <- plot_grid(NULL,p12,ncol=2,rel_widths=c(1,2),labels=c('','G'),label_size=40)
last_row <- plot_grid(p11,p13,p14,ncol=3,labels=c('H','I','J'),label_size=40)
png(paste0(sample.name,tag,'_suppfig_2.png'), width = 2000,height=2750)
plot_grid(title,title.1,top_row,mid_row1,title.2,mid_row2,mid_row3,last_row,ncol=1,
    rel_heights = c(0.1,0.05,1,1,0.05,1,1,1))
dev.off()




#### FIGURES S21-22 ####
### Epithelium Supplement
## Load samples, colors, set-up
# Load up-to-date epithelium object
load('./Epi/epi.seurat.objs.int.2023-07-27.Robj')
epi <- object ; rm(object)
epi$Dataset4 <- factor(epi$Dataset4,levels=c('Start','Mono','Co','Tri','Quad'))
load('./all.color_palettes_2.R')
sample_colors <- color_pals$sample_colors
dataset_colors <- color_pals$dataset_colors
class_colors <- color_pals$class_colors
epi_colors <- color_pals[['epi']]

## Plots
# UMAP plot of integrated archetypes - colored by Final1
if (ncol(object) < 10000){pt.size=1.5} else if (ncol(object) >= 10000){pt.size=NULL}
res <- 400
png(paste(dir_name,'_umap_1.png',sep=''), width=600*res/72, height=500*res/72,res=res)
UMAPPlot(object,group.by = 'Final1',cols = type_colors,pt.size=pt.size,label = T,repel=F,label.size = 8) + NoLegend() +
    labs(title = 'Engineered Epithelial Archetypes',caption = paste("nCells =",ncol(object))) + 
    theme(plot.title = element_text(size = 30,hjust = 0,face = "plain"),
        plot.caption = element_text(size = 18,color='red'),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'))
dev.off()

# UMAP plot of integrated archetypes - colored by Sample
if (ncol(epi) < 10000){pt.size=1.5} else if (ncol(epi) >= 10000){pt.size=NULL}
p1 <- UMAPPlot(epi,group.by = 'Sample',cols = sample_colors[levels(epi$Sample)],shuffle=T,pt.size=pt.size,label = F) + 
    labs(title = 'Integration',caption = paste("nCells =",ncol(epi))) + 
    theme(plot.title = element_text(size = 30,hjust = 0,face = "plain"),
        plot.caption = element_text(size = 14,color='red'),
        axis.text.x=element_text(size=14,color='black'),
        axis.text.y = element_text(size = 14,color='black'),
        axis.title.x = element_text(size = 16,color='black'),
        axis.title.y = element_text(size = 16,color='black'))

# Heatmap of select defining archetype markers (vertical) - AVERAGE
Idents(epi) <- epi$Final1
epi$Final1 <- factor(epi$Final1,levels=c('Cycling','Basal-like','Transitional',
        'BASC','ATI-like'))
epi$Sample <- factor(epi$Sample,levels=c("BC1P3","BC1P6","BCL5","BCEC2","BEF1",
        "BEF2","BEF3","BEF12","BEF14","BEF15","BEFM1","BEFM2","BEFM4","BEFM5","BEFM6"))
epi <- ScaleData(epi,features = rownames(epi))
Idents(epi) <- epi$Final1
genes <- c('Epcam','Wnt4','Wnt7b','Fabp5','Fzd6',  # Pan-Epithelium
        'Top2a','Mki67',  # Cycling
        'Krt5','Krt14','Tp63','Igfbp2','Anxa8','Fgfbp1','Areg',  # Basal-like
        'Krt8','Krt18','Ly6d','Il33','Lcn2',  # Transitional
        'Sftpc','Scgb1a1','Scgb3a2','Sftpb','Defb4','Lyz2','Sftpa1','Cd74','Napsa','Defb3',
        'Sftpd','Slc34a2','Pla2g1b','Lamp3','Nkx2-1',  # BASC
        'Hopx','Aqp5','Scel','Wnt11','Mal','Elf3','Clca4l','Plat','Gpha2','Krt7','Krt15','Cldn23',
        'Muc20','Pllp','Wwtr1')  # ATI-like

gene_split <- c(rep('Pan-Epi',5),rep('Cycling',2),rep('Basal-like',7),rep('Transitional',5),
    rep('BASC',15),rep('ATI-like',15))
genes.use <- intersect(rownames(GetAssayData(epi,slot = 'scale.data')),genes)
epi.avg <- AverageExpression(epi,assays = 'RNA',slot='scale.data',features=genes.use)
epi.avg <- t(scale(t(epi.avg$RNA)))
unityNormalize <- function(x){(x-min(x))/(max(x)-min(x))}
epi.avg <- t(apply(as.matrix(epi.avg), 1, unityNormalize))
epi.avg <- epi.avg[,levels(epi$Final1)]
epi.avg <- epi.avg[genes,]
dim(epi.avg)

metadata <- epi@meta.data
Pseudotime <- metadata %>%
    group_by(Final1) %>%
    summarize(average_pseudotime = mean(Lineage_avg, na.rm = TRUE))
Pseudotime <- setNames(Pseudotime$average_pseudotime, Pseudotime$Final1)[levels(epi$Final1)]
colors.spectral <- colorRamp2(breaks = c(seq(40,70,length.out=60)),colorRampPalette(brewer.pal(11,'Spectral')[-6])(60), space = "RGB")
colors.inferno <- colorRamp2(breaks = c(seq(min(epi.avg),max(epi.avg),length.out=60)), inferno(n=60), space = "RGB")

ha = HeatmapAnnotation(CellType = anno_block(gp = gpar(fill=epi_colors[levels(epi$Final1)]),show_name=T),
        Pseudotime = anno_block(gp = gpar(fill = colors.spectral(Pseudotime)),show_name=T),
        annotation_name_gp = gpar(fontsize = 14)
)

lgd1 = Legend(col_fun = colors.spectral,
    title = "Pseudotime",
    title_gp = gpar(fontsize = 14),
    labels = c("early","late"), at = c(min(Pseudotime),max(Pseudotime)),
    labels_gp = gpar(fontsize = 12),
    grid_width = unit(1, "cm"),
    legend_height = unit(3, "cm"))
lgd2 = Legend(col_fun = colors.inferno,
    title = "Gene\nexpression", 
    title_gp = gpar(fontsize = 14),
    labels = c("Low","High"), at = c(0,1),
    labels_gp = gpar(fontsize = 12),
    grid_width = unit(1, "cm"),
    legend_height = unit(3, "cm"))
pd = packLegend(lgd1,lgd2,column_gap = unit(1, "cm"),max_height = unit(30, "cm"))

epi.heatmap <- Heatmap(as.matrix(epi.avg),
        col = colors.inferno,
        column_title = levels(epi$Final1),
        column_split = factor(levels(epi$Final1),levels=levels(epi$Final1)),
        top_annotation = ha,
        row_split = factor(gene_split,levels=unique(gene_split)),
        row_title_gp = gpar(fontsize = 12),
        column_title_gp = gpar(fontsize = 14),
        row_names_gp = gpar(fontsize = 12,fontface='italic'),
        column_title_rot = 0,
        cluster_rows=F,
        cluster_columns=F,
        cluster_row_slices=F,
        show_row_names=T,
        show_column_names = F,
        show_heatmap_legend = F)
p2 = grid.grabExpr(draw(epi.heatmap,annotation_legend_list = pd,padding = unit(c(2, 2, 2, 2), "mm")))
p2 <- plot_grid(NULL,p2,ncol=2,rel_widths=c(1,21))

# Proportion bar plot by Dataset4
data0 <- as.data.frame(table(epi$Dataset4))
data2 <- prop.table(table(epi$Dataset4,epi$Final1),1)*100
data2 <- as.data.frame(data2)
p4 <- ggplot(data2, aes(x = Var1, y = Freq, fill = Var2)) + 
    geom_col(position = "fill") + scale_fill_manual(values = epi_colors[levels(epi$Final1)],name="") +
    scale_y_continuous(labels = scales::percent_format()) +
    xlab("") + ylab("% of dataset") +
    labs(title = 'Composition by Archetype') + 
    theme(plot.title = element_text(size = 30),
            legend.text = element_text(size = 14),
            axis.text.x=element_text(size=16,color='black'),  # angle = 45,vjust=0.5,
            axis.text.y = element_text(size = 14,color='black'),
            axis.title.x = element_text(size = 16,color='black'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.line = element_line(size = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
    geom_text(aes(x = Var1,y = 1,label = Freq,fill = NULL),color='red',
            data=data0,position = position_dodge(width = 0.5),vjust=-0.3,size=5)+
    scale_color_manual(values = "red")
p4 <- plot_grid(NULL,p4,ncol=2,rel_widths=c(1,21))

ps13 <- plot_grid(p1,NULL,nrow=2)
topleft <- plot_grid(ps13,p2,NULL,p4,ncol=2,rel_widths=c(2,3),rel_heights=c(2,1))

# High resolution png
res <- 400
png(paste('epi_figure_pheno_1.png'),width=1000*res/72,height=1100*res/72,res=res)
plot_grid(topleft)
dev.off()

# Gene Ontology bar plot of select epithelial archetypes
golist <- list()
golist$"Basal-like" <- c('0048732','0031099','0060429')
golist$"ATI-like" <- c('0033993','0030855','0030155','0022407','0009612','0070830','0043297')
golist$BASC <- c('0033993','0030154','0030324','0006629','0043129')
gotypes <- levels(epi$Final1)[levels(epi$Final1) != 'Cycling']
gotypes <- gotypes[gotypes != 'Transitional']
allgo <- c()
for (i in 1:length(gotypes)) {
        type <- gotypes[i]
        message(type)
        go_cat <- golist[[type]]
        gosheet <- as.data.frame(read_excel("./Epi/epi.seurat.objs.int.marks.2023-07-11.xlsm",sheet = type))
        go <- gosheet[grep(paste(go_cat,collapse="|"),gosheet$'GO biological process complete'), ]
        go <- go[order(go$'upload_1 (raw P-value)'),]
        go$cluster <- type
        go[3] <- NULL
        allgo <- rbind(allgo,go)
}
allgo$p_val <- allgo$'upload_1 (raw P-value)'
allgo$go_term <- allgo$'GO biological process complete'
allgo$go_term2 <- paste(allgo$'GO biological process complete',allgo$cluster,sep='_')
allgo$go_term2 <- factor(allgo$go_term2,levels=allgo$go_term2)
allgo$cluster <- factor(allgo$cluster,levels=gotypes)
allgo$color <- epi_colors[as.character(allgo$cluster)]
allgo$revrank <- rev(as.numeric(as.factor(seq.int(nrow(allgo)))))
ref_avg <- aggregate(revrank ~ cluster, allgo, mean)
allgo$breaks <- ref_avg[match(allgo$cluster, ref_avg$cluster), "revrank"]
allgo$breaks <- factor(allgo$breaks,levels=unique(allgo$breaks))
px <- ggplot(allgo,aes(x=factor(revrank),y=p_val,fill=go_term2,label=go_term)) +
    geom_col(position = position_dodge2(preserve='single',reverse=T),show.legend = F) + 
    scale_fill_manual(values = allgo$color) + coord_flip() +
    ylab('raw P-value (log scale)') + 
    scale_x_discrete(breaks = factor(levels(allgo$breaks)),
        labels = unique(allgo$cluster)) +
    scale_y_continuous(trans=compose_trans("log10", "reverse"),
        expand = expansion(add = c(0,14)),  # 18
        breaks = trans_breaks('log10', function(x) 10^x),
        labels = trans_format('log10', math_format(10^.x))) +
    geom_text(position = position_dodge2(0.9,preserve='single',reverse=T),hjust=-0.02,size=5) +
    geom_hline(yintercept=0.05,linetype="dashed",color = "gray") +
    theme(axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 16,color='black',angle=90,hjust = 0.5),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_blank(),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
px <- plot_grid(NULL,px,ncol=2,rel_widths=c(1,11))

# Line plot of genes of interest
# Function for line plot of two datasets, split by ident
lineplot <- function(dir_name,object.split,goi,first,last){
    type_colors <- color_pals[[dir_name]]

    myplots <- vector('list', length(goi))
    for (i in 1:length(goi)) {
        gene <- goi[i]
        message(gene)
        object.avg <- c() ; df <- data.frame()
        for (j in 1:length(object.split)){
            j <- j
            gene <- gene
            type <- names(object.split)[j]
            message(type)
            object.avg <- AverageExpression(object.split[[j]],features = gene,assays = 'RNA')
            object.avg[[j]] <- AverageExpression(object.split[[j]],features = gene,assays = 'RNA')
            object.avg[[j]] <- data.frame(t(object.avg[[j]]$RNA))
            colnames(object.avg[[j]]) <- 'Expression'
            object.avg[[j]]$Dataset4 <- rownames(object.avg[[j]])
            object.avg[[j]]$CellType <- names(object.split)[j]
            df <- rbind(df,object.avg[[j]])
        }
        df2 <- df[df$Dataset %in% c(first,last),]
        myplots[[i]] <- local({
            i <- i
            first <- first
            last <- last
            p <- ggplot(df2,aes(x = factor(Dataset4,levels = c(first,last)),
                y = Expression,group = CellType,color=CellType)) + 
                scale_color_manual(values = type_colors) +
                geom_line(size=1.2) + geom_point(size=3) + ggtitle(gene) +
                xlab('') + ylab('Gene Expression')+ NoLegend() +
                theme(plot.title = element_text(size = 24, hjust = 0.5,face='italic'),
                    axis.title.y = element_text(size = 14,color='black'),
                    axis.text.x = element_text(size=12,color='black'),
                    axis.text.y = element_text(size = 12,color='black'),
                    axis.line = element_line(size = 1),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank())
        })
    }
    s <- ggplot(df2,aes(x = factor(Dataset4,levels = c(first,last)),
        y = Expression,group = CellType,color=CellType)) + 
        scale_color_manual(values = type_colors) +
        geom_line(size=1.2) + geom_point(size=3) + 
        theme(legend.title = element_text(size = 14,color = 'black',face = 'bold'),
            legend.text = element_text(size = 12,color = 'black'),
            legend.key=element_blank(),
            legend.key.size = unit(2, "lines"),
            legend.direction="vertical")
    lgd <- cowplot::get_legend(s)
    ncols <- 3
    p1 <- plot_grid(plotlist=myplots,ncol=ncols)
    ht <- 233*ceiling(length(goi)/ncols)
    wd <- 300*ncols
    p1 <<- p1
    lgd <<- lgd
    ncols <<- ncols
    ht <<- ht
    wd <<- wd
}

goi_line <- c('Scel','Wnt11','Wnt7b','Tm4sf1','Celsr1','Fgfr2','Foxa2','Foxa1','Hes1','Lcn2','Il33','Nkx2-1')

Idents(epi) <- epi$Dataset4
epi.split <- SplitObject(epi,split.by = 'Final1')

goi <- goi_line2
end <- 'Supp'

lineplot(dir_name='epi',object.split=epi.split,goi=goi,first='Tri',last='Quad')

res <- 400
png(paste0('epi_line_',end,'2.png'),width=wd*res/72,height=ht*res/72,res=res)
plot_grid(p1,lgd,ncol=2,rel_widths=c((ncols+1),1))
dev.off()

# RNA Velocity plots - see focus.R script

# Function for Featureplots of goi
featureplot <- function(object,goi){
    genes.use <- intersect(goi,rownames(GetAssayData(object,assay = "RNA",slot = "data")))
    myplots <- vector('list', length(genes.use))
    for (i in 1:length(genes.use)) {
        gene <- genes.use[i]
        message(gene)
        p <- FeaturePlot(object,features=gene,order=T,label = F) &
            labs(title = c(gene)) & NoAxes()
        myplots[[i]] <- p
    }
    ncols <- 5
    p1 <- plot_grid(plotlist=myplots,ncol=ncols)
    ht1 <- 233*ceiling(length(genes.use)/ncols)
    wd <- 300*ncols
    p1 <<- p1
    ht1 <<- ht1
    wd <<- wd
}

# Function for Violin plots of goi
violinplot <- function(dir_name,object,goi){
    genes.use <- intersect(goi,rownames(GetAssayData(object,assay = "RNA",slot = "data")))
    myplots <- vector('list', length(genes.use))
    for (i in 1:length(genes.use)) {
        gene <- genes.use[i]
        message(gene)
        p <- VlnPlot(object,gene,group.by='Dataset4',assay='RNA',slot='data',cols=dataset_colors[levels(object$Dataset4)],pt.size=0) + 
            labs(title = gene) + NoLegend() + 
            theme(plot.title = element_text(size = 30,hjust = 0,face = "italic"),
            axis.text.x=element_text(size=14,color='black',angle=0,hjust = 0.5),
            axis.text.y = element_text(size = 14,color='black'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.title.x = element_blank()) +
            geom_vline(xintercept=1.5,linetype="longdash",color = "black")
        myplots[[i]] <- p
    }
    ncols <- 3
    p2 <- plot_grid(plotlist=myplots,ncol=ncols)
    ht2 <- 233*ceiling(length(genes.use)/ncols)
    wd <- 300*ncols
    p2 <<- p2
    ht2 <<- ht2
    wd <<- wd    
}

# Supplemental FeaturePlots & Violin Plots - Epithelium
Idents(epi) <- epi$Dataset4
epi.split <- SplitObject(epi,split.by = 'Final1')
epi.sub <- subset(epi,idents=c('Tri','Quad'))
epi.sub$Dataset4 <- droplevels(epi.sub$Dataset4)

goi <- goi_Final
end <- 'Supp'

featureplot(object=epi,goi=goi)
violinplot(dir_name='epi',object=epi.sub,goi=goi)

res <- 400
png(paste0('epi_fpvln_',end,'.png'),width=wd*res/72,height=(ht1+ht2)*res/72,res=res)
plot_grid(p1,p2,ncol=1)
dev.off()




#### FIGURES S23-24 ####
### Endothelium Supplement
## Load samples, colors, set-up
# Load up-to-date endothelium object
load('./Endo/endo.seurat.objs.int.2023-07-27.Robj')
endo <- object ; rm(object)
endo$Dataset4 <- factor(endo$Dataset4,levels=c('Start','Co','Tri','Quad'))
Idents(endo) <- endo$Final1
endo$Final1 <- factor(endo$Final1,levels=c('Progenitor','Lymphatic',
    'Cycling_Lymphatic','Cycling_Microvascular','Microvascular'))
endo <- RenameIdents(endo,
        "Cycling_Lymphatic"='Cyc_Lymph',
        "Cycling_Microvascular"='Cyc_MV')
endo$Final3 <- Idents(endo)
endo$Final3 <- factor(endo$Final3,levels=c('Progenitor','Lymphatic',
    'Cyc_Lymph','Cyc_MV','Microvascular'))
endo <- endo
endo_colors <- color_pals[['endo']]
endo_colors2 <- endo_colors
names(endo_colors2) <- c("Microvascular","Cyc_MV","Lymphatic","Cyc_Lymph","Progenitor")

## Plots
# UMAP plot of integrated archetypes - colored by Final3
if (ncol(object) < 10000){pt.size=1.5} else if (ncol(object) >= 10000){pt.size=NULL}
res <- 400
png(paste(dir_name,'_umap_1.png',sep=''), width=600*res/72, height=500*res/72,res=res)
UMAPPlot(object,group.by = 'Final3',cols = type_colors2,pt.size=pt.size,label = T,repel=F,label.size = 8) + NoLegend() +
    labs(title = 'Engineered Endothelium',caption = paste("nCells =",ncol(object))) + 
    theme(plot.title = element_text(size = 30,hjust = 0,face = "plain"),
        plot.caption = element_text(size = 18,color='red'),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'))
dev.off()

# UMAP plot of integrated archetypes - colored by Sample
if (ncol(endo) < 10000){pt.size=1.5} else if (ncol(endo) >= 10000){pt.size=NULL}
p5 <- UMAPPlot(endo,group.by = 'Sample',cols = sample_colors[levels(endo$Sample)],shuffle=T,pt.size=pt.size,label = F) + 
    labs(title = 'Integration',caption = paste("nCells =",ncol(endo))) + 
    theme(plot.title = element_text(size = 30,hjust = 0,face = "plain"),
        plot.caption = element_text(size = 14,color='red'),
        axis.text.x=element_text(size=14,color='black'),
        axis.text.y = element_text(size = 14,color='black'),
        axis.title.x = element_text(size = 16,color='black'),
        axis.title.y = element_text(size = 16,color='black'))

# Heatmap of select defining archetype markers (vertical) - AVERAGE
Idents(endo) <- endo$Final3
endo <- ScaleData(endo,features = rownames(endo))
endo$Final3 <- factor(endo$Final3,levels=c('Progenitor','Lymphatic',
    'Cyc_Lymph','Cyc_MV','Microvascular'))
genes <- c('Cdh5','Epas1','Bmpr2','Fzd4','Erg',  # Pan-Endothelial 
    'Trib3','Casp4','Gadd45a','Atf3','Vegfa','Cxcl2','Cxcl3',  # Progenitor 
    'Mmrn1','Reln','Msr1','Nr2f2','Slc14a1','Vcam1','Ccn3','Sema3a','Bmp4',  # Lymphatic 
    'Top2a','Mki67',  # Cycling
    'Aqp1','Kdr','Sox17','Vegfc','Vwf','Sox18','Plxnd1','Tspan18','Kit','Kitlg','Cd34','Ephb1',
    'Emcn','Gpihbp1','Ednrb','Lrp4','Aplnr','Cemip2','Efnb2','Fbln5','Selp','Hdac9','Ackr3')  # Microvascular
gene_split <- c(rep('Pan-Endo',5),rep('Progenitor',7),rep('Lymphatic',9),
    rep('Cycling',2),rep('Microvascular',23))
genes.use <- intersect(rownames(GetAssayData(endo,slot = 'scale.data')),genes)
endo.avg <- AverageExpression(endo,assays = 'RNA',slot='scale.data',features=genes.use)
endo.avg <- t(scale(t(endo.avg$RNA)))
unityNormalize <- function(x){(x-min(x))/(max(x)-min(x))}
endo.avg <- t(apply(as.matrix(endo.avg), 1, unityNormalize))
endo.avg <- endo.avg[,levels(endo$Final3)]
endo.avg <- endo.avg[genes,]
dim(endo.avg)

metadata <- endo@meta.data
Pseudotime <- metadata %>%
    group_by(Final3) %>%
    summarize(average_pseudotime = mean(Lineage1_ps, na.rm = TRUE))
Pseudotime <- setNames(Pseudotime$average_pseudotime, Pseudotime$Final3)[levels(endo$Final3)]
colors.spectral <- colorRamp2(breaks = c(seq(20,80,length.out=60)),colorRampPalette(brewer.pal(11,'Spectral')[-6])(60), space = "RGB")
colors.inferno <- colorRamp2(breaks = c(seq(min(endo.avg),max(endo.avg),length.out=60)), inferno(n=60), space = "RGB")

ha = HeatmapAnnotation(CellType = anno_block(gp = gpar(fill=endo_colors[levels(endo$Final1)]),show_name=T),
        Pseudotime = anno_block(gp = gpar(fill = colors.spectral(Pseudotime)),show_name=T),
        annotation_name_gp = gpar(fontsize = 14)
)

lgd1 = Legend(col_fun = colors.spectral,
    title = "Pseudotime",
    title_gp = gpar(fontsize = 14),
    labels = c("early","late"), at = c(min(Pseudotime),max(Pseudotime)),
    labels_gp = gpar(fontsize = 12),
    grid_width = unit(1, "cm"),
    legend_height = unit(3, "cm"))
lgd2 = Legend(col_fun = colors.inferno,
    title = "Gene\nexpression", 
    title_gp = gpar(fontsize = 14),
    labels = c("Low","High"), at = c(0,1),
    labels_gp = gpar(fontsize = 12),
    grid_width = unit(1, "cm"),
    legend_height = unit(3, "cm"))
pd = packLegend(lgd1,lgd2,column_gap = unit(1, "cm"),max_height = unit(30, "cm"))

endo.heatmap <- Heatmap(as.matrix(endo.avg),
        col = colors.inferno,
        column_title = levels(endo$Final3),
        column_split = factor(levels(endo$Final3),levels=levels(endo$Final3)),
        top_annotation = ha,
        row_split = factor(gene_split,levels=unique(gene_split)),
        row_title_gp = gpar(fontsize = 14),
        column_title_gp = gpar(fontsize = 14),
        row_names_gp = gpar(fontsize = 12,fontface='italic'),
        cluster_rows=F,
        cluster_columns=F,
        cluster_row_slices=F,
        show_row_names=T,
        show_column_names = F,
        show_heatmap_legend = F)
p6 = grid.grabExpr(draw(endo.heatmap,annotation_legend_list = pd,padding = unit(c(2, 2, 2, 2), "mm")))
p6 <- plot_grid(NULL,p6,ncol=2,rel_widths=c(1,21))

# Proportion bar plot by Dataset4
endo$Final3 <- factor(endo$Final3,levels=c('Microvascular','Cyc_MV','Lymphatic',
    'Cyc_Lymph','Progenitor'))
data0 <- as.data.frame(table(endo$Dataset4))
data2 <- prop.table(table(endo$Dataset4,endo$Final3),1)*100
data2 <- as.data.frame(data2)
p8 <- ggplot(data2, aes(x = Var1, y = Freq, fill = Var2)) + 
    geom_col(position = "fill") + scale_fill_manual(values = endo_colors2,name="") +
    scale_y_continuous(labels = scales::percent_format()) +
    xlab("") + ylab("% of dataset") +
    labs(title = 'Composition by Archetype') + 
    theme(plot.title = element_text(size = 30),
            legend.text = element_text(size = 14),
            axis.text.x=element_text(size=16,color='black'),  # angle = 45,vjust=0.5,
            axis.text.y = element_text(size = 14,color='black'),
            axis.title.x = element_text(size = 16,color='black'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.line = element_line(size = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
    geom_text(aes(x = Var1,y = 1,label = Freq,fill = NULL),color='red',
            data=data0,position = position_dodge(width = 0.5),vjust=-0.3,size=5)+
    scale_color_manual(values = "red")
p8 <- plot_grid(NULL,p8,ncol=2,rel_widths=c(1,21))

ps57 <- plot_grid(p5,p7,nrow=2)
topright <- plot_grid(ps57,p6,NULL,p8,ncol=2,rel_widths=c(2,3),rel_heights=c(2,1))

res <- 400
png(paste('endo_figure_pheno_1.png'),width=1000*res/72,height=1100*res/72,res=res)
plot_grid(topright)
dev.off()

# Gene Ontology bar plot of select endothelial archetypes
golist <- list()
golist$Interesting <- c('0001525','0048514','0003158','0098609','0002040','0048863','0061028')
golist$Lymphatic <- c('0033993','0030154','0002521','0048534','0050878')
golist$Homeostatic <- c('0045765','0031099','0001944')
gotypes <- levels(endo$putative_labels)[levels(endo$putative_labels) != 'Cycling2']
gotypes <- gotypes[gotypes != 'Cycling3']
allgo <- c()
for (i in 1:length(gotypes)) {
        type <- gotypes[i]
        message(type)
        go_cat <- golist[[type]]
        gosheet <- as.data.frame(read_excel("./Endo/endo.seurat.objs.int.marks.2022-10-28.xlsm",sheet = type))
        go <- gosheet[grep(paste(go_cat,collapse="|"),gosheet$'GO biological process complete'), ]
        go <- go[order(go$'upload_1 (raw P-value)'),]
        go$cluster <- type
        go[3] <- NULL
        allgo <- rbind(allgo,go)
}
allgo$p_val <- allgo$'upload_1 (raw P-value)'
allgo$go_term <- allgo$'GO biological process complete'
allgo$go_term2 <- paste(allgo$'GO biological process complete',allgo$cluster,sep='_')
allgo$go_term2 <- factor(allgo$go_term2,levels=allgo$go_term2)
allgo$cluster <- factor(allgo$cluster,levels=gotypes)
allgo$cluster2[allgo$cluster == 'Interesting'] = 'Microvascular'
allgo$cluster2[allgo$cluster == 'Lymphatic'] = 'Lymphatic'
allgo$cluster2[allgo$cluster == 'Homeostatic'] = 'Progenitor'
allgo$color <- endo_colors[allgo$cluster2]
allgo$revrank <- rev(as.numeric(as.factor(seq.int(nrow(allgo)))))
ref_avg <- aggregate(revrank ~ cluster, allgo, mean)
allgo$breaks <- ref_avg[match(allgo$cluster, ref_avg$cluster), "revrank"]
allgo$breaks <- factor(allgo$breaks,levels=unique(allgo$breaks))
px <- ggplot(allgo,aes(x=factor(revrank),y=p_val,fill=go_term2,label=go_term)) +
    geom_col(position = position_dodge2(preserve='single',reverse=T),show.legend = F) + 
    scale_fill_manual(values = allgo$color) + coord_flip() +
    ylab('raw P-value (log scale)') + 
    scale_x_discrete(breaks = factor(levels(allgo$breaks)),
        labels = unique(allgo$cluster2)) +
    scale_y_continuous(trans=compose_trans("log10", "reverse"),
        expand = expansion(add = c(0,60)),
        breaks = trans_breaks('log10', function(x) 10^x),
        labels = trans_format('log10', math_format(10^.x))) +
    geom_text(position = position_dodge2(0.9,preserve='single',reverse=T),hjust=-0.02,size=5) +
    geom_hline(yintercept=0.05,linetype="dashed",color = "gray") +
    theme(axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 16,color='black',angle=90,hjust = 0.5),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_blank(),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
px <- plot_grid(NULL,px,ncol=2,rel_widths=c(1,11))

# Line plot of genes of interest
# (see lineplot function above)
goi_line <- c('Kdr','Lrp4','Aplnr','Gpihbp1','Cemip2','Tek','Efnb2','Fbln5','Trib3','Vegfa','Bnip3','Kit')

Idents(endo) <- endo$Dataset4
endo.split <- SplitObject(endo,split.by = 'Final1')

goi <- goi_line
end <- 'Supp'

lineplot(dir_name='endo',object.split=endo.split,goi=goi,first='Tri',last='Quad')

res <- 400
png(paste0('endo_line_',end,'.png'),width=wd*res/72,height=ht*res/72,res=res)
plot_grid(p1,lgd,ncol=2,rel_widths=c((ncols+1),1))
dev.off()

# RNA Velocity plots - see focus.R script

# Supplemental FeaturePlots & Violin Plots - Endothelium
# (see featureplot & violinplot functions above)
Idents(endo) <- endo$Dataset4
endo.split <- SplitObject(endo,split.by = 'Final1')
endo.sub <- subset(endo,idents=c('Tri','Quad'))
endo.sub$Dataset4 <- droplevels(endo.sub$Dataset4)

goi <- goi_Final
end <- 'Supp'

featureplot(object=endo,goi=goi)
violinplot(dir_name='endo',object=endo.sub,goi=goi)

res <- 400
png(paste0('endo_fpvln_',end,'.png'),width=wd*res/72,height=(ht1+ht2)*res/72,res=res)
plot_grid(p1,p2,ncol=1)
dev.off()




#### FIGURES S25-26 ####
### Mesenchyme Supplement
## Load samples, colors, set-up
# Load up-to-date mesenchyme object
load('./Mes/mes.seurat.objs.int.2023-07-27.Robj')
mes <- object ; rm(object)
mes$Dataset4 <- factor(mes$Dataset4,levels=c('Start','Tri','Quad'))
mes_colors <- color_pals[['mes']]


## Plots
# UMAP plot of integrated archetypes - colored by Final1
if (ncol(mes) < 10000){pt.size=1.5} else if (ncol(mes) >= 10000){pt.size=NULL}
res <- 400
png(paste(dir_name,'_umap_1.png',sep=''), width=600*res/72, height=500*res/72,res=res)
UMAPPlot(mes,group.by = 'Final1',cols = type_colors,pt.size=pt.size,label = T,repel=F,label.size = 8) + NoLegend() +
    labs(title = 'Engineered Mesenchymal Archetypes',caption = paste("nCells =",ncol(mes))) + 
    theme(plot.title = element_text(size = 30,hjust = 0,face = "plain"),
        plot.caption = element_text(size = 18,color='red'),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'))
dev.off()

# UMAP plot of integrated archetypes - colored by Sample
if (ncol(mes) < 10000){pt.size=1.5} else if (ncol(mes) >= 10000){pt.size=NULL}
p9 <- UMAPPlot(mes,group.by = 'Sample',cols = sample_colors[levels(mes$Sample)],shuffle=T,pt.size=pt.size,label = F) + 
    labs(title = 'Integration',caption = paste("nCells =",ncol(mes))) + 
    theme(plot.title = element_text(size = 30,hjust = 0,face = "plain"),
        plot.caption = element_text(size = 14,color='red'),
        axis.text.x=element_text(size=14,color='black'),
        axis.text.y = element_text(size = 14,color='black'),
        axis.title.x = element_text(size = 16,color='black'),
        axis.title.y = element_text(size = 16,color='black'))

# Heatmap of select defining archetype markers (vertical) - AVERAGE
Idents(mes) <- mes$Final1
mes <- ScaleData(mes,features = rownames(mes))
mes$Final1 <- factor(mes$Final1,levels=c('Cycling','Alveolar','Remodeling',
        'Pericytes','Aberrant'))
genes <- c('Col1a1','Itga8','Bgn','Mgp','Pdgfra','Fzd1','Vim',  # Pan-Mesenchyme
        'Top2a','Mki67',  # Cycling
        'Meox2','Wnt2','Angpt1','Rspo3','Cldn1','Angptl4','Hp','Porcn',  # Alveolar
        'Dcn','Mmp11','Twist2','Osr2','Bmp2','Islr','Mfap4','Mmp3','Eln','Fn1','Postn',  # Remodeling 
        'Pdgfrb','Gucy1a1','Gucy1b1','Acta2','Notch3',
        'Cspg4','Ebf1','Mcam','Higd1b','Rgs5','Fgfr4',  # Pericytes 
        'Hsp90b1','Hspa5','Trib3','Ero1b','Cxcl14','Creg1')  # Aberrant
gene_split <- c(rep('Pan-Mes',7),rep('Cycling',2),rep('Alveolar',8),
    rep('Remodeling',11),rep('Pericytes',11),rep('Aberrant',6))
genes.use <- intersect(rownames(GetAssayData(mes,slot = 'scale.data')),genes)
mes.avg <- AverageExpression(mes,assays = 'RNA',slot='scale.data',features=genes.use)
mes.avg <- t(scale(t(mes.avg$RNA)))
unityNormalize <- function(x){(x-min(x))/(max(x)-min(x))}
mes.avg <- t(apply(as.matrix(mes.avg), 1, unityNormalize))
mes.avg <- mes.avg[,levels(mes$Final1)]
mes.avg <- mes.avg[genes,]
dim(mes.avg)

metadata <- mes@meta.data
Pseudotime <- metadata %>%
    group_by(Final1) %>%
    summarize(average_pseudotime = mean(Lineage_avg, na.rm = TRUE))
Pseudotime <- setNames(Pseudotime$average_pseudotime, Pseudotime$Final1)[levels(mes$Final1)]
colors.spectral <- colorRamp2(breaks = c(seq(40,90,length.out=60)),colorRampPalette(brewer.pal(11,'Spectral')[-6])(60), space = "RGB")
colors.inferno <- colorRamp2(breaks = c(seq(min(mes.avg),max(mes.avg),length.out=60)), inferno(n=60), space = "RGB")

ha = HeatmapAnnotation(CellType = anno_block(gp = gpar(fill=mes_colors[levels(mes$Final1)]),show_name=T),
        Pseudotime = anno_block(gp = gpar(fill = colors.spectral(Pseudotime)),show_name=T),
        annotation_name_gp = gpar(fontsize = 14)
)

lgd1 = Legend(col_fun = colors.spectral,
    title = "Pseudotime",
    title_gp = gpar(fontsize = 14),
    labels = c("early","late"), at = c(min(Pseudotime),max(Pseudotime)),
    labels_gp = gpar(fontsize = 12),
    grid_width = unit(1, "cm"),
    legend_height = unit(3, "cm"))
lgd2 = Legend(col_fun = colors.inferno,
    title = "Gene\nexpression", 
    title_gp = gpar(fontsize = 14),
    labels = c("Low","High"), at = c(0,1),
    labels_gp = gpar(fontsize = 12),
    grid_width = unit(1, "cm"),
    legend_height = unit(3, "cm"))
pd = packLegend(lgd1,lgd2,column_gap = unit(1, "cm"),max_height = unit(30, "cm"))

mes.heatmap <- Heatmap(as.matrix(mes.avg),
        col = colors.inferno,
        column_title = levels(mes$Final1),
        column_split = factor(levels(mes$Final1),levels=levels(mes$Final1)),
        top_annotation = ha,
        row_split = factor(gene_split,levels=unique(gene_split)),
        row_title_gp = gpar(fontsize = 14),
        column_title_gp = gpar(fontsize = 14),
        row_names_gp = gpar(fontsize = 12,fontface='italic'),
        column_title_rot = 0,
        cluster_rows=F,
        cluster_columns=F,
        cluster_column_slices=F,
        cluster_row_slices=F,
        show_row_names=T,
        show_column_names = F,
        show_heatmap_legend = F)
p10 = grid.grabExpr(draw(mes.heatmap,annotation_legend_list = pd,padding = unit(c(2, 2, 10, 2), "mm")))
p10 <- plot_grid(NULL,p10,ncol=2,rel_widths=c(1,21))

# Proportion bar plot by Dataset4
data0 <- as.data.frame(table(mes$Dataset4))
data2 <- prop.table(table(mes$Dataset4,mes$Final1),1)*100
data2 <- as.data.frame(data2)
p12 <- ggplot(data2, aes(x = Var1, y = Freq, fill = Var2)) + 
    geom_col(position = "fill") + scale_fill_manual(values = mes_colors,name="") +
    scale_y_continuous(labels = scales::percent_format()) +
    xlab("") + ylab("% of dataset") +
    labs(title = 'Composition by Archetype') + 
    theme(plot.title = element_text(size = 30),
            legend.text = element_text(size = 14),
            axis.text.x=element_text(size=16,color='black'),  # angle = 45,vjust=0.5,
            axis.text.y = element_text(size = 14,color='black'),
            axis.title.x = element_text(size = 16,color='black'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.line = element_line(size = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
    geom_text(aes(x = Var1,y = 1,label = Freq,fill = NULL),color='red',
            data=data0,position = position_dodge(width = 0.5),vjust=-0.3,size=5) +
    scale_color_manual(values = "red")
p12 <- plot_grid(NULL,p12,ncol=2,rel_widths=c(1,21))

ps911 <- plot_grid(p9,p11,nrow=2)
botleft <- plot_grid(ps911,p10,NULL,p12,ncol=2,rel_widths=c(2,3),rel_heights=c(2,1))

res <- 400
png(paste('mes_figure_pheno_1.png'),width=1000*res/72,height=1100*res/72,res=res)
plot_grid(botleft)
dev.off()

# Gene Ontology bar plot of select mesenchymal archetypes
golist <- list()
golist$General <- c('0009888','0030154','0030198')
golist$Homeostatic <- c('0033554','0001666','0042594')
golist$Mesothelium <- c('0009888','0030198','0045229','0009612','0009790')
golist$Pericytes <- c('0072359','0048514','0001525','0001944','0048762')
gotypes <- levels(object$putative_labels)[levels(object$putative_labels) != 'Cycling']
allgo <- c()
for (i in 1:length(gotypes)) {
        type <- gotypes[i]
        message(type)
        go_cat <- golist[[type]]
        gosheet <- as.data.frame(read_excel("./mes.seurat.objs.int.marks.2022-10-14.xlsm",sheet = type))
        go <- gosheet[grep(paste(go_cat,collapse="|"),gosheet$'GO biological process complete'), ]
        go <- go[order(go$'upload_1 (raw P-value)'),]
        go$cluster <- type
        go[3] <- NULL
        allgo <- rbind(allgo,go)
}
allgo$p_val <- allgo$'upload_1 (raw P-value)'
allgo$go_term <- allgo$'GO biological process complete'
allgo$go_term2 <- paste(allgo$'GO biological process complete',allgo$cluster,sep='_')
allgo$go_term2 <- factor(allgo$go_term2,levels=allgo$go_term2)
allgo$cluster <- factor(allgo$cluster,levels=gotypes)
allgo$cluster2[allgo$cluster == 'General'] = 'Alveolar'
allgo$cluster2[allgo$cluster == 'Homeostatic'] = 'Aberrant'
allgo$cluster2[allgo$cluster == 'Mesothelium'] = 'Remodeling'
allgo$cluster2[allgo$cluster == 'Pericytes'] = 'Pericytes'
allgo$color <- type_colors[allgo$cluster2]
allgo$revrank <- rev(as.numeric(as.factor(seq.int(nrow(allgo)))))
ref_avg <- aggregate(revrank ~ cluster, allgo, mean)
allgo$breaks <- ref_avg[match(allgo$cluster, ref_avg$cluster), "revrank"]
allgo$breaks <- factor(allgo$breaks,levels=unique(allgo$breaks))
px <- ggplot(allgo,aes(x=factor(revrank),y=p_val,fill=go_term2,label=go_term)) +
    geom_col(position = position_dodge2(preserve='single',reverse=T),show.legend = F) + 
    scale_fill_manual(values = allgo$color) + coord_flip() +
    ylab('raw P-value (log scale)') + 
    scale_x_discrete(breaks = factor(levels(allgo$breaks)),
        labels = unique(allgo$cluster2)) +
    scale_y_continuous(trans=compose_trans("log10", "reverse"),
        expand = expansion(add = c(0,28)),  # 20
        breaks = trans_breaks('log10', function(x) 10^x),
        labels = trans_format('log10', math_format(10^.x))) +
    geom_text(position = position_dodge2(0.9,preserve='single',reverse=T),hjust=-0.02,size=6) +
    geom_hline(yintercept=0.05,linetype="dashed",color = "gray") +
    labs(title = 'GO Terms of Engineered Mesenchymal Archetypes') +
    theme(plot.title = element_text(size = 30,hjust = 0,face = "plain"),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_blank(),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
px <- plot_grid(NULL,px,ncol=2,rel_widths=c(1,11))

# Line plot of genes of interest
# (see lineplot function above)
goi_line <- c('Pdgfra','Axin2','Wnt5a','Wnt7b','Cthrc1','Egr1','Thbs2','Plin2','Thy1','Pparg','Sfrp1','Runx1')

Idents(mes) <- mes$Dataset4
mes.split <- SplitObject(mes,split.by = 'Final1')

goi <- goi_line
end <- 'Supp'

lineplot(dir_name='mes',object.split=mes.split,goi=goi,first='Tri',last='Quad')

res <- 400
png(paste0('mes_line_',end,'.png'),width=wd*res/72,height=ht*res/72,res=res)
plot_grid(p1,lgd,ncol=2,rel_widths=c((ncols+1),1))
dev.off()

# RNA Velocity plots - see focus.R script

# Supplemental FeaturePlots & Violin Plots - Mesenchyme
# (see featureplot & violinplot functions above)
Idents(mes) <- mes$Dataset4
mes.split <- SplitObject(mes,split.by = 'Final1')
mes.sub <- subset(mes,idents=c('Tri','Quad'))
mes.sub$Dataset4 <- droplevels(mes.sub$Dataset4)

goi <- goi_line
end <- 'Supp'

featureplot(object=mes,goi=goi)
violinplot(dir_name='mes',object=mes.sub,goi=goi)

res <- 400
png(paste0('mes_fpvln_',end,'3.png'),width=wd*res/72,height=(ht2)*res/72,res=res)
plot_grid(p2,ncol=1)
dev.off()




#### FIGURES S27-30 ####
### Immune Supplement
## Load samples, colors, set-up
# Load up-to-date macrophage object
load('./Imm/mac.seurat.objs.int.2023-07-27.Robj')
imm <- object ; rm(object)
imm$Dataset4 <- factor(imm$Dataset4,levels=c('Start','Tri','Quad'))
starteng_colors <- color_pals$starteng_colors
imm_colors <- color_pals[['imm']]


## Plots
# UMAP plot of integrated archetypes - colored by Final35
if (ncol(imm) < 10000){pt.size=1.5} else if (ncol(object) >= 10000){pt.size=NULL}
res <- 400
png(paste(dir_name,'_umap_2.png',sep=''), width=600*res/72, height=500*res/72,res=res)
UMAPPlot(imm,group.by = 'Final35',cols = type_colors,pt.size=pt.size,label = T,repel=T,label.size = 8) + NoLegend() +
    labs(title = 'Engineered Macrophages',caption = paste("nCells =",ncol(imm))) + 
    theme(plot.title = element_text(size = 30,hjust = 0,face = "plain"),
        plot.caption = element_text(size = 18,color='red'),
        axis.text.x=element_text(size=18,color='black'),
        axis.text.y = element_text(size = 18,color='black'),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'))
dev.off()

# UMAP plot of integrated archetypes - colored by Sample
if (ncol(imm) < 10000){pt.size=1.5} else if (ncol(imm) >= 10000){pt.size=NULL}
p13 <- UMAPPlot(imm,group.by = 'Sample',cols = sample_colors[levels(imm$Sample)],shuffle=T,pt.size=pt.size,label = F) + 
    labs(title = 'Integration',caption = paste("nCells =",ncol(imm))) + 
    theme(plot.title = element_text(size = 30,hjust = 0,face = "plain"),
        plot.caption = element_text(size = 14,color='red'),
        axis.text.x=element_text(size=14,color='black'),
        axis.text.y = element_text(size = 14,color='black'),
        axis.title.x = element_text(size = 16,color='black'),
        axis.title.y = element_text(size = 16,color='black'))

# Heatmap of select defining archetype markers (vertical) - AVERAGE
Idents(imm) <- imm$Final35
imm$Final35 <- factor(imm$Final35,levels=c('Alveolar','Cycling','Interstitial'))
imm$Sample <- factor(imm$Sample,levels=c('FB13','FB14','BAL','BEF12','BEFM1',
        'BEFM2','BEFM4','BEFM5','BEFM6'))
imm$Dataset3 <- factor(imm$Dataset3,levels=c('Start','Eng'))
imm <- ScaleData(imm,features = rownames(imm))
genes <- c('Mcemp1','Lipa','Plet1','Plin2','Cd59','Gpd1',
    'Txnip','Fcrla','S100a6','Abcg1','Ccl6','Naaa','Cela1',  # Alveolar
    'Mki67','Top2a',  # Cycling
    'C1qb','C1qa','C1qc','Pf4','Slamf9','Cd74','Apoe','Tmem37','Pmp22')  # Interstitial
gene_split <- c(rep('Alveolar',13),rep('Cycling',2),rep('Interstitial',9))
genes.use <- intersect(rownames(GetAssayData(imm,slot = 'scale.data')),genes)
imm.avg <- AverageExpression(imm,assays = 'RNA',slot='scale.data',features=genes.use)
imm.avg <- t(scale(t(imm.avg$RNA)))
unityNormalize <- function(x){(x-min(x))/(max(x)-min(x))}
imm.avg <- t(apply(as.matrix(imm.avg), 1, unityNormalize))
imm.avg <- imm.avg[,levels(imm$Final35)]
imm.avg <- imm.avg[genes,]
dim(imm.avg)

metadata <- imm@meta.data
Pseudotime <- metadata %>%
    group_by(Final35) %>%
    summarize(average_pseudotime = mean(Lineage_avg, na.rm = TRUE))
Pseudotime <- setNames(Pseudotime$average_pseudotime, Pseudotime$Final35)[levels(imm$Final35)]
colors.spectral <- colorRamp2(breaks = c(seq(40,70,length.out=60)),colorRampPalette(brewer.pal(11,'Spectral')[-6])(60), space = "RGB")
colors.inferno <- colorRamp2(breaks = c(seq(min(imm.avg),max(imm.avg),length.out=60)), inferno(n=60), space = "RGB")

ha = HeatmapAnnotation(CellType = anno_block(gp = gpar(fill=imm_colors[levels(imm$Final35)]),show_name=T),
        Pseudotime = anno_block(gp = gpar(fill = colors.spectral(Pseudotime)),show_name=T),
        annotation_name_gp = gpar(fontsize = 14)
)

lgd1 = Legend(col_fun = colors.spectral,
    title = "Pseudotime",
    title_gp = gpar(fontsize = 14),
    labels = c("early","late"), at = c(min(Pseudotime),max(Pseudotime)),
    labels_gp = gpar(fontsize = 12),
    grid_width = unit(1, "cm"),
    legend_height = unit(3, "cm"))
lgd2 = Legend(col_fun = colors.inferno,
    title = "Gene\nexpression", 
    title_gp = gpar(fontsize = 14),
    labels = c("Low","High"), at = c(0,1),
    labels_gp = gpar(fontsize = 12),
    grid_width = unit(1, "cm"),
    legend_height = unit(3, "cm"))
pd = packLegend(lgd1,lgd2,column_gap = unit(1, "cm"),max_height = unit(30, "cm"))

imm.heatmap <- Heatmap(as.matrix(imm.avg),
        col = colors.inferno,
        column_title = levels(imm$Final35),
        column_split = factor(levels(imm$Final35),levels=levels(imm$Final35)),
        top_annotation = ha,
        row_split = factor(gene_split,levels=unique(gene_split)),
        row_title_gp = gpar(fontsize = 14),
        column_title_gp = gpar(fontsize = 14),
        row_names_gp = gpar(fontsize = 12,fontface='italic'),
        column_title_rot = 0,
        cluster_rows=F,
        cluster_columns=F,
        cluster_row_slices=F,
        show_row_names=T,
        show_column_names = F,
        show_heatmap_legend = F)
p14 = grid.grabExpr(draw(imm.heatmap,annotation_legend_list = pd,padding = unit(c(2, 2, 10, 2), "mm")))
p14 <- plot_grid(NULL,p14,ncol=2,rel_widths=c(1,21))

# Proportion bar plot by Dataset4
imm$Final35 <- factor(imm$Final35,levels=c('Alveolar','Cycling','Interstitial'))
data0 <- as.data.frame(table(imm$Dataset4))
data2 <- prop.table(table(imm$Dataset4,imm$Final35),1)*100
data2 <- as.data.frame(data2)
p16 <- ggplot(data2, aes(x = Var1, y = Freq, fill = Var2)) + 
    geom_col(position = "fill") + scale_fill_manual(values = imm_colors,name="") +
    scale_y_continuous(labels = scales::percent_format()) +
    xlab("") + ylab("% of dataset") +
    labs(title = 'Composition by Cell Type') + 
    theme(plot.title = element_text(size = 30),
            legend.text = element_text(size = 14),
            axis.text.x=element_text(size=16,color='black'),  # angle = 45,vjust=0.5,
            axis.text.y = element_text(size = 14,color='black'),
            axis.title.x = element_text(size = 16,color='black'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.line = element_line(size = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
    geom_text(aes(x = Var1,y = 1,label = Freq,fill = NULL),color='red',
            data=data0,position = position_dodge(width = 0.5),vjust=-0.3,size=5)+
    scale_color_manual(values = "red")
p16 <- plot_grid(NULL,p16,ncol=2,rel_widths=c(1,21))

ps1315 <- plot_grid(p13,p15,nrow=2)
botright <- plot_grid(ps1315,p14,NULL,p16,ncol=2,rel_widths=c(2,3),rel_heights=c(2,1))

res <- 400
png(paste('imm_figure_pheno_1.png'),width=1000*res/72,height=1100*res/72,res=res)
plot_grid(botright)
dev.off()

# Line plot of genes of interest
# (see lineplot function above)
goi_line <- c('Fabp4','Mertk','Lgmn','Il1b','Siglec5','Spp1','Lpl','Fdx1','Mt2A','Cxcl2','Lcn2','Clec10a')

Idents(imm) <- imm$Dataset4
imm.split <- SplitObject(imm,split.by = 'Final35')

goi <- goi_line
end <- 'Supp'

lineplot(dir_name='imm',object.split=imm.split,goi=goi,first='Start',last='Quad')

res <- 400
png(paste0('imm_line_',end,'.png'),width=wd*res/72,height=ht*res/72,res=res)
plot_grid(p1,lgd,ncol=2,rel_widths=c((ncols+1),1))
dev.off()

# RNA Velocity plots - see focus.R script

# Supplemental FeaturePlots & Violin Plots - Immune
# (see featureplot & violinplot functions above)
Idents(imm) <- imm$Dataset4
imm.split <- SplitObject(imm,split.by = 'Final35')
imm.sub <- subset(imm,idents=c('Start','Quad'))
imm.sub$Dataset4 <- droplevels(imm.sub$Dataset4)
goi_Final <- c('Fabp4','Mertk','Lgmn','Il1b','Siglec5','Spp1','Lpl','Fdx1')  # 8

goi <- goi_Final
end <- 'Final'

featureplot(object=imm,goi=goi)
violinplot(dir_name='imm',object=imm.sub,goi=goi)

res <- 400
png(paste0('imm_fpvln_',end,'.png'),width=wd*res/72,height=(ht1+ht2)*res/72,res=res)
plot_grid(p1,p2,ncol=1)
dev.off()

goi_act <- c('Bst2','Cd9','Xrcc5','Cdc42ep3','Pla2g2d','Slc39a2','Lpl','Rgcc','Ly6al',
        'Crip1','S100a4','Scgb1a1',  # Start_Alveolar
        'Stmn1','Tuba1b','Lgals1','Tmem176b','Tmem176a','Limd2','Ccl2','Cd14','Ccl7',
        'Ccl12','Spp1','Ms4a7',  # Start_Interstitial
        'Slc11a1','Fabp4','Mt1','Mt2A','Cxcl2','Cxcl1','Fos','Clec10a','Vsig4','Retn',
        'Clca4l','S100a9','Lcn2','Mgp',  # Eng_Alveolar
        'Ctsl','Ifitm1','Hopx','Pla2g7','Tf','Slpi','Gsta1','F13a1','Gpx3','Sparc',
        'Col1a1')  # Eng_Interstitial 

goi <- goi_act
end <- 'Act'

featureplot(object=imm,goi=goi)
violinplot(dir_name='imm',object=imm.sub,goi=goi)

res <- 400
png(paste0('imm_fpvln_',end,'1.png'),width=wd*res/72,height=ht1*res/72,res=res)
plot_grid(p1,ncol=1)
dev.off()

png(paste0('imm_fpvln_',end,'2.png'),width=wd*res/72,height=ht2*res/72,res=res)
plot_grid(p2,ncol=1)
dev.off()




#### FIGURE S31 ####
### Quad-culture comparison to native rat lung
## Load native rat lung dataset and quad-culture engineered data
load('./Native.clean.2022-08-22.Robj')
Idents(object) <- object$Chemistry
native <- subset(object,idents='v3')
native$Chemistry <- droplevels(native$Chemistry)
rm(object)
Idents(native) <- native$CellType_Final
native <- RenameIdents(native,
    'T_Killer'='NK',
    'Mac_Alv_Pro'='Cycling_Immune',
    'Smooth_Muscle'='SMC')
native$CellType_Final2 <- Idents(native)
native$CellType_Final2 <- factor(native$CellType_Final2,levels=c(
    'ATII','ATII-ATI','ATI','BASC','Tuft','Secretory','Ciliated',
    'gCap','aCap','Arterial','Venous','Lymphatic',
    'Col14_Fibroblasts','Col13_Fibroblasts','Lgr5_Fibroblasts','Mesothelium','Myofibroblasts','Pericytes','SMC',
    'Mac_Alv','Mac_Int','Mo_C','Mo_NC','Mo_Act','DC','Neutro','Baso','Eos','Mast',
    'ILC','T','NK','B','Plasma','Cycling_Immune'))

load('./bef.seurat.objs.int.2023-04-06.Robj')
Idents(object) <- object$Dataset2
quad <- subset(object,idents='Quad_E')
quad$Dataset2 <- droplevels(quad$Dataset2)
rm(object)
Idents(quad) <- quad$Final2
quad <- RenameIdents(quad,
    'Cycling_Alveolar_Macrophage'='Cycling_Mac_Alv',
    'Alveolar_Macrophage'='Mac_Alv',
    'Interstitial_Macrophage'='Mac_Int',
    'Cycling_Lymphatic'='Cycling_Endo',
    'Cycling_Microvascular'='Cycling_Endo')
quad$Final3 <- Idents(quad)
quad$Final3 <- factor(quad$Final3,levels=c(
    'Basal-like','Transitional','ATI-like','BASC','Cycling_Epi',
    'Microvascular','Endo_Progenitor','Lymphatic','Cycling_Endo',
    'Alveolar_Fibroblast','Aberrant_Fibroblast','Remodeling_Fibroblast','Pericytes','Cycling_Mes',
    'Mac_Alv','Mac_Int','Cycling_Mac_Alv'))


## Generate embeddings of each object
SeuratCluster_2 <- function(sample.name,object,assay,slot,pcs,res){
    DefaultAssay(object) <- assay
    object <- NormalizeData(object)  # comment out for engineered
    object <- FindVariableFeatures(object)  # comment out for engineered
    object <- ScaleData(object)  # comment out for engineered
    object <- RunPCA(object, npcs = 100, verbose = F)
    object <- FindNeighbors(object, dims = 1:pcs)
    object <- FindClusters(object, resolution = res)
    object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
    p <- UMAPPlot(object,group.by=slot,label=T,repel=T) + NoLegend() + 
        labs(title = paste(sample.name,'-',slot),subtitle = paste("PC's =",pcs),caption = paste("nCells =",ncol(object))) + 
        theme(plot.title = element_text(size = 24),
            plot.subtitle = element_text(size=20,color='red'),
            plot.caption = element_text(size=20,color='red'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.title.x = element_text(size = 16,color='black'),
            axis.text.x = element_text(size=12,color='black'))
    cluster.object <<- object
    p <<- p
}

# Native Rat Lung v3 (14,959 cells)
sample.name <- 'NativeRatLung'
object <- native
assay <- 'RNA'
slot <- 'CellType_Final'
pcs <- 42
res <- 0.5
tag <- paste0('.pcs.',pcs)

SeuratCluster_2(sample.name=sample.name,object=object,assay=assay,slot=slot,pcs=pcs,res=res)
native <- cluster.object

# Quad-culture engineered rat lung (10,473 cells)
sample.name <- 'EngRatLung'
object <- quad
assay <- 'integrated'
slot <- 'Final2'
pcs <- 45
res <- 0.5
tag <- paste0('.pcs.',pcs)

# modify function first
SeuratCluster_2(sample.name=sample.name,object=object,assay=assay,slot=slot,pcs=pcs,res=res)
quad <- cluster.object


## Make new metadata slots for plotting
# Native
Idents(native) <- native$CellType_Final2
native <- RenameIdents(native,
    'ATII'='1','ATII-ATI'='2','ATI'='3','BASC'='4','Tuft'='5','Secretory'='6','Ciliated'='7','gCap'='8','aCap'='9',
    'Arterial'='10','Venous'='11','Lymphatic'='12','Col14_Fibroblasts'='13','Col13_Fibroblasts'='14','Lgr5_Fibroblasts'='15',
    'Mesothelium'='16','Myofibroblasts'='17','Pericytes'='18','SMC'='19','Mac_Alv'='20','Mac_Int'='21','Mo_C'='22',
    'Mo_NC'='23','Mo_Act'='24','DC'='25','Neutro'='26','Baso'='27','Eos'='28','Mast'='29','ILC'='30','T'='31',
    'NK'='32','B'='33','Plasma'='34','Cycling_Immune'='35')
native$CellType_Finalnum <- Idents(native)
native$CellType_Finalnum <- factor(native$CellType_Finalnum,levels=levels(native$CellType_Finalnum))

Idents(native) <- native$CellType_Final2
native <- RenameIdents(native,
    'ATII'='Nat_ATII','ATII-ATI'='Nat_ATII-ATI','ATI'='Nat_ATI','BASC'='Nat_BASC','Tuft'='Nat_Tuft','Secretory'='Nat_Secretory',
    'Ciliated'='Nat_Ciliated','gCap'='Nat_gCap','aCap'='Nat_aCap','Arterial'='Nat_Arterial','Venous'='Nat_Venous',
    'Lymphatic'='Nat_Lymphatic','Col14_Fibroblasts'='Nat_Col14_Fibroblasts','Col13_Fibroblasts'='Nat_Col13_Fibroblasts',
    'Lgr5_Fibroblasts'='Nat_Lgr5_Fibroblasts','Mesothelium'='Nat_Mesothelium','Myofibroblasts'='Nat_Myofibroblasts',
    'Pericytes'='Nat_Pericytes','SMC'='Nat_SMC','Mac_Alv'='Nat_Mac_Alv','Mac_Int'='Nat_Mac_Int','Mo_C'='Nat_Mo_C',
    'Mo_NC'='Nat_Mo_NC','Mo_Act'='Nat_Mo_Act','DC'='Nat_DC','Neutro'='Nat_Neutro','Baso'='Nat_Baso','Eos'='Nat_Eos',
    'Mast'='Nat_Mast','ILC'='Nat_ILC','T'='Nat_T','NK'='Nat_NK','B'='Nat_B','Plasma'='Nat_Plasma','Cycling_Immune'='Nat_Cycling_Immune')
native$heatmapid <- Idents(native)
native$heatmapid <- factor(native$heatmapid,levels=levels(native$heatmapid))

# Quad
Idents(quad) <- quad$Final3
quad <- RenameIdents(quad,
    'Basal-like'='1','Transitional'='2','ATI-like'='3','BASC'='4','Cycling_Epi'='5','Microvascular'='6','Endo_Progenitor'='7',
    'Lymphatic'='8','Cycling_Endo'='9','Alveolar_Fibroblast'='10','Aberrant_Fibroblast'='11','Remodeling_Fibroblast'='12',
    'Pericytes'='13','Cycling_Mes'='14','Mac_Alv'='15','Mac_Int'='16','Cycling_Mac_Alv'='17')
quad$Finalnum <- Idents(quad)
quad$Finalnum <- factor(quad$Finalnum,levels=levels(quad$Finalnum))

Idents(quad) <- quad$Final3
quad <- RenameIdents(quad,
    'Basal-like'='Quad_Basal-like','Transitional'='Quad_Transitional','ATI-like'='Quad_ATI-like','BASC'='Quad_BASC',
    'Cycling_Epi'='Quad_Cycling_Epi','Microvascular'='Quad_Microvascular','Endo_Progenitor'='Quad_Endo_Progenitor',
    'Lymphatic'='Quad_Lymphatic','Cycling_Endo'='Quad_Cycling_Endo','Alveolar_Fibroblast'='Quad_Alveolar_Fibroblast',
    'Aberrant_Fibroblast'='Quad_Aberrant_Fibroblast','Remodeling_Fibroblast'='Quad_Remodeling_Fibroblast',
    'Pericytes'='Quad_Pericytes','Cycling_Mes'='Quad_Cycling_Mes','Mac_Alv'='Quad_Mac_Alv','Mac_Int'='Quad_Mac_Int',
    'Cycling_Mac_Alv'='Quad_Cycling_Mac_Alv')
quad$heatmapid <- Idents(quad)
quad$heatmapid <- factor(quad$heatmapid,levels=levels(quad$heatmapid))


## Panel A - initial umaps
class_colors <- c('Epithelium'='#75C359',
    'Endothelium'='#5B9ADA',
    'Mesenchyme'='#E64955',
    'Immune'='#8949CA')
natcell_colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(length(levels(native$CellType_Finalnum)))
quadcell_colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(length(levels(quad$Finalnum)))
native$class <- factor(native$class,levels=c('Epithelium','Endothelium','Mesenchyme','Immune'))
quad$class <- factor(quad$class,levels=c('Epithelium','Endothelium','Mesenchyme','Immune'))

p1 <- UMAPPlot(native,group.by='CellType_Finalnum',label=T,label.size=5,repel=T) + NoLegend() +
        labs(title = paste('Native Rat Lung - Cell Type'),caption = paste("nCells =",ncol(native))) + 
        theme(plot.title = element_text(size = 22,hjust = 0),
            plot.caption = element_text(size=20,color='red'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.title.x = element_text(size = 16,color='black'),
            axis.text.x = element_text(size=12,color='black'))
lgd1 <- UMAPPlot(native,group.by='CellType_Final2') +
        guides(color = guide_legend(override.aes = list(size=4), ncol=1)) +
        theme(legend.text = element_text(size = 14),legend.key.size = unit(1.1, "lines"))
lgd1 <- get_legend(lgd1)
p2 <- UMAPPlot(quad,group.by='Finalnum',label=T,label.size=5,repel=T) + NoLegend() +
        labs(title = paste('Quad Culture Engineered Rat Lung - Cell Type'),caption = paste("nCells =",ncol(quad))) + 
        theme(plot.title = element_text(size = 22,hjust = 0),
            plot.caption = element_text(size=20,color='red'),
            axis.title.y = element_text(size = 16,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.title.x = element_text(size = 16,color='black'),
            axis.text.x = element_text(size=12,color='black'))
lgd2 <- UMAPPlot(quad,group.by='Final3') +
        guides(color = guide_legend(override.aes = list(size=4), ncol=1)) +
        theme(legend.text = element_text(size = 14),legend.key.size = unit(1.1, "lines"))
lgd2 <- get_legend(lgd2)
p3 <- UMAPPlot(native,group.by='class',cols=class_colors,label=F) + NoLegend() +
        labs(title = paste('Native Rat Lung - Class')) + 
        theme(plot.title = element_text(size = 22,hjust = 0),
            axis.title.y = element_text(size = 16,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.title.x = element_text(size = 16,color='black'),
            axis.text.x = element_text(size=12,color='black'))
lgd3 <- UMAPPlot(native,group.by='class',cols=class_colors) +
        guides(color = guide_legend(override.aes = list(size=4), ncol=1)) +
        theme(legend.text = element_text(size = 14),legend.key.size = unit(1.1, "lines"))
lgd3 <- get_legend(lgd3)
p4 <- UMAPPlot(quad,group.by='class',cols=class_colors,label=F) + NoLegend() +
        labs(title = paste('Quad Culture Engineered Rat Lung - Class')) + 
        theme(plot.title = element_text(size = 22,hjust = 0),
            axis.title.y = element_text(size = 16,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.title.x = element_text(size = 16,color='black'),
            axis.text.x = element_text(size=12,color='black'))
toprow <- plot_grid(p1,lgd1,p2,lgd2,ncol=4,rel_widths=c(1,0.3,1,0.4))
botrow <- plot_grid(p3,lgd3,p4,lgd3,ncol=4,rel_widths=c(1,0.3,1,0.4))
umaps <- plot_grid(toprow,botrow,nrow=2)

res <- 400
png(paste('NatEngComp_embeddings_2.png'), width=1300*res/72,height=1200*res/72,res=res)
umaps
dev.off()


## Panel B - correlation between Native cell types and Quad culture labels
library(reshape2)
library(corrplot)
library(seriation)

reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
}

# Merge native and quad objects
# Remove irrelevant native immune populations
Idents(native) <- native$heatmapid
native <- subset(native,idents=c('Nat_Mo_C','Nat_Mo_NC','Nat_Mo_Act','Nat_DC','Nat_Neutro','Nat_Baso','Nat_Eos','Nat_Mast',
    'Nat_ILC','Nat_T','Nat_NK','Nat_B','Nat_Plasma'),invert=T)

# 2000 var genes in native (without irrelevant immune populations)
native <- FindVariableFeatures(native,nfeatures = 2000)
nat.vargenes <- VariableFeatures(native)
# Use only native variable genes that are also expressed in quad
genelist <- intersect(nat.vargenes,rownames(quad))
length(genelist)

native <- ScaleData(native,features = genelist)
quad <- ScaleData(quad,features = genelist)
native.avg <- AverageExpression(native,assays='RNA',group.by='heatmapid',slot='scale.data')
quad.avg <- AverageExpression(quad,assays='RNA',group.by='heatmapid',slot='scale.data')
native.mtx <- as.matrix(native.avg[[1]])
quad.mtx <- as.matrix(quad.avg[[1]])
merge.mtx <- cbind(native.mtx,quad.mtx)
merge.cor <- cor(merge.mtx,method = 'spearman')

# Order with optimal leaf ordering - using corrplot & seriation packages
dist2order = function(corr, method, ...) {
  d_corr = as.dist(1 - corr)
  s = seriate(d_corr, method = method, ...)
  i = get_order(s)
  return(i)
}

i = dist2order(merge.cor, 'OLO_average')

names = rbind(c(rep('Nat_Mac_Int',2),rep('Quad_Mac_Int',2)),
            c(rep('Quad_Cycling_Epi',2),rep('Nat_ATI',2)),
            c(rep('Quad_Cycling_Endo',2),rep('Nat_Lymphatic',2)),
            c(rep('Nat_Mesothelium',2),rep('Quad_Pericytes',2)))

res <- 400
png('NatEngComp_corrplot_indivnatvar_OLOavg.png', width=900*res/72,height=900*res/72,res=res)
corrplot(merge.cor[i,i],is.corr=F,method='color',col=rev(COL2('RdBu')),
    tl.col='black',cl.cex=1) %>% 
    corrRect(namesMat = names,col=c('#8949CA','#75C359','#5B9ADA','#E64955'),lwd=4)
dev.off()


## Panel C - Dot plot of canonical cell type markers in native and engineered
# Remove irrelevant native immune populations
Idents(native) <- native$heatmapid
native <- subset(native,idents=c('Nat_Mo_C','Nat_Mo_NC','Nat_Mo_Act','Nat_DC','Nat_Neutro','Nat_Baso','Nat_Eos','Nat_Mast',
    'Nat_ILC','Nat_T','Nat_NK','Nat_B','Nat_Plasma'),invert=T)
merge <- merge(x=native,y=quad)
DefaultAssay(merge) <- 'RNA'

merge$class <- factor(merge$class,levels=c('Epithelium','Endothelium','Mesenchyme','Immune'))
class_colors <- c('Epithelium'='#75C359',
    'Endothelium'='#5B9ADA',
    'Mesenchyme'='#E64955',
    'Immune'='#8949CA')
green <- colorRampPalette(brewer.pal(9,'Greens'))(2)
blue <- colorRampPalette(brewer.pal(9,'Blues'))(2)
red <- colorRampPalette(brewer.pal(9,'Reds'))(2)
purple <- colorRampPalette(brewer.pal(9,'Purples'))(2)

goi <- c(
    'Epcam','Aqp5','Pdpn','Ager','Edn3','Sftpc','Sftpb','Abca3',
    'Scgb1a1','Scgb3a2','Ccdc113','Krt5',
    'Cdh5','Ptprb','Gpihbp1','Adgre5','Kdr','Apln','Emcn',
    'Efnb2','Gja4','Fbln5','Vcam1','Mmrn1',
    'Col1a1','Pdgfra','Dcn','Itga8','Rspo1',
    'Myh11','Acta2','Pdgfrb','Gucy1b1',
    'Ptprc','Naaa','Mcemp1','C1qc','C1qb',
    'Mki67','Top2a'
)

genes <- goi

Idents(merge) <- merge$heatmapid
merge$heatmapid <- factor(merge$heatmapid,levels=c(
    'Nat_ATI','Quad_ATI-like','Nat_ATII-ATI','Nat_ATII','Nat_BASC','Quad_BASC','Nat_Tuft','Nat_Secretory','Nat_Ciliated',
    'Quad_Basal-like','Quad_Transitional','Quad_Cycling_Epi',
    'Nat_gCap','Nat_aCap','Nat_Arterial','Nat_Venous','Quad_Microvascular','Quad_Endo_Progenitor',
    'Nat_Lymphatic','Quad_Lymphatic','Quad_Cycling_Endo',
    'Nat_Col14_Fibroblasts','Nat_Lgr5_Fibroblasts','Nat_Col13_Fibroblasts','Quad_Alveolar_Fibroblast','Quad_Aberrant_Fibroblast',
    'Quad_Remodeling_Fibroblast','Nat_Mesothelium','Nat_Myofibroblasts','Nat_SMC','Nat_Pericytes',
    'Quad_Pericytes','Quad_Cycling_Mes',
    'Nat_Mac_Alv','Quad_Mac_Alv','Nat_Mac_Int','Quad_Mac_Int','Nat_Cycling_Immune','Quad_Cycling_Mac_Alv')
)
natcol <- dataset_colors[['Start']]
engcol <- dataset_colors[['Quad_E']]
yaxiscols <- rev(c(natcol,engcol,natcol,natcol,natcol,engcol,natcol,natcol,natcol,engcol,engcol,engcol,natcol,natcol,natcol,
    natcol,engcol,engcol,natcol,engcol,engcol,natcol,natcol,natcol,engcol,engcol,engcol,natcol,natcol,natcol,natcol,
    engcol,engcol,natcol,engcol,natcol,engcol,natcol,engcol))

p <- DotPlot(object=merge,features=genes,group.by='heatmapid',split.by='class',cols=class_colors,
    col.max=1,dot.scale=8,scale.min=10,scale.max=90,assay='RNA') +  
    scale_y_discrete(label=rev(levels(merge$heatmapid)),limits=rev) + NoLegend() +
    theme(axis.text.x=element_text(size = 18,angle = 45,vjust=1,hjust=1,face = 'italic',color = 'black'),
        axis.text.y = element_text(size = 18, color = yaxiscols),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(color = "black")) +
    geom_hline(yintercept=6.5,linetype="dashed",color = "darkgray") +
    geom_hline(yintercept=18.5,linetype="dashed",color = "darkgray") +
    geom_hline(yintercept=27.5,linetype="dashed",color = "darkgray") +
    geom_vline(xintercept=12.5,linetype="dashed",color = "darkgray") +
    geom_vline(xintercept=24.5,linetype="dashed",color = "darkgray") +
    geom_vline(xintercept=33.5,linetype="dashed",color = "darkgray") +
    geom_vline(xintercept=38.5,linetype="dashed",color = "darkgray")
q <- DotPlot(object=merge,features=genes,group.by='heatmapid',cols=c("lightgrey",class_colors[[1]]),
    col.max=1,dot.scale=8,scale.min=10,scale.max=90,assay='RNA') +  
    guides(color = guide_colorbar(title = 'Average Expression\n   \u2014 Epithelium',order=2)) +
    guides(size = "none") +
    theme(legend.title = element_text(size=18),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.1, "lines"))
lgd1 <- get_legend(q)
q <- DotPlot(object=merge,features=genes,group.by='heatmapid',cols=c("lightgrey",class_colors[[2]]),
    col.max=1,dot.scale=8,scale.min=10,scale.max=90,assay='RNA') +  
    guides(color = guide_colorbar(title = 'Average Expression\n   \u2014 Endothelium')) +
    guides(size = "none") +
    theme(legend.title = element_text(size=18),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.1, "lines"))
lgd2 <- get_legend(q)
q <- DotPlot(object=merge,features=genes,group.by='heatmapid',cols=c("lightgrey",class_colors[[3]]),
    col.max=1,dot.scale=8,scale.min=10,scale.max=90,assay='RNA') +  
    guides(color = guide_colorbar(title = 'Average Expression\n   \u2014 Mesenchyme')) +
    guides(size = "none") +
    theme(legend.title = element_text(size=18),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.1, "lines"))
lgd3 <- get_legend(q)
q <- DotPlot(object=merge,features=genes,group.by='heatmapid',cols=c("lightgrey",class_colors[[4]]),
    col.max=1,dot.scale=8,scale.min=10,scale.max=90,assay='RNA') +  
    guides(color = guide_colorbar(order=1,title = 'Average Expression\n   \u2014 Immune')) +
    guides(size = guide_legend(order = 2,title = 'Percent Expressed')) +
    theme(legend.title = element_text(size=18),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.1, "lines"))
lgd4 <- get_legend(q)
lgd <- plot_grid(NULL,lgd1,lgd2,lgd3,lgd4,NULL,ncol=1,rel_heights=c(1,1,1,1,1.6,1))

res <- 400
png('NatEngComp_dot_4.png', width=2100*res/72,height=1000*res/72,res=res)
plot_grid(p,NULL,lgd,ncol=3,rel_widths=c(0.97,0.02,0.097))
dev.off()





