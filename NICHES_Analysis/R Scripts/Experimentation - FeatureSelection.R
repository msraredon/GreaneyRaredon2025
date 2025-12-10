# Some experiments

# 1. What signals do Macrophages BRING TO Epi? Endo? Mes?
# We want to identify markers of these unique vectors, globally
temp <- marks$global$vector
new.mac.epi <- temp[temp$cluster=='Immune - Epithelium',]
new.mac.end <- temp[temp$cluster=='Immune - Endothelium',]
new.mac.mes <- temp[temp$cluster=='Immune - Mesenchyme',]
# Quick test
VlnPlot(global,'Tnfsf13—Tnfrsf1a',group.by = 'Vector',split.by = 'Condition',pt.size = 0) # GOOD!

# 2. What signals do Macrophages LISTEN TO from Epi? Endo? Mes?
# We want to identify markers of these unique vectors, globally
temp <- marks$global$vector
new.epi.mac <- temp[temp$cluster=='Epithelium - Immune',]
new.end.mac <- temp[temp$cluster=='Endothelium - Immune',]
new.mes.mac <- temp[temp$cluster=='Mesenchyme - Immune',]
# Quick test
VlnPlot(global,'Csf1—Csf1r',group.by = 'Vector',split.by = 'Condition',pt.size = 0) # GOOD!

# 3. How does the Epithelial niche change between condition?
# Within a given niche, what features are up or down regulated by condition?
temp <- marks$niche$condition$Epithelium
tri.epi.niche <- temp[temp$cluster=='Tri_L',]
quad.epi.niche <- temp[temp$cluster=='Quad_E',]
# Quick test
VlnPlot(global,'Thbs2—Itgb1',group.by = 'Vector',split.by = 'Condition',pt.size = 0) # GOOD!

# 4. How is Mesenchymal to Epithelial communication affected by the Addition of Macrophages?
# We are interested in the Mes-Epi vector, looking for features that are differential by condition
temp <- marks$vector$condition$`Mesenchyme - Epithelium`
tri.mes.epi <- temp[temp$cluster=='Tri_L',]
quad.mes.epi <- temp[temp$cluster=='Quad_E',]
VlnPlot(global,'Tgfb1—Itgb6',group.by = 'Vector',split.by = 'Condition',pt.size = 0) # GOOD!

# 5. How is Endothelial to Mesenchymal communication affected by the Addition of Macrophages?
# We are interested in the End-Mes vector, looking for features that are differential by condition
temp <- marks$vector$condition$`Endothelium - Mesenchyme`
tri.end.mes <- temp[temp$cluster=='Tri_L',]
quad.end.mes <- temp[temp$cluster=='Quad_E',]
VlnPlot(global,'Angptl4—Tie1',group.by = 'Vector',split.by = 'Condition',pt.size = 0) # GOOD!

# 6. What marks regenerative epithelial-autocrine signaling, globally?
# We are interested in the global perspective, and want to look at one specific archetype, epi-autocrine regenerative:
temp <- marks$global$archetype
epi.auto.regen <- temp[temp$cluster=='Epi-Auto Regenerative',]
VlnPlot(global,'Lamb3—Itga3',group.by = 'Vector',split.by = 'Condition',pt.size = 0) # GOOD!

# 7. Are there instances of markers 'moving' due to the existence of macrophages? [THIS IS COOL]
# We want features that mark specific vectors within tri that then mark *different* vectors in quad
temp1 <- marks$condition$vector$Tri_L
temp2 <- marks$condition$vector$Quad_E
# Create search index
temp1$index <- paste(temp1$cluster,temp1$gene)
temp2$index <- paste(temp2$cluster,temp2$gene)

# Optional: exclude immune vectors from this??
immune.vectors <- c('Immune - Epithelium',
                    'Immune - Endothelium',
                    'Immune - Mesenchyme',
                    'Immune - Immune',
                    'Epithelium - Immune',
                    'Endothelium - Immune',
                    'Mesenchyme - Immune')
temp1 <- temp1[!(temp1$cluster %in% immune.vectors),]
temp2 <- temp2[!(temp2$cluster %in% immune.vectors),]

# Add p threshold
p.thresh <- 0.05
temp1 <- temp1[temp1$p_val_adj<p.thresh,]
temp2 <- temp2[temp2$p_val_adj<p.thresh,]

# Indexes (vector-mechanism combinations) that are maintained, lost, or gained
things.maintained <-  temp1[temp1$index %in% temp2$index,]
things.lost <-  temp1[!(temp1$index %in% temp2$index),]
things.gained <-  temp2[!(temp2$index %in% temp1$index),]

# Features that stop marking one vector and start marking another (organized by 'markerness' in Tri)
things.moved1 <- things.lost[things.lost$gene %in% things.gained$gene,]

# Features that start marking one vector and used to mark another (organized by 'markerness' in Quad)
things.moved2 <- things.gained[things.gained$gene %in% things.lost$gene,]

# Quick test
VlnPlot(global,'Tgfb1—Tgfbr1',group.by = 'Vector',split.by = 'Condition',pt.size = 0) # GOOD!
VlnPlot(global,'Tgfb2—Tgfbr3',group.by = 'Vector',split.by = 'Condition',pt.size = 0) # GOOD!
VlnPlot(global,'Efnb1—Ephb2',group.by = 'Vector',split.by = 'Condition',pt.size = 0) # GOOD! # epithelium starts listening to mesenchyme on this pathway

