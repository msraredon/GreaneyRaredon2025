# legends for BEFM
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics")

load('Saved Objects/transcriptome.TriQuad.customEmbed.fineFilter.proofRead.2024-03-15.Robj') 

dat <- table(Idents(sub.fine),sub.fine$Condition)
dat/colSums(dat)
sqrt(dat/colSums(dat)) # sqrt transform used for size


