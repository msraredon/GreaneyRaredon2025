# This script loads previously computed marker lists and annotates mechanisms based on a set of user-defined rules

# Load marker lists
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/Saved Objects/marks.2024-03-24.Robj")

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

## Annotate each marker list (clumsy but works)

# Global (all)
iterations <- names(marks$global)
for(i in 1:length(iterations)){
temp <- marks$global[[iterations[i]]]
mechanisms.to.annotate <- temp$gene
annotations.to.add <- fantom5.rat.annotations[mechanisms.to.annotate,]
temp.out <- cbind(temp,annotations.to.add)
marks$global[[iterations[i]]] <- temp.out
}

## Condition first
# Condition.Vector
iterations <- names(marks$condition$vector)
for(i in 1:length(iterations)){
  temp <- marks$condition$vector[[iterations[i]]]
  mechanisms.to.annotate <- temp$gene
  annotations.to.add <- fantom5.rat.annotations[mechanisms.to.annotate,]
  temp.out <- cbind(temp,annotations.to.add)
  marks$condition$vector[[iterations[i]]] <- temp.out
}
# Condition.Niche
iterations <- names(marks$condition$niche)
for(i in 1:length(iterations)){
  temp <- marks$condition$niche[[iterations[i]]]
  mechanisms.to.annotate <- temp$gene
  annotations.to.add <- fantom5.rat.annotations[mechanisms.to.annotate,]
  temp.out <- cbind(temp,annotations.to.add)
  marks$condition$niche[[iterations[i]]] <- temp.out
}
# Condition.Influence
iterations <- names(marks$condition$influence)
for(i in 1:length(iterations)){
  temp <- marks$condition$influence[[iterations[i]]]
  mechanisms.to.annotate <- temp$gene
  annotations.to.add <- fantom5.rat.annotations[mechanisms.to.annotate,]
  temp.out <- cbind(temp,annotations.to.add)
  marks$condition$influence[[iterations[i]]] <- temp.out
}
# Condition.archetype
iterations <- names(marks$condition$archetype)
for(i in 1:length(iterations)){
  temp <- marks$condition$archetype[[iterations[i]]]
  mechanisms.to.annotate <- temp$gene
  annotations.to.add <- fantom5.rat.annotations[mechanisms.to.annotate,]
  temp.out <- cbind(temp,annotations.to.add)
  marks$condition$archetype[[iterations[i]]] <- temp.out
}

## Niche first
# Niche.Vector
iterations <- names(marks$niche$vector)
for(i in 1:length(iterations)){
  temp <- marks$niche$vector[[iterations[i]]]
  mechanisms.to.annotate <- temp$gene
  annotations.to.add <- fantom5.rat.annotations[mechanisms.to.annotate,]
  temp.out <- cbind(temp,annotations.to.add)
  marks$niche$vector[[iterations[i]]] <- temp.out
}
# Niche.condition
iterations <- names(marks$niche$condition)
for(i in 1:length(iterations)){
  temp <- marks$niche$condition[[iterations[i]]]
  mechanisms.to.annotate <- temp$gene
  annotations.to.add <- fantom5.rat.annotations[mechanisms.to.annotate,]
  temp.out <- cbind(temp,annotations.to.add)
  marks$niche$condition[[iterations[i]]] <- temp.out
}
# Niche.Influence
iterations <- names(marks$niche$influence)
for(i in 1:length(iterations)){
  temp <- marks$niche$influence[[iterations[i]]]
  mechanisms.to.annotate <- temp$gene
  annotations.to.add <- fantom5.rat.annotations[mechanisms.to.annotate,]
  temp.out <- cbind(temp,annotations.to.add)
  marks$niche$influence[[iterations[i]]] <- temp.out
}
# Niche.archetype
iterations <- names(marks$niche$archetype)
for(i in 1:length(iterations)){
  temp <- marks$niche$archetype[[iterations[i]]]
  mechanisms.to.annotate <- temp$gene
  annotations.to.add <- fantom5.rat.annotations[mechanisms.to.annotate,]
  temp.out <- cbind(temp,annotations.to.add)
  marks$niche$archetype[[iterations[i]]] <- temp.out
}

## Influence first
# Influence.Vector
iterations <- names(marks$influence$vector)
for(i in 1:length(iterations)){
  temp <- marks$influence$vector[[iterations[i]]]
  mechanisms.to.annotate <- temp$gene
  annotations.to.add <- fantom5.rat.annotations[mechanisms.to.annotate,]
  temp.out <- cbind(temp,annotations.to.add)
  marks$influence$vector[[iterations[i]]] <- temp.out
}
# Influence.Niche
iterations <- names(marks$influence$niche)
for(i in 1:length(iterations)){
  temp <- marks$influence$niche[[iterations[i]]]
  mechanisms.to.annotate <- temp$gene
  annotations.to.add <- fantom5.rat.annotations[mechanisms.to.annotate,]
  temp.out <- cbind(temp,annotations.to.add)
  marks$influence$niche[[iterations[i]]] <- temp.out
}
# Influence.condition
iterations <- names(marks$influence$condition)
for(i in 1:length(iterations)){
  temp <- marks$influence$condition[[iterations[i]]]
  mechanisms.to.annotate <- temp$gene
  annotations.to.add <- fantom5.rat.annotations[mechanisms.to.annotate,]
  temp.out <- cbind(temp,annotations.to.add)
  marks$influence$condition[[iterations[i]]] <- temp.out
}
# Influence.archetype
iterations <- names(marks$influence$archetype)
for(i in 1:length(iterations)){
  temp <- marks$influence$archetype[[iterations[i]]]
  mechanisms.to.annotate <- temp$gene
  annotations.to.add <- fantom5.rat.annotations[mechanisms.to.annotate,]
  temp.out <- cbind(temp,annotations.to.add)
  marks$influence$archetype[[iterations[i]]] <- temp.out
}

## Vector first
# Vector.influence
iterations <- names(marks$vector$influence)
for(i in 1:length(iterations)){
  temp <- marks$vector$influence[[iterations[i]]]
  mechanisms.to.annotate <- temp$gene
  annotations.to.add <- fantom5.rat.annotations[mechanisms.to.annotate,]
  temp.out <- cbind(temp,annotations.to.add)
  marks$vector$influence[[iterations[i]]] <- temp.out
}
# Vector.Niche
iterations <- names(marks$vector$niche)
for(i in 1:length(iterations)){
  temp <- marks$vector$niche[[iterations[i]]]
  mechanisms.to.annotate <- temp$gene
  annotations.to.add <- fantom5.rat.annotations[mechanisms.to.annotate,]
  temp.out <- cbind(temp,annotations.to.add)
  marks$vector$niche[[iterations[i]]] <- temp.out
}
# Vector.condition
iterations <- names(marks$vector$condition)
for(i in 1:length(iterations)){
  temp <- marks$vector$condition[[iterations[i]]]
  mechanisms.to.annotate <- temp$gene
  annotations.to.add <- fantom5.rat.annotations[mechanisms.to.annotate,]
  temp.out <- cbind(temp,annotations.to.add)
  marks$vector$condition[[iterations[i]]] <- temp.out
}
# Vector.archetype
iterations <- names(marks$vector$archetype)
for(i in 1:length(iterations)){
  temp <- marks$vector$archetype[[iterations[i]]]
  mechanisms.to.annotate <- temp$gene
  annotations.to.add <- fantom5.rat.annotations[mechanisms.to.annotate,]
  temp.out <- cbind(temp,annotations.to.add)
  marks$vector$archetype[[iterations[i]]] <- temp.out
}

save(marks,file = 'marks.annotated.2024-03-29.Robj')
save(fantom5.rat.annotations,file = 'fantom5.rat.annotations.Robj')
