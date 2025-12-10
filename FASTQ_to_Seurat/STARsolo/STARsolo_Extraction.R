# Define function to pull out data based on a given threshold in a specific assay
STARsolo_Extraction <- function(path, min.nUMI, threshold.assay = 'Gene'){
  
  # Define sample name based on path (will use later to tag barcodes etc.)
  sample.name <- basename(dirname(path))
  
  # Load in Data
  Gene.mtx<- Matrix::readMM(file.path(path,"Gene", "raw", "matrix.mtx" ))
  GeneFull.mtx <- Matrix::readMM(file.path(path,"GeneFull", "raw", "matrix.mtx" ))
  unspliced.mtx <- Matrix::readMM(file.path(path,"Velocyto", "raw","unspliced.mtx"))
  spliced.mtx <- Matrix::readMM(file.path(path,"Velocyto", "raw","spliced.mtx"))
  
  # Add Barcodes and Rownames
  colnames(Gene.mtx) <- read.table(file.path(path,"Gene", "raw", "barcodes.tsv" ), sep = '\t', header = FALSE)[,1]
  rownames(Gene.mtx) <- read.delim(file.path(path,"Gene", "raw", "features.tsv" ), sep = '\t', header = FALSE)[,2]
  colnames(GeneFull.mtx) <- read.table(file.path(path,"GeneFull", "raw", "barcodes.tsv" ), sep = '\t', header = FALSE)[,1]
  rownames(GeneFull.mtx) <- read.delim(file.path(path,"GeneFull", "raw", "features.tsv" ), sep = '\t', header = FALSE)[,2]
  colnames(unspliced.mtx) <- read.table(file.path(path,"Velocyto", "raw", "barcodes.tsv" ), sep = '\t', header = FALSE)[,1]
  rownames(unspliced.mtx) <- read.delim(file.path(path,"Velocyto", "raw", "features.tsv" ), sep = '\t', header = FALSE)[,2]
  colnames(spliced.mtx) <- read.table(file.path(path,"Velocyto", "raw", "barcodes.tsv" ), sep = '\t', header = FALSE)[,1]
  rownames(spliced.mtx) <- read.delim(file.path(path,"Velocyto", "raw", "features.tsv" ), sep = '\t', header = FALSE)[,2]
  
  # Add sample name to barcodes, seperated by "__"
  colnames(Gene.mtx) <- paste(sample.name, colnames(Gene.mtx), sep="__")
  colnames(GeneFull.mtx) <- paste(sample.name, colnames(GeneFull.mtx), sep="__")
  colnames(unspliced.mtx) <- paste(sample.name, colnames(unspliced.mtx), sep="__")
  colnames(spliced.mtx) <- paste(sample.name, colnames(spliced.mtx), sep="__")
  
  # Remove all barcodes, in all assays, with less then specified nUMI in the assay being used for thresholding
  if (threshold.assay == 'Gene'){
    Gene.mtx <- Gene.mtx[, Matrix::colSums(Gene.mtx)> min.nUMI]
    GeneFull.mtx <- GeneFull.mtx[, colnames(Gene.mtx)]
    unspliced.mtx <- unspliced.mtx[, colnames(Gene.mtx)]
    spliced.mtx <- spliced.mtx[, colnames(Gene.mtx)]
  }
  if (threshold.assay == 'GeneFull'){
    GeneFull.mtx <- GeneFull.mtx[, Matrix::colSums(GeneFull.mtx)> min.nUMI]
    Gene.mtx <- Gene.mtx[, colnames(GeneFull.mtx)]
    unspliced.mtx <- unspliced.mtx[, colnames(GeneFull.mtx)]
    spliced.mtx <- spliced.mtx[, colnames(GeneFull.mtx)]
  }
  
  # Create output object and save
  output <- list(Gene.mtx,GeneFull.mtx,unspliced.mtx,spliced.mtx)
  names(output) <- c('Gene','GeneFull','unspliced','spliced')
  save(output, file = paste(sample.name,threshold.assay,'Extracted.DGEs.min.nUMI',min.nUMI,'Robj',sep='.'))
  #return(output)
  
  # Message
  cat("Completed multi-assay DGE extraction for sample ", sample.name,".\n",sep="")
}

##### Pneumonectomy addenda ####
path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy_addenda/10x_v2/sample_out/TE1/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 3112
STARsolo_Extraction(path, min.nUMI)


##### Pneumonectomy v2 ####
path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/10x_v2/sample_out/mRat/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 341
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/10x_v2/sample_out/fRat/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 252
STARsolo_Extraction(path, min.nUMI)

##### Pneumonectomy v3 ####
path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/10x_v3/sample_out/P0-m/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 899
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/10x_v3/sample_out/P0-f/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 1176
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/10x_v3/sample_out/P3-m/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 1283
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/10x_v3/sample_out/P3-f/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 704
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/10x_v3/sample_out/P7-m/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 1529
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/10x_v3/sample_out/P7-f/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 1334
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/10x_v3/sample_out/P14-m/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 1032
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/10x_v3/sample_out/P14-f/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 1216
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/10x_v3/sample_out/P28-m/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 941
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/10x_v3/sample_out/P56-m/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 1109
STARsolo_Extraction(path, min.nUMI)

##### Pneumonectomy DropSeq ####
path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/DropSeq/sample_out/2-7_rat/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 333
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/DropSeq/sample_out/2-8_rat/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 149
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/DropSeq/sample_out/2-12_rat/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 160
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/DropSeq/sample_out/2-14_rat/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 200
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/DropSeq/sample_out/2-16_rat/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 154
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/DropSeq/sample_out/2-28_rat/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 131
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/DropSeq/sample_out/3-1_rat/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 139
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/DropSeq/sample_out/8-10_rat/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 146
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/DropSeq/sample_out/rat1A/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 267
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/DropSeq/sample_out/rat811/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 217
STARsolo_Extraction(path, min.nUMI)




##### BEFM #####
path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM/sample_out/BEFM4/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 6066
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM/sample_out/BEFM5/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 923
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM/sample_out/BEFM6/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 2624
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM/sample_out/TXP3-L/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 706
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM/sample_out/TXP3-R/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 853
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM/sample_out/TXP4-L/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 526
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM/sample_out/TXP4-R/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 623
STARsolo_Extraction(path, min.nUMI)

##### BEFM_addenda v2 ####
path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v2/sample_out/BC1P3/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 794
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v2/sample_out/BC1P6/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 2299
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v2/sample_out/BCL5/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 3172
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v2/sample_out/BEF1/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 2274
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v2/sample_out/BEF2/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 4062
STARsolo_Extraction(path, min.nUMI)

##### BEFM_addenda v3 ####
path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v3/sample_out/BCEC2/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 9966
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v3/sample_out/BEF3/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 4156
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v3/sample_out/BEF12/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 3191
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v3/sample_out/BEF14/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 4546
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v3/sample_out/BEF15/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 3749
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v3/sample_out/BEF16/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 4833
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v3/sample_out/BEFM1/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 7380
STARsolo_Extraction(path, min.nUMI)

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v3/sample_out/BEFM2/Solo.out'
setwd(path)
threshold.assay <- 'Gene'
min.nUMI <- 7167
STARsolo_Extraction(path, min.nUMI)

