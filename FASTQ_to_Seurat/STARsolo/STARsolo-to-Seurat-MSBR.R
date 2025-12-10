# Set directories
parent.dir <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM'
sample.names <- list.files(file.path(parent.dir, "sample_out"))

# Load Packages
require(Matrix)

# Function to read in the velocyto STAR output and save it as seperate mtxs
split.STAR.velocyto.mtx <- function(sample.names){
  STARsolo.sample.dir.path <- file.path(parent.dir, "sample_out", sample.names, "Solo.out")
  # read the velocyto raw matrix
  options(readr.show_progress = FALSE)
  velocyto.mtx <- readr::read_delim(file.path(STARsolo.sample.dir.path,"Velocyto", "raw", "matrix.mtx"),
                                    delim = ' ', skip = 3, col_names = FALSE)
  number_BC <- length(count.fields(file.path(STARsolo.sample.dir.path,"Velocyto", "raw", "barcodes.tsv")))
  number_Features <- length(count.fields(file.path(STARsolo.sample.dir.path,"Velocyto", "raw", "features.tsv")))
  spliced.mtx <- Matrix::sparseMatrix(i = velocyto.mtx$X1,j = velocyto.mtx$X2,x = velocyto.mtx$X3, dims = c(number_Features, number_BC))
  Matrix::writeMM(spliced.mtx, file.path(STARsolo.sample.dir.path,"Velocyto", "raw","spliced.mtx"))
  unspliced.mtx <- Matrix::sparseMatrix(i = velocyto.mtx$X1,j = velocyto.mtx$X2,x = velocyto.mtx$X4,dims = c(number_Features, number_BC))
  Matrix::writeMM(unspliced.mtx, file.path(STARsolo.sample.dir.path,"Velocyto", "raw","unspliced.mtx"))
  cat("Completed split of STAR velocyto output for sample ", sample.names,".\n",sep="")
}
# Run
ncores = 4
parallel::mclapply(sample.names, split.STAR.velocyto.mtx, mc.cores=ncores)


# Determine thresholds for percent mito, fraction intronic and nUMI (Taylor can show you), 


# i extract the "real" cells with this: 
# (select the pattern how mito genes are names based on your species)


min.nUMI <- 1000
min.fraction.intronic <- 0.1
max.fraction.mito <- 0.15
# pattern.mito.genes <- "^MT-"
# pattern.mito.genes <- "^mt-"
pattern.mito.genes <- "^Mt-" # Rat and Mouse


# input of this function is the sample name in the parent.dir/sample_out directory
get.cell.gene.matrix <- function(sample.names){
  STARsolo.sample.dir.path <- file.path(parent.dir, "sample_out", sample.names, "Solo.out")
  # read the GeneFull raw matrix and add the row (gene) and column (cell barcode) names
  GeneFull.mtx <- Matrix::readMM(file.path(STARsolo.sample.dir.path,"GeneFull", "raw", "matrix.mtx" ))
  # barcodes are the same for raw Genefull, raw Gene and Velocyto
  colnames(GeneFull.mtx) <- read.table(file.path(STARsolo.sample.dir.path,"GeneFull", "raw", "barcodes.tsv" ), sep = '\t', header = FALSE)[,1]
  rownames(GeneFull.mtx) <- read.delim(file.path(STARsolo.sample.dir.path,"GeneFull", "raw", "features.tsv" ), sep = '\t', header = FALSE)[,2]
  # read the velocyto raw spliced and unspliced matrix and add the row (gene) and column (cell barcode) names
  unspliced.mtx <- Matrix::readMM(file.path(STARsolo.sample.dir.path,"Velocyto", "raw","unspliced.mtx"))
  spliced.mtx <- Matrix::readMM(file.path(STARsolo.sample.dir.path,"Velocyto", "raw","spliced.mtx"))
  colnames(unspliced.mtx) <- read.table(file.path(STARsolo.sample.dir.path,"Velocyto", "raw", "barcodes.tsv" ), sep = '\t', header = FALSE)[,1]
  rownames(unspliced.mtx) <- read.delim(file.path(STARsolo.sample.dir.path,"Velocyto", "raw", "features.tsv" ), sep = '\t', header = FALSE)[,2]
  colnames(spliced.mtx) <- read.table(file.path(STARsolo.sample.dir.path,"Velocyto", "raw", "barcodes.tsv" ), sep = '\t', header = FALSE)[,1]
  rownames(spliced.mtx) <- read.delim(file.path(STARsolo.sample.dir.path,"Velocyto", "raw", "features.tsv" ), sep = '\t', header = FALSE)[,2]
  # remove all barcodes with less then specified nUMI
  GeneFull.mtx <- GeneFull.mtx[, Matrix::colSums(GeneFull.mtx)> min.nUMI]
  # remove all barcodes not in the filtered GeneFull.mtx
  unspliced.mtx <- unspliced.mtx[, colnames(GeneFull.mtx)]
  spliced.mtx <- spliced.mtx[, colnames(GeneFull.mtx)]
  # get fraction of intronic and mitochondrial genes
  fraction.intronic <- Matrix::colSums(unspliced.mtx)/(Matrix::colSums(unspliced.mtx)+Matrix::colSums(spliced.mtx))
  mito.genes <- grep(pattern = pattern.mito.genes, x = rownames(GeneFull.mtx), value = TRUE)
  fraction.mito <- Matrix::colSums(GeneFull.mtx[mito.genes,])/Matrix::colSums(GeneFull.mtx)
  # filter based on fraction.mito and fraction.intronic
  fraction.mito.keep.index <- names(fraction.mito)[fraction.mito <= max.fraction.mito]
  fraction.intronic.keep.index <- names(fraction.intronic)[fraction.intronic >= min.fraction.intronic]
  GeneFull.mtx <- GeneFull.mtx[, colnames(GeneFull.mtx) %in% intersect(fraction.mito.keep.index, fraction.intronic.keep.index)]
  # add sample name to barcodes, seperated by "__"
  colnames(GeneFull.mtx) <- paste(sample.names, colnames(GeneFull.mtx), sep="__")
  cat("Completed gene cell matrix extraction for sample ", sample.names,".\n",sep="")
  return(GeneFull.mtx)
}

sample.output <- parallel::mclapply(sample.names, get.cell.gene.matrix, mc.cores=ncores)
output.mtx <- Reduce(cbind, sample.output)
dim(output.mtx)
# remove all genes with a rowSum of 0
output.mtx <- output.mtx[! Matrix::rowSums(output.mtx) == 0,]
dim(output.mtx)
save(output.mtx, file = 'pneumonectomy.GeneFull.Robj')



# this is how I extract the velocyte stuff

velo.spliced.output.file.path <- "pneumonectomy.spliced.Robj"
velo.unspliced.output.file.path <- "pneumonectomy.unspliced.Robj"


barcodes.ids <- colnames(output.mtx)
sample.names <- list.files(file.path(parent.dir, "sample_out"))

get.velocyto.spliced.matrices <- function(sample.names){
  STARsolo.sample.dir.path <- file.path(parent.dir, "sample_out", sample.names, "Solo.out")
  barcodes.ids.df <- as.data.frame(stringr::str_split(barcodes.ids, "__", simplify = TRUE))
  colnames(barcodes.ids.df) <- c("sample.name", "barcode.id")
  # subset barcode df to barcode of that sample
  barcodes.ids.df <- barcodes.ids.df[barcodes.ids.df$sample.name==sample.names,]
  spliced.mtx <- Matrix::readMM(file.path(STARsolo.sample.dir.path,"Velocyto", "raw","spliced.mtx"))
  # filtered based on the valid barcodes
  colnames(spliced.mtx) <- read.table(file.path(STARsolo.sample.dir.path,"GeneFull", "raw", "barcodes.tsv" ), sep = '\t', header = FALSE)[,1]
  rownames(spliced.mtx) <- read.delim(file.path(STARsolo.sample.dir.path,"GeneFull", "raw", "features.tsv" ), sep = '\t', header = FALSE)[,2]
  spliced.mtx <- spliced.mtx[, colnames(spliced.mtx) %in% barcodes.ids.df$barcode.id]
  colnames(spliced.mtx) <- paste(sample.names, colnames(spliced.mtx), sep="__")
  cat("Completed gene-cell extraction for sample ", sample.names,".\n",sep="")
  return(spliced.mtx)
}

get.velocyto.unspliced.matrices <- function(sample.names){
  STARsolo.sample.dir.path <- file.path(parent.dir, "sample_out", sample.names, "Solo.out")
  barcodes.ids.df <- as.data.frame(stringr::str_split(barcodes.ids, "__", simplify = TRUE))
  colnames(barcodes.ids.df) <- c("sample.name", "barcode.id")
  # subset barcode df to barcode of that sample
  barcodes.ids.df <- barcodes.ids.df[barcodes.ids.df$sample.name==sample.names,]
  unspliced.mtx <- Matrix::readMM(file.path(STARsolo.sample.dir.path,"Velocyto", "raw","unspliced.mtx"))
  # filtered based on the valid barcodes
  colnames(unspliced.mtx) <- read.table(file.path(STARsolo.sample.dir.path,"GeneFull", "raw", "barcodes.tsv" ), sep = '\t', header = FALSE)[,1]
  rownames(unspliced.mtx) <- read.delim(file.path(STARsolo.sample.dir.path,"GeneFull", "raw", "features.tsv" ), sep = '\t', header = FALSE)[,2]
  unspliced.mtx <- unspliced.mtx[, colnames(unspliced.mtx) %in% barcodes.ids.df$barcode.id]
  colnames(unspliced.mtx) <- paste(sample.names, colnames(unspliced.mtx), sep="__")
  cat("Completed unspliced velocyto extraction for sample ", sample.names,".\n",sep="")
  return(unspliced.mtx)
}


spliced.sample.output <- parallel::mclapply(sample.names, get.velocyto.spliced.matrices, mc.cores=ncores)
spliced.output.mtx <- Reduce(cbind, spliced.sample.output)
dim(spliced.output.mtx)

unspliced.sample.output <- parallel::mclapply(sample.names, get.velocyto.unspliced.matrices, mc.cores=ncores)
unspliced.output.mtx <- Reduce(cbind, unspliced.sample.output)
dim(unspliced.output.mtx)

# remove all rows which are not in the output.mtx (there they were zero sums)
spliced.output.mtx <- spliced.output.mtx[rownames(output.mtx),]
unspliced.output.mtx <- unspliced.output.mtx[rownames(output.mtx),]

save(spliced.output.mtx, file = velo.spliced.output.file.path)
save(unspliced.output.mtx, file = velo.unspliced.output.file.path)



