
# function to read in the velocyto STAR output and save it as seperate mtxs
split.STAR.velocyto.mtx <- function(sample.names){
  STARsolo.sample.dir.path <- file.path(parent.dir, "sample_out", sample.names, "Solo.out")
  # read the velocyto raw matrix
  options(readr.show_progress = FALSE)
  velocyto.mtx <- readr::read_delim(file.path(STARsolo.sample.dir.path,"Velocyto", "raw", "matrix.mtx"),
                                    delim = ' ', skip = 3, col_names = FALSE, show_col_types=FALSE)
  number_BC <- length(count.fields(file.path(STARsolo.sample.dir.path,"Velocyto", "raw", "barcodes.tsv")))
  number_Features <- length(count.fields(file.path(STARsolo.sample.dir.path,"Velocyto", "raw", "features.tsv")))
  spliced.mtx <- Matrix::sparseMatrix(i = velocyto.mtx$X1,j = velocyto.mtx$X2,x = velocyto.mtx$X3, dims = c(number_Features, number_BC))
  Matrix::writeMM(spliced.mtx, file.path(STARsolo.sample.dir.path,"Velocyto", "raw","spliced.mtx"))
  unspliced.mtx <- Matrix::sparseMatrix(i = velocyto.mtx$X1,j = velocyto.mtx$X2,x = velocyto.mtx$X4,dims = c(number_Features, number_BC))
  Matrix::writeMM(unspliced.mtx, file.path(STARsolo.sample.dir.path,"Velocyto", "raw","unspliced.mtx"))
  cat("Completed split of STAR velocyto output for sample ", sample.names,".\n",sep="")
}
#### BEFM ####
parent.dir <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM'
sample.names <- list.files(file.path(parent.dir, "sample_out"))
ncores = 4
parallel::mclapply(sample.names, split.STAR.velocyto.mtx, mc.cores=ncores)
#### BEFM addenda v2 ####
parent.dir <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v2'
sample.names <- list.files(file.path(parent.dir, "sample_out"))
sample.names
ncores = 4
parallel::mclapply(sample.names, split.STAR.velocyto.mtx, mc.cores=ncores)
#### BEFM addenda v3 ####
parent.dir <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v3'
sample.names <- list.files(file.path(parent.dir, "sample_out"))
sample.names <- sample.names[sample.names!='BEFM2'] # Work around to deal with faulty Velocyto run (Bus error,need to redo)
sample.names
ncores = 4
parallel::mclapply(sample.names, split.STAR.velocyto.mtx, mc.cores=ncores)
#### Pneumonectomy v2 ####
parent.dir <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/10x_v2'
sample.names <- list.files(file.path(parent.dir, "sample_out"))
sample.names
ncores = 4
parallel::mclapply(sample.names, split.STAR.velocyto.mtx, mc.cores=ncores)
#### Pneumonectomy v3 ####
parent.dir <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/10x_v3'
sample.names <- list.files(file.path(parent.dir, "sample_out"))
sample.names
ncores = 4
parallel::mclapply(sample.names, split.STAR.velocyto.mtx, mc.cores=ncores)
#### Pneumonectomy DropSeq ####
parent.dir <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/DropSeq'
sample.names <- list.files(file.path(parent.dir, "sample_out"))
sample.names
ncores = 4
parallel::mclapply(sample.names, split.STAR.velocyto.mtx, mc.cores=ncores)
#### Pneumonectomy Addenda ####
parent.dir <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy_addenda/10x_v2'
sample.names <- list.files(file.path(parent.dir, "sample_out"))
sample.names
ncores = 4
parallel::mclapply(sample.names, split.STAR.velocyto.mtx, mc.cores=ncores)
#### BEFM2 second run ####
parent.dir <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM2'
sample.names <- list.files(file.path(parent.dir, "sample_out"))
sample.names
ncores = 4
parallel::mclapply(sample.names, split.STAR.velocyto.mtx, mc.cores=ncores)
#### BCL6 ####
parent.dir <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BCL6'
sample.names <- list.files(file.path(parent.dir, "sample_out"))
sample.names
ncores = 4
parallel::mclapply(sample.names, split.STAR.velocyto.mtx, mc.cores=ncores)
#### Fibroblasts ####
parent.dir <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Fibroblasts'
sample.names <- list.files(file.path(parent.dir, "sample_out"))
sample.names
ncores = 4
parallel::mclapply(sample.names, split.STAR.velocyto.mtx, mc.cores=ncores)
#### Organoids ####
parent.dir <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Organoids'
sample.names <- list.files(file.path(parent.dir, "sample_out"))
sample.names
ncores = 4
parallel::mclapply(sample.names, split.STAR.velocyto.mtx, mc.cores=ncores)