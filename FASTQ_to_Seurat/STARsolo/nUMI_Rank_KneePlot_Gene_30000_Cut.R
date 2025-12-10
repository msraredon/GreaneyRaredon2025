require(ggplot2)
Barcode_rank_plot <- function(path, nBarcodes, xmin, xmax){
  # Identify sample name
  sample.name <- basename(dirname(dirname(path)))
  # Load UMIperCellSorted
  UMIperCellSorted <- read.table(paste(path,"/UMIperCellSorted.txt",sep = ''), quote="\"", comment.char="")
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
    scale_y_continuous(name = "Droplet size", trans = "log10",
                       breaks=c(1, 10, 100, 150, 200, 500, 1000, 5000, 10000, 25000, 50000),
                       limits=c(1, 500000)) + 
    geom_vline(xintercept=nBarcodes, linetype="dashed", color="red")  +     
    geom_hline(yintercept=cutoff, linetype="dashed", color="blue")  +     
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle(sample.name)+
    annotate("text", label = paste("Rank Threshold =",nBarcodes), x = xmin, y = 4,hjust=0)+
    annotate("text", label = paste("nUMI Threshold =",cutoff), x = xmin, y = 2,hjust=0)
  
  return(p)
}



##### native v2 ####
path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/10x_v2/sample_out/mRat/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 30000 # 5500
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/10x_v2/sample_out/fRat/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 30000 # 6750
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

##### native v3 ####
path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/10x_v3/sample_out/P0-m/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 30000 # 7000
xmin <- 1
xmax <- 300000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Pneumonectomy/10x_v3/sample_out/P0-f/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 8000
xmin <- 1
xmax <- 300000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()




##### BEFM #####
path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM/sample_out/BEFM4/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 2100
xmin <- 1
xmax <- 300000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM/sample_out/BEFM5/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 7000
xmin <- 1
xmax <- 300000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM/sample_out/BEFM6/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 5000
xmin <- 1
xmax <- 300000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM/sample_out/TXP3-L/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 3800
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM/sample_out/TXP3-R/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 7300
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM/sample_out/TXP4-L/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 6000
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM/sample_out/TXP4-R/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 7300
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()


##### BEFM_addenda v2 ####
path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v2/sample_out/BC1P3/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 17000
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v2/sample_out/BC1P6/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 3000
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v2/sample_out/BCL5/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 1000
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v2/sample_out/BEF1/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 3100
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v2/sample_out/BEF2/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 1000
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

##### BEFM_addenda v3 ####
path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v3/sample_out/BCEC2/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 3000
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v3/sample_out/BEF3/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 1000
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v3/sample_out/BEF12/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 9400
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v3/sample_out/BEF14/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 8000
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v3/sample_out/BEF15/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 8000
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v3/sample_out/BEF16/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 8000
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v3/sample_out/BEFM1/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 3000
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BEFM_addenda/10x_v3/sample_out/BEFM2/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 7500
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()



##### BCL6 ####
path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/BCL6/sample_out/BCL6/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 1000
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()
##### Fibroblasts ####
path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Fibroblasts/sample_out/FB13/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 1000
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Fibroblasts/sample_out/FB14/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 1000
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()



path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Organoids/sample_out/MacAlv/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 1000
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()

path <- '/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/Organoids/sample_out/RLMVEC/Solo.out/Gene'
setwd(path)
nBarcodes <- 30000 # 1000
xmin <- 1
xmax <- 100000
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
png(paste("BarcodeKneePlot", basename(dirname(dirname(path))), nBarcodes,"png", sep="."), res=200, unit="in", height=6, width=10)
Barcode_rank_plot(path, nBarcodes, xmin, xmax)
dev.off()


