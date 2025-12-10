## argument 1: path to parent directory for the project
## eg: Rscript /home/fas/kaminski/tsa27/scratch60/10x_zUMI/Rscripts/10X_zUMIs_2.0.R /home/fas/kaminski/tsa27/scratch60/10x_zUMI
#####################################################################################################################################################
########## This Rscript is designed to use as input: the same directory labeled as "output.parent.dir" in the cutadapt Rscript
########## which contains a directory titles "sample_out" and "cutadapt_scripts"
########## The cutadapt job must be run beforehand to generate "trimmed.R1.fastq" and "trimmed.R2.fastq" in the
########## "trimmed_merged_fastq" subdirectory in each sample's output directory
##########
########## The purpose of the script is to generate a sample-specific .sh rocessing each sample through STARsolo
##########
###### Summary:
###### 1 Identify all sample folders
#loop{ 2 generate script for STAR solo} 
#####################################################################################################################################################
# 8 threads
# 60GB memory
# for V2 and V3 you need to change the length of the UMI and the whitelist

args <- commandArgs(trailingOnly=TRUE)
# args[1] is input parent directory

sample.parent.dir <- file.path(args[1], "sample_out")
scripts.parent.dir <- file.path(args[1], "scripts")
STARsolo.scripts.dir <- file.path(scripts.parent.dir, "STARsolo")
if(file.exists(STARsolo.scripts.dir)==F){
  dir.create(STARsolo.scripts.dir)
}

soup.sample.names <- list.files(file.path(sample.parent.dir))
soup.sample.dir.paths <- file.path(sample.parent.dir, soup.sample.names)
jobsub.filepath <- file.path(STARsolo.scripts.dir,"STARsolo.jobsub.bat")

STAR.path <- "/gpfs/gibbs/pi/kaminski/public/softwares/STAR-2.7.3a/bin/Linux_x86_64_static/STAR"
genome.dir <- "/gpfs/gibbs/pi/kaminski/public/Backup/Sam/genome_index_rat/Ensembl_Rnor_6.0.99/STAR_Index_noGTF"
GTF.path <-  "/gpfs/gibbs/pi/kaminski/public/Backup/Sam/genome_index_rat/Ensembl_Rnor_6.0.99/GTF/Rattus_norvegicus.Rnor_6.0.99.gtf"


CBwhitelist.path <- "/gpfs/gibbs/pi/kaminski/public/Backup/Jonas/softwares/10x_Whitelists/3M-february-2018.txt" # for 10x v3
# CBwhitelist.path <- "/gpfs/gibbs/pi/kaminski/public/Backup/Jonas/softwares/10x_Whitelists/737K-august-2016.txt" # for 10x V2
# CBwhitelist.path <- "None" # for DropSeq

# for gff files add the "--sjdbGTFtagExonParentTranscript Parent" option
STAR.options <- paste("--soloType CB_UMI_Simple --readFilesCommand zcat --soloBarcodeReadLength 0 ", 
                      "--soloFeatures Gene GeneFull Velocyto --runThreadN 8 ", 
                      "--soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --twopassMode Basic --soloStrand Forward", sep="") 
# --soloUMIlen 10 for V2, --soloUMIlen 12 for V3, --soloUMIlen 8 for DropSeq
# --soloCBlen 16 --soloUMIstart 17 for 10x
# --soloCBlen 12 --soloUMIstart 13 for DropSeq
# --soloStrand Forward for 3prime, --soloStrand Reverse for 5prime
# add this for 5prime after removel of 50bps adaptor contamination: --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3


## Wipe the old jobsub.bat file
cat("",sep="",file=jobsub.filepath, append=FALSE)
## Loop over samples to extract trimmed R1 and R2 1:length(soup.sample.dir.paths)
for(i in 1:length(soup.sample.names)){
  trimmed.folder.paths <- file.path(soup.sample.dir.paths[i], "trimmed_merged_fastq")
  #### Distinguish Read1 and Read2
  read1.fastq.path <- file.path(trimmed.folder.paths , "trimmed.R1.fastq.gz")
  read2.fastq.path <- file.path(trimmed.folder.paths , "trimmed.R2.fastq.gz")
  script.filepath <- file.path(STARsolo.scripts.dir, paste(soup.sample.names[i], "_STARsolo.sh", sep = ""))
  ##### print the script ##  -N ",num.reads.per.cell,";
  cmd.out <- NULL
  cmd.out <- paste("#!/bin/bash\n")
  cmd.out <- paste(cmd.out,"#SBATCH --time=6:00:00 -p day,pi_kaminski,bigmem --ntasks=1 --cpus-per-task=8",
                   " --mem=60000M --job-name=", paste0(soup.sample.names[i], "_STARsolo"),
                   " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n",sep="")
  cmd.out <- paste(cmd.out, STAR.path, " ", STAR.options, " --soloCBwhitelist ", CBwhitelist.path, 
                   " --genomeDir ", genome.dir, " --sjdbGTFfile ", GTF.path,   
                   " --readFilesIn ", read2.fastq.path, " ", read1.fastq.path, 
                   " --outFileNamePrefix ", soup.sample.dir.paths[i], "/", sep="")
  cat(cmd.out,file=script.filepath,append=F)
  cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
}
system(paste("chmod 700 ",jobsub.filepath,sep=""))

