#!/bin/bash
#SBATCH -p pi_kaminski,day,bigmem --time=4:00:00 --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=50000M  --job-name=star.rat.index -o /gpfs/gibbs/pi/kaminski/public/Backup/Sam/genome_index_rat/Ensembl_Rnor_6.0.99/star_genomeGenerate.rat.sh.o%J -e /gpfs/gibbs/pi/kaminski/public/Backup/Sam/genome_index_rat/Ensembl_Rnor_6.0.99/star_genomeGenerate.rat.sh.e%J

/gpfs/gibbs/pi/kaminski/public/softwares/STAR-2.7.3a/bin/Linux_x86_64_static/STAR --runMode genomeGenerate \
    --genomeDir /gpfs/gibbs/pi/kaminski/public/Backup/Sam/genome_index_rat/Ensembl_Rnor_6.0.99/STAR_Index_noGTF \
    --genomeFastaFiles /gpfs/gibbs/pi/kaminski/public/Backup/Sam/genome_index_rat/Ensembl_Rnor_6.0.99/FASTA/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
    --runThreadN 8
