#!/bin/sh

#  create_download_slurm.sh
#  
#
#  Created by Josh Kemp on 4/8/23.
#  

splitn=1 #lines to run per sbatch
acclist=SRR_Acc_list.txt


#download reference genome
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff.gz




##create function
cat << 'EOF' > download_trim_and_map.sh
#!/bin/bash
#variables
genome="GCF_000001735.4_TAIR10.1_genomic.fna"
fastq_1=${1}_1.fastq
fastq_2=${1}_2.fastq
trimmed_1=${1}_1.trimmed_fastq.gz
trimmed_2=${1}_2.trimmed_fastq.gz
out=${1}.bam
#named pipes
mkfifo ${fastq_1}
mkfifo ${fastq_2}
mkfifo ${trimmed_1}
mkfifo ${trimmed_2}
#function (requires genome index and loaded modules before running)
prefetch $1
fasterq-dump -f --threads 16 --split-3 $1 &&
bbduk.sh -Xmx100g in1=${fastq_1} in2=${fastq_2} out1=${trimmed_1} out2=${trimmed_2} ref=adapters threads=8 minlen=25 qtrim=rl trimq=10 k=25 ktrim=r hdist=1 mink=6 &&
bwa mem -M -t 16 $genome ${trimmed_1} ${trimmed_2} | samtools view -buS - > ${1}.bam &&



rm $fastq_1
rm $fastq_2
rm $trimmed_1
rm $trimmed_2

exit
EOF

###################################################################################

cat << 'EOF' > slurmheader.txt
#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --time=0-16:00:00  # max job runtime
#SBATCH --ntasks-per-node=11
#SBATCH --cpus-per-task=16  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --partition=biocrunch  # partition(s)
#SBATCH --mem=1200G  # max memory
#SBATCH -J "sorg_SAP_wholegenomes_sradownload"  # job name

EOF

cat << 'EOF' > slurmheader_index.txt
#!/bin/bash
#Submit this script with: sbatch thefilename

#SBATCH --time=0-00:15:00  # max job runtime
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --partition=biocrunch  # partition(s)
#SBATCH --mem=600G  # max memory
#SBATCH -J "build index"  # job name
EOF




cat << 'EOF' > loadmodules.txt
module load sratoolkit
module load bbmap
module load samtools
module load bwa
module load parallel
EOF

####################################################################################





## create n files dividing up the accesions into managagable chunks for each slurm batch
split -n l/${splitn} $acclist piece_acc_

## create a batch script for each file chunk
for Filechunk in piece_acc_*
    do
    cat slurmheader.txt > download_trim_and_map${Filechunk}.batch
    cat loadmodules.txt >> download_trim_and_map${Filechunk}.batch
    echo "parallel -a ${Filechunk} bash download_trim_and_map.sh" >> download_trim_and_map${Filechunk}.batch
    echo "wait" >> download_trim_and_map${Filechunk}.batch
    echo "exit" >> download_trim_and_map${Filechunk}.batch
done

## create batch script to download and create the reference genome index required for mapping
cat slurmheader_index.txt > index_genome.batch
echo "module load bwa" >> index_genome.batch
echo "module load gzip" >> index_genome.batch
echo "wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz" >> index_genome.batch &&
echo "wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff.gz" >> index_genome.batch &&
echo "genome=*GCF_*.fna.gz" >> index_genome.batch
echo "gzip -d \$genome" >> index_genome.batch
echo "genome=*GCF_*.fna" >> index_genome.batch
echo "bwa index -a bwtsw \$genome" >> index_genome.batch

#run slurm scripts. first index, and then the others once it finishes
genomeindex_jobid=$(sbatch --parsable index_genome.batch)
sbatch --dependency=afterok:$genomeindex_jobid download_trim_and_map*.batch
