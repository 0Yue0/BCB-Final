#!/bin/sh

#  sc_rna_analysis.sh
#  For BCB 546
#
#  Created by Josh Kemp on 4/20/23.


acclist=SRR_Acc_list.txt  #list of SRA accesions to download and map

reference_genome_url="https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
reference_annotation_url="https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.56.gtf.gz"
ERCC_spike_sequences_url="https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip"

splitn=2 #  number of sbatchs to split accession file into.  Makes running faster if there are many accessions in the list

mkdir sc_rnaseq_analysis #create a new directory for the analysis of the SRA dataset being downloaded and then move into it
cd sc_rnaseq_analysis


#Create Sub-directories used for this analysis
mkdir prefetch_tmp
mkdir downloaded_reads
mkdir pretrim_quality
mkdir filtered_reads
mkdir trimmed_reads
mkdir posttrim_quality
mkdir star_index
mkdir alignments
mkdir sorted_alignments
mkdir star_output


#############################

##create bash script to download and trim reads; needs index first, may need to check options settings

cat << 'EOF' > download_trim_and_map.sh  #this is a readdoc notation.  It prints everything until EOF as written to the "download_trim_and_map.sh file.  In this case I am using it to create a bash script that will take one SRA accession as an input.  $1 references the first input argument in a bash script
#!/bin/bash
set -e
set -u
set
#variables
basepath="$(pwd)"  # finds the current directory so that the program can use complete path names to specify files.  file name variables can be updated to use relative paths if needed.
genome="GCF_*.fna"  #Make sure this matches after indexing
threads=16
memory=100  #memory in gb bbduk is allowed to use
prefetch_pipe=./prefetch_tmp/${1}
fastq=${basepath}/downloaded_reads/${1}.fastq
trimmed=${basepath}/trimmed_reads/${1}.trimmed_fastq.gz
out=${basepath}/alignments/${1}.bam

prefetch $1 --max-size 80g -O prefetch_tmp  # downloads a compressed version of the raw reads.
fasterq-dump -f --threads ${threads} --split-3 --outdir downloaded_reads ${prefetch_pipe} &&  #this will decompress the reads into 3 files (split 3 option) pair 1 pair 2 and unpaired reads.
fastqc -t ${threads} -o pretrim_quality $fastq & # will show the quality data for the reads pre trimming
wait
STAR \
--outFilterMismatchNoverReadLmax .02 \
--genomeDir star_index \
--runThreadN ${threads} \
--readFilesIn $fastq \
--quantMode GeneCounts \
--outFileNamePrefix ${basepath}/star_output/${1}.genecounts \


STAR --genomeLoad Remove
wait
rm -r prefetch_tmp

    exit
EOF


###########################
#create a bash script that will take previous alignments and give us gene counts TPM using stringtie

cat << 'EOF' > alignment_readcounts.sh
#!/bin/bash
mkdir readcounts
#variables
basepath=$(pwd)  # finds the current directory so that the program can use complete path names to specify files.  file name variables can be updated to use relative paths if needed.
genome="*.fa"  #Make sure this matches after indexing
annotation="*.gtf"
threads=16
trimmed_1=${basepath}/trimmed_reads/${1}_1.trimmed_fastq.gz
trimmed_2=${basepath}/trimmed_reads/${1}_2.trimmed_fastq.gz
out=${basepath}/alignments/${1}.bam
input_sam=${basepath}/alignments/${1}.sam
sorted_sam=${basepath}/sorted_alignments/sorted${1}.sam
annotation=${basepath}/star_index/*_*.*.gt

samtools sort -o sosrted_$sorted_bam $input_sam
featureCounts -T 4 -s 2 -p -t gene -g ID -a ${annotation} -o readcounts/${1}.gene.txt
wait
exit

EOF


### Creates the Slurm header needed to run the download and mapping on HPC cluster. This will likely need to be modified.  For ISU nova users the following is a good resource to update this to reflect your specific hardware:https://www.hpc.iastate.edu/guides/nova/slurm-script-generator-for-nova

cat << 'EOF' > align_slurmheader.txt
#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --time=0-01:30:00  # max job runtime
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=16  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --partition=whatever  # partition(s)
#SBATCH --mem=200G  # max memory
#SBATCH -J "LCM_rnaSeq align"  # job name

EOF

cat << 'EOF' > readcount_slurmheader.txt
#!/bin/bash

#SBATCH --time=0-01:30:00  # max job runtime
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=16  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --partition=whatever  # partition(s)
#SBATCH --mem=200G  # max memory
#SBATCH -J "LCM_readcount"  # job name

EOF



### Creates the Slurm header needed to run on HPC cluster, this one is more specific to the slurm batch for building reference genome index, which should not take as long and is somewhat memory/disk space intensive.
cat << 'EOF' > slurmheader_index.txt
#!/bin/bash
    
#SBATCH --time=0-00:20:00  # max job runtime
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --partition=biocrunch  # partition(s)
#SBATCH --mem=200G  # max memory
#SBATCH -J "build LCM index"  # job name
EOF



## This creates a file that lists the modules we will need to load in order to run the download through alignment step.  For maximum reproducibilty this should really list the exact spack module name and not the general name.

cat << 'EOF' > align_loadmodules.txt
module load fastqc
module load sratoolkit
module load bbmap
module load samtools
module load star/2.7.6a-wq7xoea
module load parallel
EOF

### modules needed for getting readcounts form alignments
cat << 'EOF' > readcount_loadmodules.txt
module load samtools
module load featureCounts
module load parallel
EOF




## create batch script to download and create the reference genome index required for the alignment step
cat slurmheader_index.txt > index_genome.batch
echo "module load zip" >> index_genome.batch
echo "module load gzip" >> index_genome.batch
echo "module load star/2.7.6a-wq7xoea" >> index_genome.batch
echo "wget ${reference_genome_url}" >> index_genome.batch
echo "wget ${reference_annotation_url}" >> index_genome.batch
echo "wget ${ERCC_spike_sequences_url}" >> index_genome.batch
echo "ERCC=ERCC*.zip" >> index_genome.batch
echo "unzip \$ERCC" >> index_genome.batch
echo "gzip -d *.gz" >> index_genome.batch
echo "genome=*_*.*toplevel.fa" >> index_genome.batch
echo "annotation=*_*.*.gtf" >> index_genome.batch
echo "ERCC=ERCC*" >> index_genome.batch
echo "cat \${ERCC}.fa >> \$genome" >> index_genome.batch
echo "cat \${ERCC}.gtf >> \$annotation" >> index_genome.batch
echo "STAR --runThreadN 24 --runMode genomeGenerate --genomeSAindexNbases 12 --genomeDir star_index --genomeFastaFiles \$genome --sjdbGTFfile \$annotation" >> index_genome.batch




## create n files dividing up the accesions into managagable chunks for each slurm batch
split -n l/${splitn} ../$acclist piece_acc_

## create batch scripts for each file chunk
for Filechunk in piece_acc_*
do
    cat align_slurmheader.txt > download_trim_and_map${Filechunk}.batch
    cat align_loadmodules.txt >> download_trim_and_map${Filechunk}.batch
    echo "parallel -a ${Filechunk} bash download_trim_and_map.sh" >> download_trim_and_map${Filechunk}.batch
    echo "wait" >> download_trim_and_map${Filechunk}.batch
    echo "exit" >> download_trim_and_map${Filechunk}.batch
    
    cat readcount_slurmheader.txt > readcount_${Filechunk}.batch
    cat readcount_loadmodules.txt >> readcount_${Filechunk}.batch
    echo "parallel -a ${Filechunk} bash alignment_readcounts.sh" >> readcount_${Filechunk}.batch
    echo "wait" >> readcount_${Filechunk}.batch
    echo "exit" >> readcount_${Filechunk}.batch
done




#run slurm scripts. This is set up so that the downloading through mapping batches will not run until after the batch that builds the index has successfully finished.  Same applies for the next batch for the counts.

genomeindex_jobid=$(sbatch --parsable index_genome.batch)
for Chunk in piece_acc_*
    do
    rna_align_jobid=$(sbatch --dependency=afterok:$genomeindex_jobid --parsable download_trim_and_map${Chunk}.batch)
done
