#!/bin/sh

#  LCM_dataset_alignment.sh
#  For BCB 546
#
#  Created by Josh Kemp on 4/15/23.

#build index with hisat2
acclist=LCM_SRA_ACC_list.txt  #list of SRA accesions to download and map
reference_genome_url="https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
reference_annotation_url="https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.56.gtf.gz"
https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/plants/gtf/arabidopsis_thaliana
https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/plants/fasta/arabidopsis_thaliana/dna/

#Note: If using the following NCBI refences and gtf annotation, there will be issues downstream for stringtie I have no clue if the original paper used a masked or unmasked version of the genome etc
#https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz
#https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gtf.gz

splitn=2 #  number of sbatchs to split accession file into.  Makes running faster if there are many accessions in the list

mkdir LCM_data_analysis #create a new directory for the analysis of the SRA dataset being downloaded and then move into it
cd LCM_data_analysis


#Create Sub-directories used for this analysis
mkdir prefetch_tmp
mkdir downloaded_reads
mkdir pretrim_quality
mkdir trimmed_reads
mkdir posttrim_quality
mkdir alignments
mkdir sorted_alignments
mkdir stringtie


#############################

##create bash script to download and trim reads; needs index first, may need to check options settings

cat << 'EOF' > download_trim_and_map.sh  #this is a readdoc notation.  It prints everything until EOF as written to the "download_trim_and_map.sh file.  In this case I am using it to create a bash script that will take one SRA accession as an input.  $1 references the first input argument in a bash script
#!/bin/bash

#variables
basepath=$(pwd)  # finds the current directory so that the program can use complete path names to specify files.  file name variables can be updated to use relative paths if needed.
genome="GCF_*.fna"  #Make sure this matches after indexing
threads=24
memory=200  #memory in gb bbduk is allowed to use
prefetch_pipe=./prefetch_tmp/${1}
fastq_1=${basepath}/downloaded_reads/${1}_1.fastq
fastq_2=${basepath}/downloaded_reads/${1}_2.fastq
trimmed_1=${basepath}/trimmed_reads/${1}_1.trimmed_fastq.gz
trimmed_2=${basepath}/trimmed_reads/${1}_2.trimmed_fastq.gz
out=${basepath}/alignments/${1}.bam

prefetch $1 --max-size 80g -O prefetch_tmp  # downloads a compressed version of the raw reads.  Max size must be increased since the files are larger (~24gb) than the preset maximum of 20gb
fasterq-dump -f --threads ${threads} --split-3 --outdir downloaded_reads ${prefetch_pipe} &&  #this will decompress the reads into 3 files (split 3 option) pair 1 pair 2 and unpaired reads.  It appears the unpaired reads were not uploaded to the database as there are none when complete
fastqc -t ${threads} -o pretrim_quality $fastq_1 && # will show the quality data for the reads pre trimming
fastqc -t ${threads} -o pretrim_quality $fastq_2 &&
bbduk.sh -Xmx${memory}g in1=${fastq_1} in2=${fastq_2} out1=${trimmed_1} out2=${trimmed_2} threads=${threads} ref=adapters minlen=50 qtrim=r trimq=15 k=23 ktrim=r hdist=1 mink=11 tpe tbo &&  # see https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/ for details on the settings of this trimmer
fastqc -t ${threads} -o posttrim_quality $trimmed_1 &
fastqc -t ${threads} -o posttrim_quality $trimmed_2 &
wait
hisat2 -p ${threads} --dta -x ${genome} -1 ${trimmed_1} -2 ${trimmed_2} -S ${out}  #hisat2 will align the read pairs to the reference genome.  Default settings were used in the original paper.   I added --dta to so that the output is better for stringtie. Note that hisat2 is an aligner meant for rna unlike bwa-mem to allow for splicing

wait
rm -r prefetch_tmp


EOF


###########################
#create a bash script that will take previous alignments and give us gene counts TPM using stringtie

cat << 'EOF' > alignment_readcounts.sh
#!/bin/bash
#variables
basepath=$(pwd)  # finds the current directory so that the program can use complete path names to specify files.  file name variables can be updated to use relative paths if needed.
genome="GCF_*.fna"  #Make sure this matches after indexing
annotation="GCF_*.gtf"
threads=16
input_bam=${basepath}/alignments/${1}.bam
sorted_bam=${basepath}/sorted_alignments/${1}.bam
string_readcounts=${basepath}/stringtie/${1}.readcounts
    

samtools sort -o $sorted_bam $input_bam
stringtie  -G ${annotation} -eB -o $string_readcounts $sorted_bam # -G provided a annotation file to work off of.  It is not required but we have one. -e optimizes stringtie for estimating expression levels.  -B indicates that we want the output to be a table of expression values.  There is a better way to do this that would be better for estimating splicforms, but I think the paper did not do this and I'm getting tired.
wait
exit

EOF


### Creates the Slurm header needed to run the download and mapping on HPC cluster. This will likely need to be modified.  For ISU nova users the following is a good resource to update this to reflect your specific hardware:https://www.hpc.iastate.edu/guides/nova/slurm-script-generator-for-nova

cat << 'EOF' > align_slurmheader.txt
#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --time=0-4:00:00  # max job runtime
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --partition=biocrunch  # partition(s)
#SBATCH --mem=800G  # max memory
#SBATCH -J "LCM_rnaSeq align"  # job name

EOF

cat << 'EOF' > readcount_slurmheader.txt
#!/bin/bash

#SBATCH --time=0-2:00:00  # max job runtime
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --partition=biocrunch  # partition(s)
#SBATCH --mem=800G  # max memory
#SBATCH -J "LCM_readcount"  # job name

EOF



### Creates the Slurm header needed to run on HPC cluster, this one is more specific to the slurm batch for building reference genome index, which should not take as long and is somewhat memory/disk space intensive.
cat << 'EOF' > slurmheader_index.txt
#!/bin/bash
    
#SBATCH --time=0-00:15:00  # max job runtime
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --partition=biocrunch  # partition(s)
#SBATCH --mem=600G  # max memory
#SBATCH -J "build LCM index"  # job name
EOF



## This creates a file that lists the modules we will need to load in order to run the download through alignment step.  For maximum reproducibilty this should really list the exact spack module name and not the general name.

cat << 'EOF' > align_loadmodules.txt
module load fastqc
module load sratoolkit
module load bbmap
module load samtools
module load hisat2
module load parallel
EOF

### modules needed for getting readcounts form alignments
cat << 'EOF' > readcount_loadmodules.txt
module load samtools
module load stringtie
module load parallel
EOF




## create batch script to download and create the reference genome index required for the alignment step
cat slurmheader_index.txt > index_genome.batch
echo "module load bwa" >> index_genome.batch
echo "module load gzip" >> index_genome.batch
echo "module load hisat2" >> index_genome.batch
echo "wget ${reference_genome_url}" >> index_genome.batch &&
echo "wget ${reference_annotation_url}" >> index_genome.batch &&
echo "genome=*.fa.gz" >> index_genome.batch
echo "annotation=*.gtf.gz" >> index_genome.batch
echo "gzip -d \$genome" >> index_genome.batch
echo "gzip -d \${annotation}" >> index_genome.batch
echo "genome=*.fa" >> index_genome.batch
echo "annotation=*.gtf" >> index_genome.batch
echo "hisat2_extract_splice_sites.py \${annotation} > genome_splice_sites" >> index_genome.batch
echo "hisat2_extract_exons.py \${annotation} > genome_exons" >> index_genome.batch # extract exons from gtf file so hisat may add them to the index and use them to inform the alignment
echo "hisat2-build  --exon genome_exons --ss genome_splice_sites \$genome \$genome" >> index_genome.batch #extract splice sites from gtf file so hisat may add them to the index and use them to inform the alignment


#there are unlabeled features in the NCBI gtf file which cause problems for stringtie.  A side effect of how NCBI currently process certain things apparently.  The following should remove those
#There is also trans-splicing gene and at least one other transcript for which strand is not assigned, (?), which also causes problems for string tie.  The following greps everything without one of the problem question marks.  This removes almost 40000 lines and is probably a little heavy handed.  I'll test an awk script that looks just at that column or I'll just repeat everything using the ENSEMBL reference and annotation.
#cat << 'EOF' >> index_genome.batch
#grep -P '\btranscript_id\s+"[^"]+"' GCF_000001735.4_TAIR10.1_genomic.gtf > fixed_GCF_000001735.4_TAIR10.1_genomic.gtf
#grep -ve '?' fixed_GCF_000001735.4_TAIR10.1_genomic.gtf > fixed_twice_GCF_000001735.4_TAIR10.1_genomic.gtf
#EOF



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
    sbatch --dependency=afterok:$rna_align_jobid --parsable readcount_${Chunk}.batch
done
