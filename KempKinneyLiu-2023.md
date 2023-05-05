### Group members: Joshua M Kemp, Shelly Kinney, Yue Liu

&nbsp;  

# Paper Overview

***Paper:*** Single-cell RNA-seq analysis reveals ploidy-dependent and cell-specific transcriptome changes in Arabidopsis female gametophytes

***Paper link:*** https://link.springer.com/article/10.1186/s13059-020-02094-0 

***Authors:*** Qingxin Song, Atsumi Ando, Ning Jiang, Yoko Ikeda, & Z. Jeffrey Chen

***Summary:***

Basic Premise: Polyploidy and transcriptome changes

Ploidy affects cell size; however, the effects of ploidy on RNA transcripts are not well studied.

This paper uses *Arabidopsis thaliana* to compare RNA transcript level differences between diploid and tetraploid plants in female gametophytic cells.

Methods are explained; however, specific details and code about processing the RNA sequence data is not listed, nor is most of the code for the paper's figures.

The statistical analysis code is given, as well as the normalization and PCA plot code. 
Specific details about what the code is doing are limited to a broad comment about the function of the entire script.



***Goal of project:*** Recreate Figure 5 Venn Diagram and Figure S4 Heatmap

 
# Workflow Details

***Note:*** specific program details are listed within the code or in a separate file in the Code directory (Spack_Modules_used)

***Note:*** for rerunning code, make sure code includes references to correct working directories or file locations



***The following code was run in the order listed.***





## 1 Bash script: LCM_dataset_alignment.sh

This script is largely self contained, and will create then run the necessary Slurm jobs for the complete analysis. There may be some editing required, particularly on the Slurm headers, for this to be run on a particular piece of hardware.  
Within the script, with the exception of building the index, there are functions defined that will run the following steps in order for each input SRA accension.


Obtain the reference genome and its annotation.

Build an index using the reference genome, which will be used for the RNA alignment

Obtain the RNA sequence data from the SRA database

The files will be a compressed version, which must be broken up to obtain paired reads with fasterq-dump.
  
Take a FASTQC to check RNA read data

Trim the RNA reads with bbduk

Take a second FASTQC to check the trimming quality 

Align the RNA reads to the indexed genome

Sort the subsequent BAM files

Obtain the read counts from the alignment BAM files using stringtie


## 2 Bash script: sc_rnaseq_analysis.sh

This script is adapted from LCM_dataset_alignment.sh and will run in a very similar manner.


Obtain the reference genome and its annotation.

Build an index using the reference genome, reference annotation and ERCC spike-in information.

Obtain the RNA sequence data from the SRA database.
  
The files will be a compresed version, which must be broken up to obtain reads.  Fasterq-dump is used for this.

Take a FASTQC to check the RNA read data.

Use STAR Align the RNA reads to the indexed genome and obtain alignment SAM files and readcount data.



## 3 Bash script: rename_output.bash

Rename STAR output files to list the cell type information.




## 4 UNIX/R script: R_script_to_merge_readcount_files

Copy out_put files if needed.

Create new file with list of file names.

Join all read count files into a single R dataframe.

Separate out spike-ins readcount numbers from actual transcript readcounts.




## 5 R script: R_ for_Venn.R

Take in merged file of read counts.

Separate out data by cell type.

Count the number of reads per cell type.

Separate out read names.

Combine read names and counted reads per cell type.

Select only reads with counts above 0.

Take read names and compare with other cell types in venn diagrams.

Three diagrams are created: 
Diploid, 
Tetraploid,
Total Read Count

The third is presumed to be a similar diagram to Figure 5b within original paper.





## 6 UNIX code: Download_existing_LCM(ref64)expression_dataset_and_import_to_R

Download and import dataset to R.





## 7 R script: Data_visualization_heatmap.Rmd

Merge readcount files.

Select cell types for egg cell (diploid), egg cell (tetraploid), central cell (diploid), and central cell (diploid).

Reorder cell types to be grouped together.

Select specific genes, filter on gene ID, reorder rows, remove first column.

Convert data into matrix format.

Make annotations into rows and columns.

Plot data into heatmap, and export heatmap figure into a file.




# Results

## 1 Data Processing:

Read count differences; our filtering is not as strict as the paper's filtering steps.


## 2 Venn Diagram (Figure 5b):

![Venn Diagram](https://github.com/0Yue0/BCB546_Spring2023_Final/blob/main/results/%23Final_venn_diagram.png)

Overall, the central cell and egg cell have the most reads in common, followed by all three cell types, the central cell and the synergid cell, and lastly the egg cell and synergid cell have the least reads in common. These match the trends in the original paper's diagram.

Our figure has more unique reads for the egg cell and fewer unique reads for both the central cell and synergid cell than the original paper's diagram.

## 3 Heatmap (Figure S4):

![Heatmap](https://github.com/0Yue0/BCB546_Spring2023_Final/blob/main/results/Figure_S4.png)

In our analysis, reported cell type specific genes did show the specific expression pattern in certain cell types, which was the same as the original paper. 


No dosage effects due to ploidy are evident in our figure compared to the original paper. 
