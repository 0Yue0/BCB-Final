Group members: Joshua M Kemp, Shelly Kinney, Yue Liu

&nbsp;  

***Paper Overview***

Paper: Single-cell RNA-seq analysis reveals ploidy-dependent and cell-specific transcriptome changes in Arabidopsis female gametophytes

Paper link: https://link.springer.com/article/10.1186/s13059-020-02094-0 

Authors: Qingxin Song, Atsumi Ando, Ning Jiang, Yoko Ikeda, & Z. Jeffrey Chen

*Summary:*

Polyploidy and transcriptome changes

&nbsp;

Goal of project: Recreate Additionalfile1:S2; Figure 5A; Figure 3C ?; Figure 3D ?; Additionalfile1:S4; Additionalfile1:S6; Additionalfile1:S9

&nbsp;  
***Workflow Details***

slurm script: 

Obtain the reference genome and its annotation.

Build an index using the reference genome, and this will be used for the RNA alignment.

Obtain the RNA sequence data from the SRA database.

  The files will be a compressed version, which must be broken up to obtain paired reads.
  
Take a FASTQC to check RNA read data

Trim the RNA reads

Take a second FASTQC to check the trimming quality 

Align the RNA reads to the indexed genome

Sort the subsequent BAM files

Obtain the read counts from the alignment BAM files





&nbsp;  
***Results***


