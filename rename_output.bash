#!/bin/bash
mkdir ./readcounts
tail -n+2 SRR_tissue.labels | sed 's/,/ /g' > replaced_SRR_tissue.labels
while read a b; do mv ./sc_rnaseq_analysis/star_output/${a}.genecountsReadsPerGene.out.tab ./readcounts/${b}_readcounts.tab; done < replaced_SRR_tissue.labels
