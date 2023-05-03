#R version: 4.2.2
#Make sure to install necessary packages used below
#load libraries
library(VennDiagram)
library(tidyverse)
library(tidyr)
library(RColorBrewer)


#run code from "R_script_to_merge_readcount_files" before running this code
total <- read_tsv("merged_files.txt")


#separate out cell types; CC is central cell, EC is egg cell, SC is synergid cell
CCd <- total[grep("CCd", colnames(total))]
CCt <- total[grep("CCt", colnames(total))]
ECd <- total[grep("ECd", colnames(total))]
ECt <- total[grep("ECt", colnames(total))]
SC <- total[grep("DSC", colnames(total))]
CC <- total[grep("CC", colnames(total))]
EC <- total[grep("EC", colnames(total))]

#Sum the readcounts for each RNA transcript
SumCCD <- rowSums(CCd)
SumCCT <- rowSums(CCt)
SumECD <- rowSums(ECd)
SumECT <- rowSums(ECt)
SumSC <- rowSums(SC)
SumCC <- rowSums(CC)
SumEC <- rowSums(EC)

#Separate out the RNA transcript names
ReadNames <- total["Merge_Col"]

#Combine sums and read names
ReadCCD <- data_frame(ReadNames, SumCCD)
ReadCCT <- data_frame(ReadNames, SumCCT)
ReadECD <- data_frame(ReadNames, SumECD)
ReadECT <- data_frame(ReadNames, SumECT)
ReadSC <- data_frame(ReadNames, SumSC)
ReadCC <- data_frame(ReadNames, SumCC)
ReadEC <- data_frame(ReadNames, SumEC)

#Filter out all reads that are not 0
FilteredCCD <- filter(ReadCCD, SumCCD != 0)
FilteredCCT <- filter(ReadCCT, SumCCT != 0)
FilteredECD <- filter(ReadECD, SumECD != 0)
FilteredECT <- filter(ReadECT, SumECT != 0)
FilteredSC <- filter(ReadSC, SumSC != 0)
FilteredCC <- filter(ReadCC, SumCC != 0)
FilteredEC <- filter(ReadECD, SumEC != 0)

#Select out all RNA transcript names
RNCCD <- FilteredCCD["Merge_Col"]
RNCCT <- FilteredCCT["Merge_Col"]
RNECD <- FilteredECD["Merge_Col"]
RNECT <- FilteredECT["Merge_Col"]
RNSC <- FilteredSC["Merge_Col"]
RNCC <- FilteredCC["Merge_Col"]
RNEC <- FilteredEC["Merge_Col"]

#Make file into list for diagram
LRNCCD <- list(RNCCD$Merge_Col)
LRNCCT <- list(RNCCT$Merge_Col)
LRNECD <- list(RNECD$Merge_Col)
LRNECT <- list(RNECT$Merge_Col)
LRNSC <- list(RNSC$Merge_Col)
LRNCC <- list(RNCC$Merge_Col)
LRNEC <- list(RNEC$Merge_Col)

#Select colors for diagram
myCol <- brewer.pal(3, "BuPu")


#Code for Venn Diagrams came from the following website (with minor alterations):
# https://r-graph-gallery.com/14-venn-diagramm.html 

#Venn diagram for Diploid Comparisons
venn.diagram(
  x = list(LRNCCD[[1]], LRNECD[[1]], LRNSC[[1]]),
  category.names = c("Central Cell (D)" , "Egg Cell (D)" , "Synergid Cell"),
  filename = '#Diploid_venn_diagram.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)


#Venn Diagram for Tetraploid Comparisons
venn.diagram(
  x = list(LRNCCT[[1]], LRNECT[[1]], LRNSC[[1]]),
  category.names = c("Central Cell (T)" , "Egg Cell (T)" , "Synergid Cell"),
  filename = '#Tetraploid_venn_diagram.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)


#Venn diagram for Total Ploidy Comparison
venn.diagram(
  x = list(LRNCC[[1]], LRNEC[[1]], LRNSC[[1]]),
  category.names = c("Central Cell" , "Egg Cell" , "Synergid Cell"),
  filename = '#Final_venn_diagram.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

