library(Seurat)
library(BiocManager)
library(GEOquery)
library(tidyverse)

setwd("D:/Desktop/d/human_study_2")

gset<-getGEO("GSE104276", GSEMatrix=TRUE)
getGEOSuppFiles("GSE104276")

setwd("GSE104276")

# Make Seurat objects for each TXT file
file_list<-list.files()
names(file_list)<-basename(file_list)

# Files may or may not have "Gene" in column names
for (file in list.files()){
  cat("loading file: ", file, "\n")
  temp<-read.csv(file, sep = "\t")
  
  cat("made variable: ",  str_extract(file, "GW\\d+_PFC\\d_?\\d?"), "\n")
  
  # If it has "Gene" in column name, convert each gene to row name
  if ("Gene" %in% names(temp)){ 
    rownames(temp) <- temp[,1]
    temp[,1] <- NULL
  }
  
  # Change it to a matrix so it is readable to Seurat
  temp<-as.matrix(temp) 
  
  # Create Seurat object
  temp<-CreateSeuratObject(temp, project = str_extract(file, "GW\\d+_PFC\\d_?\\d?"))
  assign(str_extract(file, "GW\\d+_PFC\\d_?\\d?"), temp)
  rm(temp)
  }


