library(Seurat)
library(BiocManager)
library(GEOquery)
library(tidyverse)

setwd("D:/Desktop/d/human_study_2")

gset<-getGEO("GSE104276", GSEMatrix=TRUE)
getGEOSuppFiles("GSE104276")

setwd("GSE104276")

# Make Seurat objects for each TXT file
# Files may or may not have "Gene" in column names
for (file in list.files()){
  cat("loading file: ", file, "\n")
  temp<-read.csv(file, sep = "\t")

  # If it has "Gene" in column name, convert each gene to row name
  if ("Gene" %in% names(temp)){
    rownames(temp) <- temp[,1]
    temp[,1] <- NULL
  }

  # Change it to a matrix so it is readable to Seurat
  temp<-as.matrix(temp)

  # Create Seurat object
  temp<-CreateSeuratObject(temp, project = str_extract(file, "GW\\d+_PFC\\d_?\\d{0,2}"))
  assign(str_extract(file, "GW\\d+_PFC\\d_?\\d{0,2}"), temp)
  cat("made variable: ",  str_extract(file, "GW\\d+_PFC\\d_?\\d{0,2}"), "\n")
  rm(temp)
}

rm(file)
objects<-ls() # Save ls() to get all Seurat objects at this point. 

# Merge all the Seurat objects
human<-merge(GW8_PFC1, y = c(GW9_PFC1, GW10_PFC1_1, GW10_PFC1_2, GW10_PFC2_1, GW10_PFC2_2, GW10_PFC3,
                      GW12_PFC1, GW13_PFC1, GW16_PFC1_1, GW16_PFC1_2, GW16_PFC1_3, GW16_PFC1_4,
                      GW16_PFC1_5, GW16_PFC1_6, GW16_PFC1_7, GW16_PFC1_8, GW16_PFC1_9,
                      GW19_PFC1, GW19_PFC2, GW19_PFC3, GW23_PFC1_1, GW23_PFC1_2, GW23_PFC1_3,
                      GW23_PFC2_1, GW23_PFC2_2, GW23_PFC2_3, GW23_PFC2_4, GW23_PFC2_5,
                      GW26_PFC1_1, GW26_PFC1_2, GW26_PFC1_3, GW26_PFC1_4, GW26_PFC1_5,
                      GW26_PFC1_6, GW26_PFC1_7, GW26_PFC1_8, GW26_PFC1_9, GW26_PFC1_10),
      add.cell.ids = c("GW8_PFC1", "GW9_PFC1", "GW10_PFC1_1", "GW10_PFC1_2", "GW10_PFC2_1", "GW10_PFC2_2", "GW10_PFC3",
                       "GW12_PFC1", "GW13_PFC1", "GW16_PFC1_1", "GW16_PFC1_2", "GW16_PFC1_3", "GW16_PFC1_4",
                       "GW16_PFC1_5", "GW16_PFC1_6", "GW16_PFC1_7", "GW16_PFC1_8", "GW16_PFC1_9",
                       "GW19_PFC1", "GW19_PFC2", "GW19_PFC3", "GW23_PFC1_1", "GW23_PFC1_2", "GW23_PFC1_3",
                       "GW23_PFC2_1", "GW23_PFC2_2", "GW23_PFC2_3", "GW23_PFC2_4", "GW23_PFC2_5",
                       "GW26_PFC1_1", "GW26_PFC1_2", "GW26_PFC1_3", "GW26_PFC1_4", "GW26_PFC1_5",
                       "GW26_PFC1_6", "GW26_PFC1_7", "GW26_PFC1_8", "GW26_PFC1_9", "GW26_PFC1_10"),
      project = "HPFC")

# Remove all the Seurat objects expect the merged one
rm(list = setdiff(ls(), "human"))

# Make column of samples in meta.data
human$sample<-str_extract(rownames(human@meta.data), "^GW\\d+_PFC\\d+(_\\d+){0,1}")



order<-c()
for(i in 1:(length(objects))){
  temp<-objects[i]
  temp<-sub(pattern = "GW", replacement = "", x=temp)
  temp<-sub(pattern = "_PFC", replacement = "", x=temp)
  temp<-sub(pattern = "_", replacement = ".", x=temp)
  order[i]<-as.numeric(temp)
  rm(temp)
}


