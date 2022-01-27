library(BiocManager)
library(GEOquery)
library(Seurat)
library(SeuratDisk)
library(tidyverse)

setwd("D:/Desktop/d/human_mouse")

# Human 
# 1. From csv
human.gse<-"GSE67835_RAW"

# 2. From soft format
human.gse<-"GSE67835_family.soft"

# Set working directory from either 1 or 2
setwd(human.gse)

# Taking a look at the soft file given by GEO. The expression matrix is absent
write.table(raw.data, file = "soft.txt", sep = "\n",row.names = FALSE)


raw.data <- read.table("GSE67835_family.soft", sep="\t", header=TRUE)

# For looping to combine all the raw csv files (Try 1)
file_list<-list.files()
names(file_list)<-basename(file_list)

for (file in list.files()){
  # Create the first data if no data exists yet
  if(!exists("human")){
    human<-read_tsv(file, col_names = c("Gene", substr(file, 1, nchar(file)-4)))
  }  
  # If data already exists, then append it together
  else {temp<-read_tsv(file, col_names = c("Gene", substr(file, 1, nchar(file)-4)))
    human<-full_join(human, temp, by="Gene")
    rm(temp)
  }
}

setwd("D:/Desktop/d/human_mouse")
write_csv(human, file = "human.csv")

human<-column_to_rownames(human, "Gene")
human<-as.matrix(human)

human<-CreateSeuratObject(human, project = "Human")






