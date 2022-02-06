library(BiocManager)
library(GEOquery)
library(tidyverse)
library(mclust)
library(ggplot2)
library(ComplexHeatmap)
library(Rtsne)
library(edgeR)
library(scatterplot3d)

setwd("D:/Desktop/d/human_study_1")

# Load Human Data
## 1. From csv
human.gse<-"GSE67835_RAW"

## 2. From soft format
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
  if(!exists("human.exp")){
    human.exp<-read_tsv(file, col_names = c("Gene", substr(file, 1, nchar(file)-4)))
  }  
  # If data already exists, then append it together
  else {temp<-read_tsv(file, col_names = c("Gene", substr(file, 1, nchar(file)-4)))
    human.exp<-full_join(human.exp, temp, by="Gene")
    rm(temp)
  }
}

# Save the matrix just in case
setwd("D:/Desktop/d/human_study_1")
write_csv(human.exp, file = "human_exp.csv")

# Change data frame to expression matrix so Seurat can read it
human.exp<-read_csv("human_exp.csv") # Just read this when starting a new session

human.exp<-column_to_rownames(human.exp, "Gene")
human.exp<-as.matrix(human.exp)

# Delete the last three rows: no_feature, ambiguous, alignment_not_unique
human.exp<-human.exp[1:(nrow(human.exp)-3),]

# Histogram of number of reads per cell
hist(colSums(human.exp), xlab = "Reads", main = "Number of reads per cell")

# Normalize expression matrix
human.reduce<-log10(cpm(human.exp)+1)

# ViSNE
asinh_scale<-5
human.reduce<-asinh(human.reduce/asinh_scale)

## Check ?Rtsne and you can see that X argument needs row as obs, and col as variable
set.seed(0)
human.tsne<-Rtsne(t(human.reduce), dim = 3, check_duplicates = F) # dim = 2 makes 9 groups


# BIC
human.expBIC<-mclustBIC(human.tsne$Y)
plot(human.expBIC)

human.mclust<-Mclust(human.tsne$Y, G = 1:10)
plot(human.mclust, what = "classification")

summary(human.expBIC)
summary(human.mclust)

group<-human.mclust$classification
uncertain<-human.mclust$uncertainty
boxplot(uncertain~group)

# Plot 3D ViSNE plot
## plot() in mclust uses human.mclust$data as scatterplot
X<-human.mclust$data
color<-human.mclust$classification

plot<-scatterplot3d(X, color = factor(color), xlab = "tsne_1", ylab = "tsne_2", zlab = "tsne_3",
              main = "Unbiased Clustering")

