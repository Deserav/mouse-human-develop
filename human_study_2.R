library(Seurat)
library(BiocManager)
library(GEOquery)
library(tidyverse)
library(ggplot2)

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
objects<-ls() # (Optional) Save ls() to get all Seurat objects at this point. 

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

# Each file is a matrix of raw TPMs -> normalize by log((TPM/10)+1)
human<-NormalizeData(human, normalization.method = "LogNormalize", scale.factor = 100000)

# Cells that express at least 1000 genes used
# Genes that are expressed by at least 3 cells are used
select_c<-WhichCells(human, expression = nFeature_RNA > 1000)
select_f<-rownames(human)[Matrix::rowSums(human@assays$RNA@data != 0) > 3]

human<-subset(human, features = select_f, cells = select_c)

# Drop the two cells with no features nor genes
human<-subset(human, subset = nCount_RNA != 0)

VlnPlot(human, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(human, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Find highly variable features and plot them
human<-FindVariableFeatures(human, selection.method = "vst", nfeatures = 2000)

top10<-head(VariableFeatures(human), 10)

plot1<-VariableFeaturePlot(human)
LabelPoints(plot = plot1, points = top10, repel = T) # May gives error message that x-axis value is Inf


# Loess in FindVariableFeature performs log transformation, and the ones that become 0 are ommitted in VariableFeaturesPlot
# In this case there are none
View(subset(human@assays$RNA@meta.features, subset = vst.mean == 0))

# Scale Data
all.genes<-rownames(human)
human<-ScaleData(human, features = all.genes)

# PCA
## Option 1: Use only variable features - Use this
human<-RunPCA(human, features = VariableFeatures(object = human))

## Option 2: Use all genes - Don't use this unless necessary
human<-RunPCA(human)

# Plot PCA
DimHeatmap(human, dims = 1:5, cells = 2234, balanced = T)
ElbowPlot(human) # Choose 10 components

# Clustering
human<-FindNeighbors(human, dims = 1:10)
human<-FindClusters(human, resolution = c(0.1, 0.2, 0.3, 0.4, 0.7, 1))
human<-FindClusters(human, resolution = 0.2) # Run this so that resolution = 0.2 result is stored at seurat_clusters column

View(human@meta.data)

# We use resolution 0.2 so that 6 clusters is found
DimPlot(human, group.by = "RNA_snn_res.0.2", label = T)

set.seed(0)
human<-RunTSNE(human, dims = 1:10)
plot2<-DimPlot(human, group.by = "RNA_snn_res.0.2",reduction = "tsne")
plot3<-DimPlot(human, group.by = "orig.ident",reduction = "tsne")

plot2+plot3 # Seems like there is batch effect

# Find clusters with specific genes mentioned in the paper
plot4<-FeaturePlot(human, features = c("PTPRC", "CSF1R", "AIF1", 
                                       "HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"))
plot5<-VlnPlot(human, features = c("PTPRC", "CSF1R", "AIF1", 
                                   "HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"))

# Find markers of all clusters -> saved in to dataframe
human.markers <- FindAllMarkers(human, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Top 2 genes per cluster
human.markers %>%  group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

# All marker genes per cluster : change cluster number from 0 to 5
m0<-human.markers %>% filter(cluster == 0) %>% arrange(desc(avg_log2FC)) %>% rownames()

# Are markers in our dataframe?
c("PTPRC", "CSF1R", "AIF1", "HBA1", "HBA2","HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ") %in% m0
c("PTPRC", "CSF1R", "AIF1", "HBA1", "HBA2","HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ") %in% rownames(human.markers)

# Exclude cluster 0
human<-subset(human, idents = c("1", "2", "3", "4", "5"))

# Cluster and run tSNE from here
RunPCA(human, features = VariableFeatures(object = human))
human<-FindNeighbors(human, dims = 1:10)
human<-FindClusters(human, resolution = c(0.1, 0.2, 0.3, 0.4, 0.7, 1))

human<-RunTSNE(human, dims = 1:10)
plot8<-DimPlot(human, group.by = "RNA_snn_res.0.1",reduction = "tsne")
plot9<-DimPlot(human, group.by = "orig.ident",reduction = "tsne")

plot8 + plot9

# Find clusters with specific genes mentioned in the paper
plot10<-FeaturePlot(human, features = c("PTPRC", "CSF1R", "AIF1", 
                                       "HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"))
plot11<-VlnPlot(human, features = c("PTPRC", "CSF1R", "AIF1", 
                                   "HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"))
plot10
plot11

# Plot with markers
# PAX6: NPC, NEUROD2: excitatory neuron, GAD1: interneuron, GAD1: OPC, AQP4: astrocyte, PTPRC: microglia
plot12<-FeaturePlot(human, features = c("PAX6", "NEUROD2", "GAD1", "PDGFRA", "AQP4", "PTPRC"))
plot13<-VlnPlot(human, features =  c("PAX6", "NEUROD2", "GAD1", "PDGFRA", "AQP4", "PTPRC"))

# Data integration
obj.list<-SplitObject(human, split.by = "orig.ident")
for(i in 1:length(obj.list)){
  obj.list[[i]]<-NormalizeData(object = obj.list[[i]])
  obj.list[[i]]<-FindVariableFeatures(object = obj.list[[i]])
}

# Select integration features - Not working
features<-SelectIntegrationFeatures(object.list = obj.list)

# Find integration anchors (CCA) - Not working
anchors<-FindIntegrationAnchors(object.list = obj.list, anchor.features = features, dims = 1:10)
 
# Integrate data - Not working
human.int<- IntegrateData(anchorset = anchors)
