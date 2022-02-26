library(Seurat)
library(SeuratWrappers)
library(BiocManager)
library(GEOquery)
library(tidyverse)
library(plyr)
library(ggplot2)
library(monocle)
library(monocle3)
library(Matrix)
library(DESeq2)

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
# Genes that are expressed by at least 3 cells with greater than normalized value 1 are used
select_c<-WhichCells(human, expression = nFeature_RNA > 1000)
select_f<-rownames(human)[Matrix::rowSums(human@assays$RNA@data > 1) > 3]

human<-subset(human, features = select_f, cells = select_c)

rm(select_c, select_f)

# Drop the two cells with no features nor genes
human<-subset(human, subset = nCount_RNA != 0)

VlnPlot(human, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(human, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Find highly variable features and plot them
human<-FindVariableFeatures(human, selection.method = "vst", nfeatures = 1000)

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
plot2<-ElbowPlot(human) # Choose 5 components
plot3<-VizDimLoadings(human, dims = 1:10, reduction = "pca")

# Clustering
human<-FindNeighbors(human, dims = 1:5)
human<-FindClusters(human, resolution = c(0.1, 0.2, 0.3, 0.4, 0.7, 1))
human<-FindClusters(human, resolution = 0.2) # Run this so that resolution = 0.2 result is stored at seurat_clusters column

View(human@meta.data)

# tSNE
set.seed(0)
human<-RunTSNE(human, dims = 1:5)

# Plot Preview
# We use resolution 0.2 so that 6 clusters is found
DimPlot(human, group.by = "RNA_snn_res.0.2", label = T) # Change res.0.3 to res.0.1 or whatever
FeaturePlot(human, features = c("PTPRC", "CSF1R", "AIF1", 
                                "HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"))
VlnPlot(human, group.by = "RNA_snn_res.0.2", features = c("PTPRC", "CSF1R", "AIF1", 
                            "HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"))

plot4<-DimPlot(human, group.by = "RNA_snn_res.0.2",reduction = "tsne")
plot5<-DimPlot(human, group.by = "orig.ident",reduction = "tsne")

plot4+plot5 # Seems like there is batch effect

# Cluster 6 of res.0.2 and res.0.3 have both 16 cells which are hemoglobin cells
table(human@meta.data$RNA_snn_res.0.2)
table(human@meta.data$RNA_snn_res.0.3)
table(human@meta.data$RNA_snn_res.0.4)

# Find clusters with specific genes mentioned in the paper
plot6<-FeaturePlot(human, features = c("PTPRC", "CSF1R", "AIF1", 
                                       "HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"))
plot7<-VlnPlot(human, features = c("PTPRC", "CSF1R", "AIF1", 
                                   "HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"))

# Find markers of cluster 6 (res.0.2)
c6.markers<-FindMarkers(human, ident.1 = 6, min.pct = 0.25)
c("PTPRC", "CSF1R", "AIF1", "HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ") %in% rownames(c6.markers) # We have PTPRC, CSF1R, AIF1

# Separate cluster 6 with the rest
human.c6.excl<-subset(human, idents = c("0", "1", "2", "3", "4", "5"))
human.c6<-subset(human, idents = 6)
FeaturePlot(human.c6, features = "PTPRC")

# Standard Workflow : dims = 1:10, resolution = 0.2 looks the best
human.c6.excl <- FindVariableFeatures(object = human.c6.excl, selection.method = "vst", nfeatures = 1000)

all.genes<-rownames(human.c6.excl)
human.c6.excl <- ScaleData(object = human.c6.excl, features = all.genes)

human.c6.excl <- RunPCA(object = human.c6.excl, features = VariableFeatures(object = human.c6.excl))
plot8<-ElbowPlot(human.c6.excl) #Choose 10 components
plot9<-VizDimLoadings(human.c6.excl, dims = 1:10, reduction = "pca")

human.c6.excl <- FindNeighbors(object = human.c6.excl, dims = 1:10)
human.c6.excl <- FindClusters(object = human.c6.excl, resolution = c(0.1, 0.2, 0.3, 0.4, 0.7, 1))
human.c6.excl <- FindClusters(object = human.c6.excl, resolution = 0.2)

set.seed(0)
human.c6.excl <- RunTSNE(object = human.c6.excl, dims = 1:10)
plot10<-DimPlot(object = human.c6.excl, group.by = "RNA_snn_res.0.2", reduction = "tsne")
plot11<-FeaturePlot(human.c6.excl, features = c("PTPRC", "CSF1R", "AIF1", 
                                "HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"))
plot12<-FeaturePlot(human.c6.excl, features = c("PAX6", "NEUROD2", "GAD1", "PDGFRA", "AQP4", "PTPRC"))

# We use these results for downstream analysis
# Rename clusters by cell type
new.id<-c("Excitatory neurons", "Interneurons", "NPCs", "Astrocytes", "Microglia", "OPCs")
names(new.id)<-levels(human.c6.excl)
human.c6.excl<-RenameIdents(human.c6.excl, new.id)

plot13<-DimPlot(human.c6.excl, reduction = "tsne", label = T)
plot14<-DimPlot(human.c6.excl, reduction = "tsne", group.by = "orig.ident")

# From here we find marker of all clusters
all.markers <- FindAllMarkers(human.c6.excl, test.use = "roc", only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

t1 <- FindAllMarkers(human.c6.excl, test.use = "roc", min.pct = 0.1, logfc.threshold = 0.25)
t2 <- FindAllMarkers(human.c6.excl, test.use = "DESeq2", only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

top10.markers <- all.markers %>% group_by(cluster) %>% slice_min(n = 10, order_by = p_val)
top10.markers <- top10.markers %>% pull(gene)


# Heatmap
DoHeatmap(human.c6.excl, features = top10.markers)
DoHeatmap(human.c6.excl, features = c("PAX6", "VIM", "HES1", "EMX2",
                                      "TUBB3", "NEUROD6", "SNAP25", "SATB2",
                                      "DLX2", "DLX1", "GAD1",
                                      "OLIG1", "PDGFRA",
                                      "ATP1A2",
                                      "C1QB", "CD68"))
DoHeatmap(human.c6.excl, features = c("TBR1", "BCL11B", "POU3F2", "TUBB3", "NEUROD6", "SATB2", "SLC17A7", "FOXG1",
                                      "DLX1", "DLX2", "GAD1", "SLC32A1",
                                      "PAX6", "VIM", "HES1", "EMX1", "EMX2", "OTX2",
                                      "ATP1A2", "GFAP", "AQP4", "ALDH1L1",
                                      "ITGAM", "AIF1", "C1QB", "CD68",
                                      "PDGFRA", "OLIG1", "S100B"))

# Trajectory analysis with Monocle3
# Transferring seurat object to Monocle 3
human.pfc<-as.cell_data_set(human.c6.excl)

set.seed(0)
human.pfc <- preprocess_cds(human.pfc, num_dim = 10)

human.pfc<-reduce_dimension(human.pfc) # Default is UMAP
human.pfc <- reduce_dimension(human.pfc, reduction_method="tSNE") # Do tSNE. Both tSNE and UMAP results are saved. We won't use this here

# Preview UMAP and tSNE plot
plot_cells(human.pfc, reduction_method="tSNE", color_cells_by="ident") # plot tSNE
plot_cells(human.pfc, color_cells_by="ident") # plot UMAP


# Genrate trajectory
human.pfc<-cluster_cells(human.pfc, reduction_method = "UMAP")
human.pfc<-learn_graph(human.pfc, use_partition = T)


# Plot trajectory
plot15<-plot_cells(cds = human.pfc, color_cells_by = "ident") + theme(legend.position = "right")
plot16<-plot_cells(cds = human.pfc, color_cells_by = "orig.ident") + theme(legend.position = "right")

# Calculate pseudotime - we don't know starting cells
# human.pfc<-order_cells(human.pfc)

# View meta data
colData(human.pfc)

# Remove Monocle 3 object. We will use Monocle 2 from scratch
rm(human.pfc)

# Psuedotime analysis with Monocle2
# Create Monocle 2 object`
# Since using TPM matrix is not optimal (works but not well), rpc_matrix is made and used.
data<-as(as.matrix(human.c6.excl@assays$RNA@counts), 'sparseMatrix') # @data for normalized matrix, @counts for raw matrix. We want the normalized matrix here

pd<-new('AnnotatedDataFrame', data = human.c6.excl@meta.data)

fData<-data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd<-new('AnnotatedDataFrame', data = fData)

human.pfc<-newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = gaussianff())

# Assign cell types
pData(human.pfc)$cell_type<-revalue(as.character(pData(human.pfc)$seurat_clusters), c("0" = "excitatory neurons", "1" = "interneurons", "2" = "NPCs",
                                                                                      "3" = "astrocytes", "4" = "microglia", "5" = "OPCs"))

# Remove microglia(meso-derived cells) and interneurons (can migrate)
human.pfc<-human.pfc[, pData(human.pfc)$cell_type != "microglia"]
human.pfc<-human.pfc[, pData(human.pfc)$cell_type != "interneurons"]

## The following process is used on TPM matrix (convert TPM to reads and back) -> takes too long 
# human.pfc<-newCellDataSet(data, phenoData = pd, featureData = fd, lowerDetectionLimit = 0.1, expressionFamily = tobit(Lower = 0.1))
# 
# rpc_matrix<-relative2abs(human.pfc, method = "num_genes")
# 
# human.pfc<-newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"), phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
# 
# human.pfc<-estimateSizeFactors(human.pfc)
# human.pfc<-estimateDispersions(human.pfc)

rm(data, pd, fData, fd)

# Preplotting is not working
# plot_ordering_genes(human.pfc)
# plot_pc_variance_explained(human.pfc, norm_method = "none")

# Pseudotime analysis
set.seed(0)
human.pfc<-setOrderingFilter(human.pfc, row.names(all.markers))

human.pfc<-reduceDimension(human.pfc, norm_method = "none", max_components = 5, num_dim = 10, reduction_method = "DDRTree")

# First run orderCells
human.pfc<-orderCells(human.pfc) # State is generated here
levels(pData(human.pfc)$State)

plot_cell_trajectory(human.pfc)
plot_complex_cell_trajectory(human.pfc, color_by = "cell_type")

# Set state for root_state
GM_state<-function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    GW08_counts<-table(pData(cds)$State, pData(cds)$orig.ident)[, "GW08"]
    return(as.numeric(names(GW08_counts)[which(GW08_counts == max(GW08_counts))]))
  } else {
    return(1)
  }
}

# Second run orderCells -> doesn't change the results
human.pfc<-orderCells(human.pfc, root_state = GM_state(human.pfc))


plot17<-plot_complex_cell_trajectory(human.pfc, color_by = "cell_type", show_branch_points = T)
plot18<-plot_complex_cell_trajectory(human.pfc, color_by = "State")
plot19<-plot_complex_cell_trajectory(human.pfc, color_by = "orig.ident")

plot_cell_trajectory(human.pfc, color_by = "cell_type")

plot20<-plot_multiple_branches_pseudotime(human.pfc[c("PAX6", "SOX2", "NEUROD2", "OLIG1", "AQP4"), ],
                                  branches = c(1, 2), branches_name = c("Glial", "Neuronal"),
                                  color_by = "Branch", cell_size = 1, nrow = 1, ncol = 5)

plot_genes_branched_pseudotime(human.pfc[c("PAX6", "SOX2", "NEUROD2", "OLIG1", "AQP4"), ],
                               branch_point = 1, color_by = "cell_type", nrow = 1, ncol = 5)
