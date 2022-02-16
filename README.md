# Analyzing the Brain Development scRNA-seq Data between Mouse and Human

## 1. Introduction
This is a personal project aiming to identify the differences between the human and mouse brain development using data generated by scRNA-seq. The aim for this project is to:
- Implement scRNA-seq data from GEO to Seurat and Scanpy and process the data
- Use adult/embryo data of human/mouse, and compare
- Endeavor to figure out differences of gene expression, cell type etc

## 2. Reproducing Previous Studies
In advance, we start by reproducing existing studies to get a grip of workflow. The Gene Expression Omnibus (GEO) provides data produced by sequencing and microarray. Data deposited in the database are categorized among the following hierarchy:
- Platform: The method of generating data (ex Illumina HiSeq 4000)
- DataSets: A refined version of Series that is a set of comparable Samples
- Series: The collection of Samples that describe a whole study. Our focus.
- Samples: The sample source, method of analysis. The lowest of the hierarchy.

Here we use data from Series since it provides the raw data from a complete study. We choose a study that analyzes development differences between human and mouse each.

### 2.1 Human Study 1
#### 2.1.1 Load Data
Darmanis, et al(2015) is deposited as [GSE67835](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835). Here 466 cells were generated by scRNA-seq. We use the raw TAR file and unzip it, revealing 466 individual CSV files. Individual CSV represents one cell, and stores count for each gene, with 22,088 rows. In example is shown below.

[GSM1657871_1772078217.C03.csv](https://github.com/Deserav/mouse-human-develop/files/7948785/GSM1657871_1772078217.C03.csv)

In R, we implement all the CSV files in to one expression matrix by `for` looping all file names and `full_join` each file. Read the file again and we have the expression matrix. We found that this study was published before Seurat was developed, so we will do downstream analysis without Seurat.

[human_exp.csv](https://github.com/Deserav/mouse-human-develop/files/7955447/human_exp.csv)

Note that the last three rows are named: `no_feature`, `ambiguous` and `alignment_not_unique`. This seems to be in the data because HTSeq was used to convert reads to counts. We therefore remove the last three rows.

In case of the reference genome, the authors of the study used the hg19 with the STAR alignment method. However, the gene names that come from the file name are different to that of hg19. It is assumed that the actual names come from HGNC (HUGO Gene Nomenclature Committee).

#### 2.1.2 Data Processing
Data processing is undergone after loading dataset. The expression matrix should be normalized for downstream analysis. Matrix entries are converted to CPM(entries divided by column sum times 10^6) and then plugged into log10(1+p).

#### 2.1.3 Plotting
We now try to replicate the plots of the paper. Figure S1 shows histograms of total number of reads, fraction of mapped reads, intron/exon ratio, and gene body coverage. Only the total number of reads was replicated because all the data we have is the count matrix.
![reads_per_cell](https://user-images.githubusercontent.com/88135502/152650420-3c312f6f-fb7a-4fc9-aedd-fd89c6204a43.png)

For Figure S2 and Figure 1A, the expression matrix undergoes ViSNE by the tsne package, and model-based hierarchical clustering by mclust package. The goal is to identify 10 clusters by an unbiased method. In other words, we want 10 clusters without any genetic information. However, it is found that the tsne package does not contain a ViSNE function. Thus, tSNE is performed after arcsine transformation. Unfortunately, the best we could do was getting 8 clusters, and the plot is the following (Click plot to zoom)

| Figure 1A | Figure S2 Left | Figure S2 Right|
|-----------|----------------|----------------|
|![unbiased_clustering](https://user-images.githubusercontent.com/88135502/152650953-19434c73-1bfd-4476-832d-51e679edb6c5.png) |  ![BIC](https://user-images.githubusercontent.com/88135502/152650962-65f0f946-c141-492f-94cf-7c5daea78208.png) | ![uncertainty](https://user-images.githubusercontent.com/88135502/152650968-a7186ed1-3281-4f6a-a3e2-729c90bed994.png)|

The paper does not exactly notify how each respective graph was generated, so this was the closest we could get.

### 2.2 Human Study 2
#### 2.2.1 Load Data
Zhong, et al (2018) is deposited at [GSE104276](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104276). Download the whole data by `getGEOSuppFiles` function, and we have two TAR files and one GZ file. The TAR files store a single XLSX file that contain UMI count and UMI TPM respectively for all 2394 cells. On the other hand, the GZ file contains 39 TXT files of TPM expression matrix for each GEO Sample (GSM). Since there may be a batch effect, we decided to use the GZ file. Use a `for` loop to read each file into a `Seurat` object and `merge` it altogether.
The confounding part of this process is that the row names of each file are not uniform. 

[GSM2884059_GW8_PFC1.UMI_TPM_no_ERCC.txt](https://github.com/Deserav/mouse-human-develop/files/8010410/GSM2884059_GW8_PFC1.UMI_TPM_no_ERCC.txt)

[GSM2970391_GW23_PFC2_1.UMI_TPM_no_ERCC.txt](https://github.com/Deserav/mouse-human-develop/files/8010409/GSM2970391_GW23_PFC2_1.UMI_TPM_no_ERCC.txt)

The first file has "Gene" in the row names, while the second one does not. The `for` loop for making the `Seurat` object should consider these differences. If not the last column will be omitted.

This is a common problem when reading files, and as discussed in [this link](https://www.programmingr.com/r-error-messages/more-columns-than-column-names/), we ditch the `tidyverse` method and loop `read.csv`. Create 39 Seurat objects, and merge all of them to one and perform downstream analysis.

#### 2.2.2 Data Processing
The initial data is given by a TPM matrix. The developers of Seurat mentioned that raw expression matrix should be normalized via `NormalizeData` function, while a TPM input should not, as discussed in this [QnA](https://github.com/satijalab/seurat/issues/668). Instead, the researchers used a custom normalization method log((TM/10)+1). So we use the `NormalizeData` function setting  `LogNormallize` with scale factor 10^5.

#### 2.2.3 Plotting
In prior, we subset the data by choosing the cells that express at least 1000 genes, and genes that are expressed by at least 3 cells. The general Seurat workflow is applied with the operation of `FindVariableFeatures`, `ScaleData`, `RunPCA`, `FindNeighbors`, `FindClusters`, 'RunTSNE` `DimPlot`(this is called 'general workflow' from now on). From here the plot is generate as the following(click plot to zoom). Resolution of `FindClusters` was set to 0.2 such that the output generates six clustser, as in the paper.

|![plot2_plot3_tsne_plot_res 0 2](https://user-images.githubusercontent.com/88135502/154212301-7ee80b0c-28fc-47a4-bac1-18ea12c9a241.png)|![plot4_feature_plot_res 0 2](https://user-images.githubusercontent.com/88135502/154212320-d5c9f9ae-f888-490b-92bd-7bbbaa04933a.png)|![plot5_identity_plot_res 0 2](https://user-images.githubusercontent.com/88135502/154212348-53314632-dbc5-49a3-ab85-5f6f34ceaa4c.png)|
|-----------|----------------|----------------|

The problem is that first of all, it seemed that a batch effect occurred, because the clusters of the left plot seems similar. Second, hemoglobin genes were expressed throughout the plot, regardless of the clusters. Endeavor of solving these problems by normal means were futile.
- Integrating the data via `SelectIntegrationFeatures` did not work because some samples such as GW08 had too little cells (23 cells)
- Removing cluster 0 (expressed hemoglobin genes found by `FindAllMarkers` function) and operating the generel workflow does not show meaningful results. Plot made by this process is shown here. The plot 3 that the clustering is generally done better, but in plot 2 we can see that still a significant amount of hemoglobin genes are being expressed, which was not filtered by cluster 0. 

|1|2|3|
|-----------|----------------|----------------|
|![plot8_plot9_tsne_plot_hemo_removed_res 0 2](https://user-images.githubusercontent.com/88135502/154215691-04d996ea-19ea-4fd5-bae7-48fce2d66377.png)|![plot10_feature_plot_hemo_removed res 0 2](https://user-images.githubusercontent.com/88135502/154215723-286cad22-4ed9-4e70-8e4a-34b750a9b8e0.png)|![plot12_marker_plot_hemo_removed res 0 2](https://user-images.githubusercontent.com/88135502/154215749-2d943854-9ce2-4408-a6a4-344919a7d2d1.png)|


#### 2.2.4 Troubleshooting


