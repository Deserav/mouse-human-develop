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

The data is also subsetted by choosing the cells that express at least 1000 genes, and genes that have at least 3 cells with normalized expression greater than 1. The paper stated that this subsetting leaves 17,854 genes and 2,333 cells. On the other hand, our process leaves us with 17,962 genes and 2,344 cells.

#### 2.2.3 Plotting
 The general Seurat workflow is applied with the operation of `FindVariableFeatures`, `ScaleData`, `RunPCA`, `FindNeighbors`, `FindClusters`, `RunTSNE`, `DimPlot`(this is called 'general workflow' from now on). Our initial objective is to remove hemoglobin highlighted cells to eliminate any contamination effect by erythrocytes. Testing multiple dimensions, we found that setting `dims = 1:5` on `FindNeighbors` and `RunTSNE` captures the hemoglobin highlighted cells. It is located at cluster 6 at plot 3.

|1|2|3|
|-----------|----------------|----------------|
|![plot2](https://user-images.githubusercontent.com/88135502/154541334-3b947804-ea67-4ba5-8f8f-5fd135c96931.png)|![plot3_visdimloading](https://user-images.githubusercontent.com/88135502/154541551-6cd214cf-f22e-4bd3-a909-28bfdd9932b4.png)|![plot4_plot5_tsne res 0 2](https://user-images.githubusercontent.com/88135502/154543151-840ef9cf-e120-4876-8ba6-e0c78ddeb676.png)|

|4|5|
|-----------|----------------|
|![plot6_featureplot](https://user-images.githubusercontent.com/88135502/154543786-d7ce7f16-9cb5-4344-9cfc-b3fad40102cf.png)|![plot7_violin](https://user-images.githubusercontent.com/88135502/154543821-a1c029c0-67f0-479e-8bad-3241a90a8701.png)|

The researchers defined hemoglobin genes as HBA1, HBA2, HBB, HBD, HBE1, HBG1, HBG2, HBM, HBQ1, and HBZ. In dimensions low or high, most of the hemoglobin genes are activated by all clusters(plot 4, plot5), which made it confounding. Thus as in plot 2, most of the hemoglobin genes are found at PC5, so we use dimension 5. We use a smaller dimension as possible to pinpoint erythrocytes shown in cluster 6 of plot 5. Therefore, we remove cluster 6, which consists of 16 cells. We now have 2,328 cells left. The researchers had 2,309 cells left after this process, which means that 24 cells were removed.

We apply the general workflow to the remaining Seurat object. Our goal is to identify neural progenitor cells(NPC), excitatory neurons, interneurons, oligodendrocyte progenitor cells(OPC), astrocytes, and microglia. With a total of six clusters, PAX6, NEUROD2, GAD1, PDGFRA, AQP4, and PTPRC were used as the respective markers. Here we use `dims = 1:10` to capture as much of the markers as possible (plot 6, plot 7)

|6|7|
|-----------|----------------|
|![plot8_elbow_hg_rm](https://user-images.githubusercontent.com/88135502/154545624-d39ef2f2-2bd1-4454-a526-b9a2d6933a86.png)|![plot9_visdims_hg_rm](https://user-images.githubusercontent.com/88135502/154834888-79a10862-61a7-43ee-8cb6-bab1faab669c.png)|

The plot generated by this process finally resembles Figure 1A of our paper.
|8|9|10|
|-----------|----------------|----------------|
|![plot10_tsne_hg_rm](https://user-images.githubusercontent.com/88135502/154546177-9ab71780-1c89-4e0c-988a-0159b7e8a253.png)|![plot11_feature_hg_rm](https://user-images.githubusercontent.com/88135502/154546294-a4ba1fc4-389f-4062-995e-6306beb0eba0.png)|![plot12_marker_hg_rm](https://user-images.githubusercontent.com/88135502/154546427-26023954-a43d-41c2-b771-d4075eb9f479.png)|

The gene expression pattern is shown as:
- PAX6 : part of cluster 2 and most of cluster 3
- NEUROD2: most of cluster 0
- GAD1: most of cluster 1
- PDGFRA: most of cluster 5
- AQP4: part of cluster 3
- PTPRC: most of cluster 4.

Although there is some error, we therefore allocate our clusters as:
- 0: excitatory neurons
- 1: interneurons
- 2: NPCs
- 3: astrocyte
- 4: microglia
- 5: OPCs

|11|12|
|-----------|----------------|
|![plot13_tsn_w_names](https://user-images.githubusercontent.com/88135502/154834305-5fa56d2b-87c8-42c9-bdde-eeae484af874.png)|![plot14_tsne_w_orig ident](https://user-images.githubusercontent.com/88135502/155280958-70ba5581-4041-4b8b-9806-875d5602ba92.png)|

Compare our results with `orig.ident` and we can see that NPCs, excitatory neurons, interneurons are clustered in a development-dependently, as noted in the paper.

#### 2.2.4 Trajectory and Pseudotime analysis 
To elucidate the development order of the cells, trajectory and pseudotime analysis is performed by Monocle 2 and Monocle 3. A Seruat v4 object can be converted into a cell_data_set object of Monocle 3 by the `as.cell_data_set` function. To do a trajectory analysis, operate a sequence of the functions `preprocess_cds`, `reduce_dimension`, `cluster_cells`, `learn_graph`. The plot below is obtained. Note that `order_cells` is not applied, which implies that pseudotime is not yet calculated.

|13|14|
|-----------|----------------|
|![plot15_tractory_wo_ordercell_ident](https://user-images.githubusercontent.com/88135502/155466446-43cb59d7-3ed7-47f2-bfbb-2315774d2e8e.png)|![plot16_trajectory_wo_ordercell_orig ident](https://user-images.githubusercontent.com/88135502/155466462-1c700780-26d9-4047-aa1d-085c92ff50b2.png)|

Note that black circles are branch points and gray circles are leaf points. Starting from Branch 1 (cell type is unclear), there are three major pathways
- Branch 1 -> excitatory neurons
- Branch 1 -> NPCs -> OPCs
- Branch 1 -> Interneurons

From here we decided to operate the whole analysis from scratch by Monocle 2 because of the following
- `plot_cell_trajectory` and `plot_complex_cell_trajectory` in Monocle 2 has no equivalent in Monocle 3
- The researchers used Monocle 2
- Researchers removed microglia and intneurons during pseudotime analysis, so all calculations must be undone. 
- Hard to pinpoint starting point for psuedotime with interactive system of Monocle3 

We use the `counts` of our Seurat object and transfer it to a cell_data_set object by using the `newCellDataSet` function (it would be reasonable to use the normalzied data, but raw data shows better results). Then remove the microglial cells and interneurons. Microglial cellsare mesoderm-derived cells and interneurons are generated from ganglia, not the prefrontal cortex, so they are removed. Then use the marker genes found in Seurat, and use `setOrderingFilter`, `reduceDimension`, `orderCells`. Then plot the results by `plot_complex_cell_trajectory`.

|15|16|17|
|-----------|----------------|----------------|
|![plot17_trajectory_tree_celltype](https://user-images.githubusercontent.com/88135502/155845978-75216e89-5db1-4d96-959d-468f6495b183.png)|![plot18_trajectory_tree_states](https://user-images.githubusercontent.com/88135502/155845987-34be6ff7-cac9-495f-8c60-d35943eec0fc.png)|![plot19_trajectory_tree_orig ident](https://user-images.githubusercontent.com/88135502/155845996-66411e32-bb8a-40d1-8b42-08eae672368e.png)|

These figures are a reprodution of Figure 1 of the paper. On figure 17, we can see that the trajetory is formed according to gestation period. GW8 to GW10 are on the top, and GW23 to GW26 are on the bottom of the trajectory. However, figure  15 shows some error. The trejectory shows that OPCs and astrocytes are formed after excitatory neurons. In other words, glial cells form subsequently of neuronal cells, according to our result. Analysis of thie erroneous result will be shown next section.

Known mark genes were plotted by pseudotime to expression. The plot is shown below. It was done by the `plot_multiple_branches_pseudotime` function. 

![plot20_gene_psuedotime](https://user-images.githubusercontent.com/88135502/155853044-0465aaaf-920b-4453-805b-30e8046311dc.png)

The expression pattern shown in the figure is equivalent to figure 1e except for AQP4.
- AQP4: In our figure there is decrease of expression on neural cells, and no expression on glial cells. The original figure shows no expression on neuronal cells and increasing expresion on glial cells.
- NEUROD2: Rapid increase of expression on neural cells, and some incerease on glial cells. Similar with our original figure. 
- OLIG1: Slight increase of expression on glial cells. While the shape is similar to the original, the amount of increase is significantly lower.
- PAX6: Gradual decrease of expression on both groups. Similar with our original figure. 
- SOX2: Neuronal cells show decrease, and glial cells show increase then decrease. In the original, glial expression is relatively consistent.
