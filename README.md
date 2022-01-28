# Analyzing the Brain Development scRNA-seq Data between Mouse and Human

## 1. Introduction
This is a personal project aiming to identify the differences between the human and mouse brain development using data generated by scRNA-seq. The aim for this project is to:
- Implement scRNA-seq data from GEO to Seurat and Scanpy and process the data
- Use adult/embryo data of human/mouse, and compare
- Endeavor to figure out differences of gene expression, cell type etc

## 2. Reproducing Previous Studies
In advance, we start by reproducing existing studies to get a grip of workflow. The Gene Expression Omnibus (GEO) provides data produced by sequencing and microarray. Data deposited in the database are categorized among the following hierarchy:
- Platform: 
- DataSets:
- Series:
- Samples: 

Here we use data from Series since it provides the raw data from a complete study. We choose a study that analyzes development differences between human and mouse each.

### 2.1 Human Study
#### 2.1.1 Load Data
Darmanis, et al(2015) is deposited as [GSE67835](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835). Here 466 cells were generated by scRNA-seq. We use the raw TAR file and unzip it, revealing 466 individual CSV files. Individual CSV represents one cell, and stores count for each gene, with 22,088 rows. In example is shown below.

[GSM1657871_1772078217.C03.csv](https://github.com/Deserav/mouse-human-develop/files/7948785/GSM1657871_1772078217.C03.csv)

In R, we implement all the CSV files in to one expression matrix by for looping all file names and full join each file. Then make the expression matrix as a single Seurat object. At this point, the overall gene expression matrix is saved as a Seurat object. Below is the file of the whole experssion matrix.

[human_exp.csv](https://github.com/Deserav/mouse-human-develop/files/7955447/human_exp.csv)

Note that the last three rows are named: `no_feature`, `ambiguous` and `alignment_not_unique`. This seems to be in the data because HTSeq was used to convert reads to counts. We therefore remove the last three rows, and then implement as Seurat object.

In case of the reference genome, the authors of the study used the hg19 with the STAR alignment method. However, the gene names that come from the file name are different to that of hg19. It is assumed that the actual names come from HGNC (HUGO Gene Nomenclature Committee).

#### 2.1.2 Visualize and Data Processing
