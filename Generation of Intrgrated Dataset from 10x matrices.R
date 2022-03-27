library(Seurat)
library(dplyr)
library(tibble)
library(DoubletFinder)
library(cowplot)
library(patchwork)

#Processing of each individual 10x dataset to the integrated dataset used in Figure 1 from
#Regionally distinct trophoblast regulate barrier function and invasion in the human placenta
#https://www.biorxiv.org/content/10.1101/2022.03.21.485195v1

#Processed matrices (barcodes.tsv.gz, feature.tsv.gz, and matrix.mtx.gz) for each sample can be downloaded here: 
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198373

#Processed data as R.obj can be downloaded here:
#https://figshare.com/projects/Regionally_distinct_trophoblast_regulate_barrier_function_and_invasion_in_the_human_placenta/135191

setwd("/Volumes/KINGSTON/")
#Preprocess data before doublet detection for sample 1/8
villi_24_data <- Read10X(data.dir = "region_10x_outs/GW24.0_VC/filtered_feature_bc_matrix/")
villi_24 <- CreateSeuratObject(counts = villi_24_data, project = "villi_24", min.cells = 3, min.features = 200)
rm(villi_24_data)

villi_24[["percent.mt"]] <- PercentageFeatureSet(villi_24, pattern = "^MT-")
villi_24 <- subset(villi_24, subset = nFeature_RNA > 500 & nFeature_RNA < 6500 & percent.mt < 15)

villi_24 <- NormalizeData(villi_24, normalization.method = "LogNormalize", scale.factor = 10000)
villi_24 <- NormalizeData(villi_24)
villi_24 <- FindVariableFeatures(villi_24, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(villi_24)
villi_24 <- ScaleData(villi_24, features = all.genes)
villi_24 <- RunPCA(villi_24, verbose = FALSE)
ElbowPlot(villi_24, ndims = 50) #20 PC
villi_24 <- RunUMAP(villi_24, dims = 1:20)
villi_24 <- FindNeighbors(villi_24, reduction = "pca", dims = 1:20)
villi_24 <- FindClusters(villi_24, resolution = 0.6)
DimPlot(villi_24, reduction = "umap")

#Doublet Detection
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_villi_24 <- paramSweep_v3(villi_24, PCs = 1:20, sct = FALSE)
sweep.stats_villi_24 <- summarizeSweep(sweep.res.list_villi_24, GT = FALSE)
par(mar=c(5.1, 4.1, 4.1, 2.1))
bcmvn_villi_24 <- find.pK(sweep.stats_villi_24) #higest is pK of 0.26

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- villi_24@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- villi_24.cyto.combined@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(villi_24@meta.data))  ## Assuming 7.5% doublet formation rate - from 10x collection targetting 10000 cells
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
doublet_finder_villi_24 <- doubletFinder_v3(villi_24, PCs = 1:20, pN = 0.25, pK = 0.26, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(doublet_finder_villi_24, reduction = "umap", group.by = "DF.classifications_0.25_0.26_167")

villi_24_doublets <- as.data.frame(doublet_finder_villi_24$DF.classifications_0.25_0.26_167)
villi_24_doublets <- tibble::rownames_to_column(villi_24_doublets)
villi_24_doublets <- filter(villi_24_doublets, villi_24_doublets$`doublet_finder_villi_24$DF.classifications_0.25_0.26_167`=="Doublet")
singlets <- setdiff(colnames(doublet_finder_villi_24), villi_24_doublets$rowname)
villi_24_df <- subset(villi_24, cells = singlets)
save(villi_24_df, file = "region_processed_seurats/villi_24_df.Robj") 
rm(list=ls())

####################################################################################################

#Preprocess data before doublet detection for sample 2/8
sc_24_data <- Read10X(data.dir = "region_10x_outs/GW24.0_SC/filtered_feature_bc_matrix/")
sc_24 <- CreateSeuratObject(counts = sc_24_data, project = "sc_24", min.cells = 3, min.features = 200)
rm(sc_24_data)

sc_24[["percent.mt"]] <- PercentageFeatureSet(sc_24, pattern = "^MT-")
sc_24 <- subset(sc_24, subset = nFeature_RNA > 500 & nFeature_RNA < 6500 & percent.mt < 15)
sc_24 <- NormalizeData(sc_24, normalization.method = "LogNormalize", scale.factor = 10000)
sc_24 <- NormalizeData(sc_24)
sc_24 <- FindVariableFeatures(sc_24, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sc_24)
sc_24 <- ScaleData(sc_24, features = all.genes)
sc_24 <- RunPCA(sc_24, verbose = FALSE)
ElbowPlot(sc_24, ndims = 50) #20 PC
sc_24 <- RunUMAP(sc_24, dims = 1:20)
sc_24 <- FindNeighbors(sc_24, reduction = "pca", dims = 1:20)
sc_24 <- FindClusters(sc_24, resolution = 0.6)
DimPlot(sc_24, reduction = "umap")

#Doublet Detection
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_sc_24 <- paramSweep_v3(sc_24, PCs = 1:20, sct = FALSE)
sweep.stats_sc_24 <- summarizeSweep(sweep.res.list_sc_24, GT = FALSE)
par(mar=c(5.1, 4.1, 4.1, 2.1))
bcmvn_sc_24 <- find.pK(sweep.stats_sc_24) #higest are pK of 0.28

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- sc_24@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- sc_24.cyto.combined@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(sc_24@meta.data))  ## Assuming 7.5% doublet formation rate - from 10x collection targetting 10000 cells
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
doublet_finder_sc_24 <- doubletFinder_v3(sc_24, PCs = 1:20, pN = 0.25, pK = 0.28, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(doublet_finder_sc_24, reduction = "umap", group.by = "DF.classifications_0.25_0.28_415")

sc_24_doublets <- as.data.frame(doublet_finder_sc_24$DF.classifications_0.25_0.28_415)
sc_24_doublets <- tibble::rownames_to_column(sc_24_doublets)
sc_24_doublets <- filter(sc_24_doublets, sc_24_doublets$`doublet_finder_sc_24$DF.classifications_0.25_0.28_415`=="Doublet")
singlets <- setdiff(colnames(doublet_finder_sc_24), sc_24_doublets$rowname)
sc_24_df <- subset(sc_24, cells = singlets)
save(sc_24_df, file = "region_processed_seurats/sc_24_df.robj")
rm(list=ls())

####################################################################################################

#Preprocess data before doublet detection for sample 3/8
villi_18.2_data <- Read10X(data.dir = "region_10x_outs/GW18.2_VC/filtered_feature_bc_matrix/")
villi_18.2 <- CreateSeuratObject(counts = villi_18.2_data, project = "villi_18.2", min.cells = 3, min.features = 200)
rm(villi_18.2_data)

villi_18.2[["percent.mt"]] <- PercentageFeatureSet(villi_18.2, pattern = "^MT-")
villi_18.2 <- subset(villi_18.2, subset = nFeature_RNA > 500 & nFeature_RNA < 6500 & percent.mt < 15)
villi_18.2 <- NormalizeData(villi_18.2, normalization.method = "LogNormalize", scale.factor = 10000)
villi_18.2 <- NormalizeData(villi_18.2)
villi_18.2 <- FindVariableFeatures(villi_18.2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(villi_18.2)
villi_18.2 <- ScaleData(villi_18.2, features = all.genes)
villi_18.2 <- RunPCA(villi_18.2, verbose = FALSE)
ElbowPlot(villi_18.2, ndims = 50) #30 PC
villi_18.2 <- RunUMAP(villi_18.2, dims = 1:30)
villi_18.2 <- FindNeighbors(villi_18.2, reduction = "pca", dims = 1:30)
villi_18.2 <- FindClusters(villi_18.2, resolution = 0.5)
DimPlot(villi_18.2, reduction = "umap")

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_villi_18.2 <- paramSweep_v3(villi_18.2, PCs = 1:30, sct = FALSE)
sweep.stats_villi_18.2 <- summarizeSweep(sweep.res.list_villi_18.2, GT = FALSE)
par(mar=c(5.1, 4.1, 4.1, 2.1))
bcmvn_villi_18.2 <- find.pK(sweep.stats_villi_18.2) #higest is pK of 0.16

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- villi_18.2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- villi_18.2.cyto.combined@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(villi_18.2@meta.data))  ## Assuming 7.5% doublet formation rate - from 10x collection targetting 10000 cells
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
doublet_finder_villi_18.2 <- doubletFinder_v3(villi_18.2, PCs = 1:30, pN = 0.25, pK = 0.16, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(doublet_finder_villi_18.2, reduction = "umap", group.by = "DF.classifications_0.25_0.16_580")

villi_18.2_doublets <- as.data.frame(doublet_finder_villi_18.2$DF.classifications_0.25_0.16_580)
villi_18.2_doublets <- tibble::rownames_to_column(villi_18.2_doublets)
villi_18.2_doublets <- filter(villi_18.2_doublets, villi_18.2_doublets$`doublet_finder_villi_18.2$DF.classifications_0.25_0.16_580`=="Doublet")
singlets <- setdiff(colnames(doublet_finder_villi_18.2), villi_18.2_doublets$rowname)
villi_18.2_df <- subset(villi_18.2, cells = singlets)
save(villi_18.2_df, file = "region_processed_seurats/villi_18.2_df.robj")
rm(list=ls())

####################################################################################################

#Preprocess data before doublet detection for sample 4/8
sc_18.2_data <- Read10X(data.dir = "region_10x_outs/GW18.2_SC/filtered_feature_bc_matrix/")
sc_18.2 <- CreateSeuratObject(counts = sc_18.2_data, project = "sc_18.2", min.cells = 3, min.features = 200)
rm(sc_18.2_data)

sc_18.2[["percent.mt"]] <- PercentageFeatureSet(sc_18.2, pattern = "^MT-")
sc_18.2 <- subset(sc_18.2, subset = nFeature_RNA > 500 & nFeature_RNA < 6500 & percent.mt < 15)
sc_18.2 <- NormalizeData(sc_18.2, normalization.method = "LogNormalize", scale.factor = 10000)
sc_18.2 <- NormalizeData(sc_18.2)
sc_18.2 <- FindVariableFeatures(sc_18.2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sc_18.2)
sc_18.2 <- ScaleData(sc_18.2, features = all.genes)
sc_18.2 <- RunPCA(sc_18.2, verbose = FALSE)
ElbowPlot(sc_18.2, ndims = 50) #20 PC
sc_18.2 <- RunUMAP(sc_18.2, dims = 1:20)
sc_18.2 <- FindNeighbors(sc_18.2, reduction = "pca", dims = 1:20)
sc_18.2 <- FindClusters(sc_18.2, resolution = 0.6)
DimPlot(sc_18.2, reduction = "umap")

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_sc_18.2 <- paramSweep_v3(sc_18.2, PCs = 1:20, sct = FALSE)
sweep.stats_sc_18.2 <- summarizeSweep(sweep.res.list_sc_18.2, GT = FALSE)
par(mar=c(5.1, 4.1, 4.1, 2.1))
bcmvn_sc_18.2 <- find.pK(sweep.stats_sc_18.2) #higest are pK of 0.005

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- sc_18.2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- sc_18.2.cyto.combined@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(sc_18.2@meta.data))  ## Assuming 7.5% doublet formation rate - from 10x collection targetting 10000 cells
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
doublet_finder_sc_18.2 <- doubletFinder_v3(sc_18.2, PCs = 1:20, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(doublet_finder_sc_18.2, reduction = "umap", group.by = "DF.classifications_0.25_0.005_513")

sc_18.2_doublets <- as.data.frame(doublet_finder_sc_18.2$DF.classifications_0.25_0.005_513)
sc_18.2_doublets <- tibble::rownames_to_column(sc_18.2_doublets)
sc_18.2_doublets <- filter(sc_18.2_doublets, sc_18.2_doublets$`doublet_finder_sc_18.2$DF.classifications_0.25_0.005_513`=="Doublet")
singlets <- setdiff(colnames(doublet_finder_sc_18.2), sc_18.2_doublets$rowname)
sc_18.2_df <- subset(sc_18.2, cells = singlets)
save(sc_18.2_df, file = "region_processed_seurats/sc_18.2_df.robj")
rm(list=ls())

####################################################################################################

#Preprocess data before doublet detection for sample 5/8
villi_17.6_data <- Read10X(data.dir = "region_10x_outs/GW17.6_VC/filtered_feature_bc_matrix/")
villi_17.6 <- CreateSeuratObject(counts = villi_17.6_data, project = "villi_17.6", min.cells = 3, min.features = 200)
rm(villi_17.6_data)
villi_17.6[["percent.mt"]] <- PercentageFeatureSet(villi_17.6, pattern = "^MT-")
villi_17.6 <- subset(villi_17.6, subset = nFeature_RNA > 500 & nFeature_RNA < 6500 & percent.mt < 15)
villi_17.6 <- NormalizeData(villi_17.6, normalization.method = "LogNormalize", scale.factor = 10000)
villi_17.6 <- NormalizeData(villi_17.6)
villi_17.6 <- FindVariableFeatures(villi_17.6, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(villi_17.6)
villi_17.6 <- ScaleData(villi_17.6, features = all.genes)
villi_17.6 <- RunPCA(villi_17.6, verbose = FALSE)
ElbowPlot(villi_17.6, ndims = 50) #21 PC
villi_17.6 <- RunUMAP(villi_17.6, dims = 1:21)
villi_17.6 <- FindNeighbors(villi_17.6, reduction = "pca", dims = 1:21)
villi_17.6 <- FindClusters(villi_17.6, resolution = 0.6)
DimPlot(villi_17.6, reduction = "umap")

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_villi_17.6 <- paramSweep_v3(villi_17.6, PCs = 1:21, sct = FALSE)
sweep.stats_villi_17.6 <- summarizeSweep(sweep.res.list_villi_17.6, GT = FALSE)
par(mar=c(5.1, 4.1, 4.1, 2.1))
bcmvn_villi_17.6 <- find.pK(sweep.stats_villi_17.6) #higest are pK of 0.3 then .04

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- villi_17.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- villi_17.6.cyto.combined@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(villi_17.6@meta.data))  ## Assuming 7.5% doublet formation rate - from 10x collection targetting 10000 cells
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
doublet_finder_villi_17.6 <- doubletFinder_v3(villi_17.6, PCs = 1:21, pN = 0.25, pK = 0.04, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(doublet_finder_villi_17.6, reduction = "umap", group.by = "DF.classifications_0.25_0.04_499")

villi_17.6_doublets <- as.data.frame(doublet_finder_villi_17.6$DF.classifications_0.25_0.04_499)
villi_17.6_doublets <- tibble::rownames_to_column(villi_17.6_doublets)
villi_17.6_doublets <- filter(villi_17.6_doublets, villi_17.6_doublets$`doublet_finder_villi_17.6$DF.classifications_0.25_0.04_499`=="Doublet")
singlets <- setdiff(colnames(doublet_finder_villi_17.6), villi_17.6_doublets$rowname)
villi_17.6_df <- subset(villi_17.6, cells = singlets)
save(villi_17.6_df, file = "region_processed_seurats/villi_17.6_df.robj")
rm(list=ls())

####################################################################################################

#Preprocess data before doublet detection for sample 6/8
sc_17.6_data <- Read10X(data.dir = "region_10x_outs/GW17.6_SC/filtered_feature_bc_matrix/")
sc_17.6 <- CreateSeuratObject(counts = sc_17.6_data, project = "sc_17.6", min.cells = 3, min.features = 200)
rm(sc_17.6_data)

sc_17.6[["percent.mt"]] <- PercentageFeatureSet(sc_17.6, pattern = "^MT-")
sc_17.6 <- subset(sc_17.6, subset = nFeature_RNA > 500 & nFeature_RNA < 6500 & percent.mt < 15)
sc_17.6 <- NormalizeData(sc_17.6, normalization.method = "LogNormalize", scale.factor = 10000)
sc_17.6 <- NormalizeData(sc_17.6)
sc_17.6 <- FindVariableFeatures(sc_17.6, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sc_17.6)
sc_17.6 <- ScaleData(sc_17.6, features = all.genes)
sc_17.6 <- RunPCA(sc_17.6, verbose = FALSE)
ElbowPlot(sc_17.6, ndims = 50) #20 PC
sc_17.6 <- RunUMAP(sc_17.6, dims = 1:20)
sc_17.6 <- FindNeighbors(sc_17.6, reduction = "pca", dims = 1:20)
sc_17.6 <- FindClusters(sc_17.6, resolution = 0.6)
DimPlot(sc_17.6, reduction = "umap")

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_sc_17.6 <- paramSweep_v3(sc_17.6, PCs = 1:20, sct = FALSE)
sweep.stats_sc_17.6 <- summarizeSweep(sweep.res.list_sc_17.6, GT = FALSE)
par(mar=c(5.1, 4.1, 4.1, 2.1))
bcmvn_sc_17.6 <- find.pK(sweep.stats_sc_17.6) #higest are pK of 0.27

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- sc_17.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- sc_17.6.cyto.combined@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(sc_17.6@meta.data))  ## Assuming 7.5% doublet formation rate - from 10x collection targetting 10000 cells
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
doublet_finder_sc_17.6 <- doubletFinder_v3(sc_17.6, PCs = 1:20, pN = 0.25, pK = 0.27, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(doublet_finder_sc_17.6, reduction = "umap", group.by = "DF.classifications_0.25_0.27_485")

sc_17.6_doublets <- as.data.frame(doublet_finder_sc_17.6$DF.classifications_0.25_0.27_485)
sc_17.6_doublets <- tibble::rownames_to_column(sc_17.6_doublets)
sc_17.6_doublets <- filter(sc_17.6_doublets, sc_17.6_doublets$`doublet_finder_sc_17.6$DF.classifications_0.25_0.27_485`=="Doublet")
singlets <- setdiff(colnames(doublet_finder_sc_17.6), sc_17.6_doublets$rowname)
sc_17.6_df <- subset(sc_17.6, cells = singlets)
save(sc_17.6_df, file = "region_processed_seurats/sc_17.6_df.robj")
rm(list=ls())

####################################################################################################

#Preprocess data before doublet detection for sample 7/8
villi_23_data <- Read10X(data.dir = "region_10x_outs/GW23.0_VC/filtered_feature_bc_matrix/")
villi_23 <- CreateSeuratObject(counts = villi_23_data, project = "villi_23", min.cells = 3, min.features = 200)
rm(villi_23_data)

villi_23[["percent.mt"]] <- PercentageFeatureSet(villi_23, pattern = "^MT-")
villi_23 <- subset(villi_23, subset = nFeature_RNA > 500 & nFeature_RNA < 6500 & percent.mt < 15)
villi_23 <- NormalizeData(villi_23, normalization.method = "LogNormalize", scale.factor = 10000)
villi_23 <- NormalizeData(villi_23)
villi_23 <- FindVariableFeatures(villi_23, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(villi_23)
villi_23 <- ScaleData(villi_23, features = all.genes)
villi_23 <- RunPCA(villi_23, verbose = FALSE)
ElbowPlot(villi_23, ndims = 50) #20 PC
villi_23 <- RunUMAP(villi_23, dims = 1:20)
villi_23 <- FindNeighbors(villi_23, reduction = "pca", dims = 1:20)
villi_23 <- FindClusters(villi_23, resolution = 0.6)
DimPlot(villi_23, reduction = "umap")

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_villi_23 <- paramSweep_v3(villi_23, PCs = 1:20, sct = FALSE)
sweep.stats_villi_23 <- summarizeSweep(sweep.res.list_villi_23, GT = FALSE)
par(mar=c(5.1, 4.1, 4.1, 2.1))
bcmvn_villi_23 <- find.pK(sweep.stats_villi_23) #higest is pK of 0.005

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- villi_23@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- villi_23.cyto.combined@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(villi_23@meta.data))  ## Assuming 7.5% doublet formation rate - from 10x collection targetting 10000 cells
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
doublet_finder_villi_23 <- doubletFinder_v3(villi_23, PCs = 1:20, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(doublet_finder_villi_23, reduction = "umap", group.by = "DF.classifications_0.25_0.005_811")

villi_23_doublets <- as.data.frame(doublet_finder_villi_23$DF.classifications_0.25_0.005_811)
villi_23_doublets <- tibble::rownames_to_column(villi_23_doublets)
villi_23_doublets <- filter(villi_23_doublets, villi_23_doublets$`doublet_finder_villi_23$DF.classifications_0.25_0.005_811`=="Doublet")
singlets <- setdiff(colnames(doublet_finder_villi_23), villi_23_doublets$rowname)
villi_23_df <- subset(villi_23, cells = singlets)
save(villi_23_df, file = "region_processed_seurats/villi_23_df.robj") 
rm(list=ls())

####################################################################################################

#Preprocess data before doublet detection for sample 8/8
sc_23_data <- Read10X(data.dir = "region_10x_outs/GW23.0_SC/filtered_feature_bc_matrix/")
sc_23 <- CreateSeuratObject(counts = sc_23_data, project = "sc_23", min.cells = 3, min.features = 200)
rm(sc_23_data)

sc_23[["percent.mt"]] <- PercentageFeatureSet(sc_23, pattern = "^MT-")
sc_23 <- subset(sc_23, subset = nFeature_RNA > 500 & nFeature_RNA < 6500 & percent.mt < 15)
sc_23 <- NormalizeData(sc_23, normalization.method = "LogNormalize", scale.factor = 10000)
sc_23 <- NormalizeData(sc_23)
sc_23 <- FindVariableFeatures(sc_23, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sc_23)
sc_23 <- ScaleData(sc_23, features = all.genes)
sc_23 <- RunPCA(sc_23, verbose = FALSE)
ElbowPlot(sc_23, ndims = 50) #20 PC
sc_23 <- RunUMAP(sc_23, dims = 1:20)
sc_23 <- FindNeighbors(sc_23, reduction = "pca", dims = 1:20)
sc_23 <- FindClusters(sc_23, resolution = 0.6)
DimPlot(sc_23, reduction = "umap")

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_sc_23 <- paramSweep_v3(sc_23, PCs = 1:20, sct = FALSE)
sweep.stats_sc_23 <- summarizeSweep(sweep.res.list_sc_23, GT = FALSE)
par(mar=c(5.1, 4.1, 4.1, 2.1))
bcmvn_sc_23 <- find.pK(sweep.stats_sc_23) #higest are pK of 0.23

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- sc_23@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- sc_23.cyto.combined@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(sc_23@meta.data))  ## Assuming 7.5% doublet formation rate - from 10x collection targetting 10000 cells
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
doublet_finder_sc_23 <- doubletFinder_v3(sc_23, PCs = 1:20, pN = 0.25, pK = 0.23, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(doublet_finder_sc_23, reduction = "umap", group.by = "DF.classifications_0.25_0.23_625")

sc_23_doublets <- as.data.frame(doublet_finder_sc_23$DF.classifications_0.25_0.23_625)
sc_23_doublets <- tibble::rownames_to_column(sc_23_doublets)
sc_23_doublets <- filter(sc_23_doublets, sc_23_doublets$`doublet_finder_sc_23$DF.classifications_0.25_0.23_625`=="Doublet")
singlets <- setdiff(colnames(doublet_finder_sc_23), sc_23_doublets$rowname)
sc_23_df <- subset(sc_23, cells = singlets)
save(sc_23_df, file = "region_processed_seurats/sc_23_df.robj")


##############################
#Perform integration of all datasets using Seurat v3

#Find anchors
cyto.anchors <- FindIntegrationAnchors(object.list = list(villi_24_df, villi_18.2_df, villi_17.6_df, villi_23_df, 
                                                          sc_24_df, sc_18.2_df, sc_17.6_df, sc_23_df), dims = 1:20)

#Integrate
cyto.combined <- IntegrateData(anchorset = cyto.anchors, dims = 1:20)

DefaultAssay(cyto.combined) <- "integrated"
cyto.combined <- ScaleData(cyto.combined, verbose = FALSE)
cyto.combined <- RunPCA(cyto.combined, npcs = 50, verbose = FALSE)
tiff(file = "cyto.combined.elbow.plot.tiff", width = 5000, height = 5000, units = "px", res = 800)
ElbowPlot(cyto.combined, ndims = 50)
dev.off()

cyto.combined <- RunUMAP(cyto.combined, reduction = "pca", dims = 1:30)
cyto.combined <- FindNeighbors(cyto.combined, reduction = "pca", dims = 1:30)
cyto.combined <- FindClusters(cyto.combined, resolution = 0.6)

#Rename according to primary cell type annotations
new.cluster.ids <- c("CTB1_int", "EVT1_int", "CTB2_int", "Macro1_int", "EVT2_int",
                     "Macro2_int", "Macro3_int", "EVT3_int", "CTB_S-phase_int", "Macro4_int",  
                     "Mesenchyme1_int", "Mesenchyme2_int", "NK_Tcell_int", "CTB_G2M-phase_int", "STB_like_int",
                     "Endothelial1_int", "Macro5_int", "NK_Tcell2_int", "NK_Tcell3_MKI67_int", "19_int", 
                     "Macro6_MKI67_int", "EVT4_int", "CTB_EVT_int", "EVT5_int", "Mesenchyme3_int",
                     "Mesenchyme4_int", "26_int")

names(new.cluster.ids) <- levels(cyto.combined)
cyto.combined <- RenameIdents(cyto.combined, new.cluster.ids)

#add region metadata
Region <- factor(gsub("_.*","",cyto.combined@meta.data$orig.ident), levels = c("villi", "sc"))
cyto.combined@meta.data$Region <- Region

#add coarse cluster metadata
coarse_cluster <- c()
a <- c()
for (i in 1:nrow(cyto.combined@meta.data)) {
  cluster <- cyto.combined@meta.data[i, "seurat_clusters"]
  if(cluster %in% c(0, 2, 8, 13, 14)) {
    a <- "CTB"
    coarse_cluster <- append(coarse_cluster, a, after = length(coarse_cluster))
  } else if (cluster %in% c(1, 4, 7, 22, 23)) {
    a <- "EVT"
    coarse_cluster <- append(coarse_cluster, a, after = length(coarse_cluster))
  } else if (cluster %in% c(3, 5, 6, 9, 12, 16, 17, 18, 20)) {
    a <- "Immune"
    coarse_cluster <- append(coarse_cluster, a, after = length(coarse_cluster))
  } else if (cluster %in% c(10, 11, 15, 24, 25)) {
    a <- "Stroma"
    coarse_cluster <- append(coarse_cluster, a, after = length(coarse_cluster))
  } else if (cluster %in% c(19, 26)) {
    a <- "Uterine Epithelium"
    coarse_cluster <- append(coarse_cluster, a, after = length(coarse_cluster))
  } else {
    a <- "Doublets" #Cluster 21 - EVT4 HLA-G, HLA-A, and VIM triple positive
    coarse_cluster <- append(coarse_cluster, a, after = length(coarse_cluster))
  }
}

cyto.combined@meta.data$coarse_cluster <- coarse_cluster

#factor metadata
cyto.combined$orig.ident <- factor(cyto.combined$orig.ident, levels = c("villi_17.6", "sc_17.6", "villi_18.2", "sc_18.2",
                                                                        "villi_23", "sc_23", "villi_24", "sc_24"))

cyto.combined@active.ident <- 
  factor(cyto.combined@active.ident,
         levels = c("CTB1_int", "CTB2_int", "CTB_EVT_int", "EVT1_int", "EVT2_int",
                    "EVT3_int", "EVT4_int", "EVT5_int", "STB_like_int", "CTB_S-phase_int",
                    "CTB_G2M-phase_int", "Macro1_int", "Macro2_int", "Macro3_int", "Macro4_int", 
                    "Macro5_int", "Macro6_MKI67_int", "NK_Tcell_int", "NK_Tcell2_int", "NK_Tcell3_MKI67_int",
                    "Mesenchyme1_int", "Mesenchyme2_int", "Mesenchyme3_int", "Mesenchyme4_int", "Endothelial1_int",
                    "19_int", "26_int"))

cyto.combined@meta.data$coarse_cluster <- factor(cyto.combined@meta.data$coarse_cluster, levels = c("CTB", "EVT", "Immune", "Stroma", "Uterine Epithelium", "Doublets"))

#update cluster ids
new.cluster.ids <- c("CTB1_int", "CTB2_int", "CTB_EVT_int", "EVT1_int", "EVT2_int",           
                     "EVT3_int", "EVT4_int", "EVT5_int", "STB_like_int", "CTB_S-phase_int",    
                     "CTB_G2M-phase_int", "Macro1_int", "Macro2_int", "Macro3_int", "Macro4_int",         
                     "Macro5_int", "Macro6_MKI67_int", "NK_Tcell_int", "NK_Tcell2_int", "NK_Tcell3_MKI67_int",
                     "Mesenchyme1_int", "Mesenchyme2_int", "Mesenchyme3_int", "Mesenchyme4_int", "Endothelial1_int",   
                     "Glandular_Epithelium1_int", "Glandular_Epithelium2_int" )

names(new.cluster.ids) <- levels(cyto.combined)
cyto.combined <- RenameIdents(cyto.combined, new.cluster.ids)
cyto.combined_final_clusters <- cyto.combined
save(cyto.combined_final_clusters, file = "cyto.combined_final_clusters.Robj")

#Final file can be found here: 
#https://figshare.com/projects/Regionally_distinct_trophoblast_regulate_barrier_function_and_invasion_in_the_human_placenta/135191
#under integrated dataset


