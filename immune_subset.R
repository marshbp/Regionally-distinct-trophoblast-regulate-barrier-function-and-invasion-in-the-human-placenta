library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tibble)

#Generation of the immune cell subset of the integrated dataset used in Figure 1 from
#Regionally distinct trophoblast regulate barrier function and invasion in the human placenta
#https://www.biorxiv.org/content/10.1101/2022.03.21.485195v1

#Processed matrices (barcodes.tsv.gz, feature.tsv.gz, and matrix.mtx.gz) for each sample can be downloaded here: 
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198373

#Processed data as R.obj can be downloaded here:
#https://figshare.com/projects/Regionally_distinct_trophoblast_regulate_barrier_function_and_invasion_in_the_human_placenta/135191

setwd("/Volumes/KINGSTON/region_coarse_cluster_subsets/immune/")

#Begin with the file from either 1) Generation of Intrgrated Dataset from 10x matrices.R 
#or 2) Integrated_dataset.Robj downloaded from 
#https://figshare.com/projects/Regionally_distinct_trophoblast_regulate_barrier_function_and_invasion_in_the_human_placenta/135191
load(file = "/Volumes/KINGSTON/region_integration_final/cyto.combined_final_clusters.Robj") 

#Subset the cells (from Tim Stuart (Post-doc in Satija Lab) - https://github.com/satijalab/seurat/issues/2087)
immune_int <- SubsetData(cyto.combined_final_clusters, 
                         ident.use = c("Macro1_int", "Macro2_int", "Macro3_int", "Macro4_int", 
                                       "Macro5_int", "Macro6_MKI67_int", "NK_Tcell_int", "NK_Tcell2_int", "NK_Tcell3_MKI67_int"))

GA <- factor(gsub(".*_","",immune_int@meta.data$orig.ident), levels = c("17.6", "18.2", "23", "24"))
immune_int@meta.data$GA <- GA
#Find variable features in the subset using the RNA assay
DefaultAssay(immune_int) <- "RNA"
immune_int <- FindVariableFeatures(immune_int, selection.method = "vst", nfeatures = 2000)
#Run ScaleData on the integrated assay on the new set of variable features
DefaultAssay(immune_int) <- "integrated"
immune_int <- ScaleData(immune_int)
#Run PCA on the integrated assay using the new set of variable features
immune_int <- RunPCA(immune_int, verbose = FALSE)
#Create Elbowplot
ElbowPlot(immune_int, ndims = 50)
#Select number of PCs (https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_clustering_analysis.html)
# Determine percent of variation associated with each PC
pct <- immune_int@reductions$pca@stdev / sum(immune_int@reductions$pca@stdev) * 100
# Calculate cumulative percents for each PC
cum <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cum > 90 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1),  decreasing = T)[1] + 1 # last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
pcs <- min(co1, co2) 
pcs
#Run FindNeighbors and FindClusters using the new PC dimensions
immune_int <- RunUMAP(immune_int, dims = 1:pcs)
immune_int <- FindNeighbors(immune_int, reduction = "pca", dims = 1:pcs)
immune_int <- FindClusters(immune_int, resolution = 0.5)
DimPlot(immune_int, reduction = "umap")
DimPlot(immune_int, reduction = "umap", split.by = "Region")
DimPlot(immune_int, reduction = "umap", split.by = "orig.ident")

#Remove cluster14 (mis clustered trophoblast - CDH1+)
immune_int <- SubsetData(immune_int, 
                         ident.use = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15))
#Find variable features in the subset using the RNA assay
DefaultAssay(immune_int) <- "RNA"
immune_int <- FindVariableFeatures(immune_int, selection.method = "vst", nfeatures = 2000)

#Run ScaleData on the integrated assay on the new set of variable features
DefaultAssay(immune_int) <- "integrated"
immune_int <- ScaleData(immune_int)

#Run PCA on the integrated assay using the new set of variable features
immune_int <- RunPCA(immune_int, verbose = FALSE)

#Create Elbowplot
ElbowPlot(immune_int, ndims = 50)
#Select number of PCs (https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_clustering_analysis.html)
# Determine percent of variation associated with each PC
pct <- immune_int@reductions$pca@stdev / sum(immune_int@reductions$pca@stdev) * 100
# Calculate cumulative percents for each PC
cum <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cum > 90 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1),  decreasing = T)[1] + 1 # last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
pcs <- min(co1, co2) 
pcs
#Run FindNeighbors and FindClusters using the new PC dimensions
immune_int <- RunUMAP(immune_int, dims = 1:pcs)
immune_int <- FindNeighbors(immune_int, reduction = "pca", dims = 1:pcs)
immune_int <- FindClusters(immune_int, resolution = 0.4)
DimPlot(immune_int, reduction = "umap", label = T)
DimPlot(immune_int, reduction = "umap", split.by = "Region")
DimPlot(immune_int, reduction = "umap", split.by = "orig.ident")

#Cluster annotation (using Vento-Tormo et al. 2018 Decidual Immune cells as a reference)
new.cluster.ids <- c("Macrophage 2.1 (M)", "Macrophage 2.3 (M)", "Macrophage 2.2 (M)", "Macrophage 1 (M)", "NK (M)",
                     "Macrophage 3 (M)", "Macrophage 2.4 (M)", "T (M)", "Cycling NK (M)", "Dendritic Cell 1 (M)",
                     "Cycling M3 (M)", "Plasma (M)", "Erythrocytes")
names(new.cluster.ids) <- levels(immune_int)
immune_int <- RenameIdents(immune_int, new.cluster.ids)

immune_int@active.ident <- 
  factor(immune_int@active.ident,
         levels = c("Macrophage 1 (M)", "Macrophage 2.1 (M)", "Macrophage 2.2 (M)", "Macrophage 2.3 (M)", "Macrophage 2.4 (M)",
                    "Macrophage 3 (M)", "Cycling M3 (M)", "Dendritic Cell 1 (M)", "NK (M)", "Cycling NK (M)", "T (M)",
                    "Plasma (M)", "Erythrocytes"))
                       
DefaultAssay(immune_int) <- "RNA"
immune_int_markers <- FindAllMarkers(immune_int, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Write markers of each cluster (Figure 1 - Source Data 1 - Immune subset marker genes)
write.csv(immune_int_markers, file = "immune_int_13PC_res04_markers.csv")

#Save immune subset 
#Immune_subset.Robj - https://figshare.com/projects/Regionally_distinct_trophoblast_regulate_barrier_function_and_invasion_in_the_human_placenta/135191
save(immune_int, file = "immune_int.Robj")
