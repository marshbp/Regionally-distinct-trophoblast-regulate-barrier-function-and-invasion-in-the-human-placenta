library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tibble)

#Generation of the trophoblast subset of the integrated dataset used in Figures 2-5 from
#Regionally distinct trophoblast regulate barrier function and invasion in the human placenta
#https://www.biorxiv.org/content/10.1101/2022.03.21.485195v1

#Processed matrices (barcodes.tsv.gz, feature.tsv.gz, and matrix.mtx.gz) for each sample can be downloaded here: 
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198373

#Processed data as R.obj can be downloaded here:
#https://figshare.com/projects/Regionally_distinct_trophoblast_regulate_barrier_function_and_invasion_in_the_human_placenta/135191

setwd("/Volumes/KINGSTON/region_coarse_cluster_subsets/trophoblast/")

#Begin with the file from either 1) Generation of Intrgrated Dataset from 10x matrices.R 
#or 2) Integrated_dataset.Robj downloaded from 
#https://figshare.com/projects/Regionally_distinct_trophoblast_regulate_barrier_function_and_invasion_in_the_human_placenta/135191
load(file = "/Volumes/KINGSTON/region_integration_final/cyto.combined_final_clusters.Robj")

#Subset the cells (from Tim Stuart (Post-doc in Satija Lab) - https://github.com/satijalab/seurat/issues/2087)
trophoblast_int <- SubsetData(cyto.combined_final_clusters, 
                         ident.use = c("CTB1_int", "CTB2_int", "CTB_EVT_int", "EVT1_int", "EVT2_int",
                                       "EVT3_int", "EVT4_int", "EVT5_int", "STB_like_int", "CTB_S-phase_int",
                                       "CTB_G2M-phase_int"))
#Find variable features in the subset using the RNA assay
DefaultAssay(trophoblast_int) <- "RNA"
trophoblast_int <- FindVariableFeatures(trophoblast_int, selection.method = "vst", nfeatures = 2000)
#Run ScaleData on the integrated assay on the new set of variable features
DefaultAssay(trophoblast_int) <- "integrated"
trophoblast_int <- ScaleData(trophoblast_int)
#Run PCA on the integrated assay using the new set of variable features
trophoblast_int <- RunPCA(trophoblast_int, verbose = FALSE)
#Create Elbowplot
ElbowPlot(trophoblast_int, ndims = 50)
#Select number of PCs (https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_clustering_analysis.html)
# Determine percent of variation associated with each PC
pct <- trophoblast_int@reductions$pca@stdev / sum(trophoblast_int@reductions$pca@stdev) * 100
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
trophoblast_int <- RunUMAP(trophoblast_int, dims = 1:pcs)
trophoblast_int <- FindNeighbors(trophoblast_int, reduction = "pca", dims = 1:pcs)
trophoblast_int <- FindClusters(trophoblast_int, resolution = 0.6)
DimPlot(trophoblast_int, reduction = "umap")

#Remove contaminating clusters 11, 13, 14 (VIM+ HLA-A+) and low information cells (12)
trophoblast_int <- SubsetData(trophoblast_int, 
                              ident.use = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15", "16"))
#Find variable features in the subset using the RNA assay
DefaultAssay(trophoblast_int) <- "RNA"
trophoblast_int <- FindVariableFeatures(trophoblast_int, selection.method = "vst", nfeatures = 2000)
#Run ScaleData on the integrated assay on the new set of variable features
DefaultAssay(trophoblast_int) <- "integrated"
trophoblast_int <- ScaleData(trophoblast_int)
#Run PCA on the integrated assay using the new set of variable features
trophoblast_int <- RunPCA(trophoblast_int, verbose = FALSE)
#Create Elbowplot
ElbowPlot(trophoblast_int, ndims = 50)
#Select number of PCs (https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_clustering_analysis.html)
# Determine percent of variation associated with each PC
pct <- trophoblast_int@reductions$pca@stdev / sum(trophoblast_int@reductions$pca@stdev) * 100
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
trophoblast_int <- RunUMAP(trophoblast_int, dims = 1:pcs)
trophoblast_int <- FindNeighbors(trophoblast_int, reduction = "pca", dims = 1:pcs)
trophoblast_int <- FindClusters(trophoblast_int, resolution = 0.6)
DimPlot(trophoblast_int, reduction = "umap", label = T)
DimPlot(trophoblast_int, reduction = "umap", split.by = "Region")
DimPlot(trophoblast_int, reduction = "umap", split.by = "orig.ident")

new.cluster.ids <- c("CTB 1", "EVT 3", "EVT 2", "CTB 2", "CTB 3",
                     "CTB 4", "EVT 4", "EVT 1", "CTB S-phase", "CTB G2M-phase",  
                     "STB Precursor", "STB", "EVT Precursor")
names(new.cluster.ids) <- levels(trophoblast_int)
trophoblast_int <- RenameIdents(trophoblast_int, new.cluster.ids)
trophoblast_int@active.ident <- 
  factor(trophoblast_int@active.ident,
         levels = c("CTB 1", "CTB 2", "CTB 3", "CTB 4", "STB Precursor",
                    "STB", "CTB S-phase", "CTB G2M-phase", "EVT Precursor", "EVT 1",
                    "EVT 2", "EVT 3", "EVT 4"))

DefaultAssay(trophoblast_int) <- "RNA"
trophoblast_int_markers <- FindAllMarkers(trophoblast_int, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Write markers of each cluster (Figure 2 - Source Data 1 - Trophoblast subset marker genes)
write.csv(trophoblast_int_markers, file = "trophoblast_int_15PC_res06_markers.csv")

#Save trophoblast subset 
#Trophoblast_subset.Robj - https://figshare.com/projects/Regionally_distinct_trophoblast_regulate_barrier_function_and_invasion_in_the_human_placenta/135191
save(trophoblast_int, file = "trophoblast_int.Robj")

