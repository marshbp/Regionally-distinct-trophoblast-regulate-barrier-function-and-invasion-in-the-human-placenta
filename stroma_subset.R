library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tibble)

#Generation of the stromal cell subset of the integrated dataset used in Figure 1 from
#Regionally distinct trophoblast regulate barrier function and invasion in the human placenta
#https://www.biorxiv.org/content/10.1101/2022.03.21.485195v1

#Processed matrices (barcodes.tsv.gz, feature.tsv.gz, and matrix.mtx.gz) for each sample can be downloaded here: 
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198373

#Processed data as R.obj can be downloaded here:
#https://figshare.com/projects/Regionally_distinct_trophoblast_regulate_barrier_function_and_invasion_in_the_human_placenta/135191

setwd("/Volumes/KINGSTON/region_coarse_cluster_subsets/stroma/")

#Begin with the file from either 1) Generation of Intrgrated Dataset from 10x matrices.R 
#or 2) Integrated_dataset.Robj downloaded from 
#https://figshare.com/projects/Regionally_distinct_trophoblast_regulate_barrier_function_and_invasion_in_the_human_placenta/135191
load(file = "/Volumes/KINGSTON/region_integration_final/cyto.combined_final_clusters.Robj")

#Subset the cells (from Tim Stuart (Post-doc in Satija Lab) - https://github.com/satijalab/seurat/issues/2087)
stroma_int <- SubsetData(cyto.combined_final_clusters, 
                         ident.use = c("Mesenchyme1_int", "Mesenchyme2_int", "Mesenchyme3_int",
                                       "Mesenchyme4_int", "Endothelial1_int"))
#Add in GA information
GA <- factor(gsub(".*_","",stroma_int@meta.data$orig.ident), levels = c("17.6", "18.2", "23", "24"))
stroma_int@meta.data$GA <- GA

#Find variable features in the subset using the RNA assay
DefaultAssay(stroma_int) <- "RNA"
stroma_int <- FindVariableFeatures(stroma_int, selection.method = "vst", nfeatures = 2000)

#Run ScaleData on the integrated assay on the new set of variable features
DefaultAssay(stroma_int) <- "integrated"
stroma_int <- ScaleData(stroma_int)

#Run PCA on the integrated assay using the new set of variable features
stroma_int <- RunPCA(stroma_int, verbose = FALSE)

#Create Elbowplot
ElbowPlot(stroma_int, ndims = 50)
#Select number of PCs (https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_clustering_analysis.html)
# Determine percent of variation associated with each PC
pct <- stroma_int@reductions$pca@stdev / sum(stroma_int@reductions$pca@stdev) * 100
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
stroma_int <- RunUMAP(stroma_int, dims = 1:pcs)
stroma_int <- FindNeighbors(stroma_int, reduction = "pca", dims = 1:pcs)
stroma_int <- FindClusters(stroma_int, resolution = 0.4)

#Remove clusters 7 (mis clustered trophoblast - GATA3, CDH1), 8 (low information cells), 11 (mis clustered immune - PTPRC, CD84)
stroma_int <- SubsetData(stroma_int, 
                         ident.use = c(0, 1, 2, 3, 4, 5, 6, 9, 10))
#Find variable features in the subset using the RNA assay
DefaultAssay(stroma_int) <- "RNA"
stroma_int <- FindVariableFeatures(stroma_int, selection.method = "vst", nfeatures = 2000)

#Run ScaleData on the integrated assay on the new set of variable features
DefaultAssay(stroma_int) <- "integrated"
stroma_int <- ScaleData(stroma_int)

#Run PCA on the integrated assay using the new set of variable features
stroma_int <- RunPCA(stroma_int, verbose = FALSE)

#Create Elbowplot
ElbowPlot(stroma_int, ndims = 50)
#Select number of PCs (https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_clustering_analysis.html)
# Determine percent of variation associated with each PC
pct <- stroma_int@reductions$pca@stdev / sum(stroma_int@reductions$pca@stdev) * 100
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
stroma_int <- RunUMAP(stroma_int, dims = 1:pcs)
stroma_int <- FindNeighbors(stroma_int, reduction = "pca", dims = 1:pcs)
stroma_int <- FindClusters(stroma_int, resolution = 0.4)
DimPlot(stroma_int, reduction = "umap", label = T)
DimPlot(stroma_int, reduction = "umap", split.by = "Region")
DimPlot(stroma_int, reduction = "umap", split.by = "orig.ident")

new.cluster.ids <- c("Decidual Stroma 1 (M)", "Decidual Stroma 3 (M)", "Mesenchyme 3 (F)", "Decidual Stroma 2 (M)", "Mesenchyme 1 (F)",
                     "Endothelial Lymph (M)", "Mesenchyme 2 (F/M)", "Mesenchyme 4 (F)", "Cycling Mesenchyme (F)", "Endothelial (M)")
names(new.cluster.ids) <- levels(stroma_int)
stroma_int <- RenameIdents(stroma_int, new.cluster.ids)

stroma_int@active.ident <- 
  factor(stroma_int@active.ident,
         levels = c("Decidual Stroma 1 (M)", "Decidual Stroma 2 (M)", "Decidual Stroma 3 (M)",
                    "Mesenchyme 1 (F)", "Mesenchyme 2 (F/M)", "Mesenchyme 3 (F)", "Mesenchyme 4 (F)", "Cycling Mesenchyme (F)",
                    "Endothelial (M)", "Endothelial Lymph (M)"))

DefaultAssay(stroma_int) <- "RNA"
stroma_int_markers <- FindAllMarkers(stroma_int, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Write markers of each cluster (Figure 1 - Source Data 2 - Stroma subset marker genes)
write.csv(stroma_int_markers, file = "stroma_int_9PC_res04_markers.csv")

#Save stroma subset 
#Stroma_subset.Robj - https://figshare.com/projects/Regionally_distinct_trophoblast_regulate_barrier_function_and_invasion_in_the_human_placenta/135191
save(stroma_int, file = "stroma_int.Robj")

