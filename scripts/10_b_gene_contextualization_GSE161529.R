# --- 1. Setup ---
# Ensure all necessary packages are installed. Uncomment the following lines to install them if needed.
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("Seurat")
# BiocManager::install("tidyverse")
# BiocManager::install("patchwork")
# BiocManager::install("devtools")
# if (!requireNamespace("presto", quietly = TRUE)) devtools::install_github('immunogenomics/presto')
# install.packages("clustree")
# install.packages("RColorBrewer")
# install.packages("pheatmap")
# install.packages("extrafont")  # <-- needed for Arial on Windows

# Load libraries
library(Seurat)
library(tidyverse)
library(patchwork)
library(future)
library(clustree)
library(RColorBrewer)
library(pheatmap)
library(extrafont)

# =================================================================
# ARIAL FONT SETUP
# Run font_import() once the very first time you use this script,
# then loadfonts() on every subsequent run.
# Comment out font_import() after the first run — it is slow.
# =================================================================
# font_import()          # <-- Run ONCE only, then comment out
loadfonts(device = "win")  # Use device = "pdf" if saving PDFs too

# Set Arial as the global default font for all ggplot2 plots
theme_set(
  theme_classic(base_family = "Arial") +
    theme(
      text         = element_text(family = "Arial"),
      plot.title   = element_text(family = "Arial"),
      axis.title   = element_text(family = "Arial"),
      axis.text    = element_text(family = "Arial"),
      legend.text  = element_text(family = "Arial"),
      legend.title = element_text(family = "Arial"),
      strip.text   = element_text(family = "Arial")
    )
)

options(future.globals.maxSize = 30000 * 1024^2)
plan("multisession", workers = 1)

# --- Analysis Suffix ---
analysis_suffix <- "_NoPreneo_AllCells-FINAL"

# =================================================================
# GENE LIST — DEFINED ONCE HERE, USED EVERYWHERE
# PAK1 duplicate removed; unique() guards against any future duplicates.
# =================================================================
genes_you_want <- unique(c(
  "CENPF", "COL10A1", "COL11A1", "ENAH", "FAP", "FN1", "MMP13", "NAT1", "POSTN", "PYCARD", "RAB25", "SMC4", "SPP1", "TMEM132A",
  "ABCA8", "ADH1B", "ADIPOQ", "C7", "CFD", "CHRDL1", "CLDN5", "DCN", "DPT", "EGR1", "FABP4", "FGF2", "FHL1", "FOSB", "GAS1",
  "GRIA1", "HLF", "IRS2", "LAMA2", "MEOX1", "MFAP4", "MME", "SRPX", "ACKR1", "ACTG2", "FCGBP", "MYH11", "SELP", "TMEM132C",
  "TNN", "UBE4A", "VANGL2", "CRABP2", "NEK2", "SLC11A1", "LAD1", "CD44", "H2BC3", "MAPK11", "EPYC", "ADCY5"
))

cat(paste0("--- Total genes in analysis list: ", length(genes_you_want), " ---\n"))
cat(paste(genes_you_want, collapse = ", "), "\n")

# =================================================================
# SINGLE OUTPUT DIRECTORY — ALL PLOTS GO HERE
# =================================================================
pub_dir <- "C:/Users/BITS/Documents/27-03-26 sc rna images"
dir.create(pub_dir, recursive = TRUE, showWarnings = FALSE)
cat(paste0("--- ALL plots will be saved to: ", pub_dir, " ---\n"))

# --- Helper: saves TIFF (600 DPI, LZW compressed) with Arial font ---
save_pub_plot <- function(filename_base, plot_obj, width = 10, height = 8) {
  ggsave(
    file.path(pub_dir, paste0(filename_base, ".tiff")),
    plot       = plot_obj,
    width      = width,
    height     = height,
    dpi        = 600,
    units      = "in",
    device     = "tiff",
    compression = "lzw"   # LZW = lossless compression, keeps file sizes manageable
  )
}

# --- Helper: saves pheatmap objects as TIFF (600 DPI) with Arial font ---
save_pub_pheatmap <- function(filename_base, pheatmap_obj, width = 10, height = 8) {
  tiff(
    file.path(pub_dir, paste0(filename_base, ".tiff")),
    width       = width,
    height      = height,
    units       = "in",
    res         = 600,
    compression = "lzw",
    family      = "Arial"   # Arial font for pheatmap text
  )
  print(pheatmap_obj)
  dev.off()
}


# --- 2. Load Metadata and Define Exclusions ---
metadata <- read.csv("D:\\Documents Folder\\scrna\\metadata.csv") 
preneoplastic_samples_to_exclude <- c("GSM4909277", "GSM4909278", "GSM4909279", "GSM4909280")

# --- 3. Data Loading and Initial Object Creation ---
add_stage_metadata <- function(seurat_obj, metadata) {
  sample_name <- seurat_obj@meta.data$orig.ident[1]
  sample_info <- dplyr::filter(metadata, GSM.ID == sample_name)
  if (nrow(sample_info) > 0) {
    seurat_obj@meta.data$Tumor_Stage <- sample_info$Tumor.Stage
  } else {
    warning(paste("No metadata found for sample:", sample_name))
  }
  return(seurat_obj)
}

main_dir <- "D:\\Documents Folder\\scrna\\GSE161529" 
all_sub_dirs <- list.dirs(main_dir, recursive = FALSE)
sub_dirs <- all_sub_dirs[!basename(all_sub_dirs) %in% preneoplastic_samples_to_exclude]

if (length(sub_dirs) == 0) stop("No sample directories found after excluding preneoplastic samples.")

seurat_list <- list()
for (dir_path in sub_dirs) {
  subfolder_name <- basename(dir_path)
  cat("Processing sample:", subfolder_name, "\n")
  if (!any(file.exists(file.path(dir_path, "matrix.mtx.gz")), file.exists(file.path(dir_path, "matrix.mtx")))) {
    warning(paste("Matrix file not found in", dir_path, ". Skipping."))
    next
  }
  data <- Read10X(data.dir = dir_path)
  seurat_list[[subfolder_name]] <- CreateSeuratObject(counts = data, project = subfolder_name)
  cat("Created Seurat object for:", subfolder_name, "\n")
}

if (length(seurat_list) == 0) stop("No Seurat objects were created.")

# --- 4. Pre-processing and QC (Pre-filtering) ---
seurat_list <- lapply(seurat_list, function(obj) add_stage_metadata(obj, metadata))
seurat_list <- lapply(seurat_list, function(seurat_obj) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  return(seurat_obj)
})

lapply(names(seurat_list), function(sample_name) {
  seurat_obj <- seurat_list[[sample_name]]
  vln_plot_pre_filter <- VlnPlot(seurat_obj,
                                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                                 ncol = 3) +
    plot_annotation(title = paste("Pre-filter QC for", sample_name))
  save_pub_plot(paste0(sample_name, "_pre_filter_qc"), vln_plot_pre_filter, width = 12, height = 6)
  
  feature_scatter <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    ggtitle(paste("Feature Scatter for", sample_name))
  save_pub_plot(paste0(sample_name, "_feature_scatter"), feature_scatter, width = 8, height = 6)
})

# --- 5. Filtering ---
filter_and_subset <- function(obj, nFeature_min, nFeature_max, mt_max) {
  if (!is.null(obj)) subset(obj, subset = nFeature_RNA > nFeature_min & nFeature_RNA < nFeature_max & percent.mt < mt_max) else NULL
}

if ("GSM4909253" %in% names(seurat_list)) seurat_list[["GSM4909253"]] <- filter_and_subset(seurat_list[["GSM4909253"]], 500, 6500, 20)
if ("GSM4909254" %in% names(seurat_list)) seurat_list[["GSM4909254"]] <- filter_and_subset(seurat_list[["GSM4909254"]], 500, 7500, 20)
if ("GSM4909257" %in% names(seurat_list)) seurat_list[["GSM4909257"]] <- filter_and_subset(seurat_list[["GSM4909257"]], 500, 7000, 30)
if ("GSM4909261" %in% names(seurat_list)) seurat_list[["GSM4909261"]] <- filter_and_subset(seurat_list[["GSM4909261"]], 500, 8000, 20)
if ("GSM4909281" %in% names(seurat_list)) seurat_list[["GSM4909281"]] <- filter_and_subset(seurat_list[["GSM4909281"]], 500, 7000, 20)
if ("GSM4909282" %in% names(seurat_list)) seurat_list[["GSM4909282"]] <- filter_and_subset(seurat_list[["GSM4909282"]], 500, 5000, 20)
if ("GSM4909283" %in% names(seurat_list)) seurat_list[["GSM4909283"]] <- filter_and_subset(seurat_list[["GSM4909283"]], 500, 7000, 20)
if ("GSM4909284" %in% names(seurat_list)) seurat_list[["GSM4909284"]] <- filter_and_subset(seurat_list[["GSM4909284"]], 500, 5500, 30)
if ("GSM4909285" %in% names(seurat_list)) seurat_list[["GSM4909285"]] <- filter_and_subset(seurat_list[["GSM4909285"]], 500, 5500, 20)
if ("GSM4909286" %in% names(seurat_list)) seurat_list[["GSM4909286"]] <- filter_and_subset(seurat_list[["GSM4909286"]], 500, 5000, 30)
if ("GSM4909287" %in% names(seurat_list)) seurat_list[["GSM4909287"]] <- filter_and_subset(seurat_list[["GSM4909287"]], 500, 8000, 40)
if ("GSM4909288" %in% names(seurat_list)) seurat_list[["GSM4909288"]] <- filter_and_subset(seurat_list[["GSM4909288"]], 500, 5500, 40)
if ("GSM4909289" %in% names(seurat_list)) seurat_list[["GSM4909289"]] <- filter_and_subset(seurat_list[["GSM4909289"]], 500, 7000, 25)
if ("GSM4909290" %in% names(seurat_list)) seurat_list[["GSM4909290"]] <- filter_and_subset(seurat_list[["GSM4909290"]], 500, 6000, 25)
if ("GSM4909291" %in% names(seurat_list)) seurat_list[["GSM4909291"]] <- filter_and_subset(seurat_list[["GSM4909291"]], 500, 3300, 20)
if ("GSM4909300" %in% names(seurat_list)) seurat_list[["GSM4909300"]] <- filter_and_subset(seurat_list[["GSM4909300"]], 500, 6000, 20)
if ("GSM4909309" %in% names(seurat_list)) seurat_list[["GSM4909309"]] <- filter_and_subset(seurat_list[["GSM4909309"]], 500, 6000, 40)

seurat_list_filtered <- Filter(Negate(is.null), seurat_list)
if (length(seurat_list_filtered) < 2) stop("Fewer than two Seurat objects remain after filtering.")

# --- 6. Post-Filtering QC ---
lapply(names(seurat_list_filtered), function(sample_name) {
  seurat_obj <- seurat_list_filtered[[sample_name]]
  vln_plot_post_filter <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
    plot_annotation(title = paste("Post-filter QC Metrics for", sample_name))
  save_pub_plot(paste0(sample_name, "_post_filter_qc"), vln_plot_post_filter, width = 12, height = 6)
})

# --- 7. Merging and Integration ---
merged_seurat <- merge(x = seurat_list_filtered[[1]], y = seurat_list_filtered[2:length(seurat_list_filtered)])
cat("--- Seurat objects merged ---\n")
rm(data, seurat_list, seurat_list_filtered)

merged_seurat <- NormalizeData(merged_seurat, verbose = FALSE)
merged_seurat <- FindVariableFeatures(merged_seurat, verbose = FALSE)
merged_seurat <- ScaleData(merged_seurat, verbose = FALSE)
merged_seurat <- RunPCA(merged_seurat, npcs = 30, verbose = FALSE)

merged_seurat_obj <- IntegrateLayers(
  object = merged_seurat,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = 'integrated.rpca',
  verbose = TRUE
)
cat("--- Integration complete ---\n")
rm(merged_seurat)

# --- 8. Dimensionality Reduction and Clustering ---
elbow_plot <- ElbowPlot(merged_seurat_obj, ndims = 30)
print(elbow_plot)
save_pub_plot("ElbowPlot_PCA_dims", elbow_plot, width = 8, height = 6)

optimal_dims_overall <- 1:20
merged_seurat_obj <- RunUMAP(merged_seurat_obj, reduction = "integrated.rpca", dims = optimal_dims_overall, verbose = FALSE)
merged_seurat_obj <- FindNeighbors(merged_seurat_obj, reduction = "integrated.rpca", dims = optimal_dims_overall, graph.name = "integrated_snn", verbose = FALSE)
merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = 0.8, graph.name = "integrated_snn", verbose = FALSE)

# --- 9. UMAP Plots ---
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
color_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

plot_clusters <- DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE, label.size = 4, cols = color_vector) +
  ggtitle(paste0("UMAP by Cluster", analysis_suffix)) + NoLegend()
print(plot_clusters)
save_pub_plot(paste0("UMAP_by_Cluster", analysis_suffix), plot_clusters, width = 8, height = 7)

plot_orig_ident <- DimPlot(merged_seurat_obj, reduction = "umap", group.by = "orig.ident", cols = color_vector) +
  ggtitle(paste0("UMAP by Sample", analysis_suffix))
print(plot_orig_ident)
save_pub_plot(paste0("UMAP_by_Sample", analysis_suffix), plot_orig_ident, width = 10, height = 7)

plot_tumor_stage <- DimPlot(merged_seurat_obj, reduction = "umap", group.by = "Tumor_Stage", cols = color_vector) +
  ggtitle(paste0("UMAP by Tumor Stage", analysis_suffix))
print(plot_tumor_stage)
save_pub_plot(paste0("UMAP_by_TumorStage", analysis_suffix), plot_tumor_stage, width = 9, height = 7)

# --- 10. Cluster Markers ---
DefaultAssay(merged_seurat_obj) <- "RNA"
merged_seurat_obj <- JoinLayers(merged_seurat_obj)

cluster_markers <- FindAllMarkers(merged_seurat_obj,
                                  assay = "RNA",
                                  only.pos = TRUE,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25,
                                  verbose = FALSE)

write.csv(cluster_markers,
          file.path(pub_dir, paste0("all_cluster_markers", analysis_suffix, ".csv")),
          row.names = FALSE)

if (nrow(cluster_markers) > 0) {
  top10_markers <- cluster_markers %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC)
  write.csv(top10_markers,
            file.path(pub_dir, paste0("top10_cluster_markers", analysis_suffix, ".csv")),
            row.names = FALSE)
  
  if (nrow(top10_markers) > 0) {
    heatmap_top10 <- DoHeatmap(merged_seurat_obj, features = top10_markers$gene) +
      NoLegend() +
      scale_fill_gradientn(colors = c("blue", "white", "red"))
    print(heatmap_top10)
    save_pub_plot(paste0("Heatmap_top10_cluster_markers", analysis_suffix), heatmap_top10, width = 18, height = 12)
  }
}

# --- 11. Cell Type Annotation ---
file_path <- "C:/Users/BITS/Documents/scrna_final/integrated_seurat_obj_annotated_NoPreneo_AllCells-FINAL.rds"
merged_seurat_obj <- readRDS(file_path)

cluster_to_celltype_map <- c(
  "0"  = "luminal epithelial cell",
  "1"  = "macrophage",
  "2"  = "epithelial cell",
  "3"  = "epithelial cell",
  "4"  = "luminal epithelial cell",
  "5"  = "luminal epithelial cell",
  "6"  = "T cell",
  "7"  = "luminal epithelial cell",
  "8"  = "T cell",
  "9"  = "B cell",
  "10" = "epithelial cell",
  "11" = "endothelial cell",
  "12" = "luminal epithelial cell",
  "13" = "cancer-associated fibroblast",
  "14" = "cancer-associated fibroblast",
  "15" = "fibroblast of breast",
  "16" = "epithelial cell",
  "17" = "B cell",
  "18" = "luminal epithelial cell",
  "19" = "pericyte",
  "20" = "dendritic cell",
  "21" = "epithelial cell",
  "22" = "B cell",
  "23" = "luminal epithelial cell"
)

cluster_col <- "integrated_snn_res.0.8"
if (cluster_col %in% colnames(merged_seurat_obj@meta.data)) {
  original_clusters <- merged_seurat_obj@meta.data[[cluster_col]]
  merged_seurat_obj$cell_type_fixed <- unname(cluster_to_celltype_map[as.character(original_clusters)])
  Idents(merged_seurat_obj) <- "cell_type_fixed"
  cat("--- Cell type identities updated ---\n")
  print(table(Idents(merged_seurat_obj)))
} else {
  stop(paste("Could not find the original cluster column:", cluster_col))
}

# =================================================================
# FILTER GENE LIST AGAINST THE OBJECT — SET ONCE, NEVER CHANGED AGAIN
# =================================================================
available_genes_of_interest <- genes_you_want[genes_you_want %in% rownames(merged_seurat_obj)]
genes_not_found <- genes_you_want[!genes_you_want %in% available_genes_of_interest]

cat(paste0("\n--- Genes requested: ", length(genes_you_want), " ---\n"))
cat(paste0("--- Genes found in object and used in ALL analyses: ",
           length(available_genes_of_interest), " ---\n"))
cat(paste(available_genes_of_interest, collapse = ", "), "\n")

if (length(genes_not_found) > 0) {
  cat(paste0("--- Genes NOT found in object (check spelling/case): ",
             paste(genes_not_found, collapse = ", "), " ---\n"))
}

gene_status <- data.frame(
  gene            = genes_you_want,
  found_in_object = genes_you_want %in% rownames(merged_seurat_obj)
)
write.csv(gene_status, file.path(pub_dir, "gene_list_status.csv"), row.names = FALSE)

# =================================================================
# ANNOTATED UMAPs
# =================================================================
annotated_umap_labelled <- DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("UMAP with Corrected Cell Type Annotations")
print(annotated_umap_labelled)
save_pub_plot("UMAP_annotated_celltype_labelled", annotated_umap_labelled, width = 12, height = 9)

annotated_umap_legend <- DimPlot(merged_seurat_obj, reduction = "umap", label = FALSE) +
  ggtitle("UMAP with Corrected Cell Type Annotations")
print(annotated_umap_legend)
save_pub_plot("UMAP_annotated_celltype_legend_only", annotated_umap_legend, width = 12, height = 9)

# =================================================================
# DOT PLOT — all available genes
# =================================================================
dot_plot_ct <- DotPlot(merged_seurat_obj, features = available_genes_of_interest) +
  labs(title = "Gene Expression Across Cell Types") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial"))
print(dot_plot_ct)
save_pub_plot("DotPlot_GOI_CellType", dot_plot_ct, width = 14, height = 8)

# =================================================================
# DOHEATMAP — all available genes, by cell type
# =================================================================
heatmap_ct <- DoHeatmap(subset(merged_seurat_obj, downsample = 100),
                        features = available_genes_of_interest, label = TRUE) +
  labs(title = "Heatmap of Genes Across Cell Types")
print(heatmap_ct)
save_pub_plot("Heatmap_GOI_CellType", heatmap_ct, width = 16, height = 10)

# =================================================================
# INDIVIDUAL GENE PLOTS — Violin, Ridge, Feature (by Cell Type)
# =================================================================
cat("\n--- Generating per-gene plots (Violin, Ridge, Feature) for",
    length(available_genes_of_interest), "genes ---\n")

for (gene in available_genes_of_interest) {
  cat("--- Plotting:", gene, "---\n")
  
  vln_plot_ct <- VlnPlot(merged_seurat_obj, features = gene, pt.size = 0, assay = "RNA") + NoLegend()
  save_pub_plot(paste0("VlnPlot_", gene, "_CellType"), vln_plot_ct, width = 12, height = 6)
  
  ridge_plot_ct <- RidgePlot(merged_seurat_obj, features = gene, assay = "RNA") + NoLegend()
  save_pub_plot(paste0("RidgePlot_", gene, "_CellType"), ridge_plot_ct, width = 8, height = 7)
  
  feature_plot_ct <- FeaturePlot(merged_seurat_obj, features = gene, order = TRUE, max.cutoff = 'q95')
  save_pub_plot(paste0("FeaturePlot_", gene, "_CellType"), feature_plot_ct, width = 8, height = 7)
}

# =================================================================
# SPLIT VIOLIN PLOTS (Cell Type grouped, split by Tumor Stage)
# =================================================================
cat("\n--- Generating split violin plots for", length(available_genes_of_interest), "genes ---\n")

for (gene in available_genes_of_interest) {
  cat("--- Split violin:", gene, "---\n")
  split_vln_plot <- VlnPlot(merged_seurat_obj, features = gene,
                            split.by = "Tumor_Stage", pt.size = 0, assay = "RNA") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, family = "Arial"))
  save_pub_plot(paste0("SplitVlnPlot_", gene, "_CellType_Stage"), split_vln_plot, width = 16, height = 7)
}

# =================================================================
# GENE SIGNATURE SCORE VIOLIN
# =================================================================
if (!"GenesOfInterest_Score1" %in% colnames(merged_seurat_obj@meta.data)) {
  merged_seurat_obj <- AddModuleScore(merged_seurat_obj,
                                      features = list(available_genes_of_interest),
                                      name = "GenesOfInterest_Score", assay = "RNA")
}

signature_plot <- VlnPlot(merged_seurat_obj, features = "GenesOfInterest_Score1",
                          split.by = "Tumor_Stage", pt.size = 0) +
  labs(y = "Gene Signature Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, family = "Arial"))
print(signature_plot)
save_pub_plot("VlnPlot_GeneSignatureScore_CellType_Stage", signature_plot, width = 16, height = 7)

# =================================================================
# COMBINED IDENTITY HEATMAP (Cell Type + Stage)
# =================================================================
merged_seurat_obj$celltype_stage <- paste(merged_seurat_obj$cell_type_fixed,
                                          merged_seurat_obj$Tumor_Stage, sep = "_")
Idents(merged_seurat_obj) <- "celltype_stage"

heatmap_split <- DoHeatmap(subset(merged_seurat_obj, downsample = 50),
                           features = available_genes_of_interest, label = TRUE, size = 3) +
  labs(title = "Heatmap of GOI Expression by Cell Type and Tumor Stage")
print(heatmap_split)
save_pub_plot("Heatmap_Split_by_CellType_and_Stage", heatmap_split, width = 18, height = 10)

# =================================================================
# SINGLE-CELL HEATMAP BY TUMOR STAGE
# =================================================================
heatmap_by_stage <- DoHeatmap(subset(merged_seurat_obj, downsample = 200),
                              features = available_genes_of_interest,
                              group.by = "Tumor_Stage", label = TRUE, size = 4) +
  labs(title = "Heatmap of GOI Expression by Tumor Stage")
print(heatmap_by_stage)
save_pub_plot("Heatmap_GOI_by_TumorStage_SingleCell", heatmap_by_stage, width = 12, height = 8)

# =================================================================
# INDIVIDUAL SPLIT VIOLIN PLOTS (group by cell type, split by stage)
# =================================================================
cat("\n--- Generating individual split violin plots for", length(available_genes_of_interest), "genes ---\n")

for (gene in available_genes_of_interest) {
  cat("--- Individual split violin:", gene, "---\n")
  individual_split_vln_plot <- VlnPlot(merged_seurat_obj, features = gene,
                                       group.by = "cell_type_fixed",
                                       split.by = "Tumor_Stage", pt.size = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, family = "Arial")) +
    labs(title = paste("Expression of", gene, "by Cell Type and Stage"))
  save_pub_plot(paste0("VlnPlot_Split_", gene, "_Individual"), individual_split_vln_plot, width = 16, height = 7)
}

# =================================================================
# pheatmap — AVERAGE EXPRESSION BY CELL TYPE
# =================================================================
cat("\n--- Generating pheatmap: Average Expression by Cell Type ---\n")
Idents(merged_seurat_obj) <- "cell_type_fixed"

avg_exp_by_celltype <- AverageExpression(merged_seurat_obj,
                                         features = available_genes_of_interest,
                                         group.by = "ident", assay = "RNA")

heatmap_celltype <- pheatmap(avg_exp_by_celltype$RNA,
                             main        = "Average Gene Expression by Cell Type",
                             scale       = "row",
                             cluster_cols = TRUE,
                             cluster_rows = TRUE,
                             border_color = "grey60",
                             fontfamily  = "Arial")   # Arial for pheatmap
print(heatmap_celltype)
save_pub_pheatmap("Heatmap_AvgExp_by_CellType", heatmap_celltype, width = 12, height = 9)

# =================================================================
# pheatmap — AVERAGE EXPRESSION BY TUMOR STAGE
# =================================================================
cat("\n--- Generating pheatmap: Average Expression by Tumor Stage ---\n")
Idents(merged_seurat_obj) <- "Tumor_Stage"

avg_exp_by_stage <- AverageExpression(merged_seurat_obj,
                                      features = available_genes_of_interest,
                                      group.by = "ident", assay = "RNA")

avg_matrix <- avg_exp_by_stage$RNA

heatmap_stage <- pheatmap(avg_matrix,
                          main         = "Average Gene Expression by Tumor Stage",
                          scale        = "row",
                          cluster_cols = FALSE,
                          cluster_rows = TRUE,
                          border_color = "grey60",
                          labels_col   = c("None", "Stage 3", "Stage 2", "Stage 1"),
                          fontfamily   = "Arial")   # Arial for pheatmap
print(heatmap_stage)
save_pub_pheatmap("Heatmap_AvgExp_by_TumorStage", heatmap_stage, width = 10, height = 9)

cat("\n=================================================================\n")
cat(paste0("ALL plots saved to: ", pub_dir, "\n"))
cat(paste0("Format: TIFF, 600 DPI, LZW compression, Arial font\n"))
cat(paste0("Genes analysed: ", length(available_genes_of_interest), "/", length(genes_you_want), "\n"))
if (length(genes_not_found) > 0) {
  cat(paste0("Genes not found in object: ", paste(genes_not_found, collapse = ", "), "\n"))
}
cat("See gene_list_status.csv in the output folder for a full gene report.\n")
cat("=================================================================\n")