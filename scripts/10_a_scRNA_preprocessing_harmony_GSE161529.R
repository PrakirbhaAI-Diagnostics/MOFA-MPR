# --- 0. SCRIPT SETUP ---

# Install necessary packages if they are not already present
# install.packages(c('Seurat', 'tidyverse', 'patchwork', 'devtools', 'future'))
# devtools::install_github('immunogenomics/presto', quiet = TRUE) # For faster FindAllMarkers, optional

# Install Harmony if not already present
if (!requireNamespace("harmony", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("immunogenomics/harmony", quiet = TRUE)
}

# Load libraries
library(Seurat)
library(tidyverse)
library(patchwork)
library(future)
library(harmony) # Load harmony

# Configure parallelization and memory
options(future.globals.maxSize = 30000 * 1024^2) # 30GB, adjust as needed
plan("multisession", workers = 1) # Adjust workers based on your system's available cores for speed

# --- Analysis Suffix for Naming Outputs ---
analysis_suffix <- "_NoPreneo_AllCells_Harmony"

# --- 1. LOAD DATA, METADATA, AND DEFINE EXCLUSIONS ---
cat("--- Loading metadata and defining exclusions ---\n")
metadata <- read.csv("~/scrna/metadata.csv") # IMPORTANT: Ensure this path is correct
preneoplastic_samples_to_exclude <- c("GSM4909277", "GSM4909278", "GSM4909279", "GSM4909280")

# Function to add stage-wise metadata to Seurat objects
add_stage_metadata <- function(seurat_obj, metadata_df) { # Renamed metadata to metadata_df to avoid conflict
  sample_name <- seurat_obj@meta.data$orig.ident[1]
  sample_info <- dplyr::filter(metadata_df, GSM.ID == sample_name)
  if (nrow(sample_info) > 0) {
    seurat_obj@meta.data$Tumor_Stage <- sample_info$Tumor.Stage
  } else {
    warning(paste("No metadata found for sample:", sample_name))
    seurat_obj@meta.data$Tumor_Stage <- NA # Assign NA if no metadata found
  }
  return(seurat_obj)
}

# Load 10X data, create Seurat objects, excluding preneoplastic samples
main_dir <- "~/scrna/GSE161529" # IMPORTANT: Ensure this path is correct
all_sub_dirs <- list.dirs(main_dir, recursive = FALSE)
sub_dirs <- all_sub_dirs[!basename(all_sub_dirs) %in% preneoplastic_samples_to_exclude]

if (length(sub_dirs) == 0) {
  stop("No sample directories found after excluding preneoplastic samples. Check 'main_dir' and 'preneoplastic_samples_to_exclude'.")
}

seurat_list <- list()
for (dir_path in sub_dirs) {
  subfolder_name <- basename(dir_path)
  cat("Processing sample:", subfolder_name, "\n")
  if (!any(file.exists(file.path(dir_path, "matrix.mtx.gz")), file.exists(file.path(dir_path, "matrix.mtx")))) {
    warning(paste("Matrix file (matrix.mtx.gz or matrix.mtx) not found in", dir_path, ". Skipping this directory."))
    next
  }
  data_counts <- Read10X(data.dir = dir_path)
  seurat_list[[subfolder_name]] <- CreateSeuratObject(counts = data_counts, project = subfolder_name)
  cat("Created Seurat object for:", subfolder_name, "\n")
}

if (length(seurat_list) == 0) {
  stop("No Seurat objects were created. Check data directories and 10X output files.")
}

# Add stage-wise metadata and mitochondrial gene percentage
seurat_list <- lapply(seurat_list, function(obj) add_stage_metadata(obj, metadata))
seurat_list <- lapply(seurat_list, function(seurat_obj) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-") # Assumes human MT genes
  return(seurat_obj)
})

# --- 2. PRE-FILTERING QUALITY CONTROL ---
cat("--- Performing pre-filtering QC ---\n")
# IMPORTANT: Ensure base path "~/scrna/" is writable
qc_base_dir <- file.path("~", "scrna", "QC_Reports") # Centralized QC reports
pre_filter_dir <- file.path(qc_base_dir, paste0("PreFilter", analysis_suffix))
dir.create(pre_filter_dir, recursive = TRUE, showWarnings = FALSE)

lapply(names(seurat_list), function(sample_name) {
  seurat_obj_qc <- seurat_list[[sample_name]] # Use a temporary variable for clarity
  vln_plot_pre <- VlnPlot(seurat_obj_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1) +
    plot_annotation(title = paste("Pre-filter QC for", sample_name))
  ggsave(filename = file.path(pre_filter_dir, paste0(sample_name, "_pre_filter_qc.png")),
         plot = vln_plot_pre, width = 12, height = 6)
  
  scatter_plot_pre <- FeatureScatter(seurat_obj_qc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    ggtitle(paste("nCount_RNA vs nFeature_RNA for", sample_name))
  ggsave(filename = file.path(pre_filter_dir, paste0(sample_name, "_feature_scatter.png")),
         plot = scatter_plot_pre, width = 8, height = 6)
})
cat(paste0("Pre-filtering QC plots saved to: ", pre_filter_dir, "\n"))

# --- 3. APPLY FILTERING ---
cat("--- Applying filtering criteria to individual samples ---\n")
filter_and_subset <- function(obj, id, nF_low, nF_high, pMT_high) {
  if (is.null(obj)) return(NULL)
  original_cells <- ncol(obj)
  obj_filtered <- subset(obj, subset = nFeature_RNA > nF_low & nFeature_RNA < nF_high & percent.mt < pMT_high)
  cat(sprintf("Sample %s: Before filtering = %d cells, After filtering = %d cells\n",
              id, original_cells, ncol(obj_filtered)))
  if (ncol(obj_filtered) == 0) {
    warning(paste("Sample", id, "has 0 cells after filtering. Check filter thresholds."))
    return(NULL)
  }
  return(obj_filtered)
}

filter_params <- list(
  GSM4909253=c(500,6500,20), GSM4909254=c(500,7500,20), GSM4909257=c(500,7000,30),
  GSM4909261=c(500,8000,20), GSM4909281=c(500,7000,20), GSM4909282=c(500,5000,20),
  GSM4909283=c(500,7000,20), GSM4909284=c(500,5500,30), GSM4909285=c(500,5500,20),
  GSM4909286=c(500,5000,30), GSM4909287=c(500,8000,40), GSM4909288=c(500,5500,40),
  GSM4909289=c(500,7000,25), GSM4909290=c(500,6000,25), GSM4909291=c(500,3300,20),
  GSM4909300=c(500,6000,20), GSM4909309=c(500,6000,40)
)

seurat_list_filtered <- list()
for(sample_id in names(seurat_list)){
  if(sample_id %in% names(filter_params)){
    params <- filter_params[[sample_id]]
    seurat_list_filtered[[sample_id]] <- filter_and_subset(seurat_list[[sample_id]], sample_id, params[1], params[2], params[3])
  } else {
    warning(paste("No specific filtering parameters found for sample:", sample_id, ". It will be carried forward unfiltered by these criteria."))
    seurat_list_filtered[[sample_id]] <- seurat_list[[sample_id]] # Keep it if no specific params
  }
}
seurat_list_filtered <- Filter(Negate(is.null), seurat_list_filtered) # Remove any samples that became NULL (0 cells)

if (length(seurat_list_filtered) == 0) {
  stop("No Seurat objects remaining after applying filtering rules.")
}
if (length(seurat_list_filtered) < 2 && length(seurat_list_filtered) > 0) {
  warning("Only one Seurat object remaining. Integration will be skipped. Proceeding with single sample analysis.")
} else if (length(seurat_list_filtered) < 2) {
  stop("Not enough Seurat objects remaining for integration (less than 2).")
}

# --- 4. POST-FILTERING QUALITY CONTROL ---
cat("--- Performing post-filtering QC ---\n")
post_filter_dir <- file.path(qc_base_dir, paste0("PostFilter", analysis_suffix))
dir.create(post_filter_dir, recursive = TRUE, showWarnings = FALSE)

lapply(names(seurat_list_filtered), function(sample_name) {
  seurat_obj_qc <- seurat_list_filtered[[sample_name]]
  vln_plot_post <- VlnPlot(seurat_obj_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1) +
    plot_annotation(title = paste("Post-filter QC for", sample_name))
  ggsave(filename = file.path(post_filter_dir, paste0(sample_name, "_post_filter_qc.png")), plot = vln_plot_post, width = 12, height = 6)
})
cat(paste0("Post-filtering QC plots saved to: ", post_filter_dir, "\n"))

# --- 5. MERGE SEURAT OBJECTS ---
cat("--- Merging Seurat objects ---\n")
seurat_list_renamed <- lapply(names(seurat_list_filtered), function(sample_name) {
  seurat_obj <- seurat_list_filtered[[sample_name]]
  seurat_obj <- RenameCells(seurat_obj, add.cell.id = sample_name)
  return(seurat_obj)
})

if (length(seurat_list_renamed) == 1) {
  merged_seurat_obj <- seurat_list_renamed[[1]]
  cat("--- Only one sample object. No merge performed. 'merged_seurat_obj' is the single sample object. ---\n")
} else {
  merged_seurat_obj <- merge(seurat_list_renamed[[1]], y = seurat_list_renamed[2:length(seurat_list_renamed)], project = "MergedProject")
  cat("--- Seurat objects merged into 'merged_seurat_obj' ---\n")
}
rm(seurat_list, seurat_list_filtered, seurat_list_renamed); gc()

# --- 6. PRE-PROCESSING AND HARMONY INTEGRATION ---
cat("--- Performing standard pre-processing on merged object ---\n")
DefaultAssay(merged_seurat_obj) <- "RNA"
merged_seurat_obj <- NormalizeData(merged_seurat_obj, verbose = FALSE)
merged_seurat_obj <- FindVariableFeatures(merged_seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
merged_seurat_obj <- ScaleData(merged_seurat_obj, features = VariableFeatures(object = merged_seurat_obj), verbose = FALSE)
merged_seurat_obj <- RunPCA(merged_seurat_obj, features = VariableFeatures(object = merged_seurat_obj), npcs = 30, verbose = FALSE)

if (length(unique(merged_seurat_obj$orig.ident)) > 1) {
  cat("--- Performing Harmony integration ---\n")
  merged_seurat_obj <- RunHarmony(
    object = merged_seurat_obj,
    group.by.vars = "orig.ident", 
    reduction = "pca",
    assay.use = DefaultAssay(merged_seurat_obj), 
    project.dim = FALSE,
    reduction.save = "harmony", 
    verbose = FALSE
  )
  cat("--- Harmony integration complete. Using 'harmony' reduction for downstream analysis. ---\n")
  integration_reduction_name <- "harmony"
} else {
  cat("--- Only one unique 'orig.ident'. Skipping Harmony integration. Using 'pca' for downstream analysis. ---\n")
  integration_reduction_name <- "pca" 
}
DefaultAssay(merged_seurat_obj) <- "RNA"

# --- 7. DIMENSIONALITY REDUCTION (UMAP) & CLUSTERING ---
cat(paste0("--- Performing UMAP and clustering using '", integration_reduction_name, "' reduction ---\n"))
plots_base_dir <- file.path("~", "scrna", "Analysis_Outputs") 
plots_main_dir <- file.path(plots_base_dir, paste0("Plots", analysis_suffix))
dir.create(plots_main_dir, recursive = TRUE, showWarnings = FALSE)

dims_to_use <- 1:30 

merged_seurat_obj <- RunUMAP(
  merged_seurat_obj,
  reduction = integration_reduction_name,
  dims = dims_to_use,
  verbose = FALSE
)
merged_seurat_obj <- FindNeighbors(merged_seurat_obj, reduction = integration_reduction_name, dims = dims_to_use, verbose = FALSE)
merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = 0.5, verbose = FALSE)

# --- 8. VISUALIZE AND SAVE INTEGRATED DATA ---
cat("--- Saving UMAP plots for integrated cells ---\n")
umap_plots_dir <- file.path(plots_main_dir, paste0("UMAPs_", integration_reduction_name, "_Integrated"))
dir.create(umap_plots_dir, recursive = TRUE, showWarnings = FALSE)

plot_orig_ident <- DimPlot(merged_seurat_obj, reduction = "umap", group.by = "orig.ident") +
  ggtitle(paste0("UMAP ", casefold(integration_reduction_name, upper=TRUE), " Integrated (by Sample)", analysis_suffix))
ggsave(file.path(umap_plots_dir, paste0("umap_", integration_reduction_name, "_orig_ident", analysis_suffix, ".png")), plot = plot_orig_ident, width = 10, height = 7)

plot_tumor_stage <- DimPlot(merged_seurat_obj, reduction = "umap", group.by = "Tumor_Stage") +
  ggtitle(paste0("UMAP ", casefold(integration_reduction_name, upper=TRUE), " Integrated (by Tumor Stage)", analysis_suffix))
ggsave(file.path(umap_plots_dir, paste0("umap_", integration_reduction_name, "_tumor_stage", analysis_suffix, ".png")), plot = plot_tumor_stage, width = 10, height = 7)

plot_clusters <- DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE, group.by = "seurat_clusters") +
  ggtitle(paste0("UMAP ", casefold(integration_reduction_name, upper=TRUE), " Integrated (by Cluster)", analysis_suffix))
ggsave(file.path(umap_plots_dir, paste0("umap_", integration_reduction_name, "_clusters", analysis_suffix, ".png")), plot = plot_clusters, width = 10, height = 7)
cat(paste0("UMAP plots saved to: ", umap_plots_dir, "\n"))

rds_save_dir <- file.path("~", "scrna", "Seurat_Objects")
dir.create(rds_save_dir, recursive = TRUE, showWarnings = FALSE)
output_rds_path <- file.path(rds_save_dir, paste0("integrated_seurat_obj", analysis_suffix, ".rds"))
saveRDS(merged_seurat_obj, output_rds_path)
cat(paste0("Integrated Seurat object saved to: ", output_rds_path, "\n"))

# --- 9. FIND CLUSTER MARKERS ---
cat("--- Finding cluster markers (post-integration) ---\n")
DefaultAssay(merged_seurat_obj) <- "RNA" # Ensure RNA assay is default for marker finding
# For Seurat v5, if the object was created by merging multiple samples,
# JoinLayers might be needed before FindAllMarkers on the RNA assay.
merged_seurat_obj <- JoinLayers(merged_seurat_obj, assay = "RNA")
cat("--- Layers joined for RNA assay before FindAllMarkers ---\n")

cluster_markers <- FindAllMarkers(merged_seurat_obj,
                                  assay = "RNA",
                                  only.pos = TRUE,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25, 
                                  test.use = "wilcox", 
                                  verbose = FALSE)

markers_save_dir <- file.path(plots_base_dir, paste0("Markers", analysis_suffix))
dir.create(markers_save_dir, recursive = TRUE, showWarnings = FALSE)
output_markers_rds_path <- file.path(markers_save_dir, paste0("cluster_markers", analysis_suffix, ".rds"))
output_markers_csv_path <- file.path(markers_save_dir, paste0("cluster_markers", analysis_suffix, ".csv"))
saveRDS(cluster_markers, output_markers_rds_path)
write.csv(cluster_markers, output_markers_csv_path, row.names = FALSE)
cat(paste0("Cluster markers saved to: ", markers_save_dir, "\n"))

cat("--- Plotting top 10 cluster markers ---\n")
if (nrow(cluster_markers) > 0) {
  top10_markers <- cluster_markers %>%
    filter(p_val_adj < 0.05) %>% 
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC) %>% 
    ungroup()
  
  if (nrow(top10_markers) > 0 && length(unique(top10_markers$gene)) > 0) {
    dotplot_top10 <- DotPlot(merged_seurat_obj, features = unique(top10_markers$gene), assay = "RNA", group.by = "seurat_clusters") +
      scale_colour_gradient2(low="blue", mid="white", high="red") + 
      RotatedAxis() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title = element_blank()) + 
      ggtitle(paste0("Top 10 Markers per Cluster", analysis_suffix))
    ggsave(file.path(markers_save_dir, paste0("dotplot_top10_cluster_markers", analysis_suffix, ".png")), plot = dotplot_top10, width = 20, height = max(8, 0.5*length(unique(top10_markers$cluster)) + 3 )) 
    cat("Top 10 cluster markers dot plot saved.\n")
  } else {
    cat("No significant top 10 markers found after p-value filtering or not enough unique genes to plot.\n")
  }
} else {
  cat("No cluster markers found. Skipping top 10 marker plot.\n")
}


# --- 10. DIFFERENTIAL GENE EXPRESSION (DGE) ANALYSIS ---
cat("--- Starting DGE Analysis ---\n")
dge_base_dir <- file.path(plots_base_dir, paste0("DGE_Analysis", analysis_suffix))
dir.create(dge_base_dir, recursive = TRUE, showWarnings = FALSE)
cat(paste0("Base directory for DGE output: ", dge_base_dir, "\n"))

## DGE Comparison 1: Early Stage vs. Later/Non-Tumor
cat(paste0("\nDGE Comp1: Early Stage (1,2) vs. Later/Non-Tumor", analysis_suffix, "\n"))
current_object_for_dge <- merged_seurat_obj 
current_object_for_dge$Comparison1_Group <- NA_character_
current_object_for_dge$Comparison1_Group[which(current_object_for_dge$Tumor_Stage %in% c("1", "2"))] <- "Early_Tumor_Stage"
later_cats <- c("3", "None")
if("BRCA1" %in% levels(as.factor(current_object_for_dge$Tumor_Stage))) later_cats <- c(later_cats, "BRCA1")
current_object_for_dge$Comparison1_Group[which(current_object_for_dge$Tumor_Stage %in% later_cats)] <- "Later_Or_NonTumor"

cells_comp1 <- colnames(current_object_for_dge)[!is.na(current_object_for_dge$Comparison1_Group) & current_object_for_dge$Comparison1_Group %in% c("Early_Tumor_Stage", "Later_Or_NonTumor")]
if (length(cells_comp1) > 0) {
  dge_subset_comp1 <- subset(current_object_for_dge, cells = cells_comp1)
  dge_subset_comp1$Comparison1_Group <- factor(dge_subset_comp1$Comparison1_Group, levels = c("Early_Tumor_Stage", "Later_Or_NonTumor"))
  Idents(dge_subset_comp1) <- "Comparison1_Group"
  counts_comp1 <- table(Idents(dge_subset_comp1)); cat("Cell counts for DGE Comp1:\n"); print(counts_comp1)
  if (length(counts_comp1) == 2 && all(counts_comp1 >= 3)) {
    dge_subset_comp1 <- JoinLayers(dge_subset_comp1, assay="RNA") # JOIN LAYERS BEFORE FINDMARKERS
    cat("Layers joined for DGE Comp1 subset.\n")
    markers_c1 <- FindMarkers(dge_subset_comp1, assay="RNA", ident.1="Early_Tumor_Stage", ident.2="Later_Or_NonTumor", min.pct=0.1, logfc.threshold=0.25, only.pos=FALSE, verbose=FALSE)
    output_dir_c1 <- file.path(dge_base_dir, "EarlyStage_vs_LaterOrNonTumor")
    dir.create(output_dir_c1, recursive=TRUE, showWarnings=FALSE)
    up_c1 <- markers_c1 %>% filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC))
    if(nrow(up_c1)>0) { write.csv(up_c1, file.path(output_dir_c1, paste0("genes_up_EarlyVsLater", analysis_suffix, ".csv")), row.names=TRUE); cat("Comp1 Up genes saved.\n"); print(head(up_c1,10))} else {cat("No Comp1 Up genes.\n")}
    down_c1 <- markers_c1 %>% filter(avg_log2FC < 0) %>% arrange(avg_log2FC)
    if(nrow(down_c1)>0) { write.csv(down_c1, file.path(output_dir_c1, paste0("genes_down_EarlyVsLater", analysis_suffix, ".csv")), row.names=TRUE); cat("Comp1 Down genes saved.\n"); print(head(down_c1,10))} else {cat("No Comp1 Down genes.\n")}
  } else { warning("DGE Comp1: Groups insufficient. Skipping.") }
  if(exists("dge_subset_comp1")) rm(dge_subset_comp1)
} else { warning("DGE Comp1: No cells for groups. Skipping.") }

## DGE Comparison 2: Pairwise Between All Tumor Stages
cat(paste0("\nPairwise DGE between Tumor Stages", analysis_suffix, "\n"))
if (!is.factor(current_object_for_dge$Tumor_Stage)) current_object_for_dge$Tumor_Stage <- as.factor(current_object_for_dge$Tumor_Stage)
unique_stages_dge <- levels(current_object_for_dge$Tumor_Stage)
cat("Unique stages for pairwise DGE:", paste(unique_stages_dge, collapse=", "), "\n")
pairwise_base_dir <- file.path(dge_base_dir, "Pairwise_Stage_Comparisons")
dir.create(pairwise_base_dir, recursive=TRUE, showWarnings=FALSE)
if (length(unique_stages_dge) < 2) {
  cat("Less than 2 unique stages for pairwise DGE. Skipping.\n")
} else {
  stage_pairs <- combn(unique_stages_dge, 2, simplify = FALSE)
  for (pair in stage_pairs) {
    s1 <- pair[1]; s2 <- pair[2]
    cat(paste0("\nPairwise DGE: ", s1, " vs ", s2, "\n"))
    comp_name <- paste0(gsub("[^A-Za-z0-9_]","_",s1), "_vs_", gsub("[^A-Za-z0-9_]","_",s2))
    output_dir_pair <- file.path(pairwise_base_dir, comp_name); dir.create(output_dir_pair, recursive=TRUE, showWarnings=FALSE)
    cells_pair <- WhichCells(current_object_for_dge, expression = Tumor_Stage %in% c(s1, s2))
    if(length(cells_pair)==0){cat("No cells. Skipping.\n"); next}
    dge_subset_pair <- subset(current_object_for_dge, cells=cells_pair)
    dge_subset_pair$Tumor_Stage <- factor(dge_subset_pair$Tumor_Stage, levels=c(s1,s2))
    Idents(dge_subset_pair) <- "Tumor_Stage"
    counts_pair <- table(Idents(dge_subset_pair)); cat("Cell counts:\n"); print(counts_pair)
    if(any(counts_pair < 3)){cat("Too few cells. Skipping.\n"); if(exists("dge_subset_pair")) rm(dge_subset_pair); next}
    dge_subset_pair <- JoinLayers(dge_subset_pair, assay="RNA") # JOIN LAYERS BEFORE FINDMARKERS
    cat(paste0("Layers joined for pairwise DGE subset (", s1, " vs ", s2, ").\n"))
    markers_pair <- FindMarkers(dge_subset_pair, assay="RNA", ident.1=s1, ident.2=s2, min.pct=0.1, logfc.threshold=0.25, only.pos=FALSE, verbose=FALSE)
    if(nrow(markers_pair)==0){cat("No markers found. Skipping.\n"); if(exists("dge_subset_pair")) rm(dge_subset_pair); next}
    up_pair <- markers_pair %>% filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC))
    if(nrow(up_pair)>0){write.csv(up_pair, file.path(output_dir_pair, paste0("genes_up_",comp_name, analysis_suffix, ".csv")), row.names=TRUE); cat("Pairwise Up saved.\n")} else {cat("No Pairwise Up genes.\n")}
    down_pair <- markers_pair %>% filter(avg_log2FC < 0) %>% arrange(avg_log2FC)
    if(nrow(down_pair)>0){write.csv(down_pair, file.path(output_dir_pair, paste0("genes_down_",comp_name, analysis_suffix, ".csv")), row.names=TRUE); cat("Pairwise Down saved.\n")} else {cat("No Pairwise Down genes.\n")}
    if(exists("dge_subset_pair")) rm(dge_subset_pair)
  }
}
cat("--- DGE Analysis complete. ---\n")

# --- 11. FINALIZATION ---
cat(paste0("\n--- Analysis with suffix '", analysis_suffix, "' finished successfully. ---\n"))
plan("sequential") # Reset to default plan
gc() # Garbage collect to free up memory