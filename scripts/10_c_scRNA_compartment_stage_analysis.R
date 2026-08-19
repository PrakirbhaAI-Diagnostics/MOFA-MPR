library(Seurat)
library(tidyverse)
library(patchwork)
library(future)
library(RColorBrewer)
library(pheatmap)
library(extrafont)


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

analysis_suffix <- "_NoPreneo_AllCells-FINAL"

# =================================================================
# GENE SETS 
# =================================================================

tumour_enriched <- c(
  "CENPF", "COL10A1", "COL11A1", "ENAH", "FAP", "FN1", "MMP13",
  "NAT1", "POSTN", "PYCARD", "RAB25", "SMC4", "SPP1", "TMEM132A"
)

normal_enriched <- c(
  "ABCA8", "ADH1B", "ADIPOQ", "C7", "CFD", "CHRDL1", "CLDN5", "DCN",
  "DPT", "EGR1", "FABP4", "FGF2", "FHL1", "FOSB", "GAS1", "GRIA1",
  "HLF", "IRS2", "LAMA2", "MEOX1", "MFAP4", "MME", "SRPX"
)

early_enriched <- c(
  "ABCA8", "ACKR1", "ACTG2", "ADH1B", "C7", "CHRDL1", "FABP4", "FCGBP",
  "MEOX1", "MYH11", "SELP", "TMEM132C", "TNN", "UBE4A", "VANGL2"
)

late_enriched <- c("CENPF", "CRABP2", "NEK2", "SLC11A1")

dual_track <- c("FABP4", "ADH1B", "ABCA8", "CHRDL1", "MEOX1", "CENPF", "C7")

extra_genes <- c("LAD1", "COL11A1", "CD44", "H2BC3", "MAPK11", "EPYC")

GENE_SETS <- list(
  Tumour_Enriched = tumour_enriched,
  Normal_Enriched = normal_enriched,
  Early_Enriched  = early_enriched,
  Late_Enriched   = late_enriched,
  Dual_Track      = dual_track,
  Extra_Requested = extra_genes
)

genes_you_want <- unique(c(tumour_enriched, normal_enriched, early_enriched,
                           late_enriched, dual_track, extra_genes))

cat(paste0("--- Total unique genes in analysis list: ", length(genes_you_want), " ---\n"))
for (s in names(GENE_SETS))
  cat(sprintf("      %-16s %2d genes\n", s, length(GENE_SETS[[s]])))
cat(paste(genes_you_want, collapse = ", "), "\n")

# =================================================================
# SINGLE OUTPUT DIRECTORY 
# =================================================================
pub_dir <- file.path("results", "scRNA_validated_genes")
dir.create(pub_dir, recursive = TRUE, showWarnings = FALSE)
cat(paste0("--- ALL plots will be saved to: ", pub_dir, " ---\n"))

save_pub_plot <- function(filename_base, plot_obj, width = 10, height = 8) {
  ggsave(
    file.path(pub_dir, paste0(filename_base, ".tiff")),
    plot        = plot_obj,
    width       = width,
    height      = height,
    dpi         = 600,
    units       = "in",
    device      = "tiff",
    compression = "lzw"
  )
}

save_pub_pheatmap <- function(filename_base, pheatmap_obj, width = 10, height = 8) {
  tiff(
    file.path(pub_dir, paste0(filename_base, ".tiff")),
    width       = width,
    height      = height,
    units       = "in",
    res         = 600,
    compression = "lzw",
    family      = "Arial"
  )
  print(pheatmap_obj)
  dev.off()
}

file_path <- file.path("data", "processed", "scRNA", "integrated_seurat_obj_annotated_NoPreneo_AllCells-FINAL.rds")
merged_seurat_obj <- readRDS(file_path)
cat(sprintf("Object: %d cells x %d genes\n",
            ncol(merged_seurat_obj), nrow(merged_seurat_obj)))

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

cat("--- Tumor_Stage distribution ---\n")
print(table(merged_seurat_obj$Tumor_Stage))

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

# per-set availability, so a set that is badly covered is visible
gene_status <- do.call(rbind, lapply(names(GENE_SETS), function(s) {
  data.frame(gene = GENE_SETS[[s]], gene_set = s,
             found_in_object = GENE_SETS[[s]] %in% rownames(merged_seurat_obj))
}))
write.csv(gene_status, file.path(pub_dir, "gene_list_status_by_set.csv"), row.names = FALSE)
cat("\n--- Availability by set ---\n")
print(with(gene_status, table(gene_set, found_in_object)))

# available genes per set, used for module scores
SETS_AVAIL <- lapply(GENE_SETS, function(g) intersect(g, rownames(merged_seurat_obj)))

# =================================================================
# ANNOTATED UMAPs
# =================================================================
annotated_umap_labelled <- DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE, repel = TRUE) +
  labs(title = "UMAP — Annotated Cell Types")
print(annotated_umap_labelled)
save_pub_plot("UMAP_annotated_celltype_labelled", annotated_umap_labelled, width = 12, height = 9)

annotated_umap_legend <- DimPlot(merged_seurat_obj, reduction = "umap", label = FALSE) +
  labs(title = "UMAP — Annotated Cell Types")
print(annotated_umap_legend)
save_pub_plot("UMAP_annotated_celltype_legend_only", annotated_umap_legend, width = 12, height = 9)

# =================================================================
# DOT PLOT — all available genes
# =================================================================
Idents(merged_seurat_obj) <- "cell_type_fixed"
dot_plot_ct <- DotPlot(merged_seurat_obj, features = available_genes_of_interest) +
  labs(title = "Gene Expression Across Cell Types") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial"))
print(dot_plot_ct)
save_pub_plot("DotPlot_GOI_CellType", dot_plot_ct,
              width = max(14, length(available_genes_of_interest) * 0.32), height = 8)

# DotPlot per gene set 
for (s in names(SETS_AVAIL)) {
  g <- SETS_AVAIL[[s]]
  if (length(g) < 2) next
  dp <- DotPlot(merged_seurat_obj, features = g) +
    labs(title = paste0(gsub("_", " ", s), " — Expression Across Cell Types")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial"))
  save_pub_plot(paste0("DotPlot_", s, "_CellType"), dp,
                width = max(10, length(g) * 0.45), height = 8)
}

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
for (s in names(SETS_AVAIL)) {
  g <- SETS_AVAIL[[s]]
  if (length(g) < 2) next
  score_name <- paste0(s, "_Score")
  if (!paste0(score_name, "1") %in% colnames(merged_seurat_obj@meta.data)) {
    merged_seurat_obj <- AddModuleScore(merged_seurat_obj,
                                        features = list(g),
                                        name = score_name, assay = "RNA")
  }
  sig_plot <- VlnPlot(merged_seurat_obj, features = paste0(score_name, "1"),
                      split.by = "Tumor_Stage", pt.size = 0) +
    labs(y = paste(gsub("_", " ", s), "Score"),
         title = paste0(gsub("_", " ", s), " Signature — Cell Type x Stage")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, family = "Arial"))
  print(sig_plot)
  save_pub_plot(paste0("VlnPlot_", score_name, "_CellType_Stage"), sig_plot,
                width = 16, height = 7)
}

# combined signature score across all genes 
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
# COMPARTMENT TABLES
# =================================================================
cat("\n--- Writing compartment tables ---\n")

# (a) mean expression per cell type
Idents(merged_seurat_obj) <- "cell_type_fixed"
avg_ct <- AverageExpression(merged_seurat_obj, features = available_genes_of_interest,
                            assays = "RNA", group.by = "cell_type_fixed")$RNA
write.csv(as.data.frame(avg_ct),
          file.path(pub_dir, "table_mean_expression_by_celltype.csv"))

# (b) mean expression per cell type x stage 
avg_cs <- AverageExpression(merged_seurat_obj, features = available_genes_of_interest,
                            assays = "RNA", group.by = "celltype_stage")$RNA
write.csv(as.data.frame(avg_cs),
          file.path(pub_dir, "table_mean_expression_by_celltype_stage.csv"))

# (c) compartment PROPORTION by stage 
prop_tab <- as.data.frame.matrix(
  prop.table(table(merged_seurat_obj$cell_type_fixed,
                   merged_seurat_obj$Tumor_Stage), margin = 2) * 100
)
write.csv(prop_tab, file.path(pub_dir, "table_celltype_proportion_by_stage.csv"))
cat("Cell type proportion (%) by stage:\n")
print(round(prop_tab, 2))

# (d) argmax of mean expression
compartment_assignment <- data.frame(
  gene            = rownames(avg_ct),
  top_compartment = colnames(avg_ct)[apply(avg_ct, 1, which.max)],
  top_mean        = apply(avg_ct, 1, max),
  second_mean     = apply(avg_ct, 1, function(v) sort(v, decreasing = TRUE)[2])
)
compartment_assignment$specificity_ratio <-
  compartment_assignment$top_mean / (compartment_assignment$second_mean + 1e-9)
compartment_assignment <- merge(compartment_assignment,
                                unique(gene_status[, c("gene", "gene_set")]),
                                by = "gene", all.x = TRUE)
write.csv(compartment_assignment,
          file.path(pub_dir, "table_gene_compartment_assignment.csv"), row.names = FALSE)
cat("\n--- Compartment assignment by gene set ---\n")
print(with(compartment_assignment, table(gene_set, top_compartment)))

# (e) stacked bar of composition by stage
prop_df <- as.data.frame(table(
  Stage    = merged_seurat_obj$Tumor_Stage,
  CellType = merged_seurat_obj$cell_type_fixed
)) %>%
  group_by(Stage) %>%
  mutate(Proportion = Freq / sum(Freq) * 100) %>%
  ungroup()

p_prop <- ggplot(prop_df, aes(x = Stage, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.7, colour = "white", linewidth = 0.3) +
  labs(title = "Cell Type Composition by Tumour Stage",
       x = "Tumour Stage", y = "Proportion (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial"))
print(p_prop)
save_pub_plot("Barplot_celltype_proportion_by_stage", p_prop, width = 10, height = 8)