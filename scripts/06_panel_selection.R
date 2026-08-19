# Load required libraries
library(readr)
library(dplyr)
library(limma)
library(tidyr)

if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")

library(pROC)


# File paths
RNA_PATH <- "../data/common_samples/rna_common.csv"
CLIN_PATH <- "../data/common_samples/clin_common.csv"
CANDIDATES_PATH <- "../results/mofa_mrna_candidates.csv"
LOADINGS_PATH <- "../results/mRNA_feature_loadings.csv"
OUTPUT_DIR <- "../results/tiered_panel"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)


# Thresholds
T1_FDR_Q <- 0.05
T1_MIN_LFC <- 0.5

T2_RAW_P <- 0.05
T2_MIN_LFC <- 0.5
T2_MIN_AUC <- 0.65

T3_MIN_LOAD <- 0.3
T3_TOP_N <- 20

VALIDATED_FACTORS <- c("Factor2", "Factor5", "Factor6", "Factor7", "Factor9")

EXCLUDE <- c("EPYC", "NETO1")


# Load data
candidates <- read_csv(CANDIDATES_PATH, show_col_types = FALSE)
candidate_genes <- candidates$Gene

loadings <- read_csv(LOADINGS_PATH, show_col_types = FALSE)

clin <- read_csv(CLIN_PATH, show_col_types = FALSE)

rna <- read.csv(RNA_PATH, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)


# Prepare clinical data
if ("sample_id" %in% colnames(clin)) clin$sid <- substr(clin$sample_id, 1, 15) else clin$sid <- substr(clin[[1]], 1, 15)

clin_sub <- clin %>% select(sid, stage_classification) %>% filter(!is.na(stage_classification))


# Prepare RNA data
if ("sample_id" %in% colnames(rna)) rna$sid <- substr(rna$sample_id, 1, 15) else rna$sid <- substr(rna[[1]], 1, 15)


# Merge RNA and clinical data
merged <- inner_join(rna, clin_sub, by = "sid") %>% filter(!is.na(stage_classification))


# Define stage groups
early_mask <- merged$stage_classification == "Early"
late_mask <- merged$stage_classification == "Late"
normal_mask <- merged$stage_classification == "Normal"


# Check available genes
available <- intersect(candidate_genes, colnames(merged))

cat("Candidate genes:", length(candidate_genes), "\n")
cat("Available genes:", length(available), "\n")
cat("Stage distribution:\n")
print(table(merged$stage_classification))


# limma function
run_limma <- function(data, genes, mask1, mask2, label1, label2) {
  
  sub <- data[mask1 | mask2, ]
  sub_m1 <- mask1[mask1 | mask2]
  expr_df <- sub[, genes, drop = FALSE]
  expr_mat <- t(apply(expr_df, 1, function(x) as.numeric(as.character(x))))
  colnames(expr_mat) <- genes
  expr_mat <- t(expr_mat)
  
  grp <- factor(as.integer(sub_m1), levels = c(0, 1), labels = c(label2, label1))
  
  design <- model.matrix(~ grp)
  
  coef_name <- paste0(label1, "_vs_", label2)
  
  colnames(design) <- c("Intercept", coef_name)
  
  fit <- lmFit(expr_mat, design)
  
  fit <- eBayes(fit)
  
  results <- topTable(fit, coef = coef_name, number = Inf, adjust.method = "BH", sort.by = "P")
  
  results$Gene <- rownames(results)
  
  results %>% rename(log2FC = logFC, p_value = P.Value, FDR_q = adj.P.Val) %>% select(Gene, log2FC, p_value, FDR_q) %>% mutate(Direction = ifelse(log2FC > 0, paste0("Higher in ", label1), paste0("Higher in ", label2)))
}


# AUC function
compute_auc <- function(g1_expr, g2_expr) {
  
  y <- c(rep(1, length(g1_expr)), rep(0, length(g2_expr)))
  x <- c(g1_expr, g2_expr)
  
  valid <- !is.na(x)
  
  if (sum(valid) < 10 || length(unique(y[valid])) < 2) return(NA_real_)
  
  tryCatch({
    
    roc_obj <- roc(y[valid], x[valid], quiet = TRUE)
    
    auc_value <- as.numeric(auc(roc_obj))
    
    max(auc_value, 1 - auc_value)
    
  }, error = function(e) NA_real_)
}


# Run limma
cat("Running Early vs Late\n")

lim_EL <- run_limma(merged, available, early_mask, late_mask, "Early", "Late")

cat("Running Early vs Normal\n")

lim_EN <- run_limma(merged, available, early_mask, normal_mask, "Early", "Normal")

cat("Running Late vs Normal\n")

lim_LN <- run_limma(merged, available, late_mask, normal_mask, "Late", "Normal")


# Calculate AUC
auc_results <- lapply(available, function(gene) {
  
  expr <- as.numeric(merged[[gene]])
  
  data.frame(Gene = gene, AUC_EL = compute_auc(expr[early_mask], expr[late_mask]), AUC_EN = compute_auc(expr[early_mask], expr[normal_mask]), AUC_LN = compute_auc(expr[late_mask], expr[normal_mask]))
})

auc_df <- bind_rows(auc_results)


# Rename comparison columns
lim_EL <- lim_EL %>% rename(log2FC_EL = log2FC, p_EL = p_value, q_EL = FDR_q, dir_EL = Direction)

lim_EN <- lim_EN %>% rename(log2FC_EN = log2FC, p_EN = p_value, q_EN = FDR_q, dir_EN = Direction)

lim_LN <- lim_LN %>% rename(log2FC_LN = log2FC, p_LN = p_value, q_LN = FDR_q, dir_LN = Direction)


# Build evidence table
evidence <- lim_EL %>% left_join(lim_EN, by = "Gene") %>% left_join(lim_LN, by = "Gene") %>% left_join(auc_df, by = "Gene") %>% left_join(candidates %>% select(Gene, N_Factors, Max_Abs_Loading, Factor_Loadings), by = "Gene") %>% filter(!Gene %in% EXCLUDE)


# Tier 1
evidence <- evidence %>% mutate(
  t1_EL = !is.na(q_EL) & q_EL < T1_FDR_Q & abs(log2FC_EL) >= T1_MIN_LFC,
  t1_EN = !is.na(q_EN) & q_EN < T1_FDR_Q & abs(log2FC_EN) >= T1_MIN_LFC,
  t1_LN = !is.na(q_LN) & q_LN < T1_FDR_Q & abs(log2FC_LN) >= T1_MIN_LFC,
  Tier1 = t1_EL | t1_EN | t1_LN,
  T1_best_comparison = case_when(t1_EL ~ "Early vs Late", t1_EN ~ "Early vs Normal", t1_LN ~ "Late vs Normal", TRUE ~ NA_character_)
)

tier1 <- evidence %>% filter(Tier1)


# Tier 2
evidence <- evidence %>% mutate(
  t2_EL = !Tier1 & !is.na(p_EL) & p_EL < T2_RAW_P & abs(log2FC_EL) >= T2_MIN_LFC & !is.na(AUC_EL) & AUC_EL >= T2_MIN_AUC,
  t2_EN = !Tier1 & !is.na(p_EN) & p_EN < T2_RAW_P & abs(log2FC_EN) >= T2_MIN_LFC & !is.na(AUC_EN) & AUC_EN >= T2_MIN_AUC,
  t2_LN = !Tier1 & !is.na(p_LN) & p_LN < T2_RAW_P & abs(log2FC_LN) >= T2_MIN_LFC & !is.na(AUC_LN) & AUC_LN >= T2_MIN_AUC,
  Tier2 = !Tier1 & (t2_EL | t2_EN | t2_LN),
  T2_best_comparison = case_when(t2_EL ~ "Early vs Late", t2_EN ~ "Early vs Normal", t2_LN ~ "Late vs Normal", TRUE ~ NA_character_)
)

tier2 <- evidence %>% filter(Tier2)


# Tier 3
t3_genes <- c()

for (factor in VALIDATED_FACTORS) {
  
  factor_genes <- loadings %>% select(Gene_Symbol, all_of(factor)) %>% rename(Gene = Gene_Symbol, loading = all_of(factor)) %>% mutate(abs_load = abs(loading)) %>% arrange(desc(abs_load)) %>% head(T3_TOP_N) %>% filter(abs_load >= T3_MIN_LOAD)
  
  t3_genes <- c(t3_genes, factor_genes$Gene)
}

t3_genes <- unique(t3_genes)

t1_t2_genes <- c(tier1$Gene, tier2$Gene)

tier3_candidates <- evidence %>% filter(Gene %in% t3_genes, !Gene %in% t1_t2_genes, !Gene %in% EXCLUDE)


# Add loading information
tier3_load <- data.frame(Gene = t3_genes) %>% left_join(loadings %>% select(Gene_Symbol, all_of(VALIDATED_FACTORS)) %>% rename(Gene = Gene_Symbol), by = "Gene") %>% rowwise() %>% mutate(
  Best_factor = VALIDATED_FACTORS[which.max(abs(c_across(all_of(VALIDATED_FACTORS))))],
  Max_loading = max(abs(c_across(all_of(VALIDATED_FACTORS))), na.rm = TRUE)
) %>% ungroup()

tier3 <- tier3_candidates %>% left_join(tier3_load %>% select(Gene, Best_factor, Max_loading), by = "Gene") %>% arrange(desc(Max_loading))


# Combine tiers
tier1_out <- tier1 %>% mutate(Tier = "Tier 1", Tier_num = 1)

tier2_out <- tier2 %>% mutate(Tier = "Tier 2", Tier_num = 2)

tier3_out <- tier3 %>% mutate(Tier = "Tier 3", Tier_num = 3)


# Final panel
full_panel <- bind_rows(tier1_out, tier2_out, tier3_out) %>% arrange(Tier_num, desc(Max_Abs_Loading))


# Save outputs
write_csv(tier1_out, file.path(OUTPUT_DIR, "tier1_confirmed.csv"))

write_csv(tier2_out, file.path(OUTPUT_DIR, "tier2_plausible.csv"))

write_csv(tier3_out, file.path(OUTPUT_DIR, "tier3_mofa_nominated.csv"))

write_csv(full_panel, file.path(OUTPUT_DIR, "full_panel_all_tiers.csv"))

write_csv(evidence, file.path(OUTPUT_DIR, "evidence_all_genes.csv"))


# Print results
cat("\nTier 1:", nrow(tier1), "genes\n")

cat("Tier 2:", nrow(tier2), "genes\n")

cat("Tier 3:", nrow(tier3), "genes\n")

cat("Total panel:", nrow(full_panel), "genes\n")