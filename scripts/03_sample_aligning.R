library(data.table)
library(ggplot2)
library(tidyverse)
library(biomaRt)

# Load all processed datasets# Load all processed datasets
load("../data/processed/rna_processed_final.rda")       # RNA
load("../data/processed/mirna_processed_final.rda") # miRNA
load("../data/processed/methyl_processed_final.rda") # Methylation
load("../data/processed/CNV_processed_final.rda")   # CNV
load("../data/processed/protein_processed.rda")     # Protein

# Create standardized 15-character barcodes for each modality
standardize_barcodes <- function(barcodes) {
  substr(barcodes, 1, 15)
}

# Apply to all datasets
colnames(cnv_imputed) <- standardize_barcodes(colnames(cnv_imputed))
colnames(m_combat) <- standardize_barcodes(colnames(m_combat))
colnames(combat_counts) <- standardize_barcodes(colnames(combat_counts))
colnames(prot_clean) <- standardize_barcodes(colnames(prot_clean))
colnames(exprs_combat) <- standardize_barcodes(colnames(exprs_combat))

# Get sample lists from all modalities
sample_lists <- list(
  cnv = colnames(cnv_imputed),
  methyl = colnames(m_combat),
  mirna = colnames(combat_counts),
  prot = colnames(prot_clean),
  rna = colnames(exprs_combat)
)

# Find intersection
common_samples <- Reduce(intersect, sample_lists)
cat("Final common samples:", length(common_samples), "\n")

# Subset each modality to common samples
cnv_common <- cnv_imputed[, common_samples, drop = FALSE]
methyl_common <- m_combat[, common_samples, drop = FALSE]
mirna_common <- combat_counts[, common_samples, drop = FALSE]
prot_common <- prot_clean[, common_samples, drop = FALSE]
rna_common <- exprs_combat[, common_samples, drop = FALSE]

# Sort common_samples alphabetically (recommended)
common_samples <- sort(common_samples)

# Reorder all omics datasets to match sorted common_samples
cnv_common <- cnv_common[, common_samples, drop = FALSE]
methyl_common <- methyl_common[, common_samples, drop = FALSE]
mirna_common <- mirna_common[, common_samples, drop = FALSE]
prot_common <- prot_common[, common_samples, drop = FALSE]
rna_common <- rna_common[, common_samples, drop = FALSE]

# Check dimensions
cat("CNV dimensions:", dim(cnv_common), "\n")
cat("Methylation dimensions:", dim(methyl_common), "\n")
cat("miRNA dimensions:", dim(mirna_common), "\n")
cat("Protein dimensions:", dim(prot_common), "\n")
cat("RNA dimensions:", dim(rna_common), "\n")

# Verify sample order (should be TRUE)
all(
  identical(colnames(cnv_common), common_samples),
  identical(colnames(methyl_common), common_samples),
  identical(colnames(mirna_common), common_samples),
  identical(colnames(prot_common), common_samples),
  identical(colnames(rna_common), common_samples)
)

# Load clinical data for consistent splitting
clin_data <- read.csv("../data/metadata/common_samples_final.csv")
print(class(clin_data)) # Should output "data.frame"

# Filter clinical data to keep only rows matching the `common_samples`
clin_common <- clin_data %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, common_samples))

# Check the dimensions of the filtered clinical data
cat("Filtered clinical data dimensions:", dim(clin_common), "\n")

head(clin_common)
table(clin_common$stage_classification)

library(MOFA2)
library(data.table)
library(ggplot2)
library(tidyverse)

setdiff(common_samples, colnames(cnv_common))
setdiff(common_samples, colnames(methyl_common))
setdiff(common_samples, colnames(mirna_common))
setdiff(common_samples, colnames(prot_common))
setdiff(common_samples, colnames(rna_common))

library(ggplot2)
ggplot(clin_common, aes(x = stage_classification, fill = stage_classification)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "Stage Classification Distribution")


feature_names_list <- list(
  mRNA = rownames(rna_common),
  miRNA = rownames(mirna_common),
  Methylation = rownames(methyl_common),
  CNV = rownames(cnv_common),
  Protein = rownames(prot_common)
)

# Find duplicated features across datasets
all_features <- unlist(feature_names_list)
duplicates <- all_features[duplicated(all_features)]

cat("Number of duplicated features:", length(duplicates), "\n")
print(duplicates)




if (anyDuplicated(all_features_renamed) == 0) {
  cat("Feature names are now unique across all datasets.\n")
} else {
  cat("Duplicates still exist! Fix needed.\n")
}

setdiff(common_samples, colnames(cnv_common))
setdiff(common_samples, colnames(methyl_common))
setdiff(common_samples, colnames(mirna_common))
setdiff(common_samples, colnames(rna_common))
setdiff(common_samples, colnames(prot_common))


fwrite(as.data.table(rna_common, keep.rownames = "Feature"), "../data/rna_common.csv")
fwrite(as.data.table(mirna_common, keep.rownames = "Feature"), "../data/mirna_common.csv")
fwrite(as.data.table(methyl_common, keep.rownames = "Feature"), "../data/methyl_common.csv")
fwrite(as.data.table(cnv_common, keep.rownames = "Feature"), "../data/cnv_common.csv")
fwrite(as.data.table(prot_common, keep.rownames = "Feature"), "../data/prot_common.csv")
fwrite(as.data.table(clin_common, keep.rownames = 'sample_id'), "../data/clin_common.csv")

# Now, when creating your data_list for MOFA, use the transposed and renamed RNA data
data_list <- list(
  mRNA = rna_common,
  miRNA = mirna_common,
  Methylation = methyl_common,
  CNV = cnv_common,
  Protein = prot_common
)

