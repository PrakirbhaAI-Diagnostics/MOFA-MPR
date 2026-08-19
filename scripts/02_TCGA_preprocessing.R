load("../data/processed/RNA.rda")


# ==============================================
# Step 1: Load necessary libraries and functions
# ==============================================
library(caTools)     # For stratified train-test split
library(tidyverse)   # For data wrangling
library(DESeq2)      # For RNA processing
library(edgeR)       # For miRNA processing
library(limma)       # For voom transformation
library(minfi)       # For methylation processing
library(sva)         # For batch correction
library(wateRmelon)  # For methylation normalization
library(impute)      # For missing value imputation
library(preprocessCore) # For protein normalization
# ==============================================
# Step 2: Load clinical data and perform train-test split
# ==============================================
# Assuming clin_data has columns: "sample_id" and "stage_classification"
clin_data <- read.csv("../metadata/common_samples_final.csv")
head(clin_data)

# Standardize sample barcodes
clin_data$sample_id <- standardize_barcodes(clin_data$sample_id)


# ==============================================
# mRNA Preprocessing 
# ==============================================

# 1. Load Libraries
library(tidyverse)
library(DESeq2)
library(sva)  # For ComBat batch correction

# 2. Load Data
load("../data/processed/RNA.rda")
rna <- data
rm(data)

# 3. Extract Counts & Metadata
counts <- assay(rna, "unstranded")
row_data <- rowData(rna)
col_data <- colData(rna)

# 4. Create DESeq2 Object with Batch Info
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = col_data,
  design = ~1
)

# 5. Pre-filtering (More stringent)
keep <- rowSums(counts(dds) >= 10) >= ncol(dds)*0.2  # Keep genes with ≥10 counts in ≥20% samples
dds <- dds[keep, ]

# 6. Variance Stabilizing Transformation
dds <- estimateSizeFactors(dds)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
exprs_vst <- assay(vsd)

# 7. Batch Correction (Critical Fix!)
batch <- substr(colnames(exprs_vst), 22, 25)  # Extract plate ID from TCGA barcode
exprs_combat <- ComBat(
  dat = exprs_vst,
  batch = batch,
  par.prior = TRUE
)

# 8. Map to Gene Symbols
gene_symbols <- row_data$gene_name[match(rownames(exprs_combat), rownames(row_data))]
rownames(exprs_combat) <- ifelse(
  is.na(gene_symbols),
  rownames(exprs_combat),  # Keep Ensembl ID if no symbol
  gene_symbols
)

# 9. Remove Duplicate Genes 
exprs_combat <- exprs_combat[!duplicated(rownames(exprs_combat)), ]

# 10. Scaling 
exprs_scaled <- t(scale(t(exprs_combat)))  # Features x samples, mean=0, var=1

# 11. Save Processed Data
save(exprs_combat, file = "../data/processed/rna_processed_final.rda")

# 12. Quality Checks
cat("Final dimensions:", dim(exprs_scaled), "\n")
cat("Missing values:", sum(is.na(exprs_scaled)), "\n")
boxplot(exprs_scaled, main="Scaled mRNA Expression")

# Update boxplot to save as TIFF
boxplot(exprs_scaled, main="Scaled mRNA Expression")
tiff("plots/scaled_mrna_expression.tiff", width=6, height=6, units="in", res=600)
boxplot(exprs_scaled, main="Scaled mRNA Expression", col="lightblue", las=2, cex.axis=1.2, cex.lab=1.5, font.lab=2, family="Arial")
dev.off()

# ==============================================
# miRNA Preprocessing (Corrected for MOFA)
# ==============================================
# 1. Load Libraries
library(tidyverse)
library(edgeR)
library(limma)
library(sva) 

# 2. Load Raw Data
load("../data/processed/miRNA.rda")
mirna <- data
rm(data)


library(miRBaseConverter)
library(dplyr)

# 3. Process miRNA Data with Proper Mature miRNA Mapping
# ------------------------------------------------------
# Get miRBase v22 annotations (human mature miRNAs only)
mirbase_db <- getAllMiRNAs(version = "v22")
mature_mirnas <- mirbase_db %>%
  filter(
    grepl("^hsa", Name),  # Human miRNAs start with "hsa"
    grepl("^MIMAT", Accession)  # Mature miRNA accessions start with "MIMAT"
  ) %>%
  dplyr::select(  # Explicitly specify dplyr::select
    mimat_id = Accession,
    mature_miRNA = Name
  )

# 4. Create Expression Matrix with Strict Filtering
# ------------------------------------------------------
mirna_processed <- mirna %>%
  # Keep only mature miRNAs with valid annotations
  filter(grepl("^mature", miRNA_region),
         `cross-mapped` == "N") %>%
  separate(
    miRNA_region,
    into = c("type", "mimat_id"),
    sep = ",", 
    extra = "drop"
  ) %>%
  left_join(mature_mirnas, by = "mimat_id") %>%
  filter(!is.na(mature_miRNA)) %>%
  
  # Aggregate counts properly (Critical Fix!)
  group_by(mature_miRNA, barcode) %>%
  summarise(total_count = sum(read_count), .groups = "drop") %>%
  
  # Pivot to matrix format
  pivot_wider(
    names_from = barcode,
    values_from = total_count,
    values_fill = 0
  ) %>%
  column_to_rownames("mature_miRNA") %>%
  as.matrix()

# 5. edgeR Normalization with voom (Critical Fix!)
# ------------------------------------------------------
# Keep miRNAs expressed in ≥20% samples (CPM > 1)
dge <- DGEList(counts = mirna_processed)
keep <- rowSums(cpm(dge) > 1) >= ncol(mirna_processed)*0.2
dge <- dge[keep, , keep.lib.sizes = FALSE]

# TMM normalization + voom transformation
dge <- calcNormFactors(dge)
design <- model.matrix(~1, data = dge$samples)
v <- voom(dge, design=design, plot=FALSE)  # Proper variance stabilization

# 6. Batch Correction (Added Step!)
# ------------------------------------------------------
batch <- substr(colnames(v$E), 22, 25)  # Match mRNA batch ID positions
combat_counts <- ComBat(
  dat = v$E,
  batch = batch,
  par.prior = TRUE
)

voom_plot <- voom(dge, design=design, plot=TRUE)  

# 7. Scaling (Essential for MOFA!)
# ------------------------------------------------------
mirna_scaled <- t(scale(t(combat_counts))) 

# 8. Remove Duplicate miRNAs (Critical Fix!)
# ------------------------------------------------------
mirna_final <- mirna_scaled[!duplicated(rownames(mirna_scaled)), ]

# 9. Save Processed Data
# ------------------------------------------------------
save(combat_counts, file = "../data/processed/mirna_processed_final.rda")

# 10. Quality Checks
# ------------------------------------------------------
cat("Final dimensions:", dim(mirna_final), "\n")


# ==============================================
# Methylation Preprocessing 
# ==============================================

# 1. Load Libraries
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(wateRmelon)
library(sva)
library(tidyverse)

# 2. Load Data & Initial Processing
load("../data/processed/met.rda")
methyl <- data
rm(data)
beta_matrix <- assay(methyl, 1)

# 3. Probe Filtering (Enhanced)
# -------------------------------
# Load cross-reactive probes (Pidsley et al. 2016)
pidsley <- read.csv("../metadata/13059_2016_1066_MOESM1_ESM.csv")
cross_reactive <- pidsley$X[pidsley$Total > 0]

# Get annotation data
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Filter in one step to preserve matrix structure
bad_probes <- unique(c(
  cross_reactive,
  rownames(anno)[anno$chr %in% c("chrX", "chrY")],  # Sex chromosomes
  rownames(anno)[grepl("rs", anno$Probe_rs)],        # SNP-associated probes
  rownames(anno)[anno$Islands_Name == ""]            # Non-CGI probes
))

beta_filtered <- beta_matrix[!rownames(beta_matrix) %in% bad_probes, ]

# 4. Missing Value Handling 
# -------------------------------
# Remove probes with >10% missingness
beta_filtered <- beta_filtered[rowMeans(is.na(beta_filtered)) < 0.1, ]
set.seed(123)
# Impute remaining NAs using k-NN
beta_imputed <- impute.knn(beta_filtered, k = 5)$data

# 5. BMIQ Normalization
# -------------------------------
anno_filtered <- anno[rownames(beta_imputed), ]
design_vector <- ifelse(anno_filtered$Type == "I", 1, 2)
design_vector <- ifelse(anno_filtered$Type == "I", 1, 2)
beta_normalized <- BMIQ(beta_imputed, design.v = design_vector)$nbeta

# 6. M-value Conversion with Proper Thresholding
# -------------------------------
m_values <- log2(beta_normalized/(1 - beta_normalized))

# Winsorize extreme values (preserve biological variation)
upper <- quantile(m_values, 0.9999, na.rm = TRUE)
lower <- quantile(m_values, 0.0001, na.rm = TRUE)
m_values[m_values > upper] <- upper
m_values[m_values < lower] <- lower

# 7. Batch Correction (Enhanced)
# -------------------------------
# Extract batch info (plate + position for finer control)
batch <- paste(
  substr(colnames(m_values), 22, 25))  # Plate ID

#table(batch)
# Combat with model matrix (preserve biological variance)
mod <- model.matrix(~1, data = data.frame(colnames(m_values)))
m_combat <- ComBat(
  dat = m_values,
  batch = batch,
  mod = mod,
  par.prior = TRUE
)

# 8. Feature Scaling (Critical for MOFA!)
# -------------------------------
m_scaled <- t(scale(t(m_combat)))  

# 9. QC Checks
# -------------------------------
# Density plots
par(mfrow = c(2,2))
dev.off()
densityPlot(beta_matrix, main = "Raw Beta Values")
densityPlot(beta_normalized, main = "Normalized Beta Values")
densityPlot(m_combat, main = "Batch Corrected M-values")
densityPlot(m_scaled, main = "Scaled M-values")

pca <- prcomp(t(m_values), scale. = FALSE)
plot(pca$x[,1:2], col = as.factor(batch), 
     main = "PCA: Batch Effects before Correction")
# PCA for batch effect check
pca <- prcomp(t(m_scaled), scale. = FALSE)
plot(pca$x[,1:2], col = as.factor(batch), 
     main = "PCA: Batch Effects After Correction")

# Update density plots to save as TIFF
par(mfrow = c(2,2))
tiff("plots/density_plots.tiff", width=8, height=8, units="in", res=600)
densityPlot(beta_matrix, main = "Raw Beta Values", family="Arial", font.main=2)
densityPlot(beta_normalized, main = "Normalized Beta Values", family="Arial", font.main=2)
densityPlot(m_combat, main = "Batch Corrected M-values", family="Arial", font.main=2)
densityPlot(m_scaled, main = "Scaled M-values", family="Arial", font.main=2)
dev.off()

# Update PCA plots to save as TIFF
tiff("plots/pca_batch_effects_before.tiff", width=6, height=6, units="in", res=600)
plot(pca$x[,1:2], col = as.factor(batch), main = "PCA: Batch Effects before Correction", pch=19, cex=1.2, family="Arial", font.main=2)
dev.off()

tiff("plots/pca_batch_effects_after.tiff", width=6, height=6, units="in", res=600)
plot(pca$x[,1:2], col = as.factor(batch), main = "PCA: Batch Effects After Correction", pch=19, cex=1.2, family="Arial", font.main=2)
dev.off()

# 10. Save Processed Data
# -------------------------------
save(m_combat, file = "../data/processed/methyl_processed_final.rda")



# ==============================================
# CNV Preprocessing
# ==============================================

# 1. Load Libraries
library(SummarizedExperiment)
library(sva)
library(impute)
library(org.Hs.eg.db)
library(dplyr)
library(stringr)
# Load Data
load("../data/processed/cnv.rda")
cnv <- data

# Clean CNV Data
gistic_seg_clean <- cnv %>%
  mutate(Sample = substr(Sample, 1, 15),
         Chromosome = case_when(Chromosome == "X" ~ "chrX",
                                Chromosome == "Y" ~ "chrY",
                                TRUE ~ paste0("chr", str_remove(Chromosome, "chr")))) %>%
  dplyr::select(Sample, Chromosome, Start, End, Num_Probes, Segment_Mean)

# Convert to GenomicRanges
cnv_gr <- makeGRangesFromDataFrame(gistic_seg_clean, start.field = "Start", end.field = "End", keep.extra.columns = TRUE)

# Load genes and extract from TxDb
genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Find Overlaps between CNV regions and genes
overlaps <- findOverlaps(cnv_gr, genes)
gene_cnv <- data.frame(
  gene_id = names(genes)[subjectHits(overlaps)],
  sample = mcols(cnv_gr)$Sample[queryHits(overlaps)],
  segment_mean = mcols(cnv_gr)$Segment_Mean[queryHits(overlaps)]
) %>%
  group_by(gene_id, sample) %>%
  summarise(cnv_score = median(segment_mean), .groups = "drop")

# Use org.Hs.eg.db to map gene_id (ENSEMBL) to gene symbols
library(org.Hs.eg.db)

# Convert gene_id (ENSEMBL/Entrez) to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys = gene_cnv$gene_id, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

# Replace gene IDs with gene symbols
gene_cnv$gene_symbol <- gene_symbols

# Remove rows with missing gene symbols
gene_cnv <- gene_cnv %>% filter(!is.na(gene_symbol))

# Convert to Wide Format
cnv_matrix <- gene_cnv %>%
  pivot_wider(names_from = sample, values_from = cnv_score, values_fill = 0) %>%
  column_to_rownames("gene_symbol") %>%
  as.matrix()

# Impute Missing Values
set.seed(123)
cnv_imputed <- impute.knn(cnv_matrix, k = 5)$data

# MOFA-Specific Scaling
cnv_scaled <- t(scale(t(cnv_imputed)))

# Final Validation: Ensure no NAs remain
stopifnot(sum(is.na(cnv_scaled)) == 0)

# Check variance
gene_var <- apply(cnv_scaled, 1, var)
hist(gene_var, main = "Post-Scaling Variance")

# Save Data
save(cnv_scaled, file = "../data/processed/CNV_processed_final.rda")


# ==============================================
# PROTEIN PREPROCESSING (MOFA-OPTIMIZED)
# ==============================================

# Load Libraries
library(tidyverse)
library(impute)
library(sva)
library(preprocessCore)

# -------------------------------
# 1. Load and Prepare Data
# -------------------------------
load("../data/processed/prot.rda")
protein <- data
rm(data)

# Load Antibody Annotations
antibody_anno <- read.delim("../metadata/TCGA_antibodies_descriptions.gencode.v36.tsv") %>%
  dplyr::select(
    AGID, 
    Gene_Symbol = "gene_name",
    Peptide_Target = "peptide_target"
  )

# -------------------------------
# 2. Initial Processing
# -------------------------------
# Extract expression matrix
expr_cols <- grep('TCGA', colnames(protein), value = TRUE)
prot_matrix <- protein %>%
  dplyr::select(all_of(expr_cols)) %>%
  as.matrix()

# Create unique row identifiers
rownames(prot_matrix) <- paste(protein$AGID, protein$peptide_target, sep = "_")

# -------------------------------
# 3. Missing Value Handling
# -------------------------------
# Filter features with >10% missingness
missing_percent <- rowMeans(is.na(prot_matrix))
prot_filtered <- prot_matrix[missing_percent <= 0.1, ]

# KNN imputation (k=5)
set.seed(123)
prot_imputed <- impute.knn(prot_filtered, k = 5)$data


# -------------------------------
# 4. Normalization & Batch Correction
# -------------------------------
# Quantile normalization FIRST
prot_normalized <- normalize.quantiles(prot_imputed)
rownames(prot_normalized) <- rownames(prot_imputed)


colnames(prot_normalized) <- colnames(prot_imputed)


# -------------------------------
# 5. Annotation & Aggregation
# -------------------------------
# Map antibody IDs to gene symbols
protein_ids <- data.frame(
  Row_ID = rownames(prot_normalized),
  AGID = sapply(strsplit(rownames(prot_normalized), "_"), "[", 1)
)

annotated_proteins <- protein_ids %>% 
  left_join(antibody_anno, by = "AGID") %>%
  mutate(
    Final_Symbol = case_when(
      !is.na(Peptide_Target) ~ paste0(Gene_Symbol, "_", Peptide_Target),
      is.na(Gene_Symbol) ~ paste0("AG_", AGID),
      TRUE ~ Gene_Symbol
    )
  )

# Update rownames with enhanced symbols
rownames(prot_normalized) <- annotated_proteins$Final_Symbol

# Aggregate duplicates (handle all-NA cases)
prot_aggregated <- prot_normalized %>%
  as.data.frame() %>%
  rownames_to_column("Symbol") %>%
  group_by(Symbol) %>%
  summarise(across(everything(), ~ {
    if (all(is.na(.x))) NA else mean(.x, na.rm = TRUE)
  })) %>%
  filter(!is.na(Symbol)) %>%
  column_to_rownames("Symbol") %>%
  as.matrix()

# Remove remaining NA rows
prot_clean <- prot_aggregated[complete.cases(prot_aggregated), ]



# Scale features to mean=0, var=1
prot_scaled <- t(scale(t(prot_clean)))

# Save Data
save(prot_clean, file = "../data/processed/protein_processed.rda")