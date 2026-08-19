# Load required libraries
install.packages("devtools")
devtools::install_github("r-lib/conflicted")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("MOFA2", force = TRUE)
library(MOFA2)
library(data.table)
library(ggplot2)
library(tidyverse)

data_list <- list(
  mRNA = rna_vst,
  miRNA = combat_counts,
  Methylation = m_combat,
  CNV =  cnv_imputed,
  Protein = prot_clean
)


# Define the range of factors to run
factor_numbers <- c(10, 15, 20, 25, 30)

# Loop through the desired number of factors
for (n_factors in factor_numbers) {
  cat(paste0("Running MOFA with ", n_factors, " factors...\n"))
  
  # Create MOFA object
  MOFAobject <- create_mofa(data_list)
  
  # Data options (using defaults)
  data_opts <- get_default_data_options(MOFAobject)
  
  # Model options
  model_opts <- get_default_model_options(MOFAobject)
  model_opts$num_factors <- n_factors # Set the number of factors
  model_opts$ard_factors <- TRUE
  model_opts$spikeslab_weights <- TRUE
  
  # Training options
  train_opts <- get_default_training_options(MOFAobject)
  train_opts$convergence_mode <- "slow"
  train_opts$drop_factor_threshold <- -1
  train_opts$seed <- 42 # For reproducibility
  train_opts$freqELBO <- 1
  
  # Prepare MOFA object
  MOFAobject <- prepare_mofa(
    MOFAobject,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )
  
  # Run MOFA
  outfile <- paste0( "../results/mofa_model_", n_factors, "_factors.hdf5")
  MOFAobject <- run_mofa(
    MOFAobject,
    outfile = outfile,
    use_basilisk = TRUE
  )
  
  cat(paste0("MOFA with ", n_factors, " factors completed. Model saved to: ", outfile, "\n\n"))
  gc()
}

cat("All MOFA runs completed.\n")



model_10 <- load_model("../results/mofa_model_10_factors.hdf5")
model_15 <- load_model("../results/mofa_model_15_factors.hdf5")
model_20 <- load_model("../results/mofa_model_20_factors.hdf5")
model_25 <- load_model("../results/mofa_model_25_factors.hdf5")
model_30 <- load_model("../results/mofa_model_30_factors.hdf5")

plot_factor_cor(model_10)
plot_factor_cor(model_15)
plot_factor_cor(model_20)
plot_factor_cor(model_25)
plot_factor_cor(model_30)


plot_weights(model_10)
plot_weights(model_15)
plot_weights(model_20)
plot_weights(model_25)
plot_weights(model_30)

plot_variance_explained(model_10)
plot_variance_explained(model_15)
plot_variance_explained(model_20)
plot_variance_explained(model_25)
plot_variance_explained(model_30)

plot_convergence(model_10)
plot_convergence(model_15)
plot_convergence(model_20)
plot_convergence(model_25)
plot_convergence(model_30)

Nsamples <- sum(model_20@dimensions$N)  # Get the number of samples
colnames(clin_common)[colnames(clin_common) == "sample_id"] <- "sample"
samples_metadata(model_20) <- clin_common

plot_factors(model_20, factors = c(1,2), color_by = "stage_classification")
head(clin_common)

plot_data_heatmap(model_20,
                  view = "mRNA",
                  factor = 1,
                  features = 20,
                  cluster_rows = TRUE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE
)

plot_data_heatmap(model_20,
                  view = "miRNA",
                  factor = 1,
                  features = 20,
                  cluster_rows = TRUE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE
)

plot_data_heatmap(model_20,
                  view = "Methylation",
                  factor = 1,
                  features = 20,
                  cluster_rows = TRUE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE
)

plot_data_heatmap(model_20,
                  view = "CNV",
                  factor = 1,
                  features = 20,
                  cluster_rows = TRUE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE
)

plot_data_heatmap(model_20,
                  view = "Protein",
                  factor = 1,
                  features = 20,
                  cluster_rows = TRUE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE
)


plot_top_weights(model_10,
                 view = "mRNA",
                 factor = 1,
                 nfeatures = 10
)

plot_top_weights(model_20,
                 view = "miRNA",
                 factor = 1,
                 nfeatures = 10
)
plot_data_scatter(model_20,
                  view = "mRNA",         # view of interest
                  factor = 1,             # factor of interest
                  features = 5,           # number of features to plot (they are selected by weight)
                  add_lm = TRUE,          # add linear regression
                  color_by = "stage_classification"
)

set.seed(42)
model <- run_umap(model_20)
model <- run_tsne(model_20)

plot_dimred(model,
            method = "UMAP",  # method can be either "TSNE" or "UMAP"
            color_by = "stage_classification"
)
plot_factors(model_20, 
             factors = 1:5,
             color_by = "stage_classification"
)
plot_factors(model_20, 
             factors = 6:10,
             color_by = "stage_classification"
)
plot_factors(model_20, 
             factors = 11:15,
             color_by = "stage_classification"
)
plot_factors(model_20, 
             factors = 16:20,
             color_by = "stage_classification"
)




factors <- get_factors(model_20, factors = "all")[[1]]
head(factors)

weights <- get_weights(model_20)
head(weights)

# Save mRNA feature loadings to a CSV file
mRNA_loadings <- as.data.frame(weights$mRNA)
head(mRNA_loadings)


# Load biomaRt
library(biomaRt)

# Connect to Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")  # Human genes

# Extract Ensembl IDs from row names
ensembl_ids <- rownames(mRNA_loadings)

# Remove "_mRNA" suffix and version number (everything after ".")
ensembl_ids_clean <- gsub("_mRNA", "", ensembl_ids)  
ensembl_ids_clean <- gsub("\\..*", "", ensembl_ids_clean)  

# Retrieve gene symbols
gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"), 
  filters = "ensembl_gene_id", 
  values = ensembl_ids_clean, 
  mart = ensembl
)

# Create a mapping vector
id_to_gene <- setNames(gene_mapping$hgnc_symbol, gene_mapping$ensembl_gene_id)

# Replace Ensembl IDs with Gene Symbols (if available), else keep Ensembl ID
gene_symbols <- ifelse(id_to_gene[ensembl_ids_clean] == "" | is.na(id_to_gene[ensembl_ids_clean]), 
                       ensembl_ids_clean, id_to_gene[ensembl_ids_clean])

# Add the Gene Symbol column
mRNA_loadings$Gene_Symbol <- gene_symbols

# Remove row names
rownames(mRNA_loadings) <- NULL  

# Drop the first Gene Symbol column if duplicate exists
if ("Gene_Symbol.1" %in% colnames(mRNA_loadings)) {
  mRNA_loadings <- mRNA_loadings[, !colnames(mRNA_loadings) %in% "Gene_Symbol.1"]
}

# Move "Gene_Symbol" to the first column
mRNA_loadings <- mRNA_loadings[, c("Gene_Symbol", setdiff(colnames(mRNA_loadings), "Gene_Symbol"))]

# View updated dataframe
head(mRNA_loadings)


miRNA_loadings <- as.data.frame(weights$miRNA)
head(miRNA_loadings)
# Extract row names
rownames_fixed <- gsub("_miRNA", "", rownames(miRNA_loadings))

# Assign the cleaned row names
rownames(miRNA_loadings) <- rownames_fixed

# View updated dataframe
head(miRNA_loadings)


methyl_loadings <- as.data.frame(weights$Methylation)
head(methyl_loadings)
# Extract row names
rownames_fixed <- gsub("_Methylation", "", rownames(methyl_loadings))

# Assign the cleaned row names
rownames(methyl_loadings) <- rownames_fixed

# View updated dataframe
head(methyl_loadings)


CNV_loadings <- as.data.frame(weights$CNV)
head(CNV_loadings)
# Extract row names
rownames_fixed <- gsub("_CNV", "", rownames(CNV_loadings))

# Assign the cleaned row names
rownames(CNV_loadings) <- rownames_fixed

# View updated dataframe
head(CNV_loadings)


Prot_loadings <- as.data.frame(weights$Protein)
head(Prot_loadings)
# Extract row names
rownames_fixed <- gsub("_Protein", "", rownames(Prot_loadings))

# Assign the cleaned row names
rownames(Prot_loadings) <- rownames_fixed

head(factors)

# Save the dataframes as CSV files in R
write.csv(mRNA_loadings, "mRNA_feature_loadings.csv", row.names = FALSE)
write.csv(miRNA_loadings, "miRNA_feature_loadings.csv", row.names = FALSE)
write.csv(methyl_loadings, "methylation_feature_loadings.csv", row.names = FALSE)
write.csv(CNV_loadings, "CNV_feature_loadings.csv", row.names = FALSE)
write.csv(Prot_loadings, "protein_feature_loadings.csv", row.names = FALSE)
write.csv(factors, "latent_factors.csv", row.names = FALSE)

factors <- get_factors(model_20, factors = "all")[[1]]
head(factors)

weights <- get_weights(model_20)
head(weights)

# Save mRNA feature loadings to a CSV file
mRNA_loadings <- as.data.frame(weights$mRNA)
head(mRNA_loadings)


# Load biomaRt
library(biomaRt)

# Connect to Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")  # Human genes

# Extract Ensembl IDs from row names
ensembl_ids <- rownames(mRNA_loadings)

# Remove "_mRNA" suffix and version number (everything after ".")
ensembl_ids_clean <- gsub("_mRNA", "", ensembl_ids)  
ensembl_ids_clean <- gsub("\\..*", "", ensembl_ids_clean)  

# Retrieve gene symbols
gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"), 
  filters = "ensembl_gene_id", 
  values = ensembl_ids_clean, 
  mart = ensembl
)

# Create a mapping vector
id_to_gene <- setNames(gene_mapping$hgnc_symbol, gene_mapping$ensembl_gene_id)

# Replace Ensembl IDs with Gene Symbols (if available), else keep Ensembl ID
gene_symbols <- ifelse(id_to_gene[ensembl_ids_clean] == "" | is.na(id_to_gene[ensembl_ids_clean]), 
                       ensembl_ids_clean, id_to_gene[ensembl_ids_clean])

# Add the Gene Symbol column
mRNA_loadings$Gene_Symbol <- gene_symbols

# Remove row names
rownames(mRNA_loadings) <- NULL  

# Drop the first Gene Symbol column if duplicate exists
if ("Gene_Symbol.1" %in% colnames(mRNA_loadings)) {
  mRNA_loadings <- mRNA_loadings[, !colnames(mRNA_loadings) %in% "Gene_Symbol.1"]
}

# Move "Gene_Symbol" to the first column
mRNA_loadings <- mRNA_loadings[, c("Gene_Symbol", setdiff(colnames(mRNA_loadings), "Gene_Symbol"))]

# View updated dataframe
head(mRNA_loadings)


miRNA_loadings <- as.data.frame(weights$miRNA)
head(miRNA_loadings)
# Extract row names
rownames_fixed <- gsub("_miRNA", "", rownames(miRNA_loadings))

# Assign the cleaned row names
rownames(miRNA_loadings) <- rownames_fixed

# View updated dataframe
head(miRNA_loadings)


methyl_loadings <- as.data.frame(weights$Methylation)
head(methyl_loadings)
# Extract row names
rownames_fixed <- gsub("_Methylation", "", rownames(methyl_loadings))

# Assign the cleaned row names
rownames(methyl_loadings) <- rownames_fixed

# View updated dataframe
head(methyl_loadings)


CNV_loadings <- as.data.frame(weights$CNV)
head(CNV_loadings)
# Extract row names
rownames_fixed <- gsub("_CNV", "", rownames(CNV_loadings))

# Assign the cleaned row names
rownames(CNV_loadings) <- rownames_fixed

# View updated dataframe
head(CNV_loadings)


Prot_loadings <- as.data.frame(weights$Protein)
head(Prot_loadings)
# Extract row names
rownames_fixed <- gsub("_Protein", "", rownames(Prot_loadings))

# Assign the cleaned row names
rownames(Prot_loadings) <- rownames_fixed

head(factors)

# Save the dataframes as CSV files in R
write.csv(mRNA_loadings, "mRNA_feature_loadings.csv", row.names = TRUE)
write.csv(miRNA_loadings, "miRNA_feature_loadings.csv", row.names = TRUE)
write.csv(methyl_loadings, "methylation_feature_loadings.csv", row.names = TRUE)
write.csv(CNV_loadings, "CNV_feature_loadings.csv", row.names = TRUE)
write.csv(Prot_loadings, "protein_feature_loadings.csv", row.names = TRUE)
write.csv(factors, "latent_factors.csv", row.names = TRUE)
write.csv(clin_data, "clinical_20factor.csv", row.names = TRUE)






