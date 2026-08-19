library(TCGAbiolinks)
library(SummarizedExperiment)

# Load the CSV file using read.csv
clin_data <- read.csv("../metadata/clin_data.csv")
barcodes <-clin_data$barcode


# -----------------------------------------------
# Step 1: Query data
# -----------------------------------------------
project <- "TCGA-BRCA"

# mRNA
query_mrna <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  barcode = barcodes)

# miRNA
query_mirna <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Isoform Expression Quantification",
  workflow.type = "BCGSC miRNA Profiling")

# Methylation
query_methyl <- GDCquery(
  project = project,
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450") #parameter enforced by GDCquery

# CNV
query_cnv <- GDCquery(
  project = project,
  data.category = "Copy Number Variation",
  data.type = "Masked Copy Number Segment")


# Protein Data
query_protein <- GDCquery(
  project = project,
  data.category = "Proteome Profiling",
  data.type = "Protein Expression Quantification",
  experimental.strategy = "Reverse Phase Protein Array")  



# -----------------------------------------------
# Step 2: Find Common Samples
# -----------------------------------------------

# Define samples to exclude
unwanted_barcodes <- c(
  "TCGA-A2-A25A-01A-12R-A16F-07",
  "TCGA-LL-A440-01A-11R-A24H-07",
  "TCGA-S3-AA17-01A-11R-A41B-07",
  "TCGA-GM-A2DO-01A-11R-A18M-07",
  "TCGA-B6-A0RM-01A-11R-A084-07", 'TCGA-AR-A0U1-01A-11R-A115-07', 'TCGA-A7-A0DC-01A-11R-A00Z-07', 'TCGA-BH-A0B2-01A-11R-A10J-07'
)

# Get first 15 characters of unwanted barcodes (to match your common.samples format)
unwanted_samples <- substr(unwanted_barcodes, 1, 15)

# Find common samples
common.samples <- intersect(
  substr(getResults(query_mirna, cols = "cases"), 1, 15),
  substr(getResults(query_cnv, cols = "cases"), 1, 15)
)
common.samples <- intersect(common.samples, 
                            substr(getResults(query_mrna, cols = "cases"), 1, 15))
common.samples <- intersect(common.samples, 
                            substr(getResults(query_methyl, cols = "cases"), 1, 15))
common.samples <- intersect(common.samples,
                            substr(getResults(query_protein, cols = "cases"), 1, 15))

cat("Final common samples:", length(common.samples), "\n") #624
# -----------------------------------------------
# Step 3: Requery with common samples
# -----------------------------------------------

# mRNA
query_mrna <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  barcode = common.samples) #parameter enforced by GDCquery

# miRNA
query_mirna <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Isoform Expression Quantification",  
  barcode = common.samples,
  workflow.type = "BCGSC miRNA Profiling")

# Methylation
query_methyl <- GDCquery(
  project = project,
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450",  
  barcode = common.samples) #parameter enforced by GDCquery

# CNV
query_cnv <- GDCquery(
  project = project,
  data.category = "Copy Number Variation",
  data.type =  "Masked Copy Number Segment",
  barcode= common.samples
)

# Protein Data
query_protein <- GDCquery(
  project = project,
  data.category = "Proteome Profiling",
  data.type = "Protein Expression Quantification",
  experimental.strategy = "Reverse Phase Protein Array", 
  barcode = common.samples)  



# -----------------------------------------------
# Step 4: Download and prepare common samples
# -----------------------------------------------

GDCdownload(query = query_mrna, method = "api")
rna <- GDCprepare(query = query_mrna, save = TRUE, save.filename = "../data/processed/RNA.rda")

GDCdownload(query = query_mirna, method = "api")
mirna <- GDCprepare(query = query_mirna, save = TRUE, save.filename = "../data/processed/miRNA.rda")

GDCdownload(query = query_methyl, method = "api")
met <- GDCprepare(query = query_methyl, save = TRUE, save.filename = "../data/processed/met.rda")

GDCdownload(query = query_protein, method = "api")
prot <- GDCprepare(query = query_protein, save = TRUE, save.filename = "../data/processed/prot.rda")

GDCdownload(query = query_cnv, method = "api", files.per.chunk = 10)
cnv <- GDCprepare(query = query_cnv, save = TRUE, save.filename = "../data/processed/cnv.rda")