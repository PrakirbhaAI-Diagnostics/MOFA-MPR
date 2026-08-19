suppressPackageStartupMessages(library(data.table))

base_dir   <- "../data"
expr_src   <- file.path(base_dir, "GEO/GSE96058",
                "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv")
# Works for either the .csv or the .csv.gz - the connection is chosen below.
panel_path <- file.path(base_dir, "cv_results_corrected/tiered_panel/full_panel_all_tiers.csv")
out_path   <- file.path(base_dir, "GSE96058_panel370_expression.csv")


CHUNK <- 2000   # lines per chunk


 #destfile = expr_src, mode = "wb")  

if (!file.exists(expr_src)) stop("Expression file not found:\n  ", expr_src, call. = FALSE)
if (!file.exists(panel_path)) stop("Panel file not found:\n  ", panel_path, call. = FALSE)

is_gz <- grepl("\\.gz$", expr_src, ignore.case = TRUE)
cat("Source:", expr_src, "\n")
cat("Size  :", round(file.size(expr_src) / 1e6, 1), "MB",
    if (is_gz) "(gzipped)" else "(uncompressed)", "\n")

# PANEL LIST

panel_raw <- fread(panel_path)
stopifnot("Gene" %in% names(panel_raw))
panel <- toupper(trimws(as.character(panel_raw$Gene)))
panel <- panel[panel != "" & !is.na(panel)]

cat("\n--- Panel check ---\n")
cat("Rows in panel file :", nrow(panel_raw), "\n")
cat("Unique genes       :", length(unique(panel)), "\n")
if (anyDuplicated(panel))
  warning("Duplicate gene symbols in the panel file: ",
          paste(unique(panel[duplicated(panel)]), collapse = ", "))
panel <- unique(panel)
if (length(panel) != EXPECTED_PANEL)
  stop("Panel size ", length(panel), " != expected ", EXPECTED_PANEL,
       ". Resolve against MOFA_EML_DECISIONS_REGISTER.docx before proceeding.")
cat("Matches expected panel size of", EXPECTED_PANEL, ": PASS\n")

if ("Tier_num" %in% names(panel_raw)) {
  cat("Tier composition   : ")
  cat(paste(names(table(panel_raw$Tier_num)), table(panel_raw$Tier_num),
            sep = "=", collapse = ", "), "\n")
  cat("  (retained for annotation only - ALL 370 genes are extracted)\n")
}
cat("-------------------\n\n")


con <- if (is_gz) gzfile(expr_src, "r") else file(expr_src, "r")
on.exit(close(con), add = TRUE)

header <- readLines(con, n = 1)
n_cols <- length(strsplit(header, ",", fixed = TRUE)[[1]])
cat("Header columns:", n_cols, "\n")

kept       <- character(0)
kept_genes <- character(0)
n_read     <- 0
repeat {
  lines <- readLines(con, n = CHUNK)
  if (length(lines) == 0) break
  n_read <- n_read + length(lines)

  # first field of each line = gene symbol; strip optional surrounding quotes
  g <- sub(",.*$", "", lines)
  g <- toupper(trimws(gsub('^"|"$', "", g)))

  hit <- g %in% panel
  if (any(hit)) {
    kept       <- c(kept, lines[hit])
    kept_genes <- c(kept_genes, g[hit])
  }
  if (n_read %% 10000 == 0)
    cat("  scanned", n_read, "genes | matched", length(kept), "\r")
}
cat("\nScanned", n_read, "gene rows total.\n")


writeLines(c(header, kept), out_path)
cat("\nWrote:", out_path, "\n")
cat("Size :", round(file.size(out_path) / 1e6, 2), "MB\n\n")
cat("Next: set expr_path in gse96058_stage_I_vs_II.R to:\n  ", out_path, "\n")
