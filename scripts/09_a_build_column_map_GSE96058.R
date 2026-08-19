suppressPackageStartupMessages(library(data.table))

base_dir  <- "../data"
meta_path <- file.path(base_dir, "GSE96058 - Sheet4.csv")
expr_path <- file.path(base_dir, "GEO/GSE96058",
              "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv")
map_out   <- file.path(base_dir, "GSE96058_column_map.csv")


sm_urls <- c(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/matrix/GSE96058-GPL11154_series_matrix.txt.gz",
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/matrix/GSE96058-GPL18573_series_matrix.txt.gz"
)

read_header <- function(u, n = 100) {
  con <- tryCatch(gzcon(url(u, "rb")), error = function(e) NULL)
  if (is.null(con)) return(NULL)
  on.exit(try(close(con), silent = TRUE))
  tryCatch(readLines(con, n = n, warn = FALSE), error = function(e) NULL)
}

parse_field <- function(lines, key) {
  ln <- grep(paste0("^!", key), lines, value = TRUE)[1]
  if (is.na(ln)) return(character(0))
  v <- strsplit(ln, "\t")[[1]][-1]
  gsub('^"|"$', "", v)
}

maps <- list()
for (u in sm_urls) {
  cat("Reading:", basename(u), "\n")
  h <- read_header(u)
  if (is.null(h)) { cat("  -> could not read (skipping)\n"); next }
  ttl <- parse_field(h, "Sample_title")
  gsm <- parse_field(h, "Sample_geo_accession")
  if (!length(ttl) || length(ttl) != length(gsm)) {
    cat("  -> title/accession mismatch (skipping)\n"); next
  }
  cat("  ->", length(gsm), "samples | example title:", ttl[1], "| GSM:", gsm[1], "\n")
  maps[[u]] <- data.table(title = ttl, GSM = gsm)
}

if (!length(maps))
  stop("Could not read any series matrix header.\n",
       "Fallback: download the MINiML or SOFT family file from\n",
       "  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058\n",
       "and extract !Sample_title / !Sample_geo_accession manually.",
       call. = FALSE)

cmap <- unique(rbindlist(maps))
cat("\nTotal GSM<->title pairs:", nrow(cmap), "\n")
cat("Titles look like F-numbers?",
    all(grepl("^F[0-9]+$", head(cmap$title, 20))), "\n")

# JOIN TO THE METADATA VIA GSM 
meta <- fread(meta_path)
setnames(meta, "S ID", "SID")
cat("\nMetadata rows:", nrow(meta), "| unique GSM:", uniqueN(meta$GSM), "\n")

m <- merge(meta, cmap, by = "GSM", all.x = TRUE)
cat("Metadata rows with a title:", sum(!is.na(m$title)), "of", nrow(m), "\n")
cat("GEO samples absent from the sheet:", nrow(cmap) - sum(cmap$GSM %in% meta$GSM), "\n")
if (any(is.na(m$title)))
  warning(sum(is.na(m$title)), " metadata rows have no matching GEO title. ",
          "These will be dropped from the analysis.")

fwrite(m[, .(GSM, SID, SAMPLE, title, TStage, NStage, MStage,
             pam50_subtype, nhg, age_at_diagnosis, tumor_size)], map_out)
cat("Wrote column map:", map_out, "\n")


markers <- c("ESR1", "ERBB2", "MKI67")
con <- file(expr_path, "r")
hdr <- gsub('^"|"$', "", strsplit(readLines(con, 1), ",", fixed = TRUE)[[1]])
rows <- list()
repeat {
  ln <- readLines(con, n = 5000)
  if (!length(ln)) break
  g <- toupper(gsub('^"|"$', "", sub(",.*$", "", ln)))
  hit <- g %in% markers
  if (any(hit)) for (i in which(hit)) rows[[g[i]]] <- ln[i]
  if (length(rows) == length(markers)) break
}
close(con)
cat("\nMarker genes found:", paste(names(rows), collapse = ", "), "\n")

if (length(rows)) {
  # header has an empty first field for the gene column
  samp_cols <- hdr[-1]
  vals <- lapply(rows, function(l) {
    v <- strsplit(l, ",", fixed = TRUE)[[1]][-1]
    setNames(as.numeric(v), samp_cols)
  })
  mm <- m[!is.na(title)]
  cat("\n--- BIOLOGICAL VERIFICATION ---\n")
  if ("ESR1" %in% names(vals)) {
    e <- vals$ESR1[mm$title]
    print(round(tapply(e, mm$pam50_subtype, median, na.rm = TRUE), 2))
    cat("  ^ ESR1: expect Basal LOWEST, LumA/LumB HIGHEST\n\n")
  }
  if ("ERBB2" %in% names(vals)) {
    e <- vals$ERBB2[mm$title]
    print(round(tapply(e, mm$pam50_subtype, median, na.rm = TRUE), 2))
    cat("  ^ ERBB2: expect Her2 HIGHEST\n\n")
  }
  if ("MKI67" %in% names(vals)) {
    e <- vals$MKI67[mm$title]
    print(round(tapply(e, mm$nhg, median, na.rm = TRUE), 2))
    cat("  ^ MKI67: expect G3 > G2 > G1\n\n")
  }
  cat("If these patterns hold, the mapping is correct.\n")
  cat("If they are flat or scrambled, STOP - the mapping is wrong.\n")
}
