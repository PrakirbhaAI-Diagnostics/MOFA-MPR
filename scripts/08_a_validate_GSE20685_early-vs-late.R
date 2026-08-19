library(GEOquery)
library(dplyr)
library(readr)
options(timeout=300)

DATA_DIR   <- "../model_data"
OUTPUT_DIR <- file.path(DATA_DIR,"cv_results_corrected","geo_validation")
PANEL_PATH <- file.path(DATA_DIR,"cv_results_corrected",
                         "tiered_panel","full_panel_all_tiers.csv")
dir.create(OUTPUT_DIR, recursive=TRUE, showWarnings=FALSE)

panel       <- read_csv(PANEL_PATH, show_col_types=FALSE)
panel_genes <- panel$Gene
cat(sprintf("Panel: %d genes\n", length(panel_genes)))

# GPL570 probe mapping
cat("Loading GPL570 annotation (may take 2-3 min — large file)...\n")
gpl570     <- getGEO("GPL570")
gpl570_tbl <- Table(gpl570)
cat(sprintf("GPL570: %d probes\n", nrow(gpl570_tbl)))

# Find gene symbol column in GPL570
sym_col <- names(gpl570_tbl)[grepl("^gene.symbol$|^symbol$",
                                    names(gpl570_tbl),
                                    ignore.case=TRUE)][1]
if (is.na(sym_col))
    sym_col <- names(gpl570_tbl)[grepl("symbol",
                                        names(gpl570_tbl),
                                        ignore.case=TRUE)][1]
cat(sprintf("Gene symbol column: '%s'\n", sym_col))

probe_map <- gpl570_tbl %>%
    select(ID, all_of(sym_col)) %>%
    rename(probe_id=ID, gene_sym=all_of(sym_col)) %>%
    mutate(gene_sym=trimws(gsub(" ///.*","",
                                 gsub("///.*","", gene_sym)))) %>%
    filter(!is.na(gene_sym) & nchar(gene_sym)>0 &
           !gene_sym %in% c("---","NA",""))
cat(sprintf("Probes with symbols: %d\n", nrow(probe_map)))

map_probes_gpl <- function(expr_mat, probe_map, panel_genes) {
    matched  <- probe_map[probe_map$probe_id %in% rownames(expr_mat),]
    cat(sprintf("  Probes matched: %d / %d\n",
                nrow(matched), nrow(expr_mat)))
    if (nrow(matched)==0) { cat("  No match!\n"); return(NULL) }
    expr_sub <- expr_mat[matched$probe_id,, drop=FALSE]
    rownames(expr_sub) <- matched$gene_sym
    iqr_v    <- apply(expr_sub, 1, IQR, na.rm=TRUE)
    sym_v    <- rownames(expr_sub)
    best_idx <- tapply(seq_along(sym_v), sym_v,
                       function(i) i[which.max(iqr_v[i])])
    eg <- expr_sub[unlist(best_idx),, drop=FALSE]
    rownames(eg) <- names(best_idx)
    avail <- intersect(panel_genes, rownames(eg))
    cat(sprintf("  Unique genes: %d | Panel found: %d / %d\n",
                nrow(eg), length(avail), length(panel_genes)))
    return(eg)
}

# Download
cat("Loading GSE20685...\n")
gse   <- getGEO("GSE20685", GSEMatrix=TRUE,
                 AnnotGPL=FALSE, getGPL=FALSE)
if (is.list(gse)) gse <- gse[[1]]
expr  <- exprs(gse); pdata <- pData(gse)
cat(sprintf("Loaded: %d probes x %d samples\n", nrow(expr), ncol(expr)))


cat("\nAll pData columns:\n")
print(names(pdata))
cat("\nFirst sample characteristics:\n")
print(pdata[1,])


T_col <- "characteristics_ch1"    # ← update if needed
N_col <- "characteristics_ch1.1"  # ← update if needed
M_col <- "characteristics_ch1.2"  # ← update if needed


T_raw <- as.character(pdata[[T_col]])
N_raw <- as.character(pdata[[N_col]])
M_raw <- as.character(pdata[[M_col]])

cat(sprintf("\nT values (first 5): %s\n",
            paste(head(T_raw,5), collapse=" | ")))
cat(sprintf("N values (first 5): %s\n",
            paste(head(N_raw,5), collapse=" | ")))
cat(sprintf("M values (first 5): %s\n",
            paste(head(M_raw,5), collapse=" | ")))

# Extract T/N/M
T_stage <- regmatches(T_raw,
            regexpr("T[0-4][a-c]?", T_raw, ignore.case=TRUE))
N_stage <- regmatches(N_raw,
            regexpr("N[0-3][a-c]?(mi)?", N_raw, ignore.case=TRUE))
M_stage <- regmatches(M_raw,
            regexpr("M[01]", M_raw, ignore.case=TRUE))

# Pad to full length (regmatches drops non-matching)
T_v <- rep(NA_character_, nrow(pdata))
N_v <- rep(NA_character_, nrow(pdata))
M_v <- rep(NA_character_, nrow(pdata))
T_v[grepl("T[0-4]", T_raw, ignore.case=TRUE)] <- T_stage
N_v[grepl("N[0-3]", N_raw, ignore.case=TRUE)] <- N_stage
M_v[grepl("M[01]",  M_raw, ignore.case=TRUE)] <- M_stage

T_u <- toupper(T_v); N_u <- toupper(N_v); M_u <- toupper(M_v)

stage_num <- case_when(
    grepl("M1",  M_u) ~ 4L,
    grepl("^T1", T_u) & grepl("^N0", N_u) ~ 1L,
    grepl("^T2", T_u) & grepl("^N0", N_u) ~ 2L,
    grepl("^T[12]",T_u) & grepl("^N1",N_u) ~ 2L,
    grepl("^T3", T_u) & grepl("^N0", N_u) ~ 2L,
    grepl("^T[123]",T_u) & grepl("^N2",N_u) ~ 3L,
    grepl("^T4", T_u) ~ 3L,
    grepl("^N3", N_u) ~ 3L,
    TRUE ~ NA_integer_
)

cat("\nStage distribution after TNM conversion:\n")
print(table(stage_num, useNA="ifany"))

if (all(is.na(stage_num))) {
    cat("\nTNM extraction failed.\n")
    cat("Update T_col, N_col, M_col above and re-run.\n")
    stop("Fix column names then re-run")
}

# Map probes
expr_genes <- map_probes_gpl(expr, probe_map, panel_genes)
stopifnot(!is.null(expr_genes))

avail     <- intersect(panel_genes, rownames(expr_genes))
early_idx <- stage_num %in% c(1,2) & !is.na(stage_num)
late_idx  <- stage_num %in% c(3,4) & !is.na(stage_num)
valid_s   <- !is.na(stage_num)
cat(sprintf("Panel: %d | Early(I+II): %d | Late(III+IV): %d\n",
            length(avail), sum(early_idx), sum(late_idx)))

results <- lapply(avail, function(gene) {
    eg  <- as.numeric(expr_genes[gene,])
    kw  <- tryCatch(
        kruskal.test(eg[valid_s]~factor(stage_num[valid_s]))$p.value,
        error=function(e) NA_real_)
    sp  <- tryCatch(
        cor.test(stage_num[valid_s], eg[valid_s],
                 method="spearman", exact=FALSE),
        error=function(e) NULL)
    mwu <- tryCatch(
        wilcox.test(eg[early_idx], eg[late_idx],
                    exact=FALSE)$p.value,
        error=function(e) NA_real_)
    med <- tapply(eg, stage_num, median, na.rm=TRUE)
    data.frame(
        Gene=gene, Dataset="GSE20685",
        KW_p=kw,
        Spearman_rho=if(!is.null(sp)) sp$estimate else NA_real_,
        Spearman_p  =if(!is.null(sp)) sp$p.value  else NA_real_,
        MWU_EL_p=mwu,
        Median_S1=med["1"], Median_S2=med["2"],
        Median_S3=med["3"], Median_S4=med["4"],
        n_early=sum(early_idx), n_late=sum(late_idx),
        stringsAsFactors=FALSE)
})

res            <- bind_rows(results)
res$KW_q       <- p.adjust(res$KW_p,       method="BH")
res$Spearman_q <- p.adjust(res$Spearman_p,  method="BH")
res$MWU_EL_q   <- p.adjust(res$MWU_EL_p,   method="BH")
res$Stage_dir  <- ifelse(!is.na(res$Spearman_rho) & res$Spearman_rho>0,
                          "Increases with stage", "Decreases with stage")

out <- file.path(OUTPUT_DIR,"GSE20685_stage_results.csv")
write_csv(res, out)
cat(sprintf("Saved: %s\nKW sig: %d | Spearman sig: %d | MWU sig: %d\nGSE20685 COMPLETE\n",
            out,
            sum(res$KW_q<0.05, na.rm=TRUE),
            sum(res$Spearman_q<0.05, na.rm=TRUE),
            sum(res$MWU_EL_q<0.05, na.rm=TRUE)))
