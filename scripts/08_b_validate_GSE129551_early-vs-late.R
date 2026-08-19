library(GEOquery); library(dplyr); library(readr)
options(timeout=300)

DATA_DIR   <- "../model_data"
OUTPUT_DIR <- file.path(DATA_DIR,"cv_results_corrected","geo_validation")
PANEL_PATH <- file.path(DATA_DIR,"cv_results_corrected",
                         "tiered_panel","full_panel_all_tiers.csv")
dir.create(OUTPUT_DIR, recursive=TRUE, showWarnings=FALSE)

panel       <- read_csv(PANEL_PATH, show_col_types=FALSE)
panel_genes <- panel$Gene
cat(sprintf("Panel: %d genes\n", length(panel_genes)))

# GPL96 probe mapping
cat("Loading GPL96 annotation...\n")
gpl96     <- getGEO("GPL96")
gpl96_tbl <- Table(gpl96)
probe_map <- gpl96_tbl %>%
    select(ID, `Gene Symbol`) %>%
    rename(probe_id=ID, gene_sym=`Gene Symbol`) %>%
    mutate(gene_sym=trimws(gsub(" ///.*","",
                                 gsub("///.*","", gene_sym)))) %>%
    filter(!is.na(gene_sym) & nchar(gene_sym)>0 &
           !gene_sym %in% c("---","NA",""))
cat(sprintf("GPL96: %d probes\n", nrow(probe_map)))

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
cat("Loading GSE129551...\n")
gse   <- getGEO("GSE129551", GSEMatrix=TRUE,
                 AnnotGPL=FALSE, getGPL=FALSE)
if (is.list(gse)) gse <- gse[[1]]
expr  <- exprs(gse); pdata <- pData(gse)
cat(sprintf("Loaded: %d probes x %d samples\n", nrow(expr), ncol(expr)))

# Print characteristics to find stage
cat("\nCharacteristics (first sample):\n")
char_cols <- names(pdata)[grepl("charact|stage|title",
                                 names(pdata), ignore.case=TRUE)]
print(pdata[1, char_cols[1:min(6, length(char_cols))]])

# Find stage column
stage_num <- NULL
for (col in names(pdata)) {
    vals <- tolower(as.character(pdata[[col]]))
    if (any(grepl("stage", vals))) {
        cat(sprintf("\nStage col '%s': %s\n", col,
                    paste(unique(vals)[1:min(6,length(unique(vals)))],
                          collapse=" | ")))
        sn <- case_when(
            grepl("stage i[^iv]|stage 1$|\\bi$",  vals) ~ 1L,
            grepl("stage ii[^i]|stage 2$|\\bii$", vals) ~ 2L,
            grepl("stage iii|stage 3$|\\biii$",   vals) ~ 3L,
            grepl("stage iv|stage 4$|\\biv$",     vals) ~ 4L,
            TRUE ~ NA_integer_
        )
        if (any(!is.na(sn))) { stage_num <- sn; break }
    }
}

if (is.null(stage_num) || all(is.na(stage_num))) {
    cat("\nAuto-detection failed. Full pData row 1:\n")
    print(pdata[1,])
    stop("Update stage parsing above based on pData output")
}

cat("Stage distribution:\n")
print(table(stage_num, useNA="ifany"))

# Map probes
expr_genes <- map_probes_gpl(expr, probe_map, panel_genes)
stopifnot(!is.null(expr_genes))

avail     <- intersect(panel_genes, rownames(expr_genes))
early_idx <- stage_num %in% c(1,2) & !is.na(stage_num)
late_idx  <- stage_num %in% c(3,4) & !is.na(stage_num)
valid_s   <- !is.na(stage_num)
cat(sprintf("Panel: %d | Early: %d | Late: %d\n",
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
    data.frame(
        Gene=gene, Dataset="GSE129551",
        KW_p=kw,
        Spearman_rho=if(!is.null(sp)) sp$estimate else NA_real_,
        Spearman_p  =if(!is.null(sp)) sp$p.value  else NA_real_,
        MWU_EL_p=mwu,
        n_early=sum(early_idx), n_late=sum(late_idx),
        stringsAsFactors=FALSE)
})

res            <- bind_rows(results)
res$KW_q       <- p.adjust(res$KW_p,      method="BH")
res$Spearman_q <- p.adjust(res$Spearman_p, method="BH")
res$MWU_EL_q   <- p.adjust(res$MWU_EL_p,  method="BH")
res$Stage_dir  <- ifelse(!is.na(res$Spearman_rho) & res$Spearman_rho>0,
                          "Increases with stage", "Decreases with stage")

out <- file.path(OUTPUT_DIR,"GSE129551_stage_results.csv")
write_csv(res, out)
cat(sprintf("Saved: %s\nKW sig: %d | Spearman sig: %d | MWU sig: %d\nGSE129551 COMPLETE\n",
            out,
            sum(res$KW_q<0.05, na.rm=TRUE),
            sum(res$Spearman_q<0.05, na.rm=TRUE),
            sum(res$MWU_EL_q<0.05, na.rm=TRUE)))
