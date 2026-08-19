library(GEOquery) 
library(limma)
library(dplyr)    
library(readr)
if (!requireNamespace("pROC", quietly=TRUE)) install.packages("pROC")
library(pROC)
options(timeout=300)

DATA_DIR   <- "../model_data"
OUTPUT_DIR <- file.path(DATA_DIR,"cv_results_corrected","geo_validation")
PANEL_PATH <- file.path(DATA_DIR,"cv_results_corrected",
                         "tiered_panel","full_panel_all_tiers.csv")
dir.create(OUTPUT_DIR, recursive=TRUE, showWarnings=FALSE)

panel       <- read_csv(PANEL_PATH, show_col_types=FALSE)
panel_genes <- panel$Gene
cat(sprintf("Panel: %d genes\n", length(panel_genes)))

# GPL-based probe mapping (works for any Affymetrix platform)
get_probe_map <- function(gpl_id) {
    cat(sprintf("Loading %s annotation...\n", gpl_id))
    gpl     <- getGEO(gpl_id)
    gpl_tbl <- Table(gpl)
    # Find gene symbol column
    sym_col <- names(gpl_tbl)[grepl("^gene.symbol$|^symbol$",
                                     names(gpl_tbl), ignore.case=TRUE)][1]
    if (is.na(sym_col)) {
        sym_col <- names(gpl_tbl)[grepl("symbol",
                                         names(gpl_tbl),
                                         ignore.case=TRUE)][1]
    }
    cat(sprintf("  Gene symbol column: '%s'\n", sym_col))
    probe_map <- gpl_tbl %>%
        select(ID, all_of(sym_col)) %>%
        rename(probe_id=ID, gene_sym=all_of(sym_col)) %>%
        mutate(gene_sym=trimws(gsub(" ///.*","",
                                     gsub("///.*","", gene_sym)))) %>%
        filter(!is.na(gene_sym) & nchar(gene_sym)>0 &
               !gene_sym %in% c("---","NA",""))
    cat(sprintf("  Probes with symbols: %d\n", nrow(probe_map)))
    return(probe_map)
}

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
cat("Loading GSE109169...\n")
gse   <- getGEO("GSE109169", GSEMatrix=TRUE,
                 AnnotGPL=FALSE, getGPL=FALSE)
if (is.list(gse)) gse <- gse[[1]]
expr  <- exprs(gse); pdata <- pData(gse)
cat(sprintf("Loaded: %d probes x %d samples\n", nrow(expr), ncol(expr)))

# Identify platform
gpl_id <- annotation(gse)
cat(sprintf("Platform: %s\n", gpl_id))
probe_map <- get_probe_map(gpl_id)

# Tissue type
is_tumour <- NULL
for (col in names(pdata)) {
    vals <- tolower(as.character(pdata[[col]]))
    if (any(grepl("normal|adjacent",vals)) &&
        any(grepl("tumor|cancer|carcinoma",vals))) {
        is_tumour <- grepl("tumor|cancer|carcinoma",
                            vals, ignore.case=TRUE)
        cat(sprintf("Tissue col '%s': T=%d N=%d\n",
                    col, sum(is_tumour), sum(!is_tumour))); break
    }
}
if (is.null(is_tumour)) {
    is_tumour <- !grepl("normal|adjacent|non.tumor",
                         tolower(as.character(pdata$title)))
    cat(sprintf("Using title: T=%d N=%d\n",
                sum(is_tumour), sum(!is_tumour)))
}

# Map + analyse
expr_genes <- map_probes_gpl(expr, probe_map, panel_genes)
stopifnot(!is.null(expr_genes))

avail  <- intersect(panel_genes, rownames(expr_genes))
ep     <- expr_genes[avail,, drop=FALSE]
grp    <- factor(as.integer(is_tumour), levels=c(0,1),
                 labels=c("Normal","Tumour"))
design <- model.matrix(~grp)
colnames(design) <- c("Intercept","Tumour_vs_Normal")
fit  <- lmFit(ep, design); fit2 <- eBayes(fit)
res  <- topTable(fit2, coef="Tumour_vs_Normal",
                  number=Inf, adjust.method="BH", sort.by="P")
res$Gene <- rownames(res); res$Dataset <- "GSE109169"
res <- res %>%
    rename(log2FC=logFC, p_value=P.Value, FDR_q=adj.P.Val) %>%
    select(Gene, Dataset, log2FC, p_value, FDR_q) %>%
    mutate(sig=FDR_q<0.05,
           Direction=ifelse(log2FC>0,"Higher in Tumour","Higher in Normal"))
res$AUC <- sapply(res$Gene, function(g) {
    g1<-as.numeric(expr_genes[g,is_tumour])
    g2<-as.numeric(expr_genes[g,!is_tumour])
    y<-c(rep(1,length(g1)),rep(0,length(g2))); x<-c(g1,g2); v<-!is.na(x)
    if (sum(v)<10||length(unique(y[v]))<2) return(NA_real_)
    tryCatch({r<-roc(y[v],x[v],quiet=TRUE);a<-as.numeric(auc(r));max(a,1-a)},
             error=function(e) NA_real_)
})
out <- file.path(OUTPUT_DIR,"GSE109169_tvn_results.csv")
write_csv(res, out)
cat(sprintf("Saved: %s\nSig: %d | AUC>=0.65: %d\nGSE109169 COMPLETE\n",
            out, sum(res$sig),
            sum(!is.na(res$AUC)&res$AUC>=0.65)))
