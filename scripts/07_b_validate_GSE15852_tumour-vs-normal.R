library(GEOquery) 
library(limma)
library(dplyr)    
library(readr)
if (!requireNamespace("pROC", quietly=TRUE)) install.packages("pROC")
library(pROC)

DATA_DIR   <- "D:/Documents Folder/model_data"
OUTPUT_DIR <- file.path(DATA_DIR,"cv_results_corrected","geo_validation")
PANEL_PATH <- file.path(DATA_DIR,"cv_results_corrected",
                         "tiered_panel","full_panel_all_tiers.csv")
dir.create(OUTPUT_DIR, recursive=TRUE, showWarnings=FALSE)

options(timeout=300)
panel       <- read_csv(PANEL_PATH, show_col_types=FALSE)
panel_genes <- panel$Gene
cat(sprintf("Panel: %d genes\n", length(panel_genes)))

# Probe mapping
map_probes <- function(expr_mat, fdata) {
    sym_col <- NULL
    for (col in c("Gene Symbol","GENE_SYMBOL","gene_assignment",
                  "Gene.Symbol","Symbol","SYMBOL")) {
        if (col %in% names(fdata)) { sym_col <- col; break }
    }
    if (is.null(sym_col)) {
        cat("Symbol column not found. Available:\n")
        print(names(fdata)[1:10])
        return(NULL)
    }
    fdata$sym <- trimws(gsub("///.*","",gsub(" //.*","",
                              as.character(fdata[[sym_col]]))))
    valid     <- !is.na(fdata$sym) & nchar(fdata$sym)>0 &
                 !fdata$sym %in% c("---","NA","")
    e2        <- expr_mat[valid,,drop=FALSE]
    f2        <- fdata[valid,]
    f2$iqr    <- apply(e2, 1, IQR, na.rm=TRUE)
    f2$ridx   <- seq_len(nrow(f2))
    best      <- f2 %>%
        group_by(sym) %>%
        slice_max(iqr, n=1, with_ties=FALSE) %>% ungroup()
    eg        <- e2[best$ridx,,drop=FALSE]
    rownames(eg) <- best$sym
    cat(sprintf("  Mapped: %d unique genes\n", nrow(eg)))
    return(eg)
}

# Download
cat("Downloading GSE15852...\n")
gse   <- getGEO("GSE15852", GSEMatrix=TRUE,
                 AnnotGPL=FALSE, getGPL=FALSE)
if (is.list(gse)) gse <- gse[[1]]
expr  <- exprs(gse)
pdata <- pData(gse)
fdata <- fData(gse)
cat(sprintf("Loaded: %d probes x %d samples\n",
            nrow(expr), ncol(expr)))

# Identify tumour vs normal 
is_tumour <- NULL
for (col in names(pdata)) {
    vals <- tolower(as.character(pdata[[col]]))
    if (any(grepl("normal|reduction",vals)) &&
        any(grepl("tumor|cancer|invasive|carcinoma",vals))) {
        is_tumour <- grepl("tumor|cancer|invasive|carcinoma",
                            vals, ignore.case=TRUE)
        cat(sprintf("Tissue col '%s': T=%d N=%d\n",
                    col, sum(is_tumour), sum(!is_tumour)))
        break
    }
}
if (is.null(is_tumour)) {
    titles    <- tolower(as.character(pdata$title))
    is_tumour <- !grepl("normal|reduction|adjacent", titles)
    cat(sprintf("Using title: T=%d N=%d\n",
                sum(is_tumour), sum(!is_tumour)))
}

# Map probes 
cat("Mapping probes to genes...\n")
expr_genes <- map_probes(expr, fdata)
if (is.null(expr_genes)) {
    cat("Retrying with AnnotGPL=TRUE...\n")
    gse2       <- getGEO("GSE15852", GSEMatrix=TRUE,
                          AnnotGPL=TRUE, getGPL=TRUE)
    if (is.list(gse2)) gse2 <- gse2[[1]]
    expr_genes <- map_probes(exprs(gse2), fData(gse2))
}
stopifnot(!is.null(expr_genes))

# Analysis 
avail  <- intersect(panel_genes, rownames(expr_genes))
cat(sprintf("Panel genes found: %d / %d\n",
            length(avail), length(panel_genes)))

ep     <- expr_genes[avail,,drop=FALSE]
grp    <- factor(as.integer(is_tumour), levels=c(0,1),
                 labels=c("Normal","Tumour"))
design <- model.matrix(~grp)
colnames(design) <- c("Intercept","Tumour_vs_Normal")
fit    <- lmFit(ep, design)
fit2   <- eBayes(fit)
res    <- topTable(fit2, coef="Tumour_vs_Normal",
                   number=Inf, adjust.method="BH", sort.by="P")
res$Gene    <- rownames(res)
res$Dataset <- "GSE15852"
res <- res %>%
    rename(log2FC=logFC, p_value=P.Value, FDR_q=adj.P.Val) %>%
    select(Gene, Dataset, log2FC, p_value, FDR_q) %>%
    mutate(sig=FDR_q<0.05,
           Direction=ifelse(log2FC>0,"Higher in Tumour",
                                      "Higher in Normal"))

# AUC
res$AUC <- sapply(res$Gene, function(g) {
    g1 <- as.numeric(expr_genes[g,  is_tumour])
    g2 <- as.numeric(expr_genes[g, !is_tumour])
    y  <- c(rep(1,length(g1)), rep(0,length(g2)))
    x  <- c(g1, g2)
    v  <- !is.na(x)
    if (sum(v)<10 || length(unique(y[v]))<2) return(NA_real_)
    tryCatch({
        r <- roc(y[v], x[v], quiet=TRUE)
        a <- as.numeric(auc(r)); max(a,1-a)
    }, error=function(e) NA_real_)
})

out <- file.path(OUTPUT_DIR, "GSE15852_tvn_results.csv")
write_csv(res, out)
cat(sprintf("Saved: %s\n", out))
cat(sprintf("Sig (q<0.05): %d | AUC>=0.65: %d\n",
            sum(res$sig), sum(!is.na(res$AUC) & res$AUC>=0.65)))
cat("GSE15852 COMPLETE\n")
