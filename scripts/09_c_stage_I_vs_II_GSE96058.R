suppressPackageStartupMessages({
  library(data.table)
})


base_dir    <- "../data"

meta_path   <- file.path(base_dir, "GSE96058 - Sheet4.csv")   # CONFIRMED path
map_path    <- file.path(base_dir, "GSE96058_column_map.csv")


expr_path   <- file.path(base_dir, "GSE96058_panel370_expression.csv")


panel_path  <- file.path(base_dir, "cv_results_corrected/tiered_panel/full_panel_all_tiers.csv")

EXPECTED_PANEL <- 370

out_dir     <- file.path(base_dir, "gse96058_stage_validation")


for (p in c(meta_path, expr_path, panel_path)) {
  if (!file.exists(p)) stop("File not found:\n  ", p,
                            "\nFix the path in the CONFIG block.", call. = FALSE)
}

RECLASSIFY_T0     <- TRUE    # R2 above
EXPR_IS_LOG2      <- TRUE    # SCAN-B ships log2(FPKM+0.1)
MIN_STRATUM_N     <- 20      # min samples per stratum 
MIN_PER_GROUP     <- 5       # min Stage I and Stage II count within a stratum
RUN_EL_SUPPLEMENT <- FALSE   # Early-vs-Late
FDR_METHOD        <- "BH"
SEED              <- 1

set.seed(SEED)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
manifest <- c()
log_line <- function(...) {
  s <- paste0(...)
  message(s)
  manifest <<- c(manifest, s)
}
log_line("=== GSE96058 Stage I vs II validation ===")
log_line("Run date: ", as.character(Sys.Date()))
log_line("RECLASSIFY_T0 = ", RECLASSIFY_T0)


meta <- fread(meta_path)
setnames(meta, "S ID", "SID")
log_line("Metadata rows (GSM level): ", nrow(meta))
log_line("Unique S ID: ", uniqueN(meta$SID))

# integrity check
tnm_var <- meta[, .(nT = uniqueN(TStage), nN = uniqueN(NStage), nM = uniqueN(MStage)),
                by = SID]
stopifnot(max(tnm_var$nT, tnm_var$nN, tnm_var$nM) == 1)
log_line("TNM constant within all replicated S IDs: PASS")


cmap <- fread(map_path)
stopifnot(all(c("GSM", "title", "platform") %in% names(cmap)))
meta <- merge(meta, cmap[, .(GSM, title, platform)], by = "GSM", all.x = TRUE)
log_line("Metadata rows with an F-number title: ", sum(!is.na(meta$title)),
         " of ", nrow(meta))
if (any(is.na(meta$title)))
  warning(sum(is.na(meta$title)), " rows have no title and will be dropped. ",
          "If this is ~333, a series matrix platform failed to load - rerun 00b.")
meta[, expr_col := title]

derive_stage <- function(T, N, M, size, reclass_t0) {
  if (reclass_t0 && !is.na(T) && T == "T0" && !is.na(size)) {
    T <- if (size <= 20) "T1" else if (size <= 50) "T2" else "T3"
  }
  if (!is.na(M) && M == "M1")                       return("IV")
  if (!is.na(T) && T == "Tis")                      return("0")           # R3
  if ((!is.na(T) && T == "TX") ||
      (!is.na(N) && N == "NX"))                     return("Unknown")     # R3
  if (N %in% c("N2", "N3"))                         return("III")
  if (T == "T4")                                    return("III")
  if (T == "T3" && N == "N1")                       return("III")
  if (T == "T1")                                    return(if (N == "N0") "I" else "II")
  if (T == "T2")                                    return("II")
  if (T == "T3" && N == "N0")                       return("II")
  if (T == "T0")                                    return(if (N == "N1") "II" else "Unstageable")
  return("Unknown")
}

meta[, tumor_size := suppressWarnings(as.numeric(tumor_size))]
meta[, stage := mapply(derive_stage, TStage, NStage, MStage, tumor_size,
                       MoreArgs = list(reclass_t0 = RECLASSIFY_T0))]
meta[, stage_group := fifelse(stage %in% c("I", "II"), "Early",
                      fifelse(stage %in% c("III", "IV"), "Late", stage))]

samp <- unique(meta, by = "SID")
stage_tab <- samp[, .N, by = stage][order(-N)]
fwrite(stage_tab, file.path(out_dir, "01_stage_derivation_table.csv"))
log_line("Stage counts (unique samples):")
for (i in seq_len(nrow(stage_tab)))
  log_line("   Stage ", stage_tab$stage[i], ": ", stage_tab$N[i])

con <- file(expr_path, "r"); hdr_line <- readLines(con, n = 1); close(con)
hdr <- gsub('^"|"$', "", trimws(strsplit(hdr_line, ",", fixed = TRUE)[[1]]))
expr <- fread(expr_path, header = FALSE, skip = 1)
if (length(hdr) == ncol(expr)) {
  setnames(expr, hdr); setnames(expr, 1, "gene")
} else if (length(hdr) == ncol(expr) - 1) {
  setnames(expr, c("gene", hdr))
} else {
  stop("Header/data mismatch: ", length(hdr), " vs ", ncol(expr), call. = FALSE)
}
expr[, gene := toupper(trimws(gsub('^"|"$', "", gene)))]
expr <- expr[!duplicated(gene)]          # SCAN-B is gene-level; no probe collapse
log_line("Expression matrix: ", nrow(expr), " genes x ", ncol(expr) - 1, " columns")

panel_raw <- fread(panel_path)
stopifnot("Gene" %in% names(panel_raw))
panel <- toupper(trimws(as.character(panel_raw$Gene)))
panel <- unique(panel[panel != "" & !is.na(panel)])
log_line("Panel genes (unique): ", length(panel))
if (length(panel) != EXPECTED_PANEL)
  stop("Panel size ", length(panel), " != expected ", EXPECTED_PANEL,
       ". Resolve against the register.")
log_line("Panel size matches register value of ", EXPECTED_PANEL, ": PASS")

panel_found <- intersect(panel, expr$gene)
log_line("Panel genes present in GSE96058: ", length(panel_found),
         " of ", length(panel),
         "  (missing: ", length(setdiff(panel, expr$gene)), ")")
writeLines(setdiff(panel, expr$gene), file.path(out_dir, "panel_genes_missing.txt"))

E <- as.matrix(expr[gene %in% panel_found, -"gene"])
rownames(E) <- expr[gene %in% panel_found, gene]
if (!EXPR_IS_LOG2) E <- log2(E + 0.1)

meta_m <- meta[expr_col %in% colnames(E)]
log_line("GSM rows matched to an expression column: ", nrow(meta_m), " of ", nrow(meta))
if (nrow(meta_m) < nrow(meta) * 0.9)
  warning("Low column-match rate. Check the map has 3409 rows and that titles are F-numbers.")

E <- E[, meta_m$expr_col, drop = FALSE]
sid_vec <- as.character(meta_m$SID)
n_repl <- sum(duplicated(sid_vec))

if (n_repl > 0) {
  idx <- split(seq_along(sid_vec), sid_vec)          # named list, sorted by SID
  E_collapsed <- vapply(idx,
                        function(j) rowMeans(E[, j, drop = FALSE], na.rm = TRUE),
                        numeric(nrow(E)))
  if (nrow(E) == 1) E_collapsed <- matrix(E_collapsed, nrow = 1)
  rownames(E_collapsed) <- rownames(E)
  colnames(E_collapsed) <- names(idx)
} else {
  E_collapsed <- E
  colnames(E_collapsed) <- sid_vec
}
E_collapsed <- as.matrix(E_collapsed)
storage.mode(E_collapsed) <- "double"
# strip any stray attributes that would break identical()-style comparisons
colnames(E_collapsed) <- as.character(unname(colnames(E_collapsed)))
rownames(E_collapsed) <- as.character(unname(rownames(E_collapsed)))

log_line("Replicate columns averaged: ", n_repl,
         " -> collapsed matrix ", nrow(E_collapsed), " x ", ncol(E_collapsed))

plat_by_sid <- meta_m[, .(platform = if (uniqueN(platform) > 1) "Mixed"
                                     else as.character(platform[1])), by = SID]
log_line("Platform per sample: ",
         paste(names(table(plat_by_sid$platform)), table(plat_by_sid$platform),
               sep = "=", collapse = ", "))

pheno <- unique(meta_m, by = "SID")[, .(SID = as.character(SID), stage,
                                        stage_group, pam50 = pam50_subtype, nhg,
                                        age = age_at_diagnosis,
                                        size = tumor_size)]
pheno <- merge(pheno, plat_by_sid, by = "SID", all.x = TRUE)

cn <- as.character(colnames(E_collapsed))
if (is.null(cn) || !length(cn))
  stop("E_collapsed has no column names. The replicate collapse failed upstream.",
       call. = FALSE)

only_in_E     <- setdiff(cn, pheno$SID)
only_in_pheno <- setdiff(pheno$SID, cn)
if (length(only_in_E) || length(only_in_pheno)) {
  message("--- SAMPLE ALIGNMENT MISMATCH ---")
  message("  In expression but not metadata: ", length(only_in_E),
          if (length(only_in_E)) paste0("  e.g. ",
              paste(head(only_in_E, 5), collapse = ", ")) else "")
  message("  In metadata but not expression: ", length(only_in_pheno),
          if (length(only_in_pheno)) paste0("  e.g. ",
              paste(head(only_in_pheno, 5), collapse = ", ")) else "")
  message("  Proceeding on the intersection. If these counts are large, the ",
          "F-number map is wrong - check GSE96058_column_map.csv against the ",
          "matrix header before trusting anything downstream.")
}

common <- intersect(cn, pheno$SID)          # preserves expression column order
if (!length(common))
  stop("No samples shared between expression matrix and metadata.\n",
       "  Example expression colname: ", cn[1], "\n",
       "  Example metadata SID      : ", pheno$SID[1], "\n",
       "  The F-number map join has failed. Rerun 00b_build_column_map.R.",
       call. = FALSE)

E_collapsed <- E_collapsed[, common, drop = FALSE]
pheno       <- pheno[match(common, SID)]

stopifnot(!anyNA(pheno$SID),
          nrow(pheno) == ncol(E_collapsed),
          all(pheno$SID == colnames(E_collapsed)))
log_line("Aligned samples (expression + metadata): ", ncol(E_collapsed))

coh12 <- pheno[stage %in% c("I", "II")]
E12   <- E_collapsed[, coh12$SID, drop = FALSE]
coh12[, grp := factor(stage, levels = c("I", "II"))]
log_line("PRIMARY cohort  Stage I vs II: ",
         sum(coh12$stage == "I"), " vs ", sum(coh12$stage == "II"),
         "  (total ", nrow(coh12), ")")

# CONFOUND DIAGNOSTICS 

ct_pam <- table(coh12$stage, coh12$pam50)
ct_nhg <- table(coh12$stage, coh12$nhg, useNA = "no")
chi_pam <- chisq.test(ct_pam)
chi_nhg <- chisq.test(ct_nhg)
cramers_v <- function(tab) {
  chi <- suppressWarnings(chisq.test(tab)$statistic)
  sqrt(as.numeric(chi) / (sum(tab) * (min(dim(tab)) - 1)))
}
conf <- data.table(
  factor    = c("PAM50", "NHG"),
  chisq     = c(chi_pam$statistic, chi_nhg$statistic),
  df        = c(chi_pam$parameter, chi_nhg$parameter),
  p         = c(chi_pam$p.value, chi_nhg$p.value),
  cramers_v = c(cramers_v(ct_pam), cramers_v(ct_nhg))
)
fwrite(conf, file.path(out_dir, "02_confound_diagnostics.csv"))
write.csv(ct_pam, file.path(out_dir, "02b_stage_by_pam50_table.csv"))
write.csv(ct_nhg, file.path(out_dir, "02c_stage_by_nhg_table.csv"))
log_line("Stage~PAM50  Cramer's V = ", signif(conf$cramers_v[1], 3),
         "  p = ", signif(conf$p[1], 3))
log_line("Stage~NHG    Cramer's V = ", signif(conf$cramers_v[2], 3),
         "  p = ", signif(conf$p[2], 3))
log_line("  -> If V is non-trivial, unadjusted stage hits are confound-suspect.")

## --------------------------------------- 5. PRIMARY: UNADJUSTED MWU -------
# rank-biserial correlation r = 2*AUC - 1  (Cliff's delta, identical here)
rank_biserial <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  r  <- rank(c(x, y))
  U1 <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
  auc <- U1 / (n1 * n2)
  2 * auc - 1
}

mwu_one <- function(v, grp) {
  x <- v[grp == "I"]; y <- v[grp == "II"]
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < MIN_PER_GROUP || length(y) < MIN_PER_GROUP)
    return(list(p = NA_real_, rb = NA_real_, n1 = length(x), n2 = length(y),
                med_I = NA_real_, med_II = NA_real_))
  w <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
  list(p = w$p.value, rb = rank_biserial(x, y),
       n1 = length(x), n2 = length(y),
       med_I = median(x), med_II = median(y))
}

prim <- rbindlist(lapply(rownames(E12), function(g) {
  r <- mwu_one(E12[g, ], coh12$grp)
  data.table(gene = g, n_I = r$n1, n_II = r$n2,
             median_I = r$med_I, median_II = r$med_II,
             delta_median = r$med_II - r$med_I,
             rank_biserial = r$rb, p_MWU = r$p)
}))
prim[, q_MWU := p.adjust(p_MWU, method = FDR_METHOD)]
prim[, direction := fifelse(is.na(rank_biserial), NA_character_,
                     fifelse(rank_biserial > 0, "higher_in_StageI",
                                                "higher_in_StageII"))]
setorder(prim, q_MWU, -abs(rank_biserial))
fwrite(prim, file.path(out_dir, "03_stageI_vs_II_unadjusted.csv"))
log_line("PRIMARY MWU: ", sum(prim$q_MWU < 0.05, na.rm = TRUE),
         " of ", nrow(prim), " genes at q < 0.05")
log_line("   of which |rank-biserial| >= 0.147 (~small effect): ",
         sum(prim$q_MWU < 0.05 & abs(prim$rank_biserial) >= 0.147, na.rm = TRUE))
log_line("   max |rank-biserial| among q<0.05: ",
         signif(suppressWarnings(max(abs(prim[q_MWU < 0.05]$rank_biserial),
                                     na.rm = TRUE)), 3))

# STRATIFIED: VAN ELTEREN 

van_elteren <- function(v, grp, stratum) {
  keep <- is.finite(v) & !is.na(grp) & !is.na(stratum)
  v <- v[keep]; grp <- droplevels(factor(grp[keep])); stratum <- stratum[keep]
  Tstat <- 0; Vstat <- 0; used <- 0
  for (s in unique(stratum)) {
    idx <- stratum == s
    g <- grp[idx]; x <- v[idx]
    n1 <- sum(g == "I"); n2 <- sum(g == "II"); n <- n1 + n2
    if (n1 < MIN_PER_GROUP || n2 < MIN_PER_GROUP || n < MIN_STRATUM_N) next
    r  <- rank(x)
    W  <- sum(r[g == "I"])
    E  <- n1 * (n + 1) / 2
    # tie-corrected variance
    V  <- (n1 * n2 / (n * (n - 1))) * sum(r^2) -
          n1 * n2 * (n + 1)^2 / (4 * (n - 1))
    w  <- 1 / (n + 1)
    Tstat <- Tstat + w * (W - E)
    Vstat <- Vstat + w^2 * V
    used  <- used + 1
  }
  if (used == 0 || Vstat <= 0)
    return(list(z = NA_real_, p = NA_real_, k = used))
  z <- Tstat / sqrt(Vstat)
  list(z = z, p = 2 * pnorm(-abs(z)), k = used)
}

run_ve <- function(stratum_vec, label) {
  res <- rbindlist(lapply(rownames(E12), function(g) {
    r <- van_elteren(E12[g, ], coh12$grp, stratum_vec)
    data.table(gene = g, blocking = label, z_vanElteren = r$z,
               p_vanElteren = r$p, n_strata_used = r$k)
  }))
  res[, q_vanElteren := p.adjust(p_vanElteren, method = FDR_METHOD)]
  res
}

strat_pam  <- coh12$pam50
strat_nhg  <- coh12$nhg
strat_both <- paste(coh12$pam50, coh12$nhg, sep = "|")
strat_both[is.na(coh12$pam50) | is.na(coh12$nhg)] <- NA
strat_all  <- paste(coh12$pam50, coh12$nhg, coh12$platform, sep = "|")
strat_all[is.na(coh12$pam50) | is.na(coh12$nhg) | is.na(coh12$platform)] <- NA

ve <- rbindlist(list(
  run_ve(strat_pam,  "PAM50"),
  run_ve(strat_nhg,  "NHG"),
  run_ve(strat_both, "PAM50_x_NHG"),
  run_ve(strat_all,  "PAM50_x_NHG_x_platform")
))
setorder(ve, blocking, q_vanElteren)
fwrite(ve, file.path(out_dir, "04_stageI_vs_II_stratified.csv"))
for (b in unique(ve$blocking))
  log_line("van Elteren (block = ", b, "): ",
           sum(ve[blocking == b]$q_vanElteren < 0.05, na.rm = TRUE),
           " of ", nrow(ve[blocking == b]), " genes at q < 0.05")

# --- THE HEADLINE COMPARISON: what survives adjustment?
surv <- merge(prim[, .(gene, q_MWU, rank_biserial)],
              ve[blocking == "PAM50_x_NHG", .(gene, q_vanElteren)], by = "gene")
surv[, status := fifelse(q_MWU < 0.05 & q_vanElteren < 0.05, "robust_to_subtype_grade",
                  fifelse(q_MWU < 0.05 & q_vanElteren >= 0.05, "CONFOUNDED_lost_on_adjustment",
                  fifelse(q_MWU >= 0.05 & q_vanElteren < 0.05, "revealed_on_adjustment",
                                                               "null_both")))]
fwrite(surv, file.path(out_dir, "04b_adjustment_survival.csv"))
log_line("--- Adjustment survival (PAM50 x NHG) ---")
for (s in unique(surv$status)) log_line("   ", s, ": ", sum(surv$status == s))


lr <- rbindlist(lapply(rownames(E12), function(g) {
  d <- data.table(y = as.integer(coh12$grp == "II"),
                  x = as.numeric(scale(E12[g, ])),
                  pam50 = factor(coh12$pam50), nhg = factor(coh12$nhg),
                  plat = factor(coh12$platform))
  d <- d[complete.cases(d)]
  out <- tryCatch({
    m <- glm(y ~ x + pam50 + nhg + plat, data = d, family = binomial())
    cf <- summary(m)$coefficients
    data.table(gene = g, OR_perSD = exp(cf["x", 1]),
               p_adj_model = cf["x", 4], n_used = nrow(d))
  }, error = function(e) data.table(gene = g, OR_perSD = NA_real_,
                                    p_adj_model = NA_real_, n_used = nrow(d)))
  out
}))
lr[, q_adj_model := p.adjust(p_adj_model, method = FDR_METHOD)]
setorder(lr, q_adj_model)
fwrite(lr, file.path(out_dir, "05_stageI_vs_II_logistic.csv"))
log_line("Logistic (adj. PAM50+NHG): ", sum(lr$q_adj_model < 0.05, na.rm = TRUE),
         " of ", nrow(lr), " genes at q < 0.05")


RUN_LIMMA <- TRUE

if (RUN_LIMMA) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    warning("limma not installed. Install with:\n",
            '  if (!requireNamespace("BiocManager")) install.packages("BiocManager")\n',
            '  BiocManager::install("limma")\n',
            "Skipping limma arm.")
  } else {
    suppressPackageStartupMessages(library(limma))
    log_line("--- LIMMA ARM (limma-trend on log2 FPKM) ---")

    # GSM-level (uncollapsed) data for Stage I/II, so replicates can be blocked
    lm_meta <- meta_m[stage %in% c("I", "II")]
    lm_meta <- lm_meta[!is.na(pam50_subtype) & !is.na(nhg) & !is.na(platform)]
    E_lm <- E[, lm_meta$expr_col, drop = FALSE]
    log_line("limma cohort (library level, complete covariates): ", ncol(E_lm),
             " libraries from ", uniqueN(lm_meta$SID), " patients")
    log_line("   dropped for missing PAM50/NHG/platform: ",
             sum(meta_m$stage %in% c("I", "II")) - nrow(lm_meta), " libraries")

    grp   <- factor(lm_meta$stage, levels = c("I", "II"))
    pam50 <- factor(lm_meta$pam50_subtype)
    nhg   <- factor(lm_meta$nhg)
    plat  <- factor(lm_meta$platform)
    log_line("   Stage I: ", sum(grp == "I"), " | Stage II: ", sum(grp == "II"))
    log_line("   platform: ", paste(names(table(plat)), table(plat),
                                    sep = "=", collapse = ", "))

    fit_limma <- function(design, coef_name, label, block_by_patient = TRUE) {
      if (block_by_patient && anyDuplicated(lm_meta$SID)) {
        cf <- duplicateCorrelation(E_lm, design, block = lm_meta$SID)
        log_line("   [", label, "] duplicateCorrelation consensus = ",
                 signif(cf$consensus, 4))
        f <- lmFit(E_lm, design, block = lm_meta$SID, correlation = cf$consensus)
      } else {
        f <- lmFit(E_lm, design)
      }
      f <- eBayes(f, trend = TRUE)   # trend=TRUE : correct for log2 FPKM
      tt <- topTable(f, coef = coef_name, number = Inf, sort.by = "none")
      data.table(gene = rownames(tt), model = label,
                 logFC = tt$logFC, AveExpr = tt$AveExpr,
                 t = tt$t, p_limma = tt$P.Value, q_limma = tt$adj.P.Val)
    }

# (a) UNADJUSTED: stage only
    d_unadj <- model.matrix(~ grp)
    res_unadj <- fit_limma(d_unadj, "grpII", "unadjusted")

# (b) ADJUSTED: stage + PAM50 + NHG + platform 
    
    d_adj <- model.matrix(~ grp + pam50 + nhg + plat)
    res_adj <- fit_limma(d_adj, "grpII", "adjusted_PAM50_NHG")

    limma_all <- rbind(res_unadj, res_adj)
    fwrite(limma_all, file.path(out_dir, "09_limma_stageI_vs_II.csv"))
    for (m in unique(limma_all$model))
      log_line("limma (", m, "): ",
               sum(limma_all[model == m]$q_limma < 0.05, na.rm = TRUE),
               " of ", nrow(limma_all[model == m]), " genes at q < 0.05")

    cmp <- merge(prim[, .(gene, q_MWU, rank_biserial)],
                 res_unadj[, .(gene, q_limma_unadj = q_limma,
                               logFC_unadj = logFC)], by = "gene")
    cmp <- merge(cmp, ve[blocking == "PAM50_x_NHG", .(gene, q_vanElteren)],
                 by = "gene")
    cmp <- merge(cmp, res_adj[, .(gene, q_limma_adj = q_limma,
                                  logFC_adj = logFC)], by = "gene")

    cmp[, agree_unadjusted := (q_MWU < 0.05) == (q_limma_unadj < 0.05)]
    cmp[, agree_adjusted   := (q_vanElteren < 0.05) == (q_limma_adj < 0.05)]
    setorder(cmp, q_limma_adj)
    fwrite(cmp, file.path(out_dir, "10_framework_comparison.csv"))

    log_line("--- Framework agreement (same data, same genes) ---")
    log_line("   UNADJUSTED  MWU vs limma      : ",
             sum(cmp$agree_unadjusted), "/", nrow(cmp), " concordant calls")
    log_line("   ADJUSTED    vanElteren vs limma: ",
             sum(cmp$agree_adjusted), "/", nrow(cmp), " concordant calls")
    log_line("   Spearman of -log10 q (unadjusted): ",
             signif(cor(-log10(cmp$q_MWU + 1e-300),
                        -log10(cmp$q_limma_unadj + 1e-300),
                        method = "spearman"), 3))
    log_line("   Spearman of -log10 q (adjusted)  : ",
             signif(cor(-log10(cmp$q_vanElteren + 1e-300),
                        -log10(cmp$q_limma_adj + 1e-300),
                        method = "spearman"), 3))
    log_line("   -> High concordance = frameworks agree; a robustness sentence.")
    log_line("   -> Divergence = investigate BEFORE it reaches a reviewer.")

   
    surv_l <- cmp[, .(gene, q_limma_unadj, q_limma_adj, logFC_unadj, logFC_adj)]
    surv_l[, status := fifelse(q_limma_unadj < 0.05 & q_limma_adj < 0.05,
                               "robust_to_subtype_grade",
                        fifelse(q_limma_unadj < 0.05 & q_limma_adj >= 0.05,
                               "CONFOUNDED_lost_on_adjustment",
                        fifelse(q_limma_unadj >= 0.05 & q_limma_adj < 0.05,
                               "revealed_on_adjustment", "null_both")))]
    fwrite(surv_l, file.path(out_dir, "09b_limma_adjustment_survival.csv"))
    log_line("--- limma adjustment survival ---")
    for (s in unique(surv_l$status))
      log_line("   ", s, ": ", sum(surv_l$status == s))

   
    fit_chk <- eBayes(lmFit(E_lm, d_adj), trend = TRUE)
    png(file.path(out_dir, "11_limma_qq_residuals.png"), 1200, 500, res = 120)
    par(mfrow = c(1, 2))
    plotSA(fit_chk, main = "limma-trend: mean-variance")
    qqt(fit_chk$t[, "grpII"], df = fit_chk$df.total[1], main = "QQ: stage t-stats")
    abline(0, 1, col = "red")
    dev.off()
    log_line("Wrote 11_limma_qq_residuals.png - INSPECT before trusting limma.")
  }
}

# WITHIN-STRATUM MWU 

within <- rbindlist(lapply(c("PAM50", "NHG"), function(fac) {
  sv <- if (fac == "PAM50") coh12$pam50 else coh12$nhg
  rbindlist(lapply(setdiff(unique(sv), NA), function(s) {
    idx <- which(sv == s)
    if (length(idx) < MIN_STRATUM_N) return(NULL)
    sub_grp <- droplevels(coh12$grp[idx])
    if (nlevels(sub_grp) < 2) return(NULL)
    r <- rbindlist(lapply(rownames(E12), function(g) {
      o <- mwu_one(E12[g, idx], sub_grp)
      data.table(gene = g, factor = fac, stratum = s,
                 n_I = o$n1, n_II = o$n2,
                 rank_biserial = o$rb, p_MWU = o$p)
    }))
    r[, q_MWU := p.adjust(p_MWU, method = FDR_METHOD)]
    r
  }))
}))
fwrite(within, file.path(out_dir, "06_within_stratum_MWU.csv"))
if (nrow(within)) {
  smry <- within[, .(n_sig = sum(q_MWU < 0.05, na.rm = TRUE),
                     n_I = n_I[1], n_II = n_II[1]), by = .(factor, stratum)]
  for (i in seq_len(nrow(smry)))
    log_line("Within ", smry$factor[i], "=", smry$stratum[i],
             " (", smry$n_I[i], " vs ", smry$n_II[i], "): ",
             smry$n_sig[i], " genes q<0.05")
}

# SECONDARY: KRUSKAL-WALLIS I / II / III
coh123 <- pheno[stage %in% c("I", "II", "III")]
E123   <- E_collapsed[, coh123$SID, drop = FALSE]
log_line("Secondary KW cohort I/II/III: ",
         sum(coh123$stage == "I"), " / ", sum(coh123$stage == "II"), " / ",
         sum(coh123$stage == "III"))
sf <- factor(coh123$stage, levels = c("I", "II", "III"))
kw <- rbindlist(lapply(rownames(E123), function(g) {
  v <- E123[g, ]
  ok <- is.finite(v)
  k <- suppressWarnings(kruskal.test(v[ok] ~ sf[ok]))
 
  rho <- suppressWarnings(cor(as.integer(sf[ok]), v[ok], method = "spearman"))
  data.table(gene = g, KW_chisq = as.numeric(k$statistic), p_KW = k$p.value,
             spearman_rho_descriptive_only = rho)
}))
kw[, q_KW := p.adjust(p_KW, method = FDR_METHOD)]
setorder(kw, q_KW)
fwrite(kw, file.path(out_dir, "07_KW_stageI_II_III.csv"))
log_line("KW I/II/III: ", sum(kw$q_KW < 0.05, na.rm = TRUE), " genes q < 0.05")
log_line("   CAUTION: Stage III n is small; KW here is driven by the I/II mass.")

# SUPPLEMENTARY EARLY vs LATE 
if (RUN_EL_SUPPLEMENT) {
  cohEL <- pheno[stage_group %in% c("Early", "Late")]
  EEL   <- E_collapsed[, cohEL$SID, drop = FALSE]
  log_line("SUPPLEMENTARY Early vs Late: ", sum(cohEL$stage_group == "Early"),
           " vs ", sum(cohEL$stage_group == "Late"),
           "  <-- POWERED BY THE LATE n ONLY")
  gEL <- factor(cohEL$stage_group, levels = c("Early", "Late"))
  el <- rbindlist(lapply(rownames(EEL), function(g) {
    x <- EEL[g, gEL == "Early"]; y <- EEL[g, gEL == "Late"]
    x <- x[is.finite(x)]; y <- y[is.finite(y)]
    if (length(y) < MIN_PER_GROUP) return(NULL)
    w <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
    data.table(gene = g, n_early = length(x), n_late = length(y),
               rank_biserial = rank_biserial(x, y), p_MWU = w$p.value)
  }))
  el[, q_MWU := p.adjust(p_MWU, method = FDR_METHOD)]
  setorder(el, q_MWU)
  fwrite(el, file.path(out_dir, "08_supplementary_early_vs_late.csv"))
  log_line("Supp. EL: ", sum(el$q_MWU < 0.05, na.rm = TRUE), " genes q < 0.05",
           " (late n = ", el$n_late[1], ")")
}

## ------------------------------------------------------- 11. MANIFEST -----
writeLines(manifest, file.path(out_dir, "run_manifest.txt"))
message("\nDone. Outputs in: ", normalizePath(out_dir))
message("Paste run_manifest.txt into the decisions register before writing prose.")
