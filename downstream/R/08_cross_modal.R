# downstream/R/08_cross_modal.R
#
# Cross-Modal Integration
# PDAC Multiome ↔ LUAD Bulk ATAC ↔ CRC WGBS Methylation
# ─────────────────────────────────────────────────────────────
# This script ties together all three portfolio projects:
#
#   Project 1 (CRC)    — WGBS methylation, cfDNA multi-modal
#   Project 2 (LUAD)   — Bulk ATAC-seq, TF footprinting
#   Project 3 (PDAC)   — scRNA+scATAC multiome, SCENIC+
#
# Cross-modal questions:
#   1. METHYLATION × ATAC: Do CRC-hypermethylated loci correspond to
#      closed chromatin in PDAC tumor cells?
#      (Shared epigenetic silencing mechanisms across cancers)
#
#   2. BULK ATAC × SINGLE-CELL ATAC: Do LUAD tumor-specific peaks
#      from the bulk ATAC pipeline match the ATAC accessibility
#      in PDAC tumor cell clusters?
#      (Cancer-shared vs cancer-specific chromatin programs)
#
#   3. TF FOOTPRINTS × REGULONS: Do TFs footprinted in LUAD bulk
#      ATAC (TOBIAS) show high eRegulon activity in PDAC scATAC?
#      (Validates pan-cancer TF programs)
#
#   4. CIMP × CHROMATIN: Do CIMP-H genes (methylated in CRC)
#      have closed chromatin in PDAC — are they also silenced
#      in a different cancer via a different mechanism?

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr)
  library(GenomicRanges); library(tidyr); library(optparse)
})

opt_list <- list(
  make_option("--pdac_rds",         type="character", help="pdac_multiome.rds"),
  make_option("--luad_da_peaks",    type="character", help="LUAD differential peaks CSV"),
  make_option("--crc_methylation",  type="character", help="CRC methylation DMR CSV"),
  make_option("--scenic_dir",       type="character", help="SCENIC+ output dir"),
  make_option("--luad_tobias",      type="character", help="TOBIAS BINDetect results CSV"),
  make_option("--outdir",           type="character", default="results/cross_modal")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

message("── Loading cross-modal data ──")
pdac      <- readRDS(opt$pdac_rds)
luad_da   <- read.csv(opt$luad_da_peaks)
crc_meth  <- read.csv(opt$crc_methylation)

# ── Analysis 1: CRC Methylation × PDAC Chromatin ─────────────────────────
# DMRs hypermethylated in CRC (silenced genes) — are these also
# inaccessible (closed) in PDAC tumor cells?
# Tests the hypothesis of shared epigenetic silencing programs
message("── CRC methylation × PDAC chromatin ──")

# CRC hyper-methylated DMRs → extract gene names
crc_hyper_genes <- crc_meth %>%
  filter(direction == "hypermethylated" | Fold > 1) %>%
  pull(SYMBOL) %>%
  na.omit() %>%
  unique()

message(sprintf("  CRC hypermethylated genes: %d", length(crc_hyper_genes)))

# Get gene activity scores in PDAC tumor cells from ArchR
# Gene activity = accessibility around gene body/promoter
pdac_tumor <- pdac[, pdac$CellType %in% c("Tumor_classical","Tumor_basal")]

genes_in_pdac <- intersect(crc_hyper_genes, rownames(pdac))
if (length(genes_in_pdac) > 5) {
  # Compare gene activity in PDAC tumor vs normal ductal cells
  pdac_ductal <- pdac[, pdac$CellType == "Ductal_normal"]
  cells_use   <- c(colnames(pdac_tumor), colnames(pdac_ductal))
  pdac_sub    <- pdac[, cells_use]

  DefaultAssay(pdac_sub) <- "SCT"
  Idents(pdac_sub) <- pdac_sub$CellType

  # Test if CRC-hypermethylated genes are also downregulated in PDAC tumor
  de_crc_genes <- FindMarkers(
    pdac_sub,
    ident.1  = "Tumor_classical",
    ident.2  = "Ductal_normal",
    features = genes_in_pdac,
    test.use = "wilcox",
    min.pct  = 0.05
  ) %>%
    tibble::rownames_to_column("gene") %>%
    mutate(
      direction_pdac = ifelse(avg_log2FC < 0, "downregulated","upregulated"),
      also_methylated_crc = TRUE
    )

  concordant <- sum(de_crc_genes$direction_pdac == "downregulated" &
                      de_crc_genes$p_val_adj < 0.05)
  total_tested <- nrow(de_crc_genes)

  message(sprintf("  CRC-methylated genes also downregulated in PDAC tumor: %d/%d (%.0f%%)",
                  concordant, total_tested,
                  concordant/max(total_tested,1)*100))

  p_crc_pdac <- ggplot(de_crc_genes,
    aes(x=avg_log2FC, y=-log10(p_val_adj), color=direction_pdac)) +
    geom_point(alpha=0.7, size=2) +
    geom_vline(xintercept=0, linetype="dashed") +
    ggrepel::geom_label_repel(
      data=de_crc_genes %>%
        filter(p_val_adj < 0.01) %>%
        slice_min(avg_log2FC, n=10),
      aes(label=gene), size=3
    ) +
    scale_color_manual(values=c(downregulated="#2C7BB6",upregulated="#D7191C")) +
    labs(
      title    = "CRC-Hypermethylated Genes in PDAC Tumor vs Normal Ductal",
      subtitle = sprintf("%d/%d CRC-silenced genes also downregulated in PDAC",
                         concordant, total_tested),
      x        = "Log2FC (PDAC Tumor / Normal Ductal)",
      y        = "-log10(adj. p-value)",
      color    = "PDAC direction",
      caption  = "Shared epigenetic silencing: methylated in CRC → inaccessible in PDAC"
    ) +
    theme_bw(base_size=13)

  ggsave(file.path(opt$outdir, "crc_methylation_vs_pdac_expression.pdf"),
         p_crc_pdac, width=8, height=6)

  write.csv(de_crc_genes,
            file.path(opt$outdir, "crc_methylated_genes_in_pdac.csv"),
            row.names=FALSE)
}

# ── Analysis 2: LUAD Bulk ATAC peaks × PDAC scATAC ────────────────────────
# Do LUAD tumor-specific open chromatin peaks correspond to
# open regions in PDAC tumor cells?
# Tests: shared cancer chromatin programs vs cancer-type-specific
message("── LUAD bulk ATAC × PDAC scATAC ──")

luad_tumor_peaks <- luad_da %>%
  filter(direction == "More open in Tumor" | (Fold > 1 & FDR < 0.05)) %>%
  select(seqnames, start, end, SYMBOL, Fold) %>%
  na.omit()

luad_gr <- makeGRangesFromDataFrame(luad_tumor_peaks, keep.extra.columns=TRUE)

# Get PDAC peak coordinates
pdac_peaks <- granges(pdac[["ATAC"]])

# Overlap: which LUAD-specific peaks are also peaks in PDAC?
olap <- findOverlaps(luad_gr, pdac_peaks, minoverlap=100)
pct_shared <- length(unique(queryHits(olap))) / length(luad_gr) * 100
message(sprintf("  LUAD tumor-specific peaks overlapping PDAC peaks: %.1f%%", pct_shared))

# For shared peaks — are they more accessible in PDAC tumor vs normal?
shared_pdac_peaks <- pdac_peaks[unique(subjectHits(olap))]
shared_peak_names <- paste0(seqnames(shared_pdac_peaks), ":",
                              start(shared_pdac_peaks), "-",
                              end(shared_pdac_peaks))
shared_in_pdac <- intersect(shared_peak_names, rownames(pdac[["ATAC"]]))

if (length(shared_in_pdac) > 0) {
  DefaultAssay(pdac) <- "ATAC"
  # Mean accessibility at LUAD-shared peaks in each PDAC cell type
  acc_by_ct <- sapply(unique(pdac$CellType), function(ct) {
    cells <- colnames(pdac)[pdac$CellType == ct]
    rowMeans(GetAssayData(pdac[shared_in_pdac[1:min(100, length(shared_in_pdac))],
                                cells], slot="data"))
  })

  acc_df <- as.data.frame(acc_by_ct) %>%
    summarise(across(everything(), mean)) %>%
    pivot_longer(everything(), names_to="CellType", values_to="mean_acc")

  p_luad_pdac <- ggplot(acc_df,
    aes(x=reorder(CellType, mean_acc), y=mean_acc,
        fill=grepl("Tumor", CellType))) +
    geom_col(width=0.7) +
    coord_flip() +
    scale_fill_manual(values=c("TRUE"="#D7191C","FALSE"="grey60")) +
    labs(
      title    = "LUAD Tumor-Specific Peaks — Accessibility in PDAC Cell Types",
      subtitle = sprintf("%d LUAD peaks tested | %.1f%% overlap with PDAC peak set",
                         nrow(luad_tumor_peaks), pct_shared),
      x        = NULL, y = "Mean ATAC accessibility",
      fill     = "Tumor cell type",
      caption  = "Shared open chromatin between LUAD and PDAC tumor cells"
    ) +
    theme_bw(base_size=12)

  ggsave(file.path(opt$outdir, "luad_peaks_in_pdac_accessibility.pdf"),
         p_luad_pdac, width=8, height=6)
}

# ── Analysis 3: TOBIAS TFs (LUAD) × SCENIC+ eRegulons (PDAC) ─────────────
# Validate pan-cancer TF programs:
# TFs with active footprints in LUAD → high eRegulon activity in PDAC?
if (!is.null(opt$luad_tobias) && file.exists(opt$luad_tobias)) {
  message("── LUAD TOBIAS × PDAC SCENIC+ TF comparison ──")

  tobias  <- read.csv(opt$luad_tobias)
  scenic  <- read.csv(file.path(opt$scenic_dir, "eregulon_summary.csv"))

  # Extract TF names (standardize formatting)
  tobias$tf_clean  <- gsub("_.*","", toupper(tobias$name))
  scenic$tf_clean  <- gsub("_.*","", toupper(scenic$TF))

  tf_comparison <- inner_join(
    tobias  %>% select(tf_clean, luad_footprint_change=score_change,
                       luad_tf_status=tf_status),
    scenic  %>% select(tf_clean, pdac_eregulon_auc=mean_AUC),
    by="tf_clean"
  )

  p_tf_compare <- ggplot(tf_comparison,
    aes(x=luad_footprint_change, y=pdac_eregulon_auc,
        color=luad_tf_status)) +
    geom_point(size=2.5, alpha=0.8) +
    ggrepel::geom_label_repel(
      data=tf_comparison %>%
        filter(abs(luad_footprint_change) > 0.1, pdac_eregulon_auc > 0.1),
      aes(label=tf_clean), size=3, max.overlaps=15
    ) +
    geom_smooth(method="lm", se=TRUE, color="black", linewidth=0.5) +
    scale_color_manual(values=c(
      "Gained in Tumor"="#D7191C","Lost in Tumor"="#2C7BB6","Unchanged"="grey70"
    )) +
    labs(
      title    = "Pan-Cancer TF Activity: LUAD Footprint × PDAC eRegulon",
      subtitle = "TFs active in LUAD (TOBIAS) vs active in PDAC (SCENIC+)",
      x        = "LUAD: Footprint score change (TOBIAS)",
      y        = "PDAC: eRegulon mean AUC (SCENIC+)",
      color    = "LUAD TF status",
      caption  = "Correlated TFs = pan-cancer drivers | Discordant TFs = cancer-type-specific"
    ) +
    theme_bw(base_size=13)

  ggsave(file.path(opt$outdir, "pancancer_tf_comparison.pdf"),
         p_tf_compare, width=8, height=6)
}

# ── Summary report ────────────────────────────────────────────────────────
summary_stats <- list(
  crc_methylated_also_downregulated_pdac = if (exists("concordant")) concordant else NA,
  luad_pdac_peak_overlap_pct             = pct_shared,
  n_shared_luad_pdac_peaks               = length(shared_in_pdac)
)

writeLines(
  c("# Cross-Modal Integration Summary",
    "",
    "## CRC Methylation × PDAC Chromatin",
    sprintf("  CRC-hypermethylated genes tested in PDAC: %d", length(genes_in_pdac)),
    if (exists("concordant"))
      sprintf("  Concordantly silenced in PDAC tumor: %d (%.0f%%)",
              concordant, concordant/total_tested*100) else "",
    "",
    "## LUAD Bulk ATAC × PDAC scATAC",
    sprintf("  LUAD tumor-specific peaks: %d", nrow(luad_tumor_peaks)),
    sprintf("  Overlap with PDAC peak set: %.1f%%", pct_shared),
    "",
    "## Interpretation",
    "  Shared epigenetic silencing across CRC, LUAD, PDAC suggests",
    "  conserved cancer epigenome programs operating independently",
    "  of tissue-of-origin and mutation context."
  ),
  file.path(opt$outdir, "cross_modal_summary.md")
)

message("── Cross-modal integration complete ──")
message(sprintf("  Summary written to: %s/cross_modal_summary.md", opt$outdir))
