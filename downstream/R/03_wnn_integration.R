# downstream/R/03_wnn_integration.R
#
# Weighted Nearest Neighbor (WNN) Integration
# RNA (Seurat) + ATAC (ArchR/Signac) → Joint multiome embedding
# ─────────────────────────────────────────────────────────────
# WNN Algorithm (Hao et al. 2021 Cell):
#
# For each cell, WNN computes:
#   1. k-nearest neighbors in RNA space (using harmony PCA)
#   2. k-nearest neighbors in ATAC space (using harmony LSI)
#   3. A per-cell WEIGHT for each modality based on how well
#      that modality predicts the cell's neighbors in the other
#
# The weight reflects INFORMATION CONTENT per cell:
#   - Cell with high-quality RNA + poor ATAC → high RNA weight
#   - Cell with high-quality ATAC + poor RNA → high ATAC weight
#   - Most cells → intermediate weight (~0.5 each)
#
# This is more principled than simple concatenation or averaging
# because it adapts to data quality per cell rather than assuming
# both modalities are equally informative everywhere.
#
# The result: a WNN graph where edges reflect similarity in the
# COMBINED epigenomic + transcriptomic space → better cell type
# separation than either modality alone.
# ─────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ArchR)
  library(ggplot2)
  library(dplyr)
  library(optparse)
})
set.seed(42)

opt_list <- list(
  make_option("--seurat_rds",   type="character", help="pdac_rna.rds from 02_seurat_rna.R"),
  make_option("--archr_dir",    type="character", help="ArchR_PDAC/ directory"),
  make_option("--atac_lsi",     type="character", help="atac_lsi_harmony.rds"),
  make_option("--atac_peaks",   type="character", help="atac_peak_matrix.rds"),
  make_option("--outdir",       type="character", default="results/wnn")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

# ── Load RNA Seurat object ─────────────────────────────────────────────────
message("── Loading RNA Seurat object ──")
pdac <- readRDS(opt$seurat_rds)

# ── Load ATAC data from ArchR exports ─────────────────────────────────────
message("── Loading ATAC data from ArchR ──")
lsi_matrix  <- readRDS(opt$atac_lsi)   # cells × LSI components
peak_matrix <- readRDS(opt$atac_peaks)  # cells × peaks sparse matrix

# ── Align barcodes between modalities ─────────────────────────────────────
# Critical step: ensure exactly the same cells are in both modalities.
# In multiome, barcodes match but some cells may be filtered in one
# modality but not the other during QC.
shared_barcodes <- intersect(colnames(pdac),
                              rownames(lsi_matrix))

message(sprintf("  RNA cells: %d", ncol(pdac)))
message(sprintf("  ATAC cells: %d", nrow(lsi_matrix)))
message(sprintf("  Shared cells: %d", length(shared_barcodes)))

pdac <- pdac[, shared_barcodes]
lsi_matrix  <- lsi_matrix[shared_barcodes, ]

# ── Add ATAC assay to Seurat object (Signac) ──────────────────────────────
# Signac extends Seurat with a ChromatinAssay that stores:
#   - Peak counts matrix
#   - Fragment file path (for computing signal tracks, footprinting)
#   - Genome annotation (for peak-to-gene distance calculations)
message("── Creating ChromatinAssay ──")

# Subset peak matrix to shared barcodes
peak_counts <- assay(peak_matrix)[, shared_barcodes]

# Create Signac ChromatinAssay
chrom_assay <- CreateChromatinAssay(
  counts     = peak_counts,
  sep        = c(":", "-"),       # peak name format: "chr1:1000-2000"
  genome     = "hg38",
  fragments  = opt$fragments_path %||% NULL,
  min.cells  = 10,
  min.features = 200
)

pdac[["ATAC"]] <- chrom_assay

# ── Add LSI reduction ──────────────────────────────────────────────────────
# Transfer ArchR's Harmony-corrected LSI into Seurat
pdac[["lsi"]] <- CreateDimReducObject(
  embeddings = lsi_matrix[shared_barcodes, 1:30],
  key        = "LSI_",
  assay      = "ATAC"
)

# ── Normalize ATAC assay ──────────────────────────────────────────────────
# TF-IDF normalization (standard for chromatin accessibility data)
DefaultAssay(pdac) <- "ATAC"
pdac <- RunTFIDF(pdac)
pdac <- FindTopFeatures(pdac, min.cutoff="q0")

# ── WNN: Find multimodal neighbors ────────────────────────────────────────
# This is the key step — computes per-cell modality weights
# and builds the joint KNN graph
message("── Computing WNN ──")

DefaultAssay(pdac) <- "SCT"   # use SCT for RNA reduction

pdac <- FindMultiModalNeighbors(
  object         = pdac,
  reduction.list = list("harmony", "lsi"),  # RNA PCA, ATAC LSI
  dims.list      = list(1:30, 2:30),        # dim 1 of LSI = depth (skip it)
  modality.weight.name = c("RNA.weight","ATAC.weight"),
  verbose        = TRUE
)

# ── WNN UMAP ──────────────────────────────────────────────────────────────
pdac <- RunUMAP(
  pdac,
  nn.name        = "weighted.nn",
  reduction.name = "UMAP_WNN",
  reduction.key  = "wnnUMAP_",
  verbose        = TRUE
)

# ── WNN clustering ────────────────────────────────────────────────────────
pdac <- FindClusters(
  pdac,
  graph.name = "wsnn",
  algorithm  = 3,        # SLM algorithm (better for multimodal)
  resolution = 0.6,
  verbose    = TRUE
)
pdac$WNN_clusters <- Idents(pdac)

# ── Compare RNA, ATAC, and WNN UMAPs ──────────────────────────────────────
# This 3-panel figure is the key visualisation for a WNN paper/portfolio:
# Shows that WNN captures structure from BOTH modalities
p_compare <- cowplot::plot_grid(
  DimPlot(pdac, reduction="UMAP_RNA",  group.by="WNN_clusters",
          label=TRUE, pt.size=0.1) + ggtitle("RNA UMAP") + NoLegend(),
  DimPlot(pdac, reduction="UMAP_ATAC", group.by="WNN_clusters",
          label=TRUE, pt.size=0.1) + ggtitle("ATAC UMAP") + NoLegend(),
  DimPlot(pdac, reduction="UMAP_WNN",  group.by="WNN_clusters",
          label=TRUE, pt.size=0.1) + ggtitle("WNN UMAP (joint)"),
  ncol=3
)
ggsave(file.path(opt$outdir, "umap_rna_atac_wnn_comparison.pdf"),
       p_compare, width=18, height=6)

# ── Visualize modality weights ─────────────────────────────────────────────
# This is a unique plot only possible with WNN:
# Shows which cells rely more on RNA vs ATAC for their neighborhood
p_weights <- cowplot::plot_grid(
  FeaturePlot(pdac, reduction="UMAP_WNN",
              features="RNA.weight", pt.size=0.1) +
    scale_color_gradient2(low="#2C7BB6", mid="white", high="#D7191C",
                          midpoint=0.5) +
    ggtitle("RNA modality weight per cell"),
  FeaturePlot(pdac, reduction="UMAP_WNN",
              features="ATAC.weight", pt.size=0.1) +
    scale_color_gradient2(low="#2C7BB6", mid="white", high="#D7191C",
                          midpoint=0.5) +
    ggtitle("ATAC modality weight per cell"),
  ncol=2
)
ggsave(file.path(opt$outdir, "wnn_modality_weights.pdf"),
       p_weights, width=12, height=5)

# ── Cell type annotation on WNN clusters ──────────────────────────────────
# Map cluster IDs to PDAC cell types using RNA markers from 02_seurat_rna.R
# (adjusted based on WNN cluster composition)
cell_type_annotations <- c(
  "0"="Tumor_classical", "1"="myCAF",       "2"="CD8_T_exhausted",
  "3"="Macrophage_M2",   "4"="Tumor_basal", "5"="iCAF",
  "6"="Ductal_normal",   "7"="CD4_Treg",    "8"="Endothelial",
  "9"="B_cell",          "10"="NK_cell",    "11"="apCAF",
  "12"="Acinar",         "13"="Macrophage_M1"
)

pdac$CellType <- recode(as.character(pdac$WNN_clusters),
                         !!!cell_type_annotations)

# Cell type colour palette
cell_type_colors <- c(
  Tumor_classical="#D7191C", Tumor_basal="#E85D04",
  myCAF="#F4A261",  iCAF="#E9C46A",  apCAF="#2A9D8F",
  Ductal_normal="#264653",
  Macrophage_M2="#6A4C93", Macrophage_M1="#9B5DE5",
  CD8_T_exhausted="#0077B6", CD4_Treg="#00B4D8",
  B_cell="#48CAE4",  NK_cell="#90E0EF",
  Endothelial="#ADE8F4", Acinar="#F8F9FA"
)

p_celltypes <- DimPlot(pdac, reduction="UMAP_WNN",
                        group.by="CellType",
                        cols=cell_type_colors,
                        label=TRUE, label.size=3,
                        pt.size=0.3) +
  labs(title="PDAC Multiome — WNN Cell Types",
       subtitle="Joint RNA + ATAC embedding") +
  theme(legend.text=element_text(size=8))

ggsave(file.path(opt$outdir, "umap_wnn_celltypes.pdf"),
       p_celltypes, width=10, height=8)

# ── Save integrated object ────────────────────────────────────────────────
saveRDS(pdac, file.path(opt$outdir, "pdac_multiome.rds"))
message("── WNN integration complete ──")
message("  Saved: pdac_multiome.rds")
message("  Next steps:")
message("    04_peak_to_gene.R   — cis-regulatory element linkage")
message("    05_scenic_plus.R    — TF regulon inference (run Python first)")
message("    06_trajectory.R     — pseudotime along CAF axis")
message("    07_cellchat.R       — cell-cell communication")
