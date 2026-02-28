# downstream/R/06_trajectory.R
#
# Trajectory / Pseudotime Analysis — Monocle3
# Biological focus: Pancreatic Stellate Cell → CAF activation axis
# ─────────────────────────────────────────────────────────────
# PDAC creates a unique desmoplastic stroma driven by the activation
# of pancreatic stellate cells (PSC) into cancer-associated fibroblasts.
# This is not a binary switch — it is a continuous epigenetic transition:
#
#   Quiescent PSC → activated PSC → iCAF → myCAF
#                                 ↘ apCAF
#
# Monocle3 learns this trajectory using UMAP graph abstraction.
# With multiome data, we can ask:
#   - Which ATAC peaks open/close along the trajectory?
#   - Which TFs (chromVAR) change activity along the trajectory?
#   - Which genes drive stellate→CAF transition?
#
# Method:
#   Monocle3 builds a 'principal graph' through the WNN UMAP,
#   then assigns pseudotime as distance along that graph from
#   a user-specified root (quiescent stellate cells).

suppressPackageStartupMessages({
  library(monocle3)
  library(Seurat)
  library(SeuratWrappers)   # as.cell_data_set()
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(optparse)
})
set.seed(42)

opt_list <- list(
  make_option("--multiome_rds", type="character"),
  make_option("--archr_dir",    type="character"),
  make_option("--outdir",       type="character", default="results/trajectory")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

message("── Loading multiome object ──")
pdac <- readRDS(opt$multiome_rds)

# ── Subset to stromal / fibroblast lineage ────────────────────────────────
# Trajectory analysis works best within a lineage.
# Include stellate cells and all CAF subtypes + some tumor for context.
stromal_types <- c("Stellate_cells","myCAF","iCAF","apCAF","Tumor_classical")
pdac_stroma   <- pdac[, pdac$CellType %in% stromal_types]
message(sprintf("Stromal + tumor cells: %d", ncol(pdac_stroma)))

# ── Convert to CellDataSet (Monocle3 format) ──────────────────────────────
cds <- as.cell_data_set(pdac_stroma, assay="SCT")

# Transfer WNN UMAP embedding to Monocle3
reducedDims(cds)[["UMAP"]] <- Embeddings(pdac_stroma, reduction="UMAP_WNN")

# Monocle3 requires cluster assignments
cds@clusters[["UMAP"]]$clusters <- factor(pdac_stroma$CellType)

# ── Learn principal graph ──────────────────────────────────────────────────
# The principal graph is a set of connected nodes fit to the UMAP
# topology. It captures branching trajectories (one lineage can
# split into iCAF and myCAF from a common activated stellate cell state).
message("── Learning principal graph ──")
cds <- learn_graph(
  cds,
  use_partition = FALSE,  # use single connected graph
  close_loop    = FALSE,  # trajectories don't loop
  learn_graph_control = list(
    minimal_branch_len = 5,   # minimum cells per branch
    ncenter            = 100  # graph resolution
  )
)

# ── Assign root node (quiescent stellate cells) ───────────────────────────
# Root = cells with lowest ACTA2 expression (quiescent PSC marker)
# and highest LRAT expression (stellate quiescence marker)
root_cells <- colnames(pdac_stroma)[
  pdac_stroma$CellType == "Stellate_cells"
][1:min(20, sum(pdac_stroma$CellType == "Stellate_cells"))]

cds <- order_cells(cds, root_cells=root_cells)

# ── Plot trajectory ────────────────────────────────────────────────────────
p_traj_celltype <- plot_cells(
  cds,
  color_cells_by   = "CellType",
  label_groups_by_cluster = FALSE,
  label_leaves     = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 3
) +
  labs(title="Stellate Cell → CAF Trajectory (Monocle3)",
       subtitle="Arrows indicate direction of pseudotime progression") +
  theme(legend.text=element_text(size=9))

ggsave(file.path(opt$outdir, "trajectory_celltype.pdf"),
       p_traj_celltype, width=8, height=7)

p_traj_pseudo <- plot_cells(
  cds,
  color_cells_by      = "pseudotime",
  label_cell_groups   = FALSE,
  label_leaves        = FALSE,
  label_branch_points = FALSE
) +
  scale_color_viridis_c(option="magma") +
  labs(title="Pseudotime — Stellate → CAF Activation Axis",
       subtitle="Root = quiescent pancreatic stellate cells")

ggsave(file.path(opt$outdir, "trajectory_pseudotime.pdf"),
       p_traj_pseudo, width=7, height=6)

# ── Find trajectory-variable genes ────────────────────────────────────────
# graph_test identifies genes whose expression changes significantly
# along the principal graph (Moran's I spatial autocorrelation test)
message("── Testing trajectory-variable genes (Moran's I) ──")
pr_test_res <- graph_test(
  cds,
  neighbor_graph = "principal_graph",
  cores          = 4
)

pr_deg <- pr_test_res %>%
  filter(q_value < 0.05, morans_I > 0.1) %>%
  arrange(desc(morans_I))

write.csv(pr_deg, file.path(opt$outdir, "trajectory_variable_genes.csv"),
          row.names=FALSE)

message(sprintf("  Trajectory-variable genes: %d", nrow(pr_deg)))

# ── Plot gene expression along pseudotime ─────────────────────────────────
# Key genes for stellate → CAF transition:
#   Early (quiescent): LRAT, RBP1, PDGFRA (stellate markers)
#   Mid (activation):  VIM, ACTA2, FAP (early activation)
#   Late myCAF:        MMP11, POSTN, HOPX
#   Late iCAF:         IL6, CXCL12, HAS2
traj_genes <- c("LRAT","PDGFRA","VIM","ACTA2","FAP","MMP11","POSTN","IL6","CXCL12")
traj_genes <- intersect(traj_genes, rownames(cds))

if (length(traj_genes) > 0) {
  p_genes <- plot_genes_in_pseudotime(
    cds[traj_genes, ],
    color_cells_by = "CellType",
    min_expr       = 0.5,
    ncol           = 3
  )
  ggsave(file.path(opt$outdir, "trajectory_key_genes.pdf"),
         p_genes, width=12, height=8)
}

# ── ATAC: chromatin changes along trajectory ──────────────────────────────
# Transfer pseudotime back to Seurat and correlate with ATAC peaks
pdac_stroma$pseudotime <- pseudotime(cds)[colnames(pdac_stroma)]

# Find peaks that correlate with pseudotime
# (proxy for peaks that open during stellate→CAF activation)
pseudo_vals  <- pdac_stroma$pseudotime
peak_mat     <- GetAssayData(pdac_stroma, assay="ATAC")
finite_cells <- !is.na(pseudo_vals) & is.finite(pseudo_vals)

message("── Correlating ATAC peaks with pseudotime ──")
peak_pseudo_cor <- apply(peak_mat[, finite_cells], 1, function(x) {
  tryCatch(
    cor(x, pseudo_vals[finite_cells], method="spearman"),
    error = function(e) NA
  )
})

pseudo_peaks <- data.frame(
  peak        = names(peak_pseudo_cor),
  spearman_r  = peak_pseudo_cor
) %>%
  filter(!is.na(spearman_r)) %>%
  arrange(desc(abs(spearman_r)))

write.csv(pseudo_peaks,
          file.path(opt$outdir, "pseudotime_correlated_peaks.csv"),
          row.names=FALSE)

# Peaks opening along CAF activation (positive correlation with pseudotime)
top_opening <- head(pseudo_peaks %>% filter(spearman_r > 0.3), 20)
message(sprintf("  Peaks opening during CAF activation: %d",
                sum(pseudo_peaks$spearman_r > 0.3, na.rm=TRUE)))

# ── Plot: pseudotime vs key TF accessibility ──────────────────────────────
pdac_stroma$pseudo_bin <- cut(pdac_stroma$pseudotime,
                               breaks=10, labels=FALSE)

pseudo_tf_df <- data.frame(
  pseudotime = pdac_stroma$pseudotime,
  pseudo_bin = pdac_stroma$pseudo_bin,
  CellType   = pdac_stroma$CellType
) %>%
  filter(!is.na(pseudotime))

p_pseudo_density <- ggplot(pseudo_tf_df,
  aes(x=pseudotime, fill=CellType)) +
  geom_density(alpha=0.6) +
  scale_fill_manual(values=c(
    Stellate_cells="#264653", myCAF="#F4A261",
    iCAF="#E9C46A", apCAF="#2A9D8F", Tumor_classical="#D7191C"
  )) +
  labs(title="Cell Type Distribution Along CAF Pseudotime",
       subtitle="Stellate cells at low pseudotime → CAF subtypes diverge at high pseudotime",
       x="Pseudotime", y="Density", fill="Cell Type") +
  theme_bw(base_size=13)

ggsave(file.path(opt$outdir, "pseudotime_celltype_density.pdf"),
       p_pseudo_density, width=8, height=5)

# Save pseudotime metadata for cross-modal script
pdac_meta <- data.frame(
  cell      = colnames(pdac_stroma),
  pseudotime = pdac_stroma$pseudotime,
  CellType  = pdac_stroma$CellType
)
write.csv(pdac_meta, file.path(opt$outdir, "pseudotime_metadata.csv"),
          row.names=FALSE)

message("── Trajectory analysis complete ──")
