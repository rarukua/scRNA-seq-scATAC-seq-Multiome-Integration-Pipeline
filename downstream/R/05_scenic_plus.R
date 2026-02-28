# downstream/R/05_scenic_plus.R
#
# SCENIC+ eRegulon visualization
# Reads Python SCENIC+ outputs → publication-quality figures

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr)
  library(ComplexHeatmap); library(circlize); library(optparse)
})

opt_list <- list(
  make_option("--multiome_rds",    type="character"),
  make_option("--scenic_plus_dir", type="character"),
  make_option("--outdir",          type="character", default="results/scenic_plus_viz")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

pdac     <- readRDS(opt$multiome_rds)
auc_df   <- read.csv(file.path(opt$scenic_plus_dir, "eRegulon_AUC.csv"),
                      row.names=1, check.names=FALSE)
reg_sum  <- read.csv(file.path(opt$scenic_plus_dir, "eregulon_summary.csv"))

# Align cells
shared <- intersect(colnames(pdac), rownames(auc_df))
auc_mat <- t(as.matrix(auc_df[shared, ]))

# ── PDAC TF groups ─────────────────────────────────────────────────────────
pdac_tf_groups <- list(
  "KRAS-AP1"     = c("FOSL2","JUNB","FOSB","FOS","JUN","JUND"),
  "Lineage-Cls"  = c("FOXA1","FOXA2","GATA6","HNF1B"),
  "Basal-EMT"    = c("TP63","ZEB1","ZEB2","SNAI2","VIM"),
  "Immunosuppres"= c("STAT3","IRF4","NFKB1","RELA"),
  "CAF-Stroma"   = c("SMAD3","TGFB1","TWIST1","ACTA2")
)

# ── Figure 1: eRegulon heatmap across cell types ──────────────────────────
mean_auc_by_ct <- do.call(cbind, lapply(unique(pdac$CellType), function(ct) {
  cells <- colnames(pdac)[pdac$CellType == ct]
  cells <- intersect(cells, colnames(auc_mat))
  if (length(cells) < 5) return(NULL)
  rowMeans(auc_mat[, cells, drop=FALSE])
}))
colnames(mean_auc_by_ct) <- unique(pdac$CellType)
mean_auc_by_ct <- mean_auc_by_ct[, !is.null(colnames(mean_auc_by_ct))]

# Keep regulons with high variance across cell types (informative ones)
row_var <- apply(mean_auc_by_ct, 1, var)
top_regs <- names(sort(row_var, decreasing=TRUE))[1:min(50, nrow(mean_auc_by_ct))]
mat_plot  <- t(scale(t(mean_auc_by_ct[top_regs, ])))

# Row annotation: TF group
tf_group_vec <- sapply(rownames(mat_plot), function(tf) {
  grp <- names(which(sapply(pdac_tf_groups, function(g) any(grepl(tf, g, ignore.case=TRUE)))))
  if (length(grp)) grp[1] else "Other"
})
row_anno <- rowAnnotation(
  TF_Group = tf_group_vec,
  col      = list(TF_Group = c(
    "KRAS-AP1"="darkred","Lineage-Cls"="steelblue",
    "Basal-EMT"="darkorange","Immunosuppres"="purple",
    "CAF-Stroma"="darkgreen","Other"="grey80"
  )),
  annotation_legend_param = list(title="TF Group")
)

pdf(file.path(opt$outdir, "scenic_eregulon_heatmap.pdf"), width=12, height=14)
Heatmap(mat_plot,
        name               = "Scaled\neRegulon AUC",
        col                = colorRamp2(c(-2,0,2), c("#2C7BB6","white","#D7191C")),
        right_annotation   = row_anno,
        column_names_rot   = 45,
        row_names_gp       = gpar(fontsize=7),
        column_names_gp    = gpar(fontsize=9),
        clustering_method_rows    = "ward.D2",
        clustering_method_columns = "ward.D2",
        column_title       = "PDAC Cell Types",
        row_title          = "Top 50 eRegulons",
        heatmap_legend_param = list(title="Scaled AUC"))
dev.off()

# ── Figure 2: AP-1 regulon detail (KRAS → AP-1 → target genes) ────────────
# This figure tells the biological story in PDAC:
# KRAS mutation → sustained ERK → AP-1 TF activation →
# AP-1 opens enhancers → activates proliferation/survival genes
fosl2_col <- grep("FOSL2", colnames(auc_df), value=TRUE, ignore.case=TRUE)[1]

if (!is.na(fosl2_col)) {
  pdac$FOSL2_AUC <- auc_df[colnames(pdac), fosl2_col]
  p_fosl2 <- FeaturePlot(pdac, reduction="UMAP_WNN",
                          features="FOSL2_AUC", pt.size=0.3) +
    scale_color_gradientn(colors=c("grey90","#FDAE61","#D7191C")) +
    labs(title="FOSL2 (AP-1) eRegulon Activity",
         subtitle="KRAS-driven AP-1 activation — highest in tumor cells",
         color="eRegulon\nAUC")
  ggsave(file.path(opt$outdir, "fosl2_eregulon_umap.pdf"),
         p_fosl2, width=7, height=6)
}

# ── Figure 3: Regulon specificity — which TFs define each cell type? ───────
specificity <- as.data.frame(mean_auc_by_ct) %>%
  tibble::rownames_to_column("TF") %>%
  tidyr::pivot_longer(-TF, names_to="CellType", values_to="AUC") %>%
  group_by(TF) %>%
  mutate(specificity_score = (AUC - mean(AUC)) / (sd(AUC) + 1e-6)) %>%
  ungroup() %>%
  group_by(CellType) %>%
  slice_max(specificity_score, n=5) %>%
  ungroup()

p_spec <- ggplot(specificity,
  aes(x=CellType, y=reorder(TF, specificity_score),
      fill=specificity_score, size=AUC)) +
  geom_point(shape=21, stroke=0.3) +
  scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C",
                       midpoint=0, name="Specificity") +
  scale_size_continuous(range=c(2,8), name="Mean AUC") +
  labs(title="Cell-Type-Specific TF Regulons (SCENIC+)",
       subtitle="Top 5 specific TFs per cell type",
       x=NULL, y="TF Regulon") +
  theme_bw(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(file.path(opt$outdir, "scenic_tf_specificity_dotplot.pdf"),
       p_spec, width=11, height=8)

message("── SCENIC+ visualization complete ──")
