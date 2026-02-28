# downstream/R/07_cellchat.R
#
# Cell-Cell Communication Analysis — CellChat
# ─────────────────────────────────────────────────────────────
# Reference: Jin et al. 2021 Nature Communications
#
# CellChat infers intercellular communication from scRNA-seq by:
#   1. Using a curated database of ligand-receptor (L-R) pairs
#   2. Computing communication probability from L-R expression levels
#      weighted by cell fraction (more cells expressing → stronger signal)
#   3. Modeling multi-subunit complexes (e.g. TGF-β receptor complex)
#   4. Identifying signaling pathways from groups of L-R pairs
#
# Why multiome data improves CellChat:
#   - Cell type annotation is more accurate (WNN uses both RNA + ATAC)
#   - Rarer cell populations are better resolved → better communication inference
#   - ATAC confirms that receptor genes have accessible promoters
#     (a ligand-receptor pair where the receptor's promoter is closed
#     cannot realistically mediate communication)
#
# PDAC biological focus:
#   - Tumor → Macrophage: CSF1 → CSF1R (macrophage recruitment)
#   - Tumor → CAF: TGF-β → TGFBR (CAF activation)
#   - CAF → T cell: CXCL12 → CXCR4 (T cell exclusion)
#   - Macrophage → T cell: PD-L1 → PD-1 (exhaustion induction)
#   - Tumor → T cell: all immunosuppressive signals

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(optparse)
})
options(stringsAsFactors=FALSE)

opt_list <- list(
  make_option("--multiome_rds", type="character"),
  make_option("--outdir",       type="character", default="results/cellchat")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

message("── Loading multiome object ──")
pdac <- readRDS(opt$multiome_rds)

# ── Prepare CellChat input ────────────────────────────────────────────────
# CellChat uses normalized RNA expression (not raw counts)
DefaultAssay(pdac) <- "SCT"
expr_mat  <- GetAssayData(pdac, slot="data")
meta      <- data.frame(
  cell_type = pdac$CellType,
  row.names = colnames(pdac)
)

# ── Create CellChat object ────────────────────────────────────────────────
cellchat <- createCellChat(
  object   = expr_mat,
  meta     = meta,
  group.by = "cell_type"
)

# Use human L-R database (CellChatDB.human contains 2,021 validated L-R pairs)
CellChatDB  <- CellChatDB.human
# Focus on Secreted Signaling (most relevant for TME)
# Options: "Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor"
CellChatDB_use <- subsetDB(CellChatDB,
                             search=c("Secreted Signaling","Cell-Cell Contact"))
cellchat@DB    <- CellChatDB_use

# ── Preprocessing ─────────────────────────────────────────────────────────
# Subset genes to those in the L-R database (speeds up computation)
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# ── Compute communication probabilities ───────────────────────────────────
# population.size=TRUE: weight by cell fraction
# (important for PDAC where tumor cells dominate)
message("── Computing communication probabilities ──")
cellchat <- computeCommunProb(
  cellchat,
  type           = "triMean",  # robust aggregation
  population.size = TRUE
)

# Filter communications with few cells
cellchat <- filterCommunication(cellchat, min.cells=10)

# Infer pathway-level communication (groups L-R pairs)
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate cell-level → network-level statistics
cellchat <- aggregateNet(cellchat)

# ── Summary statistics ────────────────────────────────────────────────────
n_interactions <- nrow(cellchat@net$count)
message(sprintf("  Inferred interactions: %d", sum(cellchat@net$count)))
message(sprintf("  Signaling pathways: %d",
                length(cellchat@netP$pathways)))

# ── Figure 1: Global interaction network ──────────────────────────────────
# Chord diagram showing total communication strength between cell types
pdf(file.path(opt$outdir, "cellchat_global_chord.pdf"), width=10, height=10)
netVisual_circle(
  cellchat@net$count,
  vertex.weight    = as.numeric(table(cellchat@idents)),
  weight.scale     = TRUE,
  label.edge.weight= FALSE,
  title.name       = "Number of Interactions — PDAC TME"
)
dev.off()

# ── Figure 2: Heatmap — interaction strength ──────────────────────────────
pdf(file.path(opt$outdir, "cellchat_interaction_heatmap.pdf"), width=9, height=8)
netVisual_heatmap(
  cellchat,
  color.heatmap = "Reds",
  title.name    = "Interaction Strength (PDAC TME)"
)
dev.off()

# ── Figure 3: PDAC-relevant signaling pathways ────────────────────────────
# Pathways of biological interest in PDAC immunosuppression
pdac_pathways <- c("TGFb","CXCL","PD-L1","CSF","MHC-II","VEGF",
                   "IL6","BAFF","PDGF","FGF")
present_pathways <- intersect(pdac_pathways, cellchat@netP$pathways)

for (pw in present_pathways) {
  tryCatch({
    pdf(file.path(opt$outdir, sprintf("pathway_%s.pdf", pw)),
        width=8, height=8)
    netVisual_aggregate(cellchat, signaling=pw,
                        layout="chord",
                        vertex.receiver=c("CD8_T_exhausted","CD4_Treg"),
                        title.name=sprintf("%s Signaling", pw))
    dev.off()
  }, error=function(e) message(sprintf("  Skipped %s: %s", pw, e$message)))
}

# ── Figure 4: Bubble plot — specific L-R pairs ────────────────────────────
# Focus on tumor → immune communication
p_bubble <- netVisual_bubble(
  cellchat,
  sources.use = which(levels(cellchat@idents) %in%
                        c("Tumor_classical","Tumor_basal","myCAF","iCAF")),
  targets.use = which(levels(cellchat@idents) %in%
                        c("CD8_T_exhausted","CD4_Treg","Macrophage_M2")),
  remove.isolate = TRUE,
  max.dataset    = 1
)
ggsave(file.path(opt$outdir, "cellchat_tumor_immune_bubble.pdf"),
       p_bubble, width=12, height=10)

# ── Figure 5: Signaling role analysis ────────────────────────────────────
# Rank cell types by their signaling role:
#   Sender vs Receiver vs Mediator vs Influencer
p_roles <- netAnalysis_signalingRole_scatter(
  cellchat,
  title = "Signaling Role in PDAC TME\n(Sender vs Receiver)"
) +
  labs(caption = "Tumor cells: dominant senders\nExhausted T cells: dominant receivers")
ggsave(file.path(opt$outdir, "cellchat_signaling_roles.pdf"),
       p_roles, width=7, height=6)

# ── Figure 6: PD-1/PD-L1 axis specifically ───────────────────────────────
# The key immunosuppression axis — which cells express PD-L1?
if ("PD-L1" %in% cellchat@netP$pathways) {
  pdf(file.path(opt$outdir, "pdl1_pd1_network.pdf"), width=9, height=9)
  netVisual_aggregate(cellchat, signaling="PD-L1",
                      layout="hierarchy",
                      vertex.receiver=c("CD8_T_exhausted"),
                      title.name="PD-L1 → PD-1 Immunosuppression Axis")
  dev.off()
}

# ── Export full L-R pair table ─────────────────────────────────────────────
lr_pairs <- subsetCommunication(cellchat) %>%
  arrange(desc(prob))

write.csv(lr_pairs, file.path(opt$outdir, "all_lr_pairs.csv"),
          row.names=FALSE)

# Save CellChat object
saveRDS(cellchat, file.path(opt$outdir, "cellchat_pdac.rds"))
message("── CellChat analysis complete ──")
message(sprintf("  %d L-R pairs across %d pathways",
                nrow(lr_pairs), length(cellchat@netP$pathways)))
