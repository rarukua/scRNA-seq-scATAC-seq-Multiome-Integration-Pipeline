# downstream/R/02_seurat_rna.R
#
# Seurat: scRNA-seq processing from Multiome gene expression matrix
# ─────────────────────────────────────────────────────────────
# The RNA component of multiome data is processed exactly like
# standard scRNA-seq — the key difference is that the BARCODE
# whitelist must match the ATAC whitelist for WNN integration.
#
# Cell Ranger ARC guarantees this — it uses the same barcode
# space for both modalities. So every cell in the RNA matrix
# has a corresponding entry in the ATAC fragments file.
#
# Processing strategy:
#   - SCTransform normalization (preferred over log-normalization
#     for integration across multiple patients with varying depths)
#   - Harmony batch correction (patient effects in RNA)
#   - Standard UMAP + clustering
#   - Cell type annotation using known PDAC markers
#   - Export Seurat object for WNN integration

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(optparse)
})

set.seed(42)

opt_list <- list(
  make_option("--matrix_dirs", type="character",
              help="Comma-separated paths to filtered_feature_bc_matrix/ dirs"),
  make_option("--sample_names",type="character"),
  make_option("--outdir",      type="character", default="results/seurat"),
  make_option("--min_genes",   type="integer",   default=500),
  make_option("--max_pct_mt",  type="double",    default=20.0),
  make_option("--min_cells",   type="integer",   default=3)
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

matrix_dirs  <- strsplit(opt$matrix_dirs,  ",")[[1]]
sample_names <- strsplit(opt$sample_names, ",")[[1]]

# ── Step 1: Load and merge all samples ────────────────────────────────────
message("── Loading RNA matrices ──")

seurat_list <- mapply(function(dir, sname) {
  counts <- Read10X(data.dir=dir)
  # Multiome matrix has two assays: "Gene Expression" and "Peaks"
  # Keep only Gene Expression for RNA processing
  if (is.list(counts)) counts <- counts[["Gene Expression"]]

  so <- CreateSeuratObject(
    counts     = counts,
    project    = sname,
    min.cells  = opt$min_cells,
    min.features = opt$min_genes
  )
  so$sample <- sname
  return(so)
}, matrix_dirs, sample_names, SIMPLIFY=FALSE)

# Merge all samples
pdac <- merge(seurat_list[[1]],
              y          = seurat_list[-1],
              add.cell.ids = sample_names,
              project    = "PDAC_Multiome")

message(sprintf("Total cells pre-QC: %d", ncol(pdac)))

# ── Step 2: QC filtering ──────────────────────────────────────────────────
pdac[["pct_mt"]] <- PercentageFeatureSet(pdac, pattern="^MT-")
pdac[["pct_ribo"]]<- PercentageFeatureSet(pdac, pattern="^RP[SL]")

# QC violin plots
p_qc <- VlnPlot(pdac,
                features = c("nFeature_RNA","nCount_RNA","pct_mt"),
                group.by = "sample",
                pt.size  = 0,
                ncol     = 3)
ggsave(file.path(opt$outdir, "qc_violin.pdf"), p_qc, width=14, height=5)

pdac <- subset(pdac,
               subset = nFeature_RNA  > opt$min_genes &
                        nFeature_RNA  < 8000          &  # doublet proxy
                        pct_mt        < opt$max_pct_mt &
                        nCount_RNA    > 1000)

message(sprintf("Cells post-QC: %d", ncol(pdac)))

# ── Step 3: SCTransform normalization ─────────────────────────────────────
# SCTransform is preferred over log-normalization for multi-sample
# multiome data because:
#   - Regresses out sequencing depth variation
#   - Models negative binomial distribution (more accurate for scRNA)
#   - Better handles cross-sample technical variation
#   - vars.to.regress: remove cell cycle and mitochondrial effects
pdac <- SCTransform(pdac,
                    vars.to.regress  = c("pct_mt","nCount_RNA"),
                    verbose          = TRUE,
                    return.only.var.genes = FALSE)

# ── Step 4: PCA ───────────────────────────────────────────────────────────
pdac <- RunPCA(pdac, npcs=50, verbose=FALSE)

# Elbow plot to choose optimal PC number
ElbowPlot(pdac, ndims=50)
ggsave(file.path(opt$outdir, "pca_elbowplot.pdf"), width=6, height=4)

# ── Step 5: Harmony batch correction ──────────────────────────────────────
# PDAC multiome from 8 patients — patient effects dominate without correction
# Harmony corrects PCA embeddings by iteratively projecting cells into
# a shared space, then updating cluster centroids
pdac <- RunHarmony(
  object     = pdac,
  group.by.vars = "sample",
  reduction  = "pca",
  dims.use   = 1:30,
  assay.use  = "SCT",
  verbose    = TRUE
)

# ── Step 6: UMAP + clustering ─────────────────────────────────────────────
pdac <- RunUMAP(pdac,
                reduction  = "harmony",
                dims       = 1:30,
                reduction.name = "UMAP_RNA")

pdac <- FindNeighbors(pdac, reduction="harmony", dims=1:30)
pdac <- FindClusters(pdac,  resolution=0.6,
                     algorithm=1)   # Louvain

p_umap_rna <- DimPlot(pdac, reduction="UMAP_RNA",
                       group.by=c("seurat_clusters","sample"),
                       label=TRUE, ncol=2)
ggsave(file.path(opt$outdir, "umap_rna_clusters.pdf"), p_umap_rna,
       width=14, height=6)

# ── Step 7: Cell type annotation ──────────────────────────────────────────
# PDAC marker genes from Raghavan et al. 2021 Cell and Moffitt et al. 2019
pdac_markers <- list(
  Tumor_classical   = c("KRT19","KRT8","EPCAM","TFF1","MUC5AC"),
  Tumor_basal       = c("KRT5","KRT6A","TP63","S100A2"),
  Ductal_normal     = c("KRT19","CFTR","SLC4A4","SPP1"),
  iCAF              = c("IL6","CXCL12","HAS2","CFD","PDGFRA"),    # inflammatory CAF
  myCAF             = c("ACTA2","MMP11","POSTN","HOPX","TAGLN"),  # myofibroblastic CAF
  apCAF             = c("SLPI","HLA-DRA","CD74","PDPN"),          # antigen-presenting CAF
  Macrophage_M1     = c("IL1B","CXCL10","CD80","TNF"),
  Macrophage_M2     = c("MRC1","CD163","MARCO","IL10"),
  CD8_T_cell        = c("CD8A","CD8B","GZMB","PRF1"),
  CD4_T_cell        = c("CD4","IL7R","TCF7"),
  Treg              = c("FOXP3","IL2RA","CTLA4"),
  Exhausted_T       = c("PDCD1","HAVCR2","LAG3","TIGIT"),
  B_cell            = c("CD79A","MS4A1","IGHM"),
  NK_cell           = c("NCAM1","KLRB1","NKG7"),
  Endothelial       = c("PECAM1","VWF","CDH5","PLVAP"),
  Acinar            = c("PRSS1","CPA1","AMY2A")
)

# Compute module scores for each cell type
for (ct in names(pdac_markers)) {
  genes_present <- intersect(pdac_markers[[ct]], rownames(pdac))
  if (length(genes_present) < 2) next
  pdac <- AddModuleScore(pdac,
                         features = list(genes_present),
                         name     = paste0(ct, "_score"))
}

# Feature plots for key markers
p_features <- FeaturePlot(pdac,
                           features  = c("KRT19","ACTA2","CD68","CD8A","FOXP3","PDCD1"),
                           reduction = "UMAP_RNA",
                           ncol      = 3)
ggsave(file.path(opt$outdir, "umap_rna_key_markers.pdf"),
       p_features, width=12, height=8)

# ── Step 8: Find cluster markers ──────────────────────────────────────────
message("── Finding cluster markers ──")
Idents(pdac) <- "seurat_clusters"
cluster_markers <- FindAllMarkers(pdac,
                                   only.pos  = TRUE,
                                   min.pct   = 0.25,
                                   logfc.threshold = 0.5,
                                   test.use  = "wilcox")

write.csv(cluster_markers,
          file.path(opt$outdir, "cluster_markers.csv"),
          row.names=FALSE)

# ── Save Seurat object for WNN integration ─────────────────────────────────
saveRDS(pdac, file.path(opt$outdir, "pdac_rna.rds"))
message("── RNA processing complete ──")
message("  Saved: pdac_rna.rds")
message("  Next: 03_wnn_integration.R")
