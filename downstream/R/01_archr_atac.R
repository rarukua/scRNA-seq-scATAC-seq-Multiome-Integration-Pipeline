# downstream/R/01_archr_atac.R
#
# ArchR: scATAC-seq processing from Multiome fragments
# ─────────────────────────────────────────────────────────────
# In multiome data, the ATAC component is processed first with ArchR
# to get a clean LSI embedding and peak set. This is then handed to
# Signac/Seurat for WNN integration with the RNA component.
#
# Key multiome-specific considerations vs standalone scATAC:
#   - Cell barcodes are SHARED between RNA and ATAC — use the
#     SAME barcode whitelist for both to ensure 1:1 cell matching
#   - TSS enrichment thresholds can be slightly relaxed because
#     the RNA modality provides a second quality filter
#   - Peak calling benefits from all cells (not just per-cluster)
#     because multiome gives better cell type annotation upfront

suppressPackageStartupMessages({
  library(ArchR)
  library(parallel)
  library(optparse)
})

set.seed(42)
addArchRThreads(threads=detectCores()-2)
addArchRGenome("hg38")

opt_list <- list(
  make_option("--fragments_dir", type="character"),
  make_option("--sample_names",  type="character", help="Comma-separated sample names"),
  make_option("--outdir",        type="character", default="results/archr"),
  make_option("--min_tss",       type="double",    default=3.0,
              help="Min TSS enrichment (relaxed for multiome — RNA provides 2nd filter)"),
  make_option("--min_frags",     type="integer",   default=1000)
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
setwd(opt$outdir)

sample_names <- strsplit(opt$sample_names, ",")[[1]]

# ── Step 1: Create Arrow files ────────────────────────────────────────────
# One Arrow file per sample — disk-based HDF5 storage
message("── Creating Arrow files ──")

fragment_files <- file.path(opt$fragments_dir,
                             paste0(sample_names, "/outs/atac_fragments.tsv.gz"))
names(fragment_files) <- sample_names

arrow_files <- createArrowFiles(
  inputFiles      = fragment_files,
  sampleNames     = sample_names,
  minTSS          = opt$min_tss,
  minFrags        = opt$min_frags,
  addTileMat      = TRUE,
  addGeneScoreMat = TRUE,
  # Multiome-specific: force barcode match to 10x whitelist
  # (ensures same barcodes as RNA matrix)
  force           = TRUE
)

# ── Step 2: Doublet detection ─────────────────────────────────────────────
doublet_scores <- addDoubletScores(
  input     = arrow_files,
  k         = 10,
  knnMethod = "UMAP",
  LSIMethod = 1
)

# ── Step 3: Create ArchR project ──────────────────────────────────────────
proj <- ArchRProject(
  ArrowFiles      = arrow_files,
  outputDirectory = "ArchR_PDAC",
  copyArrows      = TRUE
)
proj <- filterDoublets(proj)

message(sprintf("Cells after doublet filtering: %d", nCells(proj)))

# ── Step 4: LSI dimensionality reduction ──────────────────────────────────
# For multiome, 2-iteration LSI on the TileMatrix is standard.
# Iteration 1 removes depth-driven variation
# Iteration 2 captures biological signal
proj <- addIterativeLSI(
  ArchRProj    = proj,
  useMatrix    = "TileMatrix",
  name         = "IterativeLSI",
  iterations   = 2,
  clusterParams = list(
    resolution  = c(0.2),
    sampleCells = 10000,
    n.start     = 10
  ),
  varFeatures  = 25000,
  dims         = 1:30
)

# ── Step 5: Harmony batch correction ──────────────────────────────────────
# Multiple PDAC patients → must correct for patient batch effects
# Harmony operates on LSI components, preserving biological variation
# while removing patient-specific technical/genetic variation
proj <- addHarmony(
  ArchRProj    = proj,
  reducedDims  = "IterativeLSI",
  name         = "Harmony",
  groupBy      = "Sample"   # correct by patient/sample
)

# ── Step 6: Clustering ────────────────────────────────────────────────────
proj <- addClusters(
  input       = proj,
  reducedDims = "Harmony",
  method      = "Seurat",
  name        = "ATAC_Clusters",
  resolution  = 0.6
)

# ── Step 7: UMAP ──────────────────────────────────────────────────────────
proj <- addUMAP(
  ArchRProj   = proj,
  reducedDims = "Harmony",
  name        = "UMAP_ATAC",
  nNeighbors  = 30,
  minDist     = 0.5,
  metric      = "cosine"
)

# ── Step 8: Peak calling per cluster (pseudobulk MACS2) ───────────────────
# Pseudobulk peak calling: aggregate all cells in each cluster
# before calling peaks — more sensitive than per-cell calling
# and more biologically meaningful than joint peak calling
proj <- addGroupCoverages(
  ArchRProj     = proj,
  groupBy       = "ATAC_Clusters",
  minCells      = 50,
  maxCells      = 500
)

proj <- addReproduciblePeakSet(
  ArchRProj     = proj,
  groupBy       = "ATAC_Clusters",
  pathToMacs2   = findMacs2(),
  reproducibility = "2",     # peak must be reproducible in ≥2 pseudobulks
  peaksPerCell  = 500
)

proj <- addPeakMatrix(proj)

# ── Step 9: chromVAR TF deviation scores ──────────────────────────────────
# Compute per-cell TF activity from motif accessibility
# PDAC-relevant motifs: KRAS downstream TFs (AP-1, ETS), lineage TFs
proj <- addMotifAnnotations(
  ArchRProj = proj,
  motifSet  = "JASPAR2022",
  name      = "Motif"
)

proj <- addBgdPeaks(proj)

proj <- addDeviationsMatrix(
  ArchRProj      = proj,
  peakAnnotation = "Motif",
  force          = TRUE
)

# ── Step 10: Gene score matrix QC ─────────────────────────────────────────
# PDAC cell type markers for annotation reference
pdac_markers <- list(
  Tumor_cells    = c("KRT19","KRT8","EPCAM","MUC5AC","CA19-9","KRAS"),
  Ductal_normal  = c("KRT19","CFTR","SLC4A4"),
  Stellate_cells = c("ACTA2","FAP","COL1A1","PDGFRA","THY1"),
  CAF            = c("ACTA2","FAP","POSTN","PDPN","MMP11"),
  Macrophage     = c("CD68","LYZ","MARCO","MRC1","CD163"),
  T_cell         = c("CD3D","CD8A","CD4","FOXP3","PDCD1"),
  B_cell         = c("CD79A","MS4A1","CD19"),
  Endothelial    = c("PECAM1","VWF","CDH5"),
  NK_cell        = c("NCAM1","KLRB1","NKG7","GNLY")
)

p_markers_atac <- plotEmbedding(
  ArchRProj    = proj,
  colorBy      = "GeneScoreMatrix",
  name         = unlist(pdac_markers),
  embedding    = "UMAP_ATAC",
  quantCut     = c(0.01, 0.95)
)
plotPDF(p_markers_atac,
        name     = "UMAP_ATAC_marker_gene_scores",
        width    = 5, height = 5,
        addDOC   = FALSE,
        ArchRProj = proj)

# ── Save project for WNN integration ──────────────────────────────────────
# Export key objects needed by Signac/Seurat WNN step:
#   1. LSI matrix (for WNN weight computation)
#   2. Peak matrix (for Signac)
#   3. Cell metadata (cluster assignments, QC metrics)
message("── Exporting for WNN integration ──")

# LSI matrix for Seurat
lsi_matrix <- getReducedDims(proj, reduceDims="Harmony")
saveRDS(lsi_matrix,
        file.path(opt$outdir, "atac_lsi_harmony.rds"))

# Peak matrix as Seurat-compatible sparse matrix
peak_matrix <- getMatrixFromProject(proj, useMatrix="PeakMatrix")
saveRDS(peak_matrix,
        file.path(opt$outdir, "atac_peak_matrix.rds"))

# Cell metadata
cell_meta <- as.data.frame(getCellColData(proj))
write.csv(cell_meta,
          file.path(opt$outdir, "atac_cell_metadata.csv"))

# TF deviation matrix
tf_devs <- getMatrixFromProject(proj, useMatrix="MotifMatrix")
saveRDS(tf_devs,
        file.path(opt$outdir, "atac_tf_deviations.rds"))

saveArchRProject(ArchRProj=proj,
                 outputDirectory="ArchR_PDAC",
                 load=FALSE)

message("── ArchR ATAC processing complete ──")
message("  Exported: atac_lsi_harmony.rds, atac_peak_matrix.rds")
message("  Next: run 02_seurat_rna.R, then 03_wnn_integration.R")
