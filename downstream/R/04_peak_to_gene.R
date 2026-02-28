# downstream/R/04_peak_to_gene.R
#
# Peak-to-Gene Linkage — Cis-Regulatory Element Analysis
# ─────────────────────────────────────────────────────────────
# This analysis links ATAC peaks to the genes they regulate,
# using co-variation across cells in the multiome dataset.
#
# Method (Signac LinkPeaks / ArchR addPeak2GeneLinks):
#   For each gene:
#     1. Identify all peaks within 500kb of the gene TSS
#     2. Compute correlation between peak accessibility and
#        gene expression ACROSS CELLS (using WNN neighbors
#        to smooth the sparse single-cell signal)
#     3. Test significance vs background correlation
#        (permuted peak-gene pairs at same distance)
#
# Why this requires MULTIOME (not separate RNA + ATAC):
#   - Same cell ID links RNA and ATAC measurements directly
#   - No assumption about cell type matching between datasets
#   - Captures cell-state-specific regulatory relationships
#   - Distance-matched background is accurate at single-cell resolution
#
# Biological output for PDAC:
#   - Which enhancers regulate KRAS, MYC, TP53 in tumor cells?
#   - Do CAF subtypes use different enhancers for ACTA2/FAP?
#   - Which peaks near immune checkpoint genes are accessible
#     in exhausted T cells?
# ─────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
  library(optparse)
})
set.seed(42)

opt_list <- list(
  make_option("--multiome_rds", type="character", help="pdac_multiome.rds"),
  make_option("--outdir",       type="character", default="results/peak2gene")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

message("── Loading multiome object ──")
pdac <- readRDS(opt$multiome_rds)
DefaultAssay(pdac) <- "ATAC"

# ── Compute peak-to-gene links ────────────────────────────────────────────
# LinkPeaks tests correlation between peak accessibility and
# gene expression, using distance-matched background peaks
# to control for genomic distance effects
message("── Computing peak-to-gene links (this may take 30-60 mins) ──")

pdac <- LinkPeaks(
  object          = pdac,
  peak.assay      = "ATAC",
  expression.assay= "SCT",
  genes.use       = NULL,      # test all genes
  distance        = 5e5,       # 500kb window around TSS
  min.cells       = 10,
  score_cutoff    = 0.05,
  verbose         = TRUE
)

# Extract link results
links <- Links(pdac)
links_df <- as.data.frame(links) %>%
  filter(!is.na(score)) %>%
  arrange(pvalue)

message(sprintf("  Total significant links: %d", sum(links_df$pvalue < 0.05)))

write.csv(links_df,
          file.path(opt$outdir, "peak_gene_links.csv"),
          row.names=FALSE)

# ── PDAC-relevant gene regulatory landscapes ──────────────────────────────
# Focus on key PDAC oncogenes and tumor suppressor genes
pdac_key_genes <- c(
  # Oncogenes
  "KRAS","MYC","CCND1","CDK6","ERBB2",
  # Stroma / CAF
  "ACTA2","FAP","POSTN","MMP11","IL6",
  # Immune
  "PDCD1","CD274","HAVCR2","CTLA4","TIGIT",
  # Tumor suppressors
  "TP53","CDKN2A","SMAD4","BRCA2"
)

# ── Plot 1: Coverage + links at key loci ──────────────────────────────────
# CoveragePlot shows:
#   Top track: aggregated ATAC signal per cell type at the locus
#   Middle: linked peaks (arcs connecting peak to gene TSS)
#   Bottom: gene model
message("── Plotting coverage tracks ──")

for (gene in c("KRAS","ACTA2","PDCD1","MYC")) {
  tryCatch({
    p_cov <- CoveragePlot(
      object    = pdac,
      region    = gene,
      features  = gene,
      group.by  = "CellType",
      extend.upstream   = 20000,
      extend.downstream = 5000,
      links     = TRUE,       # show peak-gene links as arcs
      peaks     = TRUE        # show called peaks
    )
    ggsave(file.path(opt$outdir, sprintf("coverage_links_%s.pdf", gene)),
           p_cov, width=10, height=12)
  }, error=function(e) message(sprintf("  Skipped %s: %s", gene, e$message)))
}

# ── Plot 2: Links summary — distribution of link distances ────────────────
p_dist <- ggplot(links_df %>% filter(pvalue < 0.05),
                  aes(x=abs(distance)/1000, fill=score > 0)) +
  geom_histogram(bins=50, alpha=0.7) +
  scale_fill_manual(values=c("TRUE"="#D7191C","FALSE"="#2C7BB6"),
                    labels=c("Positive link","Negative link")) +
  labs(title="Distance Distribution of Significant Peak-Gene Links",
       x="Distance from TSS (kb)", y="Number of links",
       fill=NULL) +
  theme_bw(base_size=13)

ggsave(file.path(opt$outdir, "link_distance_histogram.pdf"),
       p_dist, width=7, height=5)

# ── Plot 3: Genes with most regulatory links ───────────────────────────────
# Genes with many linked peaks are likely under complex regulatory control
# In PDAC, MYC and KRAS are known to have many distal enhancers
top_regulated <- links_df %>%
  filter(pvalue < 0.01) %>%
  count(gene, name="n_links") %>%
  slice_max(n_links, n=20) %>%
  mutate(is_pdac_gene = gene %in% pdac_key_genes)

p_top_genes <- ggplot(top_regulated,
  aes(x=reorder(gene, n_links), y=n_links, fill=is_pdac_gene)) +
  geom_col(width=0.7) +
  coord_flip() +
  scale_fill_manual(values=c("TRUE"="#D7191C","FALSE"="grey60")) +
  labs(title="Genes with Most Significant cis-Regulatory Links",
       subtitle="Peak-gene correlation across multiome cells",
       x=NULL, y="Number of linked peaks (p<0.01)",
       fill="PDAC key gene") +
  theme_bw(base_size=12)

ggsave(file.path(opt$outdir, "top_regulated_genes.pdf"),
       p_top_genes, width=7, height=7)

# ── Plot 4: Cell-type-specific links at ACTA2 (CAF marker) ────────────────
# ACTA2 is strongly expressed in myCAF — do its regulatory peaks show
# myCAF-specific accessibility?
acta2_links <- links_df %>%
  filter(gene == "ACTA2", pvalue < 0.05) %>%
  arrange(pvalue)

if (nrow(acta2_links) > 0) {
  # Get accessibility at ACTA2-linked peaks per cell type
  acta2_peaks  <- acta2_links$peak
  acta2_mat    <- GetAssayData(pdac, assay="ATAC")[
                     acta2_peaks[acta2_peaks %in% rownames(pdac)], ]
  acta2_by_ct  <- as.data.frame(t(as.matrix(acta2_mat))) %>%
    mutate(CellType = pdac$CellType) %>%
    group_by(CellType) %>%
    summarise(across(everything(), mean)) %>%
    tidyr::pivot_longer(-CellType, names_to="peak", values_to="accessibility")

  p_acta2 <- ggplot(acta2_by_ct,
    aes(x=CellType, y=peak, fill=accessibility)) +
    geom_tile() +
    scale_fill_viridis_c(option="magma") +
    labs(title="Accessibility at ACTA2-Linked Peaks by Cell Type",
         subtitle="myCAF-specific peaks drive ACTA2 expression in stroma",
         x=NULL, y="Linked peak", fill="Mean\naccessibility") +
    theme_bw(base_size=10) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          axis.text.y=element_text(size=6))

  ggsave(file.path(opt$outdir, "acta2_linked_peaks_heatmap.pdf"),
         p_acta2, width=10, height=6)
}

# Save updated object with links
saveRDS(pdac, file.path(opt$outdir, "pdac_multiome_linked.rds"))
message("── Peak-to-gene linkage complete ──")
