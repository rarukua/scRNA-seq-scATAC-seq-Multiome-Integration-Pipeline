# scRNA-seq + scATAC-seq Multiome Integration Pipeline
## Pancreatic Ductal Adenocarcinoma (PDAC)

> **Portfolio project** | 10x Genomics Multiome (simultaneous RNA + ATAC from same cell)  
> Framework stack: ArchR · Signac · Seurat WNN · SCENIC+ · Monocle3 · CellChat  
> Cross-modal integration with bulk ATAC (LUAD pipeline) and WGBS methylation (CRC pipeline)

---

## Why Multiome?

Standard single-cell experiments measure either RNA **or** chromatin accessibility.
10x Genomics Multiome measures **both simultaneously from the same cell**, enabling:

```
Single cell nucleus
        │
        ├── RNA library  ──► gene expression (what the cell is doing)
        └── ATAC library ──► chromatin accessibility (what the cell CAN do)

The same cell ID links both measurements → no integration uncertainty
```

This is fundamentally more powerful than integrating separate RNA and ATAC
experiments because:
- No batch effects between modalities (same cell, same timepoint)
- Direct peak-to-gene linkage within individual cells
- TF regulons can be validated: is the TF expressed AND is its target accessible?
- Epigenetic priming visible: accessible promoters with no RNA yet = poised state

---

## Biological Question — PDAC

Pancreatic ductal adenocarcinoma has a uniquely complex tumor microenvironment (TME):
- **Cancer cells**: KRAS-driven, highly plastic epigenome
- **Cancer-associated fibroblasts (CAF)**: activated from pancreatic stellate cells
- **Macrophages**: polarized toward immunosuppressive M2 state
- **T cells**: exhausted, excluded from tumor core
- **Ductal cells**: normal epithelium adjacent to tumor

This pipeline asks:

1. What **chromatin accessibility programs** distinguish PDAC tumor cells from normal ductal epithelium?
2. Which **transcription factor regulons** (TF + target genes + accessible peaks) drive each cell state?
3. How do pancreatic stellate cells **transition to CAFs** — what epigenetic changes accompany this?
4. Which **cis-regulatory elements** (enhancers/promoters) control key PDAC oncogenes?
5. What **signals** do tumor cells send to suppress T cell function in the TME?
6. How does the PDAC chromatin landscape compare to the CRC methylation landscape at shared oncogenic loci?

---

## Analytical Architecture

```
┌──────────────────────────────────────────────────────────────────────┐
│  UPSTREAM — Nextflow                                                   │
│  10x Multiome FASTQ → Cell Ranger ARC → fragments.tsv.gz + matrix     │
│  Per-sample, parallelized, standardized outputs                        │
└───────────────────────────────┬──────────────────────────────────────┘
                                │
              ┌─────────────────┼──────────────────┐
              ▼                 ▼                  ▼
        ArchR (ATAC)      Seurat (RNA)        Both together
        LSI + peaks       SCTransform         WNN joint embedding
        chromVAR TF       clustering          unified cell types
              │                 │                  │
              └─────────────────┼──────────────────┘
                                ▼
              ┌──────────────────────────────────────┐
              │  Downstream R / Python               │
              │                                      │
              │  01 ArchR QC + LSI + chromVAR        │
              │  02 Seurat RNA preprocessing         │
              │  03 WNN joint embedding (Signac)     │
              │  04 Peak-to-gene linkage             │
              │  05 SCENIC+ TF regulon inference     │
              │  06 Monocle3 trajectory (CAF axis)   │
              │  07 CellChat communication           │
              │  08 Cross-modal integration          │
              └──────────────────────────────────────┘
```

---

## Repository Structure

```
multiome-pdac-pipeline/
│
├── README.md
│
├── upstream/
│   ├── workflows/
│   │   └── multiome.nf              Cell Ranger ARC → ArchR inputs
│   └── modules/
│       ├── cellranger_arc.nf        10x Multiome preprocessing
│       ├── fastqc.nf
│       └── multiqc.nf
│
├── downstream/
│   ├── R/
│   │   ├── 01_archr_atac.R          ArchR: QC, LSI, peaks, chromVAR
│   │   ├── 02_seurat_rna.R          Seurat: SCTransform, clustering
│   │   ├── 03_wnn_integration.R     Signac WNN: joint embedding
│   │   ├── 04_peak_to_gene.R        Cis-regulatory element linkage
│   │   ├── 05_scenic_plus.R         R wrapper for SCENIC+ output
│   │   ├── 06_trajectory.R          Monocle3: stellate → CAF axis
│   │   ├── 07_cellchat.R            Cell-cell communication
│   │   └── 08_cross_modal.R         Link to bulk ATAC + methylation
│   ├── python/
│   │   └── scenic_plus.py           SCENIC+ GRN inference (Python)
│   └── notebooks/
│       ├── Multiome_ATAC_Report.Rmd
│       ├── Multiome_RNA_Report.Rmd
│       ├── WNN_Integration_Report.Rmd
│       └── Cross_Modal_Report.Rmd
│
├── docs/
│   ├── multiome_vs_separate.md      Why multiome > integrated separate assays
│   ├── wnn_explained.md             WNN algorithm walkthrough
│   ├── scenic_plus_explained.md     GRN inference method
│   └── tool_rationale.md
│
├── samplesheets/
│   └── multiome_samples.csv
│
└── assets/
    ├── blacklist_hg38.bed
    └── pdac_marker_genes.csv
```

---

## Dataset

| Dataset | Source | Cells | Description |
|---------|--------|-------|-------------|
| PDAC 10x Multiome | GEO GSE131886 | ~30,000 |  |
| Bulk ATAC (cross-modal) | This repo: atac-luad-pipeline | — | Used for cross-cancer chromatin comparison |
| WGBS cfDNA (cross-modal) | This repo: cfdna-multimodal-pipeline | — | CRC methylation at shared loci |

---

## Quickstart

```bash
# Step 1: Upstream (Nextflow)
nextflow run upstream/workflows/multiome.nf \
  --input samplesheets/multiome_samples.csv \
  --genome hg38 \
  --outdir results/upstream \
  -profile docker

# Step 2: ArchR ATAC processing
Rscript downstream/R/01_archr_atac.R \
  --fragments_dir results/upstream/fragments \
  --outdir results/archr

# Step 3: Seurat RNA processing
Rscript downstream/R/02_seurat_rna.R \
  --matrix_dir results/upstream/filtered_feature_bc_matrix \
  --outdir results/seurat

# Step 4: WNN integration
Rscript downstream/R/03_wnn_integration.R \
  --archr_dir results/archr \
  --seurat_rds results/seurat/pdac_rna.rds \
  --outdir results/wnn

# Step 5: SCENIC+ (Python)
python downstream/python/scenic_plus.py \
  --multiome_dir results/wnn \
  --genome hg38 \
  --outdir results/scenic_plus

# Steps 6-8: Trajectory, CellChat, cross-modal
Rscript downstream/R/06_trajectory.R --wnn_rds results/wnn/pdac_multiome.rds
Rscript downstream/R/07_cellchat.R   --wnn_rds results/wnn/pdac_multiome.rds
Rscript downstream/R/08_cross_modal.R
```

---

## References

1. Raghavan S et al. (2021) Microenvironment drives cell state, plasticity and drug response in PDAC. *Cell* 
2. Granja JM et al. (2021) ArchR. *Nat Genet* — scATAC framework
3. Hao Y et al. (2021) Seurat v4/v5. *Cell* — RNA framework + WNN
4. Stuart T et al. (2021) Signac. *Nat Methods* — scATAC + WNN in R
5. Bravo González-Blas C et al. (2023) SCENIC+. *Nat Methods* — TF regulon GRN
6. Trapnell C et al. (2014) Monocle. *Nat Biotechnol* — Trajectory analysis
7. Jin S et al. (2021) CellChat. *Nat Commun* — Cell-cell communication
