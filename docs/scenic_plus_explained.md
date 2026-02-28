# SCENIC+ — eGRN Inference Explained

## What SCENIC+ Infers

SCENIC+ outputs **eRegulons** — enhancer-driven Gene Regulatory Networks.

Each eRegulon is a tuple: **(TF, binding peaks, target genes)**

```
Example eRegulon: FOSL2 (AP-1)
  ├── TF: FOSL2 (expressed in tumor cells, high chromVAR score)
  ├── Binding peaks: [chr1:123000-123200, chr7:456000-456300, ...]
  │     (peaks containing AP-1 motif + footprint evidence from ATAC)
  └── Target genes: [MMP7, VEGFA, MYC, CCND1, ...]
        (genes co-expressed with FOSL2 in RNA space)

Interpretation: FOSL2 binds these enhancer peaks and activates these genes.
This is the mechanism by which KRAS → ERK → AP-1 drives oncogenesis.
```

---

## pySCENIC vs SCENIC+

| Feature | pySCENIC | SCENIC+ |
|---------|----------|---------|
| Input | scRNA-seq only | scRNA-seq + scATAC-seq |
| TF-target evidence | Co-expression only | Co-expression + motif + footprint |
| Peak evidence | None | Each regulon has ATAC peaks as mechanism |
| Directionality | Activating regulons | Activating + repressing |
| False positive rate | Higher | Lower (requires epigenetic evidence) |
| Requires multiome | No | Yes (or matched separate RNA+ATAC) |

---

## The Three Pillars of SCENIC+ Evidence

### 1. cisTopic (chromatin accessibility programs)
LDA topic modelling on the peak matrix. Each topic = a group of co-accessible
peaks = a regulatory program. Topics are denoised (not sparse like raw peaks)
and biologically interpretable.

### 2. pyCisTarget (motif enrichment in topics)
For each cisTopic topic: which TF motifs are enriched in the topic's peaks?
Uses the cisTarget databases (500Mb of pre-computed motif scores).
Output: TF → accessible region assignments.

### 3. GRNBoost2 (TF-gene co-expression)
Gradient boosted trees identify which TFs best predict each gene's expression
across cells. Fast, handles the non-linear relationships in scRNA-seq.
Output: TF → target gene weights.

### Integration
eRegulons are inferred where a TF has:
- Motif evidence in accessible peaks (steps 1+2)
- Co-expression with target genes (step 3)
- The peaks are within 500kb of the target gene TSS

---

## SCENIC+ in the Context of This Pipeline

The unique value of running SCENIC+ on PDAC multiome data:

```
KRAS mutation
     ↓
ERK/MAPK activation
     ↓
AP-1 TF family (FOSL2, JUNB) — phosphorylated, nuclear
     ↓ (inferred from eRegulon + ATAC footprint at AP-1 motifs)
FOSL2 eRegulon peaks — enhancers opened by AP-1 binding
     ↓ (peak-to-gene links from 04_peak_to_gene.R)
Target genes: MMP7, VEGFA, MYC, CCND1... — oncogenic program
     ↓
Tumor proliferation, invasion, immune evasion
```

Without ATAC (pySCENIC only): you see FOSL2 co-expressed with MYC, but
you cannot tell if FOSL2 is directly regulating MYC or if both are
downstream of a third factor.

With SCENIC+ (multiome): FOSL2 footprints are detected in the MYC enhancer
region at chr8q24. The epigenetic evidence makes the regulatory claim
mechanistically grounded.

---

## AUCell Scores

SCENIC+ outputs per-cell AUC scores for each eRegulon.
AUC = Area Under the Recovery Curve = how enriched are the eRegulon's
target genes in the top-expressed genes of that cell.

High AUC = the cell is using this regulatory program.
Plotting AUC on UMAP shows which cell types express which regulons.

The AUC score is the same as in pySCENIC — it's the way SCENIC
quantifies "is this regulon active in this cell?"

---

## Common Interview Questions

**Q: "How is SCENIC+ different from just running motif enrichment?"**
A: Motif enrichment (HOMER/TOBIAS) tells you which motifs are over-represented
in accessible peaks. SCENIC+ additionally requires that the TF bearing that
motif is co-expressed with target genes near those peaks. This adds two
additional layers of evidence — expression context and genomic proximity —
dramatically reducing false positives.

**Q: "What's the limitation of SCENIC+?"**
A: Three main ones:
1. Requires high-quality ATAC data (>2000 fragments/cell recommended)
2. cisTopic training time (~2-4 hours on 30k cells, 16 CPUs)
3. The motif database may miss non-canonical or context-specific TF binding
   (solved partially by using TOBIAS footprints as additional evidence)

**Q: "Why cisTopic rather than just using the peak matrix directly?"**
A: The peak matrix is extremely sparse (most peaks have 0 or 1 reads per cell).
LDA denoising through cisTopic aggregates signal across co-accessible peaks
into topics, making the per-cell chromatin state estimate much more robust.
It's analogous to why scRNA-seq tools use PCA/embedding rather than raw counts.
