# WNN (Weighted Nearest Neighbor) — Algorithm Explained

## The Problem WNN Solves

When you have two measurements per cell (RNA + ATAC), you could:

1. **Concatenate** both feature spaces → same weight to both, ignores quality variation
2. **Average** both UMAPs → same problem, loses local structure
3. **Use only RNA** → ignores the epigenome
4. **Use only ATAC** → ignores transcriptome

None of these is optimal. The fundamental insight of WNN is:
**different cells have different data quality in each modality, and the
algorithm should adapt to this.**

---

## The WNN Algorithm (Hao et al. 2021 Cell)

For each cell `i`:

### Step 1: Compute within-modality neighbors
- Find k-nearest neighbors of cell `i` in RNA space → `KNN_RNA(i)`
- Find k-nearest neighbors of cell `i` in ATAC space → `KNN_ATAC(i)`

### Step 2: Compute cross-modality prediction accuracy
- Use `KNN_RNA(i)` to predict cell `i`'s ATAC profile
  → prediction error = `err_predict_ATAC_from_RNA(i)`
- Use `KNN_ATAC(i)` to predict cell `i`'s RNA profile
  → prediction error = `err_predict_RNA_from_ATAC(i)`

### Step 3: Compute per-cell modality weights
```
w_RNA(i)  = 1 / (1 + exp( err_predict_RNA_from_ATAC(i) ))
w_ATAC(i) = 1 - w_RNA(i)
```

A cell where RNA predicts ATAC neighbors well → ATAC is redundant given RNA →
**lower ATAC weight** for that cell.

A cell where ATAC provides unique information not captured by RNA →
**higher ATAC weight** for that cell.

### Step 4: Build WNN graph
The combined KNN graph uses distance:
```
d_WNN(i,j) = w_RNA(i) × d_RNA(i,j) + w_ATAC(i) × d_ATAC(i,j)
```

### Step 5: UMAP + clustering on WNN graph
Same as standard scRNA-seq, but using the weighted graph.

---

## What the Modality Weights Tell You Biologically

When you plot modality weights on the UMAP (as in `03_wnn_integration.R`):

- **High RNA weight regions**: cell types well-defined by gene expression
  (T cells, B cells — strong transcriptional identity)

- **High ATAC weight regions**: cell types where epigenome adds information
  beyond RNA (transitional states, progenitors — chromatin predicts future state)

- **Equal weights**: cells where both modalities contribute equally
  (most mature, stable cell types)

In PDAC, stellate cells transitioning to CAF often show elevated ATAC weight
because the chromatin remodeling PRECEDES gene expression changes — the
epigenome is "ahead" of the transcriptome during cell state transitions.

---

## WNN vs Simple Integration Methods

| Method | Handles quality variation | Per-cell adaptation | Computationally feasible |
|--------|--------------------------|---------------------|--------------------------|
| Concatenation | ✗ | ✗ | ✓ |
| MOFA+ | Partial | ✗ | ✓ |
| Harmony (per modality) | ✓ | ✗ | ✓ |
| **WNN** | **✓** | **✓** | **✓** |
| totalVI | ✓ | ✓ | ✗ (GPU) |

---

## Common Interview Questions About WNN

**Q: "Why dim 2:30 for LSI instead of 1:30?"**
A: LSI dimension 1 is dominated by sequencing depth (like PC1 in bulk RNA-seq
being library size). Including it would give RNA cells with high ATAC depth
artificially high ATAC weight, not reflecting biology. Dimensions 2:30 capture
biological variation.

**Q: "Why SCTransform for the RNA modality?"**
A: Standard log-normalization inflates variance at high-expression genes and
creates false structure when integrating across patients with different sequencing
depths. SCTransform models the technical variance explicitly and regresses it out,
giving cleaner PCA embeddings for WNN.

**Q: "How do you validate the WNN clustering is better?"**
A: Three approaches:
1. Cell type marker gene specificity improves (silhouette score on known markers)
2. WNN UMAP separates known co-localizing populations that single-modality fails to resolve
3. Cluster-level DA between WNN clusters is more statistically clean (fewer mixed cells)
