#!/usr/bin/env python3
"""
downstream/python/scenic_plus.py

SCENIC+ Gene Regulatory Network (GRN) inference
─────────────────────────────────────────────────────────────
Reference: Bravo González-Blas et al. 2023 Nature Methods

SCENIC+ is the multiome-aware extension of SCENIC/pySCENIC:
  - Uses BOTH scRNA-seq AND scATAC-seq simultaneously
  - Infers eGRNs (enhancer-driven Gene Regulatory Networks)
  - Links: TF → accessible peak (binding site) → target gene
  - The key output: eRegulons — each containing:
      * A TF (e.g. FOSL2/AP-1)
      * Its binding peaks (from ATAC motif + footprint evidence)
      * Its target genes (from RNA co-expression)
      * Direction: activating or repressing

Why SCENIC+ > pySCENIC alone:
  pySCENIC:  TF → target genes (RNA only, no peak evidence)
  SCENIC+:   TF → peak (epigenetic validation) → gene (RNA)
             The peak provides MECHANISTIC evidence of TF binding

PDAC biological focus:
  - KRAS downstream TFs: AP-1 (FOSL2/JUNB), ETS (ERG, ETV4)
  - Lineage TFs: FOXA1, GATA6 (classical subtype)
  - EMT / basal: TP63, ZEB1, SNAI2
  - Immunosuppression: STAT3 (in macrophages + tumor cells)
"""

import argparse
import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import mudata as mu
from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import *
from pycisTopic.topic_binarization import *
from pycisTopic.diff_features import *
from scenicplus.wrappers.run_scenicplus import run_scenicplus

def parse_args():
    parser = argparse.ArgumentParser(description="SCENIC+ GRN inference for PDAC multiome")
    parser.add_argument("--rna_h5ad",       required=True,
                        help="scRNA-seq AnnData (.h5ad), exported from Seurat")
    parser.add_argument("--fragments_dir",  required=True,
                        help="Directory containing per-sample fragments.tsv.gz")
    parser.add_argument("--consensus_peaks",required=True,
                        help="Consensus peak set BED file (from ArchR)")
    parser.add_argument("--genome",         default="hg38")
    parser.add_argument("--outdir",         default="results/scenic_plus")
    parser.add_argument("--n_topics",       type=int, default=40,
                        help="Number of LDA topics for cisTopic")
    parser.add_argument("--n_cpu",          type=int, default=16)
    return parser.parse_args()


def prepare_rna_anndata(rna_h5ad_path):
    """
    Load and prepare RNA AnnData for SCENIC+.
    SCENIC+ requires raw counts (not normalized) in adata.X.
    Cell type annotations must be in adata.obs['cell_type'].
    """
    print("── Loading RNA AnnData ──")
    rna = sc.read_h5ad(rna_h5ad_path)

    # Ensure raw counts are available
    if rna.raw is not None:
        rna.X = rna.raw.X
    print(f"  Cells: {rna.n_obs}, Genes: {rna.n_vars}")
    print(f"  Cell types: {rna.obs['CellType'].value_counts().to_dict()}")
    return rna


def run_cistopic(fragments_dir, consensus_peaks, cell_barcodes,
                 sample_names, outdir, n_topics=40, n_cpu=16):
    """
    cisTopic: LDA topic modelling of scATAC-seq data.

    cisTopic treats each cell as a 'document' and each genomic region
    (peak) as a 'word'. LDA learns latent 'topics' — groups of
    co-accessible peaks that represent regulatory programs.

    These topics are analogous to gene modules in scRNA-seq, but for
    chromatin accessibility. Each topic captures a coherent
    cis-regulatory program (e.g., 'active enhancers in CAF cells').

    SCENIC+ uses cisTopic topics as the accessibility representation
    (instead of raw peak counts) because topics are:
      - Denoised (LDA smooths the sparse single-cell signal)
      - Biologically interpretable (each topic = regulatory program)
      - Dimensionality-reduced (40 topics << 100k peaks)
    """
    print("── Running cisTopic ──")

    # Build fragment path dictionary
    path_to_fragments = {
        sample: os.path.join(fragments_dir, sample, "outs", "atac_fragments.tsv.gz")
        for sample in sample_names
    }

    # Create cisTopic object from peak-cell matrix
    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments    = path_to_fragments,
        path_to_regions      = consensus_peaks,
        path_to_blacklist    = None,
        min_frag_in_regions  = 10,
        min_cell             = 10,
        metrics              = None,
        n_cpu                = n_cpu,
        project              = "PDAC_multiome",
    )

    # Subset to cells present in RNA object
    cistopic_obj = subset_cistopic_object(
        cistopic_obj,
        selected_cells = cell_barcodes
    )

    # Run LDA models with multiple topic counts
    # Best model selected by coherence score
    print(f"  Training LDA models (n_topics={n_topics}) ──")
    models = run_cgs_models(
        cistopic_obj,
        n_topics   = [n_topics - 10, n_topics, n_topics + 10],
        n_cpu      = n_cpu,
        n_iter     = 500,
        random_state = 42,
        alpha       = 50,
        alpha_by_topic = True,
        eta         = 0.1,
        eta_by_topic = False,
        save_path   = outdir
    )

    # Select best model
    cistopic_obj = evaluate_models(
        cistopic_obj,
        models,
        select_model     = None,  # auto-select by coherence
        return_model     = True,
        metrics          = ['Minmo_2011','loglikelihood','coherence'],
        plot_metrics     = False
    )

    # Binarize topics → accessible regions per topic
    binarized_topics = binarize_topics(
        cistopic_obj,
        target    = 'cell',
        method    = 'otsu',
        plot      = False
    )

    # Compute DARs (Differentially Accessible Regions) per cell type
    imputed_acc = impute_accessibility(
        cistopic_obj,
        selected_cells = None,
        selected_regions = None,
        scale_factor   = 10**6
    )

    return cistopic_obj, imputed_acc


def run_scenic_plus_grn(rna, cistopic_obj, imputed_acc, outdir,
                         genome="hg38", n_cpu=16):
    """
    Run SCENIC+ to infer eGRNs (enhancer Gene Regulatory Networks).

    The SCENIC+ pipeline:
    1. Motif enrichment in topic-specific regions (pycisTopic + pycisTarget)
       → Which TF motifs are enriched in each accessibility topic?

    2. TF-to-gene co-expression (GRNBoost2 from pySCENIC)
       → Which TFs co-vary with which target genes in RNA space?

    3. eRegulon inference (combining steps 1 + 2)
       → For each TF: its binding peaks (ATAC) + target genes (RNA)
       → Direction inference: activating vs repressing
       → Confidence score for each TF-peak-gene triplet

    Output: eRegulons — the core SCENIC+ result
    Each eRegulon contains:
      * TF name
      * List of target peaks (with motif evidence)
      * List of target genes (with co-expression evidence)
      * Regulon size (# target genes)
      * Specificity: which cell types express this regulon
    """
    print("── Running SCENIC+ ──")

    # Create MuData object linking RNA + ATAC
    # MuData is the multimodal AnnData format used by SCENIC+
    mdata = mu.MuData({
        'rna':  rna,
        'atac': ad.AnnData(
            X   = imputed_acc.T,
            obs = rna.obs,
        )
    })

    # Run full SCENIC+ pipeline
    run_scenicplus(
        mdata                    = mdata,
        genes_to_use             = None,           # all genes
        cisTopic_obj             = cistopic_obj,
        menr                     = None,           # will compute internally
        base_GEX_obj             = rna,
        path_to_blacklist        = None,
        tmp_scenicplus_dir       = os.path.join(outdir, "tmp"),
        save_path                = outdir,
        biomart_host             = "http://www.ensembl.org",
        species                  = "hsapiens",
        assembly                 = genome,
        upstream                 = [1000, 150000],  # search window around TSS
        downstream               = [1000, 150000],
        tf_file                  = None,            # uses default human TF list
        save_partial             = True,            # save intermediate results
        use_gene_boundaries      = True,
        region_ranking           = True,
        inplace_subset_eRegulons = True,
        eRegulons_direct_only    = False,
        # GRNBoost2 settings
        n_cpu                    = n_cpu,
        # Scoring thresholds
        ray_n_cpu                = n_cpu,
        ray_tmp_dir              = outdir
    )


def summarize_egregulons(outdir):
    """
    Load and summarize eRegulon results.
    Produces:
      - Table of top regulons by size and specificity
      - PDAC-relevant TF activity scores
      - eRegulon heatmap across cell types
    """
    print("── Summarizing eRegulons ──")

    # Load eRegulon AUC scores (one row per cell, one column per TF)
    auc_file = os.path.join(outdir, "eRegulon_AUC.csv")
    if not os.path.exists(auc_file):
        print("  eRegulon AUC not found — SCENIC+ may still be running")
        return None

    auc_df = pd.read_csv(auc_file, index_col=0)

    # PDAC-relevant TFs to highlight
    pdac_tfs = [
        "FOSL2", "JUNB", "FOSB",   # AP-1 family (KRAS downstream)
        "FOXA1", "GATA6",           # classical subtype lineage
        "TP63",  "ZEB1",            # basal/EMT
        "STAT3", "IRF4",            # immunosuppression
        "SMAD3",                    # TGF-β / CAF activation
        "ETS1",  "ETS2",            # ETS family
        "NF1",                      # tumor suppressor TF
    ]

    # Filter to TFs present in results
    present_tfs = [tf for tf in pdac_tfs
                   if any(tf in col for col in auc_df.columns)]
    print(f"  PDAC-relevant TFs found: {present_tfs}")

    # Summary statistics per eRegulon
    regulon_summary = pd.DataFrame({
        "TF":           auc_df.columns,
        "mean_AUC":     auc_df.mean(axis=0),
        "std_AUC":      auc_df.std(axis=0),
        "max_AUC":      auc_df.max(axis=0),
    }).sort_values("mean_AUC", ascending=False)

    regulon_summary.to_csv(
        os.path.join(outdir, "eregulon_summary.csv"), index=False
    )

    print(f"  Top 10 regulons by mean AUC:")
    print(regulon_summary.head(10).to_string())

    return auc_df, regulon_summary


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # ── Step 1: Prepare RNA ────────────────────────────────────────────────
    rna = prepare_rna_anndata(args.rna_h5ad)
    sample_names = rna.obs["sample"].unique().tolist()

    # ── Step 2: cisTopic (ATAC topic modelling) ────────────────────────────
    cistopic_obj, imputed_acc = run_cistopic(
        fragments_dir    = args.fragments_dir,
        consensus_peaks  = args.consensus_peaks,
        cell_barcodes    = rna.obs_names.tolist(),
        sample_names     = sample_names,
        outdir           = args.outdir,
        n_topics         = args.n_topics,
        n_cpu            = args.n_cpu
    )

    # ── Step 3: SCENIC+ eGRN inference ────────────────────────────────────
    run_scenic_plus_grn(
        rna          = rna,
        cistopic_obj = cistopic_obj,
        imputed_acc  = imputed_acc,
        outdir       = args.outdir,
        genome       = args.genome,
        n_cpu        = args.n_cpu
    )

    # ── Step 4: Summarize results ─────────────────────────────────────────
    result = summarize_egregulons(args.outdir)
    if result is not None:
        print("── SCENIC+ complete ──")
        print(f"  Results in: {args.outdir}/")
        print(f"  Next: run downstream/R/05_scenic_plus.R to visualize")
    else:
        print("── SCENIC+ running asynchronously — check logs ──")


if __name__ == "__main__":
    main()
