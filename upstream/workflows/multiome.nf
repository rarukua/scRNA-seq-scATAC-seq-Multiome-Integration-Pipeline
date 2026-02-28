// upstream/workflows/multiome.nf
//
// 10x Genomics Multiome upstream pipeline
// ─────────────────────────────────────────────────────────────
// 10x Multiome ATAC + Gene Expression generates THREE FASTQ sets
// per sample from a single library:
//
//   RNA FASTQ:  R1=barcode+UMI (28bp), R2=cDNA (90bp), I1=sample index
//   ATAC FASTQ: R1=read1 (50bp), R2=barcode (24bp), R3=read2 (49bp)
//
// Cell Ranger ARC processes BOTH simultaneously:
//   - Aligns RNA reads → gene expression matrix
//   - Aligns ATAC reads with chromatin-aware aligner → fragments.tsv.gz
//   - Links cell barcodes between modalities (same whitelist)
//   - Calls peaks jointly (better power than single modality)
//   - Outputs per-barcode QC for both RNA and ATAC
//
// WHY Nextflow is appropriate here:
//   - Cell Ranger ARC is computationally intensive (32+ CPUs, 64GB RAM)
//   - Multiple patient samples must be run independently in parallel
//   - Standardized outputs feed into ALL downstream tools
//   - No biological decisions made at this stage — pure preprocessing
// ─────────────────────────────────────────────────────────────

nextflow.enable.dsl=2

include { CELLRANGER_ARC } from '../modules/cellranger_arc'
include { FASTQC         } from '../modules/fastqc'
include { MULTIQC        } from '../modules/multiqc'

workflow MULTIOME {
    take:
    ch_samples  // [meta, rna_fastqs, atac_fastqs]

    main:

    // ── Step 1: Raw FASTQ QC ──────────────────────────────────────────────
    // Run FastQC on both RNA and ATAC FASTQ sets separately
    ch_all_fastqs = ch_samples
        .flatMap { meta, rna_fqs, atac_fqs ->
            rna_fqs.collect  { fq -> [[id: "${meta.id}_RNA",  type:"rna"],  fq] } +
            atac_fqs.collect { fq -> [[id: "${meta.id}_ATAC", type:"atac"], fq] }
        }
    FASTQC(ch_all_fastqs)

    // ── Step 2: Cell Ranger ARC ───────────────────────────────────────────
    // Single command processes RNA + ATAC together.
    // Key outputs:
    //   filtered_feature_bc_matrix/    RNA count matrix (gene × cell)
    //   atac_fragments.tsv.gz          ATAC fragment file (cell × position)
    //   per_barcode_metrics.csv        Per-cell QC: RNA UMIs, ATAC frags, TSS
    //   summary.csv                    Run-level QC statistics
    //   atac_peaks.bed                 Peaks called jointly on all cells
    //   atac_peak_annotation.tsv       Peak annotation (nearest gene, feature)
    CELLRANGER_ARC(
        ch_samples,
        params.cellranger_arc_ref  // pre-built reference for hg38
    )

    // ── MultiQC ──────────────────────────────────────────────────────────
    ch_qc = Channel.empty()
        .mix(FASTQC.out.zip)
        .mix(CELLRANGER_ARC.out.summary)
        .collect()
    MULTIQC(ch_qc)

    emit:
    // All outputs needed by downstream R/Python tools
    filtered_matrix  = CELLRANGER_ARC.out.filtered_matrix  // → Seurat
    fragments        = CELLRANGER_ARC.out.fragments         // → ArchR + Signac
    per_barcode_qc   = CELLRANGER_ARC.out.per_barcode_qc   // → QC filtering
    peaks            = CELLRANGER_ARC.out.peaks             // → reference peak set
    summary          = CELLRANGER_ARC.out.summary           // → MultiQC
    multiqc_html     = MULTIQC.out.report
}
