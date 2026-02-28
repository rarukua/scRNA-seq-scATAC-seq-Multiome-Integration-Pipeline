// upstream/modules/cellranger_arc.nf
//
// Cell Ranger ARC — 10x Multiome preprocessing
// Handles simultaneous RNA + ATAC alignment and demultiplexing

process CELLRANGER_ARC {
    tag "${meta.id}"
    publishDir "${params.outdir}/cellranger_arc/${meta.id}", mode: 'copy'

    // Cell Ranger ARC is resource-intensive — needs dedicated node
    cpus   32
    memory '128 GB'
    time   '24h'

    input:
    tuple val(meta), path(rna_fastqs), path(atac_fastqs)
    path  reference   // Cell Ranger ARC genome reference

    output:
    tuple val(meta), path("${meta.id}/outs/filtered_feature_bc_matrix/"),  emit: filtered_matrix
    tuple val(meta), path("${meta.id}/outs/atac_fragments.tsv.gz"),         emit: fragments
    tuple val(meta), path("${meta.id}/outs/atac_fragments.tsv.gz.tbi"),     emit: fragments_idx
    tuple val(meta), path("${meta.id}/outs/per_barcode_metrics.csv"),       emit: per_barcode_qc
    tuple val(meta), path("${meta.id}/outs/atac_peaks.bed"),                emit: peaks
    tuple val(meta), path("${meta.id}/outs/summary.csv"),                   emit: summary
    tuple val(meta), path("${meta.id}/outs/web_summary.html"),              emit: web_summary

    script:
    // Build the libraries CSV that Cell Ranger ARC requires
    // Format: fastqs,sample,library_type
    // library_type must be exactly "Gene Expression" or "Chromatin Accessibility"
    """
    # Create libraries CSV linking FASTQ paths to modality
    cat > libraries.csv << EOF
fastqs,sample,library_type
${rna_fastqs.join(',')},${meta.id}_RNA,Gene Expression
${atac_fastqs.join(',')},${meta.id}_ATAC,Chromatin Accessibility
EOF

    cellranger-arc count \\
        --id            ${meta.id} \\
        --reference     ${reference} \\
        --libraries     libraries.csv \\
        --localcores    ${task.cpus} \\
        --localmem      120 \\
        --description   "PDAC multiome ${meta.id}"

    # Verify critical outputs exist
    [ -f "${meta.id}/outs/atac_fragments.tsv.gz" ] || \\
        { echo "ERROR: fragments file missing"; exit 1; }
    [ -f "${meta.id}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" ] || \\
        { echo "ERROR: matrix barcodes missing"; exit 1; }
    """
}
