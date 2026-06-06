/*
 * ============================================================================
 *  Module: dark_gene.nf
 *  Description: Dark Gene CNV Analysis for Single Gene NIPT (cfDNA)
 *
 *  Analyses copy-number dosage for pseudogene/paralog-heavy regions:
 *    SMN1/SMN2 (chr5)  → SMA screening
 *    HBA1/HBA2 (chr16) → Alpha thalassemia screening
 *    GBA1/GBAP1 (chr1) → Gaucher disease (dosage only)
 *
 *  Approach (cfDNA-adapted from gx-exome dark gene analysis):
 *    1. samtools bedcov — regional read depth per dark gene window
 *    2. samtools mpileup at SMN1 c.840 discriminating position — SMN1 fraction
 *    3. Autosome background normalization (idxstats)
 *    4. FF-aware dosage interpretation (dark_gene_nipt.py)
 *       observed_ratio = maternal_ratio × (1 − FF) + fetal_ratio × FF
 *
 *  Uses the pre-target-extract dedup BAM (ch_bam_full) so reads outside
 *  the panel target footprint remain visible — necessary for SMN1/2 region
 *  if not fully covered by probes, and consistent with UPD analysis.
 *
 *  Tools required: samtools, bcftools, python3 (all in docker_internal image)
 *  No extra containers needed.
 *
 *  Outputs:
 *    dark_gene/<sample_id>.dark_gene_report.json
 *    dark_gene/<sample_id>.dark_gene_depth.tsv
 * ============================================================================
 */

process DARK_GENE_CNV_ANALYSIS {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/dark_gene", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    tuple val(sample_id), path(ff_json)
    path(dark_gene_bed)
    path(reference_fasta)
    path(reference_fai)

    output:
    tuple val(sample_id), path("${sample_id}.dark_gene_report.json"), emit: dark_gene_json
    tuple val(sample_id), path("${sample_id}.dark_gene_depth.tsv"),   emit: dark_gene_tsv

    script:
    def ref_arg = reference_fasta.name != 'NO_REF' ? "--ref-fasta ${reference_fasta}" : ""
    """
    export HOME=\$PWD

    python3 ${projectDir}/scripts/dark_gene_nipt.py \\
        --bam           ${bam} \\
        --dark-gene-bed ${dark_gene_bed} \\
        --ff-json       ${ff_json} \\
        --sample-id     ${sample_id} \\
        ${ref_arg} \\
        --output-json   ${sample_id}.dark_gene_report.json \\
        --output-tsv    ${sample_id}.dark_gene_depth.tsv
    """
}
