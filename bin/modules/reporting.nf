/*
 * ============================================================================
 *  Module: reporting.nf
 *  Description: Final report generation aggregating all pipeline results
 * ============================================================================
 */

process GENERATE_REPORT {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_qc_json)
    tuple val(sample_id), path(bam_qc_json)
    tuple val(sample_id), path(ff_json)
    tuple val(sample_id), path(variant_report_json)
    tuple val(sample_id),      path(annotated_vcf)

    output:
    tuple val(sample_id), path("${sample_id}_final_report.json"), emit: report_json
    tuple val(sample_id), path("${sample_id}_final_report.html"), emit: report_html

    script:
    def ann_arg = (annotated_vcf.name != 'NO_ANNOTATED_VCF') ? "--annotated_vcf ${annotated_vcf}" : ""
    """
    export HOME=\$PWD
    python3 ${projectDir}/scripts/generate_report.py \\
        --sample_id ${sample_id} \\
        --fastq_qc ${fastq_qc_json} \\
        --bam_qc ${bam_qc_json} \\
        --fetal_fraction ${ff_json} \\
        --variant_report ${variant_report_json} \\
        ${ann_arg} \\
        --output_dir . \\
        --output_prefix ${sample_id}
    """
}
