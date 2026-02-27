/*
 * ============================================================================
 *  Module: reporting.nf
 *  Description: Final report generation aggregating all pipeline results
 * ============================================================================
 */

process GENERATE_REPORT {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}/report", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_qc_json)
    tuple val(sample_id), path(bam_qc_json)
    tuple val(sample_id), path(ff_json)
    tuple val(sample_id), path(variant_report_json)

    output:
    tuple val(sample_id), path("${sample_id}.final_report.json"), emit: report_json
    tuple val(sample_id), path("${sample_id}.final_report.html"), emit: report_html

    script:
    """
    python3 ${projectDir}/scripts/generate_report.py \\
        --sample-id ${sample_id} \\
        --fastq-qc ${fastq_qc_json} \\
        --bam-qc ${bam_qc_json} \\
        --fetal-fraction ${ff_json} \\
        --variant-report ${variant_report_json} \\
        --output-json ${sample_id}.final_report.json \\
        --output-html ${sample_id}.final_report.html
    """
}
