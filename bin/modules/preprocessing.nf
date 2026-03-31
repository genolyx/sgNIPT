/*
 * ============================================================================
 *  Module: preprocessing.nf
 *  Description: FASTQ preprocessing with fastp (adapter trimming, QC filtering)
 * ============================================================================
 */

process FASTP {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R{1,2}.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}.fastp.json"),              emit: fastp_json
    tuple val(sample_id), path("${sample_id}.fastp.html"),              emit: fastp_html

    script:
    if (reads instanceof List && reads.size() == 2)
        """
        fastp \\
            --in1 ${reads[0]} \\
            --in2 ${reads[1]} \\
            --out1 ${sample_id}_trimmed_R1.fastq.gz \\
            --out2 ${sample_id}_trimmed_R2.fastq.gz \\
            --json ${sample_id}.fastp.json \\
            --html ${sample_id}.fastp.html \\
            --thread ${task.cpus} \\
            --qualified_quality_phred ${params.fastp_qualified_quality} \\
            --length_required ${params.fastp_min_length} \\
            --detect_adapter_for_pe \\
            --correction \\
            --overrepresentation_analysis \\
            ${params.fastp_extra_args ?: ''}
        """
    else
        """
        fastp \\
            --in1 ${reads} \\
            --out1 ${sample_id}_trimmed_R1.fastq.gz \\
            --json ${sample_id}.fastp.json \\
            --html ${sample_id}.fastp.html \\
            --thread ${task.cpus} \\
            --qualified_quality_phred ${params.fastp_qualified_quality} \\
            --length_required ${params.fastp_min_length} \\
            ${params.fastp_extra_args ?: ''}
        """
}

process FASTQ_QC {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(fastp_json)

    output:
    tuple val(sample_id), path("${sample_id}.fastq_qc.json"), emit: fastq_qc_json

    script:
    def threshold_arg = params.fastq_qc_thresholds ? "--thresholds ${params.fastq_qc_thresholds}" : ""
    """
    export HOME=\$PWD
    python3 ${projectDir}/scripts/fastq_qc.py \\
        --fastp-json ${fastp_json} \\
        --sample-id ${sample_id} \\
        ${threshold_arg} \\
        --output ${sample_id}.fastq_qc.json
    """
}

// Used only in BAM input mode (--input_bam): generates a placeholder FASTQ QC JSON
// so that GENERATE_REPORT receives a consistent input tuple even without FASTQ QC data.
process MOCK_FASTQ_QC {
    tag "${sample_id}"
    label 'process_low'

    input:
    val(sample_id)

    output:
    tuple val(sample_id), path("${sample_id}.fastq_qc.json"), emit: fastq_qc_json

    script:
    """
    python3 -c "
import json
json.dump({
    'sample_id': '${sample_id}',
    'source': 'bam_input_mode',
    'note': 'FASTQ QC not available (BAM input mode — alignment skipped)',
    'qc_passed': None
}, open('${sample_id}.fastq_qc.json', 'w'), indent=2)
"
    """
}
