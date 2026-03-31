/*
 * ============================================================================
 *  Module: bam_qc.nf
 *  Description: BAM quality control, coverage analysis, and MultiQC
 *
 *  Panel: Twist UCL_SingleGeneNIPT TE-96276661 (hg38)
 *  - Probe coverage QC using probes_bed
 *  - Zero-probe region flagging using zero_probe_bed
 *  - Panel-tuned thresholds for ~780 kb probe footprint
 * ============================================================================
 */

process SAMTOOLS_FLAGSTAT {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.flagstat.txt"), emit: flagstat

    script:
    """
    samtools flagstat -@ ${task.cpus} ${bam} > ${sample_id}.flagstat.txt
    """
}

process SAMTOOLS_STATS {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.stats.txt"), emit: stats

    script:
    """
    samtools stats -@ ${task.cpus} ${bam} > ${sample_id}.stats.txt
    """
}

process MOSDEPTH {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path(target_bed)

    output:
    tuple val(sample_id), path("${sample_id}.mosdepth.summary.txt"),          emit: summary
    tuple val(sample_id), path("${sample_id}.regions.bed.gz"),                emit: regions
    tuple val(sample_id), path("${sample_id}.mosdepth.global.dist.txt"),      emit: global_dist
    tuple val(sample_id), path("${sample_id}.mosdepth.region.dist.txt"),      emit: region_dist

    script:
    """
    mosdepth \\
        --by ${target_bed} \\
        --thresholds 1,10,20,50,100,200,500,1000 \\
        --mapq ${params.min_mapping_quality} \\
        --no-per-base \\
        -t ${task.cpus} \\
        ${sample_id} \\
        ${bam}
    """
}

process ON_TARGET_COUNT {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path(target_bed)

    output:
    tuple val(sample_id), path("${sample_id}.on_target_count.txt"), emit: on_target_count

    script:
    """
    samtools view -c -L ${target_bed} -q ${params.min_mapping_quality} ${bam} \
        > ${sample_id}.on_target_count.txt
    """
}

process BAM_QC_EVALUATE {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(flagstat)
    tuple val(sample_id), path(stats)
    tuple val(sample_id), path(mosdepth_summary)
    tuple val(sample_id), path(mosdepth_regions)
    tuple val(sample_id), path(on_target_count)
    path(probes_bed)
    path(zero_probe_bed)

    output:
    tuple val(sample_id), path("${sample_id}.bam_qc.json"), emit: bam_qc_json

    script:
    def threshold_arg = params.bam_qc_thresholds ? "--thresholds ${params.bam_qc_thresholds}" : ""
    def probes_arg    = probes_bed.name != 'NO_PROBES_BED' ? "--probes-bed ${probes_bed}" : ""
    def zero_arg      = zero_probe_bed.name != 'NO_ZERO_PROBE_BED' ? "--zero-probe-bed ${zero_probe_bed}" : ""

    """
    export HOME=\$PWD
    # Decompress mosdepth regions BED if needed
    if [[ "${mosdepth_regions}" == *.gz ]]; then
        gunzip -c ${mosdepth_regions} > regions_decompressed.bed
        REGIONS_FILE="regions_decompressed.bed"
    else
        REGIONS_FILE="${mosdepth_regions}"
    fi

    python3 ${projectDir}/scripts/bam_qc.py \\
        --sample-id ${sample_id} \\
        --flagstat ${flagstat} \\
        --stats ${stats} \\
        --mosdepth-summary ${mosdepth_summary} \\
        --mosdepth-regions \${REGIONS_FILE} \\
        --on-target-count ${on_target_count} \\
        ${probes_arg} \\
        ${zero_arg} \\
        ${threshold_arg} \\
        --output ${sample_id}.bam_qc.json
    """
}

process MULTIQC {
    tag "multiqc"
    label 'process_low'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path('*')

    output:
    path("multiqc_report.html"),     emit: report
    path("multiqc_report_data"),     emit: data

    script:
    """
    multiqc . \\
        --force \\
        --filename multiqc_report \\
        --title "Single Gene NIPT Pipeline QC Report" \\
        --comment "Twist UCL_SingleGeneNIPT TE-96276661 (hg38) - Aggregated QC" \\
        -o .
    """
}
