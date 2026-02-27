/*
 * ============================================================================
 *  Module: bam_qc.nf
 *  Description: BAM quality control, coverage analysis, and MultiQC
 * ============================================================================
 */

process SAMTOOLS_FLAGSTAT {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    tuple val(sample_id), path(bai)

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
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    tuple val(sample_id), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.stats.txt"), emit: stats

    script:
    """
    samtools stats -@ ${task.cpus} ${bam} > ${sample_id}.stats.txt
    """
}

process SAMTOOLS_IDXSTATS {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    tuple val(sample_id), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.idxstats.txt"), emit: idxstats

    script:
    """
    samtools idxstats ${bam} > ${sample_id}.idxstats.txt
    """
}

process MOSDEPTH {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    tuple val(sample_id), path(bai)
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

process BAM_QC_EVALUATE {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(flagstat)
    tuple val(sample_id), path(stats)
    tuple val(sample_id), path(mosdepth_summary)
    tuple val(sample_id), path(mosdepth_regions)
    tuple val(sample_id), path(on_target_count)

    output:
    tuple val(sample_id), path("${sample_id}.bam_qc.json"), emit: bam_qc_json

    script:
    def threshold_arg = params.bam_qc_thresholds ? "--thresholds ${params.bam_qc_thresholds}" : ""

    // Decompress mosdepth regions if gzipped
    """
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
    path("multiqc_data"),            emit: data

    script:
    """
    multiqc . \\
        --force \\
        --title "Single Gene NIPT Pipeline QC Report" \\
        --comment "Aggregated QC metrics for all samples" \\
        -o .
    """
}
