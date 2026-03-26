/*
 * ============================================================================
 *  Module: alignment.nf
 *  Description: Read alignment with BWA-MEM2 and BAM post-processing
 * ============================================================================
 */

process BWA_MEM2_ALIGN {
    tag "${sample_id}"
    label 'process_high'
    publishDir "${params.outdir}/${sample_id}/alignment", mode: 'copy', pattern: '*.flagstat'

    input:
    tuple val(sample_id), path(reads)
    path(bwa_index)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"),     emit: sorted_bam
    tuple val(sample_id), path("${sample_id}.sorted.bam.bai"), emit: sorted_bai

    script:
    def ref_prefix = bwa_index.find { it.name.endsWith('.fa') || it.name.endsWith('.fasta') || it.name.endsWith('.fa.gz') }?.name?.replaceAll(/\.(fa|fasta)(\.gz)?$/, '') ?: params.genome_prefix
    def rg = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}_lib"
    // BWA uses index prefix from bwa_index_dir; BWA-MEM2 uses reference_fasta path
    def bwa_idx = "${params.bwa_index_dir}/${file(params.reference_fasta).name}"
    def align_cmd = params.aligner == 'bwa'
        ? "bwa mem -t ${task.cpus} -R '${rg}' ${params.bwa_extra_args ?: ''} ${bwa_idx}"
        : "bwa-mem2 mem -t ${task.cpus} -R '${rg}' ${params.bwa_extra_args ?: ''} ${params.reference_fasta}"

    if (reads instanceof List && reads.size() == 2)
        """
        ${align_cmd} \\
            ${reads[0]} ${reads[1]} \\
        | samtools sort \\
            -@ ${task.cpus} \\
            -m ${params.samtools_sort_mem} \\
            -o ${sample_id}.sorted.bam \\
            -

        samtools index -@ ${task.cpus} ${sample_id}.sorted.bam
        """
    else
        """
        ${align_cmd} \\
            ${reads} \\
        | samtools sort \\
            -@ ${task.cpus} \\
            -m ${params.samtools_sort_mem} \\
            -o ${sample_id}.sorted.bam \\
            -

        samtools index -@ ${task.cpus} ${sample_id}.sorted.bam
        """
}

process MARK_DUPLICATES {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/alignment", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    tuple val(sample_id), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"),          emit: dedup_bam
    tuple val(sample_id), path("${sample_id}.dedup.bam.bai"),      emit: dedup_bai
    tuple val(sample_id), path("${sample_id}.dedup_metrics.txt"),  emit: dedup_metrics

    script:
    if (params.use_picard_markdup)
        """
        picard MarkDuplicates \\
            INPUT=${bam} \\
            OUTPUT=${sample_id}.dedup.bam \\
            METRICS_FILE=${sample_id}.dedup_metrics.txt \\
            REMOVE_DUPLICATES=false \\
            ASSUME_SORTED=true \\
            VALIDATION_STRINGENCY=LENIENT \\
            CREATE_INDEX=true

        mv ${sample_id}.dedup.bai ${sample_id}.dedup.bam.bai || true
        if [ ! -f ${sample_id}.dedup.bam.bai ]; then
            samtools index ${sample_id}.dedup.bam
        fi
        """
    else
        """
        # samtools markdup requires fixmate (MS tag) first
        samtools sort -n -@ ${task.cpus} -o ${sample_id}.nmsrt.bam ${bam}
        samtools fixmate -m ${sample_id}.nmsrt.bam ${sample_id}.fixmate.bam
        samtools sort -@ ${task.cpus} -o ${sample_id}.fixmate_srt.bam ${sample_id}.fixmate.bam

        samtools markdup \\
            -@ ${task.cpus} \\
            --write-index \\
            -s \\
            ${sample_id}.fixmate_srt.bam \\
            ${sample_id}.dedup.bam

        # Capture markdup stats
        samtools markdup -s ${sample_id}.fixmate_srt.bam /dev/null 2> ${sample_id}.dedup_metrics.txt || true

        if [ ! -f ${sample_id}.dedup.bam.bai ]; then
            samtools index -@ ${task.cpus} ${sample_id}.dedup.bam
        fi
        """
}

process EXTRACT_TARGET_READS {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/alignment", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    tuple val(sample_id), path(bai)
    path(target_bed)

    output:
    tuple val(sample_id), path("${sample_id}.target.bam"),     emit: target_bam
    tuple val(sample_id), path("${sample_id}.target.bam.bai"), emit: target_bai
    tuple val(sample_id), path("${sample_id}.on_target_count.txt"), emit: on_target_count

    script:
    """
    # Extract reads overlapping target regions
    samtools view \\
        -@ ${task.cpus} \\
        -b \\
        -L ${target_bed} \\
        -q ${params.min_mapping_quality} \\
        -F 3852 \\
        ${bam} \\
        -o ${sample_id}.target.bam

    samtools index -@ ${task.cpus} ${sample_id}.target.bam

    # Count on-target reads
    samtools view -c ${sample_id}.target.bam > ${sample_id}.on_target_count.txt
    """
}
