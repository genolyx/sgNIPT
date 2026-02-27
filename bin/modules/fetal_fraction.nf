/*
 * ============================================================================
 *  Module: fetal_fraction.nf
 *  Description: Fetal fraction estimation from cfDNA
 * ============================================================================
 */

process COLLECT_FRAGMENT_SIZES {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/fetal_fraction", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    tuple val(sample_id), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.fragment_sizes.tsv"),  emit: fragment_sizes
    tuple val(sample_id), path("${sample_id}.fragment_stats.json"), emit: fragment_stats

    script:
    """
    python3 ${projectDir}/scripts/collect_fragment_sizes.py \\
        --bam ${bam} \\
        --sample-id ${sample_id} \\
        --output-sizes ${sample_id}.fragment_sizes.tsv \\
        --output-stats ${sample_id}.fragment_stats.json \\
        --max-reads ${params.ff_max_reads} \\
        --count-chroms
    """
}

process MPILEUP_AT_SNPS {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/fetal_fraction", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    tuple val(sample_id), path(bai)
    path(target_bed)
    path(reference_fasta)
    path(reference_fai)

    output:
    tuple val(sample_id), path("${sample_id}.target_pileup.txt"), emit: pileup

    script:
    """
    samtools mpileup \\
        -f ${reference_fasta} \\
        -l ${target_bed} \\
        -q ${params.min_mapping_quality} \\
        -Q ${params.min_base_quality} \\
        --max-depth ${params.mpileup_max_depth} \\
        ${bam} \\
        > ${sample_id}.target_pileup.txt
    """
}

process VARIANT_CALL_FOR_FF {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/fetal_fraction", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    tuple val(sample_id), path(bai)
    path(target_bed)
    path(reference_fasta)
    path(reference_fai)

    output:
    tuple val(sample_id), path("${sample_id}.ff_variants.vcf"), emit: ff_vcf

    script:
    if (params.ff_variant_caller == 'freebayes')
        """
        freebayes \\
            -f ${reference_fasta} \\
            -t ${target_bed} \\
            --min-mapping-quality ${params.min_mapping_quality} \\
            --min-base-quality ${params.min_base_quality} \\
            --min-alternate-fraction ${params.ff_min_alt_fraction} \\
            --min-alternate-count ${params.ff_min_alt_count} \\
            --pooled-continuous \\
            --report-monomorphic \\
            ${bam} \\
            > ${sample_id}.ff_variants.raw.vcf

        # Filter for SNPs only
        bcftools view \\
            -v snps \\
            ${sample_id}.ff_variants.raw.vcf \\
            > ${sample_id}.ff_variants.vcf
        """
    else
        """
        # Use bcftools mpileup + call as default
        bcftools mpileup \\
            -f ${reference_fasta} \\
            -R ${target_bed} \\
            -q ${params.min_mapping_quality} \\
            -Q ${params.min_base_quality} \\
            -d ${params.mpileup_max_depth} \\
            -a FORMAT/AD,FORMAT/DP,FORMAT/SP \\
            ${bam} \\
        | bcftools call \\
            -m \\
            -v \\
            --ploidy 2 \\
            -Ov \\
        | bcftools view \\
            -v snps \\
            > ${sample_id}.ff_variants.vcf
        """
}

process ESTIMATE_FETAL_FRACTION {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}/fetal_fraction", mode: 'copy'

    input:
    tuple val(sample_id), path(ff_vcf)
    tuple val(sample_id), path(fragment_stats)
    tuple val(sample_id), path(fragment_sizes)

    output:
    tuple val(sample_id), path("${sample_id}.fetal_fraction.json"),        emit: ff_json
    tuple val(sample_id), path("${sample_id}.fetal_fraction.detail.json"), emit: ff_detail

    script:
    def config_arg = params.ff_config ? "--config ${params.ff_config}" : ""

    // Extract Y and autosome reads from fragment stats
    """
    # Extract chromosome read counts from fragment stats
    Y_READS=\$(python3 -c "import json; d=json.load(open('${fragment_stats}')); print(d.get('y_chromosome_reads', 0))")
    AUTO_READS=\$(python3 -c "import json; d=json.load(open('${fragment_stats}')); print(d.get('autosome_reads', 0))")

    python3 ${projectDir}/scripts/fetal_fraction.py \\
        --sample-id ${sample_id} \\
        --vcf ${ff_vcf} \\
        --y-reads \${Y_READS} \\
        --autosome-reads \${AUTO_READS} \\
        --fragment-sizes ${fragment_sizes} \\
        ${config_arg} \\
        --output ${sample_id}.fetal_fraction.json
    """
}
