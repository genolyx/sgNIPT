/*
 * ============================================================================
 *  Module: variant_calling.nf
 *  Description: Variant calling at target loci and fetal variant analysis
 * ============================================================================
 */

process VARIANT_CALL_TARGET {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    tuple val(sample_id), path(bai)
    path(target_bed)
    path(reference_fasta)
    path(reference_fai)

    output:
    tuple val(sample_id), path("${sample_id}.target_variants.vcf.gz"),     emit: target_vcf
    tuple val(sample_id), path("${sample_id}.target_variants.vcf.gz.tbi"), emit: target_vcf_idx

    script:
    if (params.variant_caller == 'freebayes')
        """
        freebayes \\
            -f ${reference_fasta} \\
            -t ${target_bed} \\
            --min-mapping-quality ${params.min_mapping_quality} \\
            --min-base-quality ${params.min_base_quality} \\
            --min-alternate-fraction ${params.vc_min_alt_fraction} \\
            --min-alternate-count ${params.vc_min_alt_count} \\
            --pooled-continuous \\
            ${bam} \\
        | bcftools sort \\
            -Oz \\
            -o ${sample_id}.target_variants.vcf.gz

        bcftools index -t ${sample_id}.target_variants.vcf.gz
        """
    else if (params.variant_caller == 'gatk')
        """
        gatk Mutect2 \\
            -R ${reference_fasta} \\
            -I ${bam} \\
            -L ${target_bed} \\
            --min-base-quality-score ${params.min_base_quality} \\
            --callable-depth ${params.vc_min_depth} \\
            -O ${sample_id}.target_variants.raw.vcf.gz

        # Filter variants
        gatk FilterMutectCalls \\
            -R ${reference_fasta} \\
            -V ${sample_id}.target_variants.raw.vcf.gz \\
            -O ${sample_id}.target_variants.vcf.gz

        bcftools index -t ${sample_id}.target_variants.vcf.gz
        """
    else
        """
        # Default: bcftools mpileup + call
        bcftools mpileup \\
            -f ${reference_fasta} \\
            -R ${target_bed} \\
            -q ${params.min_mapping_quality} \\
            -Q ${params.min_base_quality} \\
            -d ${params.mpileup_max_depth} \\
            -a FORMAT/AD,FORMAT/DP,FORMAT/SP,FORMAT/ADF,FORMAT/ADR \\
            ${bam} \\
        | bcftools call \\
            -m \\
            -v \\
            --ploidy 2 \\
            -Oz \\
            -o ${sample_id}.target_variants.vcf.gz

        bcftools index -t ${sample_id}.target_variants.vcf.gz
        """
}

process FILTER_VARIANTS {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf)
    tuple val(sample_id), path(vcf_idx)

    output:
    tuple val(sample_id), path("${sample_id}.filtered_variants.vcf"),  emit: filtered_vcf

    script:
    """
    bcftools view \\
        -i 'FORMAT/DP >= ${params.vc_min_depth}' \\
        ${vcf} \\
    | bcftools norm \\
        -m -both \\
        -Ov \\
        -o ${sample_id}.filtered_variants.vcf
    """
}

process ANNOTATE_VARIANTS {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}/variants", mode: 'copy'

    when:
    params.run_annotation

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("${sample_id}.annotated_variants.vcf"), emit: annotated_vcf

    script:
    if (params.annotation_tool == 'snpsift' && params.clinvar_vcf)
        """
        SnpSift annotate \\
            ${params.clinvar_vcf} \\
            ${vcf} \\
            > ${sample_id}.annotated_variants.vcf
        """
    else if (params.annotation_tool == 'bcftools' && params.annotation_db)
        """
        bcftools annotate \\
            -a ${params.annotation_db} \\
            -c CHROM,POS,REF,ALT,INFO \\
            ${vcf} \\
            -o ${sample_id}.annotated_variants.vcf
        """
    else
        """
        # No annotation database configured; pass through
        cp ${vcf} ${sample_id}.annotated_variants.vcf
        """
}

process VARIANT_ANALYSIS {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf)
    tuple val(sample_id), path(ff_json)
    path(target_bed)

    output:
    tuple val(sample_id), path("${sample_id}.variant_report.json"), emit: variant_report

    script:
    def config_arg = params.variant_analysis_config ? "--config ${params.variant_analysis_config}" : ""
    """
    # Extract fetal fraction value from JSON
    FF=\$(python3 -c "import json; d=json.load(open('${ff_json}')); print(d.get('primary_fetal_fraction', 0.05))")

    python3 ${projectDir}/scripts/variant_analysis.py \\
        --sample-id ${sample_id} \\
        --vcf ${vcf} \\
        --target-bed ${target_bed} \\
        --fetal-fraction \${FF} \\
        ${config_arg} \\
        --output ${sample_id}.variant_report.json
    """
}
