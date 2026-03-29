/*
 * ============================================================================
 *  Module: variant_calling.nf
 *  Description: Variant calling at target loci and fetal variant analysis
 *
 *  Panel: Twist UCL_SingleGeneNIPT TE-96276661 (hg38)
 *  - Uses All_target_regions BED for variant calling
 *  - Zero-probe region flagging for 29 uncovered target regions
 *  - Gene list validation against singleGene_NIPT gene list
 *
 *  Caller: bcftools mpileup + call (fixed)
 *    cfDNA NIPT requires detection of fetal-specific variants at VAF ≈ FF/2
 *    (typically 3–10%). Germline callers (GATK HC, DeepVariant, Strelka2)
 *    treat sub-20% VAF as noise and filter them out. bcftools reports allele
 *    counts without diploid genotype assumptions, preserving the fetal signal.
 * ============================================================================
 */

process VARIANT_CALL_TARGET {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path(target_bed)
    path(reference_fasta)
    path(reference_fai)

    output:
    tuple val(sample_id), path("${sample_id}.target_variants.vcf.gz"),     emit: target_vcf
    tuple val(sample_id), path("${sample_id}.target_variants.vcf.gz.tbi"), emit: target_vcf_idx

    script:
    """
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


process VARIANT_ANALYSIS {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf)
    tuple val(sample_id), path(ff_json)
    path(target_bed)
    path(zero_probe_bed)
    path(gene_list)

    output:
    tuple val(sample_id), path("${sample_id}.variant_report.json"), emit: variant_report
    tuple val(sample_id), path("${sample_id}.fetal_variants.vcf"),  emit: fetal_vcf

    script:
    def config_arg    = params.variant_analysis_config ? "--config ${params.variant_analysis_config}" : ""
    def zero_arg      = zero_probe_bed.name != 'NO_ZERO_PROBE_BED' ? "--zero-probe-bed ${zero_probe_bed}" : ""
    def gene_list_arg = gene_list.name != 'NO_GENE_LIST' ? "--gene-list ${gene_list}" : ""

    """
    export HOME=\$PWD
    # Extract fetal fraction value from JSON (handle None when FF estimation failed)
    FF=\$(python3 -c "import json; d=json.load(open('${ff_json}')); v=d.get('primary_fetal_fraction'); print(v if v is not None else 0.05)")

    python3 ${projectDir}/scripts/variant_analysis.py \\
        --sample-id ${sample_id} \\
        --vcf ${vcf} \\
        --target-bed ${target_bed} \\
        --fetal-fraction \${FF} \\
        ${zero_arg} \\
        ${gene_list_arg} \\
        ${config_arg} \\
        --output ${sample_id}.variant_report.json \\
        --output-vcf ${sample_id}.fetal_variants.vcf
    """
}
