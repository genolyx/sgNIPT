/*
 * ============================================================================
 *  Module: upd_detection.nf
 *  Description: Uniparental Disomy (UPD) detection from maternal cfDNA
 *
 *  Detects UPD on the seven imprinting chromosomes
 *  (chr6, chr7, chr11, chr14, chr15, chr16, chr20) using ONLY the maternal
 *  cfDNA mixture — no paternal sample required.
 *
 *  Strategy
 *  ────────
 *  1. UPD_VARIANT_CALL   : bcftools mpileup + call over an extended UPD BED
 *                          (gene loci + flanking on the 7 UPD chroms +
 *                          background autosomes for an empirical null).
 *                          We omit "-v": ALL sites with AD/DP are kept,
 *                          including mat-hom sites (needed for the patUPD
 *                          Mendelian-error signal).
 *  2. UPD_DETECT         : python upd_detection.py
 *                          Per-site mat-genotype classification + per-chrom
 *                          BAF likelihoods + Mendelian-error excess + final
 *                          per-chrom calls.
 *
 *  See bin/scripts/upd_detection.py for the full signal model.
 * ============================================================================
 */

process UPD_VARIANT_CALL {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/upd", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path(upd_bed)
    path(reference_fasta)
    path(reference_fai)

    output:
    tuple val(sample_id), path("${sample_id}.upd_variants.vcf.gz"),     emit: upd_vcf
    tuple val(sample_id), path("${sample_id}.upd_variants.vcf.gz.tbi"), emit: upd_vcf_idx

    script:
    """
    bcftools mpileup \\
        -f ${reference_fasta} \\
        -R ${upd_bed} \\
        -q ${params.min_mapping_quality} \\
        -Q ${params.min_base_quality} \\
        -d ${params.upd_mpileup_max_depth} \\
        -a FORMAT/AD,FORMAT/DP \\
        ${bam} \\
    | bcftools call \\
        -m \\
        --ploidy 2 \\
        -Ou \\
    | bcftools view \\
        -v snps \\
        -m2 -M2 \\
        -Ou \\
    | bcftools sort \\
        -Oz \\
        -o ${sample_id}.upd_variants.vcf.gz

    bcftools index -t ${sample_id}.upd_variants.vcf.gz
    """
}


process UPD_DETECT {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/upd", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf)
    tuple val(sample_id), path(ff_json)
    path(upd_bed)
    path(common_snps_vcf)

    output:
    tuple val(sample_id), path("${sample_id}.upd_report.json"), emit: upd_report
    tuple val(sample_id), path("${sample_id}.upd_per_snp.tsv"), emit: upd_per_snp

    script:
    def common_arg = (common_snps_vcf.name != 'NO_COMMON_SNPS') \
                     ? "--common-snps ${common_snps_vcf}" : ""
    def cfg_arg    = params.upd_config ? "--config ${params.upd_config}" : ""
    """
    export HOME=\$PWD
    FF=\$(python3 -c "import json; d=json.load(open('${ff_json}')); v=d.get('primary_fetal_fraction'); print(v if v is not None else 0.05)")

    python3 ${projectDir}/scripts/upd_detection.py \\
        --sample-id ${sample_id} \\
        --vcf ${vcf} \\
        --fetal-fraction \${FF} \\
        --upd-bed ${upd_bed} \\
        ${common_arg} \\
        ${cfg_arg} \\
        --output ${sample_id}.upd_report.json \\
        --per-snp-tsv ${sample_id}.upd_per_snp.tsv
    """
}
