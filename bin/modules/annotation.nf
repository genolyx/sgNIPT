/*
 * ============================================================================
 *  Module: annotation.nf
 *  Description: VEP annotation for sgNIPT variant calls
 *
 *  Mirrors carrier-screening/modules/annotation.nf so both pipelines
 *  share the same VEP invocation pattern and CSQ field set.
 *
 *  VEP fields written to INFO/CSQ:
 *    Gene, SYMBOL, HGVSc, HGVSp, Consequence, IMPACT,
 *    gnomADe_AF, gnomADg_AF, dbSNP (Existing_variation),
 *    SIFT, PolyPhen, MANE_SELECT, CANONICAL
 *
 *  Requires:
 *    - params.vep_cache_dir  : Pre-downloaded Ensembl VEP cache (release 111+)
 *    - params.reference_fasta: GRCh38 FASTA (same as alignment step)
 *    - Docker image          : ensemblorg/ensembl-vep:release_111
 * ============================================================================
 */

process VEP_ANNOTATION {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/variants", mode: 'copy'

    // Use the official Ensembl VEP Docker image
    container 'ensemblorg/ensembl-vep:release_111'

    input:
    tuple val(sample_id), path(vcf)
    path(reference_fasta)
    path(reference_fai)

    output:
    tuple val(sample_id), path("${sample_id}.vep.vcf.gz"),     emit: vep_vcf
    tuple val(sample_id), path("${sample_id}.vep.vcf.gz.tbi"), emit: vep_vcf_idx
    tuple val(sample_id), path("${sample_id}.vep_summary.html"), emit: vep_summary, optional: true

    script:
    def cache_dir = params.vep_cache_dir ?: '/opt/vep/.vep'
    def assembly  = 'GRCh38'
    def species   = 'homo_sapiens'
    // Input may be .vcf or .vcf.gz — VEP handles both
    def input_vcf = vcf.name.endsWith('.gz') ? vcf : vcf
    """
    vep \\
        --input_file  ${input_vcf} \\
        --output_file ${sample_id}.vep.vcf \\
        --format vcf \\
        --vcf \\
        --species ${species} \\
        --assembly ${assembly} \\
        --cache \\
        --dir_cache ${cache_dir} \\
        --offline \\
        --fasta ${reference_fasta} \\
        --everything \\
        --canonical \\
        --mane \\
        --hgvs \\
        --hgvsg \\
        --sift b \\
        --polyphen b \\
        --af_gnomade \\
        --af_gnomadg \\
        --check_existing \\
        --symbol \\
        --gene_phenotype \\
        --biotype \\
        --numbers \\
        --total_length \\
        --no_progress \\
        --fork ${task.cpus} \\
        --stats_file ${sample_id}.vep_summary.html \\
        --warning_file ${sample_id}.vep_warnings.txt

    # Compress and index
    bgzip -@ ${task.cpus} ${sample_id}.vep.vcf
    tabix -p vcf ${sample_id}.vep.vcf.gz
    """
}
