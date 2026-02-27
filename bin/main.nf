#!/usr/bin/env nextflow

/*
 * ============================================================================
 *  Single Gene NIPT Pipeline
 * ============================================================================
 *  Non-Invasive Prenatal Testing pipeline for single gene disorder screening
 *  using cell-free DNA (cfDNA) from maternal plasma.
 *
 *  Pipeline flow:
 *    FASTQ → Preprocessing → Alignment → BAM QC → Fetal Fraction Estimation
 *    → Variant Calling → Variant Analysis → Report Generation
 *
 *  Author  : Single Gene NIPT Team
 *  Version : 1.0.0
 * ============================================================================
 */

nextflow.enable.dsl = 2

// ── Print pipeline header ───────────────────────────────────────────────────
log.info """
╔══════════════════════════════════════════════════════════════════╗
║           Single Gene NIPT Pipeline v${params.pipeline_version}                  ║
║  Non-Invasive Prenatal Testing for Single Gene Disorders        ║
╚══════════════════════════════════════════════════════════════════╝

  Input samplesheet : ${params.input}
  Target BED        : ${params.target_bed}
  Reference         : ${params.reference_fasta}
  Output directory  : ${params.outdir}
  Variant caller    : ${params.variant_caller}
  FF variant caller : ${params.ff_variant_caller}
  ─────────────────────────────────────────────────────────────────
"""

// ── Include modules ─────────────────────────────────────────────────────────
include { FASTP; FASTQ_QC }                                            from './modules/preprocessing'
include { BWA_MEM2_ALIGN; MARK_DUPLICATES; EXTRACT_TARGET_READS }      from './modules/alignment'
include { SAMTOOLS_FLAGSTAT; SAMTOOLS_STATS; SAMTOOLS_IDXSTATS;
          MOSDEPTH; BAM_QC_EVALUATE; MULTIQC }                         from './modules/bam_qc'
include { COLLECT_FRAGMENT_SIZES; MPILEUP_AT_SNPS;
          VARIANT_CALL_FOR_FF; ESTIMATE_FETAL_FRACTION }               from './modules/fetal_fraction'
include { VARIANT_CALL_TARGET; FILTER_VARIANTS;
          ANNOTATE_VARIANTS; VARIANT_ANALYSIS }                        from './modules/variant_calling'
include { GENERATE_REPORT }                                            from './modules/reporting'


// ── Helper: Parse samplesheet ───────────────────────────────────────────────
def parse_samplesheet(samplesheet_path) {
    Channel
        .fromPath(samplesheet_path)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def sample_id = row.sample_id
            def reads = []

            if (row.fastq_2 && row.fastq_2.trim()) {
                reads = [file(row.fastq_1), file(row.fastq_2)]
            } else {
                reads = file(row.fastq_1)
            }

            return tuple(sample_id, reads)
        }
}


// ── Main workflow ───────────────────────────────────────────────────────────
workflow {

    // ── Input channels ──────────────────────────────────────────────────
    ch_reads       = parse_samplesheet(params.input)
    ch_target_bed  = Channel.fromPath(params.target_bed)
    ch_ref_fasta   = Channel.fromPath(params.reference_fasta)
    ch_ref_fai     = Channel.fromPath("${params.reference_fasta}.fai")
    ch_bwa_index   = Channel.fromPath("${params.bwa_index_dir}/*").collect()

    // ══════════════════════════════════════════════════════════════════════
    // STEP 1: FASTQ Preprocessing
    // ══════════════════════════════════════════════════════════════════════
    FASTP(ch_reads)

    FASTQ_QC(FASTP.out.fastp_json)

    // ══════════════════════════════════════════════════════════════════════
    // STEP 2: Alignment
    // ══════════════════════════════════════════════════════════════════════
    BWA_MEM2_ALIGN(
        FASTP.out.trimmed_reads,
        ch_bwa_index,
    )

    // ══════════════════════════════════════════════════════════════════════
    // STEP 3: Mark Duplicates
    // ══════════════════════════════════════════════════════════════════════
    MARK_DUPLICATES(
        BWA_MEM2_ALIGN.out.sorted_bam,
        BWA_MEM2_ALIGN.out.sorted_bai,
    )

    // ══════════════════════════════════════════════════════════════════════
    // STEP 4: Extract Target Region Reads
    // ══════════════════════════════════════════════════════════════════════
    EXTRACT_TARGET_READS(
        MARK_DUPLICATES.out.dedup_bam,
        MARK_DUPLICATES.out.dedup_bai,
        ch_target_bed.first(),
    )

    // ══════════════════════════════════════════════════════════════════════
    // STEP 5: BAM Quality Control
    // ══════════════════════════════════════════════════════════════════════
    SAMTOOLS_FLAGSTAT(
        MARK_DUPLICATES.out.dedup_bam,
        MARK_DUPLICATES.out.dedup_bai,
    )

    SAMTOOLS_STATS(
        MARK_DUPLICATES.out.dedup_bam,
        MARK_DUPLICATES.out.dedup_bai,
    )

    SAMTOOLS_IDXSTATS(
        MARK_DUPLICATES.out.dedup_bam,
        MARK_DUPLICATES.out.dedup_bai,
    )

    MOSDEPTH(
        MARK_DUPLICATES.out.dedup_bam,
        MARK_DUPLICATES.out.dedup_bai,
        ch_target_bed.first(),
    )

    BAM_QC_EVALUATE(
        SAMTOOLS_FLAGSTAT.out.flagstat,
        SAMTOOLS_STATS.out.stats,
        MOSDEPTH.out.summary,
        MOSDEPTH.out.regions,
        EXTRACT_TARGET_READS.out.on_target_count,
    )

    // ══════════════════════════════════════════════════════════════════════
    // STEP 6: Fetal Fraction Estimation
    // ══════════════════════════════════════════════════════════════════════

    // 6a. Collect fragment size distribution (for supplementary FF estimation)
    COLLECT_FRAGMENT_SIZES(
        MARK_DUPLICATES.out.dedup_bam,
        MARK_DUPLICATES.out.dedup_bai,
    )

    // 6b. Variant calling at target SNP sites for FF estimation
    VARIANT_CALL_FOR_FF(
        MARK_DUPLICATES.out.dedup_bam,
        MARK_DUPLICATES.out.dedup_bai,
        ch_target_bed.first(),
        ch_ref_fasta.first(),
        ch_ref_fai.first(),
    )

    // 6c. Estimate fetal fraction
    ESTIMATE_FETAL_FRACTION(
        VARIANT_CALL_FOR_FF.out.ff_vcf,
        COLLECT_FRAGMENT_SIZES.out.fragment_stats,
        COLLECT_FRAGMENT_SIZES.out.fragment_sizes,
    )

    // ══════════════════════════════════════════════════════════════════════
    // STEP 7: Variant Calling at Target Loci
    // ══════════════════════════════════════════════════════════════════════
    VARIANT_CALL_TARGET(
        EXTRACT_TARGET_READS.out.target_bam,
        EXTRACT_TARGET_READS.out.target_bai,
        ch_target_bed.first(),
        ch_ref_fasta.first(),
        ch_ref_fai.first(),
    )

    FILTER_VARIANTS(
        VARIANT_CALL_TARGET.out.target_vcf,
        VARIANT_CALL_TARGET.out.target_vcf_idx,
    )

    // Optional annotation
    if (params.run_annotation) {
        ANNOTATE_VARIANTS(FILTER_VARIANTS.out.filtered_vcf)
        ch_analysis_vcf = ANNOTATE_VARIANTS.out.annotated_vcf
    } else {
        ch_analysis_vcf = FILTER_VARIANTS.out.filtered_vcf
    }

    // ══════════════════════════════════════════════════════════════════════
    // STEP 8: Variant Analysis (Fetal genotype determination)
    // ══════════════════════════════════════════════════════════════════════
    VARIANT_ANALYSIS(
        ch_analysis_vcf,
        ESTIMATE_FETAL_FRACTION.out.ff_json,
        ch_target_bed.first(),
    )

    // ══════════════════════════════════════════════════════════════════════
    // STEP 9: Report Generation
    // ══════════════════════════════════════════════════════════════════════
    GENERATE_REPORT(
        FASTQ_QC.out.fastq_qc_json,
        BAM_QC_EVALUATE.out.bam_qc_json,
        ESTIMATE_FETAL_FRACTION.out.ff_json,
        VARIANT_ANALYSIS.out.variant_report,
    )

    // ══════════════════════════════════════════════════════════════════════
    // STEP 10: MultiQC (aggregate all QC outputs)
    // ══════════════════════════════════════════════════════════════════════
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files
        .mix(FASTP.out.fastp_json.map { it[1] })
        .mix(SAMTOOLS_FLAGSTAT.out.flagstat.map { it[1] })
        .mix(SAMTOOLS_STATS.out.stats.map { it[1] })
        .mix(SAMTOOLS_IDXSTATS.out.idxstats.map { it[1] })
        .mix(MOSDEPTH.out.global_dist.map { it[1] })
        .mix(MOSDEPTH.out.region_dist.map { it[1] })
        .collect()

    MULTIQC(ch_multiqc_files)
}


// ── Workflow completion handler ─────────────────────────────────────────────
workflow.onComplete {
    log.info """
    ══════════════════════════════════════════════════════════════════
    Pipeline execution summary
    ══════════════════════════════════════════════════════════════════
      Completed at : ${workflow.complete}
      Duration     : ${workflow.duration}
      Success      : ${workflow.success}
      Exit status  : ${workflow.exitStatus}
      Work dir     : ${workflow.workDir}
      Output dir   : ${params.outdir}
    ══════════════════════════════════════════════════════════════════
    """.stripIndent()
}

workflow.onError {
    log.error """
    ══════════════════════════════════════════════════════════════════
    Pipeline execution FAILED
    ══════════════════════════════════════════════════════════════════
      Error message: ${workflow.errorMessage}
    ══════════════════════════════════════════════════════════════════
    """.stripIndent()
}
