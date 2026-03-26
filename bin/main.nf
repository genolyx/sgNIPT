#!/usr/bin/env nextflow

/*
 * ============================================================================
 *  Single Gene NIPT Pipeline v1.1.0
 * ============================================================================
 *  Non-Invasive Prenatal Testing pipeline for single gene disorder screening
 *  using cell-free DNA (cfDNA) from maternal plasma.
 *
 *  Twist Custom Panel: UCL_SingleGeneNIPT_TE-96276661_hg38
 *    - 5,138 target regions (~1.07 Mb), 183 genes, 22 chromosomes
 *    - No Y-chromosome targets → ChrX-based supplementary FF
 *    - 29 zero-probe regions flagged in QC
 *
 *  Pipeline flow:
 *    FASTQ → Preprocessing → Alignment → BAM QC → Fetal Fraction Estimation
 *    → Variant Calling → Variant Analysis → Report Generation
 *
 *  Author  : Single Gene NIPT Team
 *  Version : 1.1.0
 * ============================================================================
 */

nextflow.enable.dsl = 2

// ── Print pipeline header ───────────────────────────────────────────────────
log.info """
╔══════════════════════════════════════════════════════════════════╗
║           Single Gene NIPT Pipeline v${params.pipeline_version}                  ║
║  Non-Invasive Prenatal Testing for Single Gene Disorders        ║
║  Panel: Twist UCL_SingleGeneNIPT_TE-96276661 (hg38)            ║
╚══════════════════════════════════════════════════════════════════╝

  Input samplesheet : ${params.input}
  Target BED        : ${params.target_bed}
  Probes BED        : ${params.probes_bed ?: 'not provided'}
  FF SNPs BED       : ${params.ff_snps_bed ?: 'using target_bed'}
  Zero-probe BED    : ${params.zero_probe_bed ?: 'not provided'}
  Gene list         : ${params.gene_list ?: 'not provided'}
  Reference         : ${params.reference_fasta}
  Aligner           : ${params.aligner}
  Output directory  : ${params.outdir}
  Variant caller    : ${params.variant_caller}
  FF variant caller : ${params.ff_variant_caller}
  ─────────────────────────────────────────────────────────────────
"""

// ── Include modules ─────────────────────────────────────────────────────────
include { FASTP; FASTQ_QC }                                            from './modules/preprocessing'
include { BWA_MEM2_ALIGN; MARK_DUPLICATES; EXTRACT_TARGET_READS }      from './modules/alignment'
include { SAMTOOLS_FLAGSTAT; SAMTOOLS_STATS;
          MOSDEPTH; ON_TARGET_COUNT; BAM_QC_EVALUATE; MULTIQC }        from './modules/bam_qc'
include { COLLECT_FRAGMENT_SIZES; SAMTOOLS_IDXSTATS;
          MPILEUP_AT_SNPS; VARIANT_CALL_FOR_FF;
          ESTIMATE_FETAL_FRACTION }                                     from './modules/fetal_fraction'
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

// ── Helper: Join bam + bai into tuple(sample_id, bam, bai) ─────────────────
def join_bam_bai(ch_bam, ch_bai) {
    ch_bam
        .join(ch_bai, by: 0)
        .map { sample_id, bam, bai -> tuple(sample_id, bam, bai) }
}


// ── Main workflow ───────────────────────────────────────────────────────────
workflow {

    // ── Input channels ──────────────────────────────────────────────────
    ch_reads       = parse_samplesheet(params.input)
    ch_target_bed  = Channel.fromPath(params.target_bed)
    ch_ref_fasta   = Channel.fromPath(params.reference_fasta)
    ch_ref_fai     = Channel.fromPath("${params.reference_fasta}.fai")
    ch_bwa_index   = Channel.fromPath("${params.bwa_index_dir}/*").collect()

    // FF SNPs BED: use dedicated ff_snps_bed if provided, otherwise fall back to target_bed
    ch_ff_snps_bed = params.ff_snps_bed
        ? Channel.fromPath(params.ff_snps_bed)
        : Channel.fromPath(params.target_bed)

    // Optional BED / gene list channels (provide dummy file name when absent)
    ch_probes_bed     = params.probes_bed     ? Channel.fromPath(params.probes_bed)     : Channel.of(file('NO_PROBES_BED'))
    ch_zero_probe_bed = params.zero_probe_bed ? Channel.fromPath(params.zero_probe_bed) : Channel.of(file('NO_ZERO_PROBE_BED'))
    ch_gene_list      = params.gene_list      ? Channel.fromPath(params.gene_list)      : Channel.of(file('NO_GENE_LIST'))

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

    // Create combined bam+bai channel for downstream processes
    ch_dedup = join_bam_bai(
        MARK_DUPLICATES.out.dedup_bam,
        MARK_DUPLICATES.out.dedup_bai,
    )

    // ══════════════════════════════════════════════════════════════════════
    // STEP 4: Extract Target Region Reads
    // ══════════════════════════════════════════════════════════════════════
    EXTRACT_TARGET_READS(
        MARK_DUPLICATES.out.dedup_bam,
        MARK_DUPLICATES.out.dedup_bai,
        ch_target_bed.first(),
    )

    ch_target = join_bam_bai(
        EXTRACT_TARGET_READS.out.target_bam,
        EXTRACT_TARGET_READS.out.target_bai,
    )

    // ══════════════════════════════════════════════════════════════════════
    // STEP 5: BAM Quality Control
    // ══════════════════════════════════════════════════════════════════════
    SAMTOOLS_FLAGSTAT(ch_dedup)

    SAMTOOLS_STATS(ch_dedup)

    MOSDEPTH(
        ch_dedup,
        ch_target_bed.first(),
    )

    ON_TARGET_COUNT(
        ch_dedup,
        ch_target_bed.first(),
    )

    BAM_QC_EVALUATE(
        SAMTOOLS_FLAGSTAT.out.flagstat,
        SAMTOOLS_STATS.out.stats,
        MOSDEPTH.out.summary,
        MOSDEPTH.out.regions,
        ON_TARGET_COUNT.out.on_target_count,
        ch_probes_bed.first(),
        ch_zero_probe_bed.first(),
    )

    // ══════════════════════════════════════════════════════════════════════
    // STEP 6: Fetal Fraction Estimation
    // ══════════════════════════════════════════════════════════════════════

    // 6a. Collect fragment size distribution + chromosome read counts
    //     (ChrX-based FF supplementary since no Y-chr targets in panel)
    COLLECT_FRAGMENT_SIZES(ch_dedup)

    // 6b. samtools idxstats for ChrX-based FF estimation
    SAMTOOLS_IDXSTATS(ch_dedup)

    // 6c. Variant calling at FF SNP sites
    //     Uses ff_snps_bed (dedicated common SNP list) for better FF estimation
    VARIANT_CALL_FOR_FF(
        ch_dedup,
        ch_ff_snps_bed.first(),
        ch_ref_fasta.first(),
        ch_ref_fai.first(),
    )

    // 6d. Estimate fetal fraction (SNP-based primary + ChrX/fragment supplementary)
    ESTIMATE_FETAL_FRACTION(
        VARIANT_CALL_FOR_FF.out.ff_vcf,
        SAMTOOLS_IDXSTATS.out.idxstats,
        COLLECT_FRAGMENT_SIZES.out.fragment_stats,
        COLLECT_FRAGMENT_SIZES.out.fragment_sizes,
    )

    // ══════════════════════════════════════════════════════════════════════
    // STEP 7: Variant Calling at Target Loci
    // ══════════════════════════════════════════════════════════════════════
    VARIANT_CALL_TARGET(
        ch_target,
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
        ch_zero_probe_bed.first(),
        ch_gene_list.first(),
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
