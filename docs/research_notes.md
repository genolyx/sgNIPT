# Single Gene NIPT Pipeline - Research Notes

## Fetal Fraction Calculation Methods (SNP-based)

### Method 1: Zhang et al. - Dual VAF approach
- Mother homozygous ref (AA), fetus heterozygous (AB): VAFp ≈ FF/2
- Mother homozygous alt (BB), fetus heterozygous (AB): VAFm ≈ 1 - FF/2
- FF = 2 * mean(VAFp) or FF = 2 * (1 - mean(VAFm))

### Method 2: Jiang et al. - Informative SNP approach  
- Select SNPs where mother is homozygous ref (AA) and fetus is heterozygous (AB)
- For each SNP: q = ref reads, p = alt reads
- FFi = 2p / (p + q)
- FF = mean of all FFi values

### Key Considerations:
- Read depth and number of SNPs are critical factors
- Deeper coverage + more SNPs = better results
- SNPs should have common MAF (>5%) for reliability
- Non-maternal alleles could be from fetus OR sequencing errors
- Somatic variant calling approach needed (not germline)

## Pipeline Architecture for Single Gene NIPT

### Input
- Fastq files from Illumina sequencer
- Target BED file (~200-300 loci)
- hg38 reference genome

### Steps
1. **FASTQ QC**: FastQC, Trimmomatic/fastp
2. **Alignment**: BWA-MEM2 to hg38
3. **BAM Processing**: Sort, MarkDuplicates, BQSR (optional for targeted)
4. **Target Region Extraction**: samtools view with BED file
5. **Fetal Fraction Estimation**: SNP-based method using informative SNPs
6. **Variant Calling**: For target loci - identify fetal-specific variants
7. **Variant Annotation**: Annotate pathogenic variants
8. **QC Metrics**: Coverage, mapping rate, duplicate rate, on-target rate
9. **Report Generation**: Summary report

### Fetal DNA Extraction Logic
- In cfDNA, fetal fraction is typically 5-20%
- Fetal-specific alleles: alleles present in cfDNA but NOT in maternal genome
- Informative SNPs: positions where mother is homozygous, fetus carries different allele
- Need to distinguish fetal variants from sequencing noise using statistical thresholds

### Tools
- fastp: adapter trimming, quality filtering
- BWA-MEM2: alignment
- samtools: BAM processing
- picard: MarkDuplicates
- GATK/freebayes: variant calling
- bcftools: variant filtering
- MultiQC: QC aggregation
