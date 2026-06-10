#!/usr/bin/env python3
"""Check coverage at variant positions across multiple BAMs."""
import sys
import pysam

def check_coverage(tsv_path, bam_paths):
    variants = []
    with open(tsv_path) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            chrom, pos, ref, alt, gene, label = parts[:6]
            variants.append((chrom, int(pos), ref, alt, gene, label))

    bam_names = [p.split('/')[-1].replace('.target.bam', '') for p in bam_paths]

    header = ['gene', 'label', 'chrom', 'pos', 'ref', 'alt'] + bam_names
    print('\t'.join(header))

    for chrom, pos, ref, alt, gene, label in variants:
        depths = []
        for bam_path in bam_paths:
            try:
                bam = pysam.AlignmentFile(bam_path, 'rb')
                depth = 0
                for col in bam.pileup(chrom, pos - 1, pos, truncate=True,
                                      min_base_quality=0, min_mapping_quality=0):
                    if col.reference_pos == pos - 1:
                        depth = col.nsegments
                bam.close()
            except Exception:
                depth = 0
            depths.append(str(depth))
        print('\t'.join([gene, label, chrom, str(pos), ref, alt] + depths))

if __name__ == '__main__':
    tsv = sys.argv[1]
    bams = sys.argv[2:]
    check_coverage(tsv, bams)
