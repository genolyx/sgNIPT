#!/usr/bin/env python3
"""
spike_indel_bam.py - Spike SNVs and indels into an existing BAM file.

Supports:
  - SNV: simple base substitution
  - DEL: deletion (modifies CIGAR to add D operation, removes query bases)
  - INS: insertion (modifies CIGAR to add I operation, inserts bases into query)

Each variant is spiked at a target VAF (default 0.50).
"""

import argparse
import pysam
import random
import sys
import os
from typing import Optional, List, Tuple, Dict

# ─────────────────────────────────────────────
# CIGAR operation constants
OP_M = 0   # match/mismatch (consumes ref+query)
OP_I = 1   # insertion (consumes query)
OP_D = 2   # deletion  (consumes ref)
OP_N = 3   # intron/skip (consumes ref)
OP_S = 4   # soft clip (consumes query)
OP_H = 5   # hard clip (consumes neither)
OP_P = 6   # padding
OP_EQ = 7  # sequence match
OP_X = 8   # sequence mismatch

CONSUMES_REF   = {OP_M, OP_D, OP_N, OP_EQ, OP_X}
CONSUMES_QUERY = {OP_M, OP_I, OP_S, OP_EQ, OP_X}


def walk_cigar(read: pysam.AlignedSegment):
    """Yield (ref_pos, query_pos, op, op_len) for each CIGAR block."""
    ref_pos   = read.reference_start
    query_pos = 0
    for op, length in read.cigartuples:
        yield ref_pos, query_pos, op, length
        if op in CONSUMES_REF:
            ref_pos += length
        if op in CONSUMES_QUERY:
            query_pos += length


def find_query_pos_for_ref(read: pysam.AlignedSegment, target_ref: int) -> Optional[int]:
    """Return the query position corresponding to reference position `target_ref`."""
    ref_pos = read.reference_start
    query_pos = 0
    for op, length in read.cigartuples:
        if op in (OP_M, OP_EQ, OP_X):
            for i in range(length):
                if ref_pos == target_ref:
                    return query_pos
                ref_pos += 1
                query_pos += 1
        elif op == OP_I:
            query_pos += length
        elif op == OP_D:
            for i in range(length):
                if ref_pos == target_ref:
                    return None  # ref pos is in a deletion
                ref_pos += 1
        elif op == OP_N:
            ref_pos += length
        elif op in (OP_S,):
            query_pos += length
        elif op in (OP_H, OP_P):
            pass
    return None


def modify_read_snv(read: pysam.AlignedSegment, ref_pos: int, alt_base: str,
                    rng: random.Random) -> Optional[pysam.AlignedSegment]:
    """Replace the base at ref_pos with alt_base (SNV)."""
    q = find_query_pos_for_ref(read, ref_pos)
    if q is None:
        return None
    seq = list(read.query_sequence)
    if seq[q].upper() == alt_base.upper():
        return None  # already ALT
    seq[q] = alt_base
    # Save qualities BEFORE copy+set_sequence (pysam resets query_qualities
    # to None whenever query_sequence is assigned)
    saved_quals = read.query_qualities
    new_read = read.__copy__()
    new_read.query_sequence = ''.join(seq)   # ← wipes query_qualities
    new_read.query_qualities = saved_quals   # ← restore
    return new_read


def modify_read_deletion(read: pysam.AlignedSegment, del_ref_start: int,
                         del_len: int, rng: random.Random) -> Optional[pysam.AlignedSegment]:
    """
    Spike a deletion: remove `del_len` query bases starting at ref position `del_ref_start`
    and insert a D CIGAR operation.
    """
    if read.reference_start > del_ref_start:
        return None
    if read.reference_end < del_ref_start + del_len:
        return None

    seq = list(read.query_sequence)
    quals = list(read.query_qualities) if read.query_qualities else None
    new_cigar: List[Tuple[int, int]] = []
    new_seq: List[str] = []
    new_quals: List[int] = [] if quals else None

    ref_pos = read.reference_start
    query_pos = 0
    deletion_done = False

    for op, length in read.cigartuples:
        if deletion_done:
            # Copy remaining CIGAR unchanged
            new_cigar.append((op, length))
            if op in CONSUMES_QUERY:
                chunk = seq[query_pos:query_pos + length]
                new_seq.extend(chunk)
                if quals:
                    new_quals.extend(quals[query_pos:query_pos + length])
                query_pos += length
            if op in CONSUMES_REF:
                ref_pos += length
            continue

        if op in (OP_M, OP_EQ, OP_X):
            ref_end = ref_pos + length
            if del_ref_start >= ref_pos and del_ref_start < ref_end:
                # Deletion starts within this M block
                if del_ref_start + del_len > ref_end:
                    # Deletion spans beyond this block — too complex, skip
                    return None

                pre_len = del_ref_start - ref_pos  # bases before deletion in ref
                post_len = ref_end - del_ref_start - del_len  # bases after deletion in ref

                # pre-deletion bases
                if pre_len > 0:
                    new_cigar.append((op, pre_len))
                    new_seq.extend(seq[query_pos:query_pos + pre_len])
                    if quals:
                        new_quals.extend(quals[query_pos:query_pos + pre_len])
                    query_pos += pre_len

                # skip del_len query bases (they are being deleted)
                query_pos += del_len

                # add D operation
                new_cigar.append((OP_D, del_len))

                # post-deletion bases
                if post_len > 0:
                    new_cigar.append((op, post_len))
                    new_seq.extend(seq[query_pos:query_pos + post_len])
                    if quals:
                        new_quals.extend(quals[query_pos:query_pos + post_len])
                    query_pos += post_len

                ref_pos = ref_end
                deletion_done = True
            else:
                # No overlap with deletion
                new_cigar.append((op, length))
                new_seq.extend(seq[query_pos:query_pos + length])
                if quals:
                    new_quals.extend(quals[query_pos:query_pos + length])
                query_pos += length
                ref_pos += length

        elif op == OP_I:
            new_cigar.append((op, length))
            new_seq.extend(seq[query_pos:query_pos + length])
            if quals:
                new_quals.extend(quals[query_pos:query_pos + length])
            query_pos += length

        elif op == OP_D:
            if ref_pos <= del_ref_start < ref_pos + length:
                # Requested deletion overlaps existing deletion — skip
                return None
            new_cigar.append((op, length))
            ref_pos += length

        elif op == OP_N:
            new_cigar.append((op, length))
            ref_pos += length

        elif op == OP_S:
            new_cigar.append((op, length))
            new_seq.extend(seq[query_pos:query_pos + length])
            if quals:
                new_quals.extend(quals[query_pos:query_pos + length])
            query_pos += length

        elif op in (OP_H, OP_P):
            new_cigar.append((op, length))

    if not deletion_done:
        return None

    new_read = read.__copy__()
    new_read.query_sequence = ''.join(new_seq)   # ← wipes query_qualities
    new_read.cigartuples = new_cigar
    # Restore qualities; new_quals has del_len fewer entries than original
    if new_quals is not None:
        new_read.query_qualities = new_quals
    # Remove stale MD/NM tags (they reference old sequence; set to 0 as placeholder)
    for tag in ('MD', 'NM'):
        if new_read.has_tag(tag):
            try:
                new_read.set_tag('NM', 0) if tag == 'NM' else new_read.set_tag('MD', '')
            except Exception:
                pass
    return new_read


def modify_read_insertion(read: pysam.AlignedSegment, ins_ref_pos: int,
                          ins_seq: str, rng: random.Random) -> Optional[pysam.AlignedSegment]:
    """
    Spike an insertion: insert `ins_seq` after the base at ref position `ins_ref_pos`.
    ins_ref_pos is the last base BEFORE the insertion (the anchor base).
    """
    if read.reference_start > ins_ref_pos:
        return None
    if read.reference_end <= ins_ref_pos:
        return None

    seq = list(read.query_sequence)
    quals = list(read.query_qualities) if read.query_qualities else None
    ins_len = len(ins_seq)
    # Choose a quality score for inserted bases (use median read quality or 30)
    if quals:
        ins_qual = int(sorted(quals)[len(quals)//2])
    else:
        ins_qual = 30

    new_cigar: List[Tuple[int, int]] = []
    new_seq: List[str] = []
    new_quals: List[int] = [] if quals else None

    ref_pos = read.reference_start
    query_pos = 0
    insertion_done = False

    for op, length in read.cigartuples:
        if insertion_done:
            new_cigar.append((op, length))
            if op in CONSUMES_QUERY:
                new_seq.extend(seq[query_pos:query_pos + length])
                if quals:
                    new_quals.extend(quals[query_pos:query_pos + length])
                query_pos += length
            if op in CONSUMES_REF:
                ref_pos += length
            continue

        if op in (OP_M, OP_EQ, OP_X):
            ref_end = ref_pos + length
            if ins_ref_pos >= ref_pos and ins_ref_pos < ref_end:
                # Insertion anchor is within this M block
                pre_len = ins_ref_pos - ref_pos + 1  # include anchor base
                post_len = ref_end - ins_ref_pos - 1

                # pre-insertion bases (including anchor)
                if pre_len > 0:
                    new_cigar.append((op, pre_len))
                    new_seq.extend(seq[query_pos:query_pos + pre_len])
                    if quals:
                        new_quals.extend(quals[query_pos:query_pos + pre_len])
                    query_pos += pre_len

                # inserted bases
                new_cigar.append((OP_I, ins_len))
                new_seq.extend(list(ins_seq))
                if quals:
                    new_quals.extend([ins_qual] * ins_len)

                # post-insertion bases
                if post_len > 0:
                    new_cigar.append((op, post_len))
                    new_seq.extend(seq[query_pos:query_pos + post_len])
                    if quals:
                        new_quals.extend(quals[query_pos:query_pos + post_len])
                    query_pos += post_len

                ref_pos = ref_end
                insertion_done = True
            else:
                new_cigar.append((op, length))
                new_seq.extend(seq[query_pos:query_pos + length])
                if quals:
                    new_quals.extend(quals[query_pos:query_pos + length])
                query_pos += length
                ref_pos += length

        elif op == OP_I:
            new_cigar.append((op, length))
            new_seq.extend(seq[query_pos:query_pos + length])
            if quals:
                new_quals.extend(quals[query_pos:query_pos + length])
            query_pos += length

        elif op == OP_D:
            new_cigar.append((op, length))
            ref_pos += length

        elif op == OP_N:
            new_cigar.append((op, length))
            ref_pos += length

        elif op == OP_S:
            new_cigar.append((op, length))
            new_seq.extend(seq[query_pos:query_pos + length])
            if quals:
                new_quals.extend(quals[query_pos:query_pos + length])
            query_pos += length

        elif op in (OP_H, OP_P):
            new_cigar.append((op, length))

    if not insertion_done:
        return None

    new_read = read.__copy__()
    new_read.query_sequence = ''.join(new_seq)   # ← wipes query_qualities
    new_read.cigartuples = new_cigar
    # Restore qualities; new_quals has ins_len more entries than original
    if new_quals is not None:
        new_read.query_qualities = new_quals
    for tag in ('MD', 'NM'):
        if new_read.has_tag(tag):
            try:
                new_read.set_tag('NM', 0) if tag == 'NM' else new_read.set_tag('MD', '')
            except Exception:
                pass
    return new_read


def spike_variant(bam_in: pysam.AlignmentFile, chrom: str, variant: dict,
                  target_vaf: float, rng: random.Random) -> Tuple[List, int, int]:
    """
    Return (list_of_modified_reads, n_attempted, n_succeeded).
    `variant` keys: type (SNV/DEL/INS), ref_pos, ref_allele, alt_allele
    For DEL: del_ref_start = first deleted ref base, del_len = number of deleted bases
    For INS: ins_ref_pos = anchor ref base, ins_seq = inserted sequence
    """
    vtype = variant['type']
    modified_reads: List[pysam.AlignedSegment] = []
    n_ok = 0
    n_fail = 0

    if vtype == 'SNV':
        ref_pos = variant['ref_pos']
        alt_base = variant['alt_allele']
        for read in bam_in.fetch(chrom, ref_pos, ref_pos + 1):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if not read.cigartuples:
                continue
            if rng.random() >= target_vaf:
                continue
            new_read = modify_read_snv(read, ref_pos, alt_base, rng)
            if new_read is not None:
                modified_reads.append(new_read)
                n_ok += 1
            else:
                n_fail += 1

    elif vtype == 'DEL':
        del_start = variant['del_ref_start']
        del_len = variant['del_len']
        for read in bam_in.fetch(chrom, del_start, del_start + del_len):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if not read.cigartuples:
                continue
            if rng.random() >= target_vaf:
                continue
            new_read = modify_read_deletion(read, del_start, del_len, rng)
            if new_read is not None:
                modified_reads.append(new_read)
                n_ok += 1
            else:
                n_fail += 1

    elif vtype == 'INS':
        ins_ref_pos = variant['ins_ref_pos']
        ins_seq = variant['ins_seq']
        for read in bam_in.fetch(chrom, ins_ref_pos, ins_ref_pos + 1):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if not read.cigartuples:
                continue
            if rng.random() >= target_vaf:
                continue
            new_read = modify_read_insertion(read, ins_ref_pos, ins_seq, rng)
            if new_read is not None:
                modified_reads.append(new_read)
                n_ok += 1
            else:
                n_fail += 1

    return modified_reads, n_ok + n_fail, n_ok


def main():
    parser = argparse.ArgumentParser(description='Spike SNVs and indels into a BAM file.')
    parser.add_argument('--bam-in',   required=True, help='Input BAM (indexed)')
    parser.add_argument('--bam-out',  required=True, help='Output BAM path')
    parser.add_argument('--vaf',      type=float, default=0.5, help='Target VAF for all variants (default 0.5)')
    parser.add_argument('--seed',     type=int, default=42, help='Random seed')
    args = parser.parse_args()

    rng = random.Random(args.seed)

    # ─── Variant definitions (hg38) ────────────────────────────────────────────
    # Format for DEL: del_ref_start = first ref base being deleted (0-based)
    #                 del_len = number of bases deleted
    # Format for INS: ins_ref_pos = anchor (last matching ref base before insertion, 0-based)
    #                 ins_seq = bases to insert
    # Format for SNV: ref_pos = 0-based ref position
    #                 alt_allele = single character
    VARIANTS = [
        # BRCA1 c.68_69delAG (185delAG)
        # ClinVar: POS=43124027 REF=ACT ALT=A  → delete CT at 43124028-43124029
        {
            'label': 'BRCA1 c.68_69delAG (185delAG)',
            'chrom': 'chr17',
            'type':  'DEL',
            'del_ref_start': 43124027,  # 0-based: del CT starting here (was 43124028 VCF=1-based → 43124027 0-based)
            'del_len': 2,
        },
        # BRCA1 c.5266dupC (5382insC)
        # ClinVar: POS=43057062 REF=T ALT=TG → insert G after 43057062 (0-based: 43057061)
        {
            'label': 'BRCA1 c.5266dupC (5382insC)',
            'chrom': 'chr17',
            'type':  'INS',
            'ins_ref_pos': 43057061,  # 0-based anchor position (VCF POS 43057062 - 1)
            'ins_seq': 'G',
        },
        # BRCA2 c.5946delT (6174delT)
        # ClinVar ID=9325: POS=32340300 REF=GT ALT=G → delete T at 32340301 (0-based: 32340300)
        {
            'label': 'BRCA2 c.5946delT (6174delT)',
            'chrom': 'chr13',
            'type':  'DEL',
            'del_ref_start': 32340300,  # 0-based (VCF POS 32340300, REF=GT, ALT=G → del 32340301 VCF = 32340300 0-based)
            'del_len': 1,
        },
        # BRCA1 c.2612C>T
        # BRCA1 minus strand: coding C>T → genomic G>A at chr17:43092919
        {
            'label': 'BRCA1 c.2612C>T',
            'chrom': 'chr17',
            'type':  'SNV',
            'ref_pos': 43092918,  # 0-based (VCF POS 43092919 - 1)
            'alt_allele': 'A',   # G→A on + strand
        },
        # BRCA1 c.3548A>G
        # BRCA1 minus strand: coding A>G → genomic T>C at chr17:43091983
        {
            'label': 'BRCA1 c.3548A>G',
            'chrom': 'chr17',
            'type':  'SNV',
            'ref_pos': 43091982,  # 0-based
            'alt_allele': 'C',   # T→C on + strand
        },
        # BRCA2 c.2806_2809delAAAC
        # ClinVar ID=9322: POS=32337160 REF=TAAAC ALT=T → del AAAC at 32337161-32337164
        {
            'label': 'BRCA2 c.2806_2809delAAAC',
            'chrom': 'chr13',
            'type':  'DEL',
            'del_ref_start': 32337160,  # 0-based (VCF POS 32337160, anchor T, then AAAC at 32337161 VCF = 32337160 0-based)
            'del_len': 4,
        },
        # BRCA2 c.5164_5165delAG
        # ClinVar ID=51791: POS=32339518 REF=CAG ALT=C → del AG at 32339519-32339520
        {
            'label': 'BRCA2 c.5164_5165delAG',
            'chrom': 'chr13',
            'type':  'DEL',
            'del_ref_start': 32339518,  # 0-based (VCF POS 32339518, anchor C, del AG at 32339519 VCF = 32339518 0-based)
            'del_len': 2,
        },
    ]

    # ─── Process BAM ────────────────────────────────────────────────────────────
    tmp_out = args.bam_out + '.tmp_unsorted.bam'
    bam_in = pysam.AlignmentFile(args.bam_in, 'rb')
    bam_out = pysam.AlignmentFile(tmp_out, 'wb', header=bam_in.header)

    # Build a set of modified read names per variant to avoid double-writing
    # Strategy: collect all modified reads per variant, write them, skip originals for modified ones

    print("=== Collecting modified reads for each variant ===", file=sys.stderr)
    # Map: (chrom, read_name, is_read1) → modified read
    modified_map: Dict[Tuple, pysam.AlignedSegment] = {}

    for v in VARIANTS:
        chrom = v['chrom']
        modified_reads, n_total, n_ok = spike_variant(bam_in, chrom, v, args.vaf, rng)
        print(f"  [{v['label']}]  attempted={n_total}  spiked={n_ok}  "
              f"actual_vaf={n_ok/max(n_total,1):.2%}", file=sys.stderr)
        for r in modified_reads:
            key = (r.reference_name, r.query_name, r.is_read1)
            modified_map[key] = r

    print(f"\nTotal modified reads: {len(modified_map)}", file=sys.stderr)

    # ─── Write output ────────────────────────────────────────────────────────────
    print("Writing output BAM ...", file=sys.stderr)
    bam_in.reset()
    written = 0
    modified_written = 0
    for read in bam_in.fetch(until_eof=True):
        key = (read.reference_name, read.query_name, read.is_read1)
        if key in modified_map:
            bam_out.write(modified_map[key])
            del modified_map[key]
            modified_written += 1
        else:
            bam_out.write(read)
        written += 1

    bam_in.close()
    bam_out.close()

    print(f"  Written: {written} reads total, {modified_written} modified", file=sys.stderr)

    # ─── Sort and index ──────────────────────────────────────────────────────────
    print(f"Sorting → {args.bam_out} ...", file=sys.stderr)
    pysam.sort('-o', args.bam_out, tmp_out)
    os.remove(tmp_out)
    print(f"Indexing ...", file=sys.stderr)
    pysam.index(args.bam_out)
    print(f"Done. Output: {args.bam_out}", file=sys.stderr)


if __name__ == '__main__':
    main()
