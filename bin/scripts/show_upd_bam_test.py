#!/usr/bin/env python3
"""Print a summary table of UPD BAM end-to-end test results.

Usage:
  show_upd_bam_test.py <out_dir> [<ff_labels>...]
  e.g.  show_upd_bam_test.py /tmp/v2 05 10
"""
import json, os, sys, glob

def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else "/tmp/updtest6"
    ff_labels = sys.argv[2:] if len(sys.argv) > 2 else None
    scenarios = ["normal", "matUPD_iso", "matUPD_het", "patUPD_iso"]

    hdr = ("%-26s  %5s  %-28s  %5s  %5s  %10s  %10s  %8s"
           % ("scenario", "FF", "chr7_call",
              "n_het", "n_hom", "llr_matiso", "llr_patiso", "patupd_z"))
    print()
    print(hdr)
    print("-" * len(hdr))

    def _print_row(key, ff_str, path):
        if not os.path.exists(path):
            print("%-26s  NO_OUTPUT_FILE" % key)
            return
        r = json.load(open(path))
        c = next((x for x in r.get("per_chrom_upd", []) if x.get("chrom") == "chr7"), {})
        llr = c.get("llr_vs_normal", {})
        row = ("%-26s  %5s  %-28s  %5d  %5d  %10.2f  %10.2f  %8.2f" % (
            key, ff_str,
            str(c.get("call", "?")),
            c.get("n_mat_het", 0), c.get("n_mat_hom", 0),
            llr.get("matupd_iso", 0), llr.get("patupd_iso", 0),
            c.get("patupd_z", 0),
        ))
        print(row)

    if ff_labels:
        for fl in ff_labels:
            for s in scenarios:
                path = os.path.join(out_dir, "%s_ff%s_upd.json" % (s, fl))
                _print_row("%s_FF%s" % (s, fl), "0." + fl, path)
            print()
    else:
        for s in scenarios:
            path = os.path.join(out_dir, "%s_upd.json" % s)
            _print_row(s, "?", path)
        print()

if __name__ == "__main__":
    main()

