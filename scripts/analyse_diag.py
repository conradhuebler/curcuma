#!/usr/bin/env python3
"""
Quick reader for curcuma MD diagnostics JSONL files.

Format (one JSON object per line):
  {step, time_fs, energy{...}, charges[N], cn[N], gradient_norm[N], hb_count, xb_count}

Usage:
  python3 analyse_diag.py <basename>.diag.jsonl
"""
import json
import sys


def main(path):
    print(f"{'step':>8} {'time_fs':>10} {'E_total':>14} {'max|q|':>8} {'max(CN)':>8} {'hb':>3} {'xb':>3}")
    for line in open(path):
        rec = json.loads(line)
        qs = rec.get("charges", [])
        cn = rec.get("cn", [])
        qmax = max(abs(q) for q in qs) if qs else 0.0
        cnmax = max(cn) if cn else 0.0
        e_total = sum(rec.get("energy", {}).values())
        print(f"{rec['step']:>8d} {rec['time_fs']:>10.2f} {e_total:>14.6f} "
              f"{qmax:>8.4f} {cnmax:>8.3f} {rec.get('hb_count', 0):>3d} "
              f"{rec.get('xb_count', 0):>3d}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1])
