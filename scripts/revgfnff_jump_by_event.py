import re, statistics, sys
from pathlib import Path
TERMS = ["bond","angle","tors","inv","brep","nbrep","coul","disp","hb","over","batm"]
for r in sys.argv[1:]:
    txt = Path(r).read_text(); groups = {}; prom_s = []; scan = {}
    for line in txt.splitlines():
        m = re.search(r"REACT bond (formed|broken): (\w+)-(\w+) \(blend .*r_scan ([\d.]+) w_scan ([\d.]+)", line)
        if m: scan[(m.group(2), m.group(3))] = (float(m.group(4)), float(m.group(5)))
        m = re.search(r"blend promoted: .* at s = ([\d.]+)", line)
        if m: prom_s.append(float(m.group(1)))
        m = re.search(r"REACT jump terms \[kJ/mol\] \((\w+)\): (.*?) \| s ([-\d.]+) swap ([-+\d.]+) w_ws ([\d.]+) (?:c_ws [\d.]+ )?r_ws ([\d.]+) ij (\d+) (\d+)", line)
        if m:
            d = {t: float(v) for t, v in re.findall(r"(\w+) ([-+\d.]+)", m.group(2)) if t in TERMS}
            d["_s"] = float(m.group(3)); d["_w"] = float(m.group(5)); d["_r"] = float(m.group(6))
            groups.setdefault(m.group(1), []).append(d)
    med = lambda v: statistics.median([abs(x) for x in v]) if v else 0.0
    print(f"== {Path(r).parent.name}: promotes n={len(prom_s)} s<0.5: {sum(s<0.5 for s in prom_s)} median s {statistics.median(prom_s) if prom_s else '-'}")
    for k, rows in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        bonded = [d["bond"]+d["angle"]+d["tors"]+d["inv"] for d in rows]
        rep = [d["brep"]+d["nbrep"] for d in rows]
        onb = [d["disp"]+d["hb"]+d["over"]+d["batm"] for d in rows]
        coul = [d["coul"] for d in rows]
        tot = [b+p+o+c for b,p,o,c in zip(bonded,rep,onb,coul)]
        s0 = sum(1 for d in rows if d["_s"] < 0.05 or d["_s"] < 0); s1 = sum(1 for d in rows if d["_s"] >= 0.95)
        print(f"  {k:12s} n={len(rows):3d} |total| {med(tot):6.1f} | bonded {med(bonded):6.1f} rep {med(rep):5.1f} otherNB {med(onb):4.1f} coul {med(coul):5.1f} | s<0.05: {s0} s>=0.95: {s1}")
    # w_scan vs w_ws check
    ws = [(d["_w"], d["_r"]) for d in groups.get("begin_form", []) + groups.get("begin_break", [])]
    sc = list(scan.values())
    print(f"  scan lines with r/w: {len(sc)}; ws begin lines: {len(ws)}; sample scan {sc[:2]} ws {ws[:2]}")
