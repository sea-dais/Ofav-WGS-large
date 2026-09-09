#!/usr/bin/env python3
"""
check_site_consistency.py

For one pool, join informative sites to observed ALT-read counts and report:
  - the distribution of observed ALT frequency, split by expected class (0.25 / 0.75)
  - how many sites fall on the "wrong side" of 0.5 (inconsistent with the
    assumed parental genotype), which forces the block fitter to rail.
  - the same for hom_x_hom_diff sites (expected 0.5, FIXED) as a mapping-bias
    calibration: any spread there is purely technical.

Usage:
  python3 check_site_consistency.py parental_table.tsv poolsANGSD PB-11x7-PL 5
"""
import sys, gzip
BASE = {'A':0,'C':1,'G':2,'T':3}

parental, pooldir, pool, minreads = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])

# load sites we care about: informative het x hom (exp 0.25/0.75) and hom x hom diff (0.5)
want = {}
for i, line in enumerate(open(parental)):
    f = line.rstrip("\n").split("\t")
    if i == 0:
        H = {n:j for j,n in enumerate(f)}; continue
    klass = f[H["klass"]]
    chrom = f[H["CHROM"]]; pos = int(f[H["POS"]])
    ref = f[H["REF"]]; alt = f[H["ALT"]]
    exp = f[H["exp_alt_freq"]]
    if klass in ("het_x_hom","hom_x_het") and f[H["block_id"]] != ".":
        want[(chrom,pos)] = (ref, alt, float(exp), klass)
    elif klass == "hom_x_hom_diff":
        want[(chrom,pos)] = (ref, alt, 0.5, "homdiff")

# join to pool counts
posf = f"{pooldir}/freq_{pool}.pos.gz"
cntf = f"{pooldir}/freq_{pool}.counts.gz"
buckets = {"exp0.25":[], "exp0.75":[], "homdiff":[]}
with gzip.open(posf,"rt") as pf, gzip.open(cntf,"rt") as cf:
    pf.readline(); cf.readline()
    for pl, cl in zip(pf, cf):
        pc = pl.split(); chrom, pos = pc[0], int(pc[1])
        w = want.get((chrom,pos))
        if w is None: continue
        ref, alt, exp, klass = w
        cnt = [int(x) for x in cl.split()]
        refn = cnt[BASE[ref]]; altn = cnt[BASE[alt]]
        n = refn + altn
        if n < minreads: continue
        freq = altn / n
        if klass == "homdiff":
            buckets["homdiff"].append(freq)
        elif abs(exp-0.25) < 1e-6:
            buckets["exp0.25"].append(freq)
        elif abs(exp-0.75) < 1e-6:
            buckets["exp0.75"].append(freq)

def summarize(name, vals, expected):
    if not vals:
        print(f"{name}: no sites"); return
    vals.sort()
    n = len(vals)
    med = vals[n//2]
    mean = sum(vals)/n
    wrong = sum(1 for v in vals if (expected==0.25 and v>0.5) or (expected==0.75 and v<0.5))
    print(f"{name}: n={n}  mean={mean:.3f}  median={med:.3f}  expected={expected}")
    if expected in (0.25, 0.75):
        print(f"        on WRONG side of 0.5: {wrong} ({100*wrong/n:.1f}%)")
    # crude histogram
    h=[0]*10
    for v in vals:
        b=min(int(v*10),9); h[b]+=1
    print("        hist:", " ".join(f"{c}" for c in h))

print(f"=== pool {pool}, minreads {minreads} ===")
summarize("het x hom-ref (exp 0.25)", buckets["exp0.25"], 0.25)
summarize("het x hom-alt (exp 0.75)", buckets["exp0.75"], 0.75)
summarize("hom x hom-diff (exp 0.50, FIXED = mapping-bias calibration)", buckets["homdiff"], 0.5)


