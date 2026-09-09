#!/usr/bin/env python3
"""
fit_block_t.py  --  Stage A of the block-level TRD test.

For each phased het x hom block, in each pool, fit a single transmission ratio
t = fraction of the segregating parent's gametes carrying haplotype A, by
beta-binomial maximum likelihood over the block's SNP read counts. Test t != 0.5
with a likelihood-ratio test.

INPUTS
------
--parental parental_table.tsv        (from build_parental_table.py)
--pooldir  poolsANGSD                 (holds freq_<POOL>.pos.gz and .counts.gz)
--pools    P1-7x11-PL,P2-7x11-PL,...  (comma list; pool label taken from name)
--minreads 10                         (per-site floor on REFn+ALTn)
--maxreads 100                        (per-site upper cap; drops repeat-collapse sites)
--min-snps 5                          (blocks with fewer usable SNPs are skipped)
--out      block_t_perpool.tsv

Only het_x_hom / hom_x_het blocks with a real block_id are used (the clean,
single-segregating-parent case). het_x_het is skipped in this stage.

Blocks are keyed on CHROM + seg_parent + PS, because WhatsHap PS ids are only
unique within a contig; the same PS number recurs across scaffolds.

OUTPUT (one row per block x pool with >= --min-snps usable SNPs)
---------------------------------------------------------------
pool direction seg_parent chrom block_id n_snps n_reads t_hat rho_hat
ll_full ll_null LRT p_value sign_changes boundary

Requires: numpy, scipy.
"""

import argparse, gzip, sys, math
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import chi2
from scipy.special import betaln

BASE = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
BOUNDARY_EPS = 1e-3

# ------------------------------------------------------------------ IO helpers
def read_parental(path):
    """(chrom,pos) -> SNP dict, restricted to phased het_x_hom / hom_x_het sites."""
    sites = {}
    with open(path) as fh:
        for i, line in enumerate(fh):
            f = line.rstrip("\n").split("\t")
            if i == 0:
                H = {name: j for j, name in enumerate(f)}
                continue
            klass = f[H["klass"]]
            if klass not in ("het_x_hom", "hom_x_het"):
                continue
            block_id = f[H["block_id"]]
            if block_id == ".":
                continue
            chrom = f[H["CHROM"]]; pos = int(f[H["POS"]])
            ref = f[H["REF"]]; alt = f[H["ALT"]]
            seg = f[H["seg_parent"]]
            hapA = f[H["hapA_allele"]]           # allele on seg parent's haplotype 1
            exp = float(f[H["exp_alt_freq"]])    # 0.25 or 0.75
            if hapA not in (ref, alt):
                continue                          # phase missing/odd; skip
            # hom parent's ALT dosage:
            #   exp=0.75 <=> hom parent is hom-ALT (contributes ALT)
            #   exp=0.25 <=> hom parent is hom-REF (contributes REF)
            hom_alt = 1 if abs(exp - 0.75) < 1e-6 else 0
            hapA_is_alt = 1 if hapA == alt else 0
            hapB_is_alt = 1 - hapA_is_alt         # biallelic: hapB is the other allele
            sites[(chrom, pos)] = dict(
                chrom=chrom, pos=pos, ref=ref, alt=alt, seg=seg,
                # CHROM in the key: PS ids are only unique within a contig.
                block=f"{chrom}:{seg}:{block_id}",
                hom_alt=hom_alt, hapA_is_alt=hapA_is_alt, hapB_is_alt=hapB_is_alt,
            )
    return sites


def _count_lines_gz(path):
    n = 0
    with gzip.open(path, "rt") as fh:
        for _ in fh:
            n += 1
    return n


def read_pool_counts(pooldir, pool, sites, minreads, maxreads):
    """Attach (const, coefA, coefB, altn, n) to sites present in this pool.

       Expected ALT freq for a SNP is  mu(t) = const + coefA*t + coefB*(1-t).
       Returns dict block_key -> {"seg":..., "snps":[(const,coefA,coefB,altn,n), ...]}.
    """
    posf = f"{pooldir}/freq_{pool}.pos.gz"
    cntf = f"{pooldir}/freq_{pool}.counts.gz"

    # Guard against silent desync between the two ANGSD outputs.
    npos = _count_lines_gz(posf)
    ncnt = _count_lines_gz(cntf)
    if npos != ncnt:
        sys.exit(f"[error] {pool}: pos.gz ({npos} lines) and counts.gz ({ncnt} lines) "
                 f"differ; files are out of sync, refusing to zip.")

    blocks = {}
    with gzip.open(posf, "rt") as pf, gzip.open(cntf, "rt") as cf:
        pf.readline(); cf.readline()                       # headers
        for pline, cline in zip(pf, cf):
            pc = pline.split()
            chrom, pos = pc[0], int(pc[1])
            s = sites.get((chrom, pos))
            if s is None:
                continue
            cc = cline.split()
            cnt = [int(x) for x in cc]                      # totA totC totG totT
            refn = cnt[BASE[s["ref"]]]
            altn = cnt[BASE[s["alt"]]]
            n = refn + altn
            if n < minreads or n > maxreads:                # lower floor + upper cap
                continue
            const = 0.5 * s["hom_alt"]
            coefA = 0.5 * s["hapA_is_alt"]
            coefB = 0.5 * s["hapB_is_alt"]
            b = blocks.setdefault(s["block"], {"seg": s["seg"], "snps": []})
            b["snps"].append((const, coefA, coefB, altn, n))
    return blocks

# ------------------------------------------------------- beta-binomial likelihood
def _bb_logpmf(k, n, mu, rho):
    """Beta-binomial log pmf. mean mu in (0,1), overdispersion rho in (0,1).
       rho->0 reduces to the binomial. alpha/beta parameterization."""
    mu = min(max(mu, 1e-6), 1 - 1e-6)
    if rho <= 1e-9:
        return (k * math.log(mu) + (n - k) * math.log(1 - mu)
                + math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1))
    s = (1 - rho) / rho
    a = mu * s
    b = (1 - mu) * s
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
            + betaln(k + a, n - k + b) - betaln(a, b))


def _neg_ll(t, rho, snps):
    ll = 0.0
    for const, coefA, coefB, k, n in snps:
        mu = const + coefA * t + coefB * (1 - t)
        ll += _bb_logpmf(k, n, mu, rho)
    return -ll


def _best_rho(t, snps, rho0=0.05):
    """Optimize rho at fixed t over [1e-6, 0.5]. Returns (rho_hat, neg_ll)."""
    res = minimize_scalar(lambda r: _neg_ll(t, min(max(r, 1e-6), 0.5), snps),
                          bounds=(1e-6, 0.5), method="bounded",
                          options={"xatol": 1e-5})
    return min(max(res.x, 1e-6), 0.5), res.fun


def fit_block(snps, grid=101):
    """Fit (t, rho) by profiling t on a grid and optimizing rho at each point.
       More robust than joint Nelder-Mead for small blocks on a bounded surface.
       Returns t_hat, rho_hat, ll_full, ll_null(t=0.5), LRT, p, boundary_flag."""
    ts = np.linspace(0.0, 1.0, grid)
    best_t, best_rho, best_negll = 0.5, 0.05, float("inf")
    for t in ts:
        rho, negll = _best_rho(t, snps)
        if negll < best_negll:
            best_negll, best_t, best_rho = negll, t, rho

    # local refinement of t around the grid optimum
    step = ts[1] - ts[0]
    lo, hi = max(0.0, best_t - step), min(1.0, best_t + step)
    res = minimize_scalar(lambda t: _best_rho(t, snps)[1],
                          bounds=(lo, hi), method="bounded",
                          options={"xatol": 1e-5})
    if res.fun < best_negll:
        best_t = min(max(res.x, 0.0), 1.0)
        best_rho, best_negll = _best_rho(best_t, snps)

    t_hat = best_t
    rho_hat = best_rho
    ll_full = -best_negll

    # null model: t fixed at 0.5, optimize rho only
    _, negll0 = _best_rho(0.5, snps)
    ll_null = -negll0

    lrt = max(2.0 * (ll_full - ll_null), 0.0)

    boundary = 1 if (t_hat < BOUNDARY_EPS or t_hat > 1 - BOUNDARY_EPS) else 0
    # Standard df=1 chi-square. NOTE: when t_hat is on the {0,1} boundary the
    # asymptotic null is a 50:50 mixture of chi2_0 and chi2_1, so p is only
    # approximate there; the boundary flag marks those rows for cautious use
    # (and appropriate handling) in the downstream aggregation.
    p = chi2.sf(lrt, df=1)
    return t_hat, rho_hat, ll_full, ll_null, lrt, p, boundary


def switch_flag(snps, t_hat):
    """Crude switch-error signature: number of sign changes in per-SNP residuals
       (obs-exp) along the position-ordered SNP list. A clean block has few; a
       single WhatsHap switch tends to give one dominant flip. Reported as a
       diagnostic count, thresholded downstream."""
    signs = []
    for const, coefA, coefB, k, n in snps:
        mu = const + coefA * t_hat + coefB * (1 - t_hat)
        obs = k / n if n else mu
        signs.append(1 if obs > mu else -1)
    return sum(1 for a, b in zip(signs, signs[1:]) if a != b)

# ------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parental", required=True)
    ap.add_argument("--pooldir", required=True)
    ap.add_argument("--pools", required=True,
                    help="comma list, e.g. "
                         "P1-7x11-PL,P2-7x11-PL,P3-7x11-PL,PA-11x7-PL,PB-11x7-PL,PC-11x7-PL")
    ap.add_argument("--minreads", type=int, default=10)
    ap.add_argument("--maxreads", type=int, default=100)
    ap.add_argument("--min-snps", type=int, default=5)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sites = read_parental(args.parental)
    sys.stderr.write(f"[info] block-usable parental sites: {len(sites)}\n")

    with open(args.out, "w") as out:
        out.write("\t".join([
            "pool", "direction", "seg_parent", "chrom", "block_id",
            "n_snps", "n_reads", "t_hat", "rho_hat",
            "ll_full", "ll_null", "LRT", "p_value", "sign_changes", "boundary"
        ]) + "\n")

        for pool in args.pools.split(","):
            direction = "7x11" if "7x11" in pool else ("11x7" if "11x7" in pool else "NA")
            blocks = read_pool_counts(args.pooldir, pool, sites,
                                      args.minreads, args.maxreads)
            n_fit = 0
            for block_key, d in blocks.items():
                snps = d["snps"]
                if len(snps) < args.min_snps:
                    continue
                n_reads = sum(n for *_, n in snps)
                t_hat, rho_hat, llf, lln, lrt, p, boundary = fit_block(snps)
                sc = switch_flag(snps, t_hat)
                seg = d["seg"]
                chrom_out, seg_out, bid = block_key.rsplit(":", 2)
                out.write("\t".join(str(x) for x in [
                    pool, direction, seg, chrom_out, bid, len(snps), n_reads,
                    f"{t_hat:.4f}", f"{rho_hat:.4f}",
                    f"{llf:.3f}", f"{lln:.3f}", f"{lrt:.3f}", f"{p:.3e}", sc, boundary
                ]) + "\n")
                n_fit += 1
            sys.stderr.write(f"[info] {pool}: fit {n_fit} blocks\n")


if __name__ == "__main__":
    main()