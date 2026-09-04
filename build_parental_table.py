#!/usr/bin/env python3
"""
Build the parental expected-frequency + phase + block table.

INPUT (stdin): output of
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT:%PS]\n' phased.vcf.gz
  (optionally pre-filtered to informative sites; this script re-derives the class anyway)

Sample order in the query is whatever the VCF has; we locate parents by NAME,
so set the two parent names below.

OUTPUT (stdout): TSV with header
  CHROM POS REF ALT gt_p1 gt_p2 klass seg_parent block_id hapA_allele exp_alt_freq informative

Conventions
-----------
- "haplotype A" of the segregating parent = its haplotype-1 (first slot of GT).
  WhatsHap keeps the same physical haplotype in slot 1 across a phase set (PS),
  so hapA is consistent within a block.
- hapA_allele = allele on slot 1: for "0|1" -> REF, for "1|0" -> ALT.
- exp_alt_freq = expected ALT-allele frequency in the offspring pool under
  Mendelian transmission, given both parental genotypes.

Class / expected ALT frequency table (P(offspring allele = ALT)):
  het x hom_ref  (0/1 x 0/0) -> 0.25   (one parent segregates ALT at 1/2, other gives REF)
  het x hom_alt  (0/1 x 1/1) -> 0.75
  het x het      (0/1 x 0/1) -> 0.50   (both segregate; flagged, handle separately)
  hom_ref x hom_alt          -> 0.50   but FIXED (all offspring het): zero segregation
                                        variance -> informative=0 (QC only)
  hom x hom (same)           -> uninformative -> dropped
"""

import sys

# ---- SET YOUR PARENT SAMPLE NAMES ----
P1_NAME = "11-AD"
P2_NAME = "7-AD"
# --------------------------------------

def parse_field(field):
    """'11-AD=0|1:17712' -> ('11-AD', '0|1', '17712')  (PS may be '.')"""
    name, rest = field.split("=", 1)
    gt, ps = rest.split(":", 1)
    return name, gt, ps

def gt_kind(gt):
    """Classify a genotype string. Returns one of:
       'homref','homalt','het', or 'other' (missing/multiallelic-ish)."""
    if gt in ("0|0", "0/0"):
        return "homref"
    if gt in ("1|1", "1/1"):
        return "homalt"
    if gt in ("0|1", "1|0", "0/1", "1/0"):
        return "het"
    return "other"

def is_phased(gt):
    return "|" in gt

def hapA_allele(gt, ref, alt):
    """Allele on haplotype-1 (slot 1) for a phased het GT."""
    a1 = gt.split("|")[0]
    return ref if a1 == "0" else alt

def main():
    out = sys.stdout
    out.write("\t".join([
        "CHROM","POS","REF","ALT","gt_p1","gt_p2","klass",
        "seg_parent","block_id","hapA_allele","exp_alt_freq","informative"
    ]) + "\n")

    for line in sys.stdin:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
        sample_fields = parts[4:]

        gts = {}
        pss = {}
        for f in sample_fields:
            name, gt, ps = parse_field(f)
            gts[name] = gt
            pss[name] = ps

        if P1_NAME not in gts or P2_NAME not in gts:
            continue  # parents not found in this record

        g1, g2 = gts[P1_NAME], gts[P2_NAME]
        k1, k2 = gt_kind(g1), gt_kind(g2)

        # skip anything with an unparseable/missing parent genotype
        if k1 == "other" or k2 == "other":
            continue

        klass = None
        seg_parent = "."
        block_id = "."
        hapA = "."
        exp = "."
        informative = 0

        # --- classify ---
        homs = {"homref", "homalt"}

        if k1 == "het" and k2 in homs:
            klass = "het_x_hom"
            seg_parent = P1_NAME
            exp = 0.25 if k2 == "homref" else 0.75
            informative = 1
            if is_phased(g1):
                block_id = pss[P1_NAME]
                hapA = hapA_allele(g1, ref, alt)
        elif k2 == "het" and k1 in homs:
            klass = "hom_x_het"
            seg_parent = P2_NAME
            exp = 0.25 if k1 == "homref" else 0.75
            informative = 1
            if is_phased(g2):
                block_id = pss[P2_NAME]
                hapA = hapA_allele(g2, ref, alt)
        elif k1 == "het" and k2 == "het":
            klass = "het_x_het"
            seg_parent = "both"
            exp = 0.50
            informative = 1
            # record P1's block/phase in the main columns; P2 handled downstream
            if is_phased(g1):
                block_id = pss[P1_NAME]
                hapA = hapA_allele(g1, ref, alt)
        elif k1 in homs and k2 in homs:
            if k1 != k2:
                klass = "hom_x_hom_diff"   # all offspring het; fixed, no segregation
                exp = 0.50
                informative = 0            # QC only, not a distortion test
            else:
                klass = "hom_x_hom_same"   # uninformative
                exp = 0.0 if k1 == "homref" else 1.0
                informative = 0
        else:
            continue

        out.write("\t".join([
            chrom, pos, ref, alt, g1, g2, klass,
            seg_parent, str(block_id), hapA, str(exp), str(informative)
        ]) + "\n")

if __name__ == "__main__":
    main()