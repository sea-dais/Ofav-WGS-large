# Segregation-distortion pipeline — working record

Reciprocal *Orbicella faveolata* cross. Goal: detect segregation distortion at
the level of phased haplotype blocks, and use the reciprocal design to split
autosomal distortion (same in both cross directions) from parent-of-origin
effects (flips with direction).

- **Parents:** 11-AD, 7-AD
- **Pools:** 7×11 → P1, P2, P3 · 11×7 → PA, PB, PC (larval pools, ~100 larvae each)
- **Reference:** GCF_002042975.1_ofav_dov_v1 (RefSeq names: NW_ scaffolds, NC_ mito)
- **Phasing:** WhatsHap `phase --ped` (read-backed + pedigree), ~9,980 phased hets, NG50 ~595 kb

This file is the *why*. The clean re-runnable version is `ofav_segdist_pipeline.sh`.

---

## The idea (so the steps make sense)

Each parental genotype pair sets an expected offspring ALT frequency under
Mendelian transmission:

| parents (seg × other) | expected ALT freq | use |
|---|---|---|
| het × hom-ref | 0.25 | informative |
| het × hom-alt | 0.75 | informative |
| het × het     | 0.50 | informative but messy (both segregate) — defer |
| hom-ref × hom-alt | 0.50 but FIXED | QC only (mapping bias / error check) |
| hom × hom same | — | uninformative, dropped |

For a **het × hom block**, the segregating parent's phased SNPs all track one
haplotype. Fit a single transmission ratio *t* (fraction of gametes carrying
haplotype A) across all the block's SNPs; test *t* ≠ 0.5. Coherence across SNPs
in a block is where the power comes from at our depth, and it doubles as a
switch-error check (a block splitting into two opposite-deviation runs = phasing
artifact, not two distorters).

Prior art doing essentially this: Corbett-Detig et al. 2017 (Drosophila pooled
embryos, expected 1:3 ratio, ~30×); Seymour et al. 2019 Heredity (Arabidopsis,
>500 F2 pools, **beta-binomial + FDR** — the statistical template). Our additions:
the reciprocal cross (parent-of-origin) and the phased-block-level fit.

---

## Step 1 — informative sites file

Built from the phased VCF: biallelic SNPs, ≥1 parent het, both parents
confidently genotyped, **REF/ALT fixed as major/minor** so ANGSD doesn't
re-polarize per pool, mitochondrion (NC_*) removed (maternally inherited,
non-Mendelian — a "transmission ratio" there is meaningless).

- Hardened parental depth floor to DP≥10 (up from the DP≥6 used pre-phasing):
  a het miscalled as hom flips expected freq 0.5→0.25/0.75 and manufactures
  fake distortion, so parents specifically get a stricter floor.
- Result: 6,529 nuclear informative sites in the ANGSD file.

## Step 2 — per-pool ANGSD (the frequencies)

One run per pool, six total, shared `-sites` polarization. Key flags and why:
- `-doMajorMinor 3` — major/minor from the sites file. **Without this the six
  pools can pick different minor alleles and the frequencies aren't comparable.**
- `-doMaf 1` — knownEM = ALT-allele frequency (since minor = our ALT).
- `-dumpCounts 3` — exact per-allele A,C,G,T counts (row-aligned to `.pos.gz`),
  for an honest beta-binomial on read counts. (`-dumpCounts 2` gives only total
  depth — not enough; we re-ran with 3.)
- **No `-SNP_pval`** — we're not discovering SNPs. A SNP filter would drop
  near-fixed sites, which are exactly the distortion signal. Never add it here.
- Depths came out ~23× (five pools) / 28× (PB). Lower than the 50× target, so
  the block-level aggregation is doing real work; single-SNP tests would be hopeless.

## Step 3 — parental expected/phase/block table

`build_parental_table.py`: reads `bcftools query '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT:%PS]\n'`
and emits per site: class, segregating parent, block (PS) id, hapA allele,
expected ALT freq, informative flag.
- Phase convention: haplotype A = the parent's slot-1 allele. WhatsHap keeps the
  same physical haplotype in slot 1 across a PS, so hapA is consistent in a block.
  `0|1`→hapA=REF, `1|0`→hapA=ALT.
- Block key is `seg_parent:block_id` — the two parents' phase sets are independent
  even when they share a PS integer.
- Counts: 18,799 informative sites; 14,465 block-usable het×hom across 1,493 blocks.
  Intersection with the (harder-filtered) ANGSD sites = **5,216 measured & block-usable**.
  Power tiers (measured sites): ≥3 SNPs → 586 blocks, ≥5 → 383, ≥10 → 156.
- TODO (deferred): the 14,465 vs 5,216 gap is the DP≥10 hardening vs the DP≥6
  ANGSD path. Could re-run ANGSD on the fuller block-usable set to ~double markers,
  at the cost of some lower-confidence parental calls. Decide after first-pass results.

---

## Debugging traps hit (so we don't repeat them)

1. **`head` of the phased VCF looked all-homozygous / unphased.** False alarm —
   the first records are hom sites at a contig start. `grep -c '|'` showed 17,252
   phased records. Don't judge a phased VCF by its first 10 lines.

2. **`$REF` unset in a fresh shell** → `-ref -doMajorMinor` → "fastafile
   '-doMajorMinor' doesn't exist, will exit". Interactive vars don't survive across
   shells/nodes. In the pylauncher `.cmds`, paths must be **expanded at build time**
   (echoed into the file), not left as `$REF`.

3. **The big one: ANGSD 0.935 `-sites` was broken.** Symptom: runs complete, but
   "retained after filtering: 0" — every site dropped, `.mafs.gz` empty — while
   `samtools depth` showed real coverage at those exact positions and a no-`-sites`
   ANGSD run counted depth fine. Ruled out, in order: wrong output dir; proper-pair
   flags; quality filters; empty regions; contig-name mismatch (BAM vs ref were
   **identical**, 1933 contigs, `comm -3` empty); sites-file delimiter (clean tabs,
   4 cols); position sort; stale index; major/minor parsing (2-col + `-doMajorMinor 1`
   still zero); single known-covered site (depth 93, still zero). Conclusion by
   elimination: the `-sites` implementation in the 0.935 pre-release binary.
   **Fix: conda env with angsd=0.940, rebuild the sites index under it.** Worked.
   - Lesson: ANGSD `.idx`/`.bin` are version-keyed; always rebuild after a version change.

4. **Mito leaked into the sites file** and sorted to the top (NC_ < NW_), which is
   what made the VCF *look* mitochondrial. Strip NC_ — and it shouldn't be in a
   nuclear Mendelian analysis anyway.

5. **`.cmds` builder wrote only 1 line.** Paste mangled the `for … done` loop.
   Fix: build via a heredoc script (`cat > build_cmds.sh <<'EOF' … EOF; bash it`)
   so the terminal can't eat the loop. Always `wc -l angsd.cmds` (must = 6).

6. **Login-node etiquette.** bcftools/awk/python streaming = fine on login node.
   ANGSD, mapping, anything multi-core/long = idev or sbatch.