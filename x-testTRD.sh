idev -p skx-dev -N 1 -n 1 -t 02:00:00 -A IBN21018

python test_fit_block_t.py

# Temporary directory: /tmp/fit_block_t_test_gt0van89
Random seed:         12345

Running:
python fit_block_t.py --parental /tmp/fit_block_t_test_gt0van89/parental_table.tsv --pooldir /tmp/fit_block_t_test_gt0van89/poolsANGSD --pools TEST-7x11-PL --minreads 10 --maxreads 100 --min-snps 5 --out /tmp/fit_block_t_test_gt0van89/block_t_perpool.tsv

Production-script STDERR:
[info] block-usable parental sites: 80
[info] TEST-7x11-PL: fit 4 blocks

========================================================================
RESULTS
========================================================================

Block       True t     t_hat   rho_hat       LRT       p-value    Status
------------------------------------------------------------------------
B1            0.50    0.4230    0.0655     1.712     1.908e-01      PASS
B2            0.75    0.8417    0.0675    21.681     3.220e-06      PASS
B3            0.90    0.8948    0.0508    20.444     6.141e-06      PASS
B4            0.25    0.2077    0.0472    17.057     3.628e-05      PASS

========================================================================
SANITY CHECKS
========================================================================

PASS: four blocks were recovered
PASS: all blocks contain 20 usable SNPs
PASS: B1 recovered t ≈ 0.50
PASS: B2 recovered direction t > 0.50
PASS: B3 recovered direction t > 0.50
PASS: B4 recovered direction t < 0.50

========================================================================
OVERALL: PASS
========================================================================
