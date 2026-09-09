#!/usr/bin/env python3

"""
Null-calibration test for fit_block_t.py

Purpose
-------
Test whether the LRT p-values from fit_block_t.py are approximately
calibrated when the true transmission ratio is t = 0.5.

This does NOT use real ANGSD files. Instead, it:

  1. Simulates beta-binomial SNP counts under t = 0.5.
  2. Fits each simulated block using the actual fit_block() function
     from fit_block_t.py.
  3. Collects the resulting LRT p-values.
  4. Checks the fraction below 0.05, 0.01, and 0.001.
  5. Also reports the distribution of t_hat and rho_hat.

Run:

    python test_null_calibration.py

Requires:
    numpy
    scipy

Important:
-----------
This imports fit_block() directly from fit_block_t.py, so make sure
the two files are in the same directory.
"""

import numpy as np

from fit_block_t import fit_block


# ================================================================
# CONFIGURATION
# ================================================================

SEED = 12345

N_BLOCKS = 1000

# Number of SNPs per simulated block.
#
# We deliberately vary this to make the test more realistic.
MIN_SNPS = 5
MAX_SNPS = 25

# Read depth range.
MIN_DEPTH = 10
MAX_DEPTH = 100

# Simulated beta-binomial overdispersion.
#
# We'll vary this somewhat between blocks.
MIN_RHO = 0.01
MAX_RHO = 0.10

# Tolerance used when judging calibration.
#
# With 1000 simulations:
#
#   expected 5% = 50
#
# Standard binomial sampling error is approximately +/- 1.4 percentage
# points, so don't expect exactly 5.000%.
#
# We flag values outside roughly a broad 95% range.
CALIBRATION_TOLERANCE = 0.015


rng = np.random.default_rng(SEED)


# ================================================================
# SIMULATION
# ================================================================

def beta_binomial_count(mu, n, rho):
    """
    Draw k from a beta-binomial distribution.

    Same mean/overdispersion parameterization as fit_block_t.py.
    """

    if rho <= 0:
        return rng.binomial(n, mu)

    s = (1.0 - rho) / rho

    alpha = mu * s
    beta = (1.0 - mu) * s

    p = rng.beta(alpha, beta)

    return rng.binomial(n, p)


def simulate_block(n_snps, rho):
    """
    Simulate one block under:

        TRUE t = 0.5

    using several parental configurations.

    At t=0.5 all configurations should have the expected Mendelian
    allele frequency:

        mu = 0.25   for REF/REF homozygous parent
        mu = 0.75   for ALT/ALT homozygous parent

    We randomly mix the four phase/orientation configurations.
    """

    snps = []

    for _ in range(n_snps):

        config = rng.integers(0, 4)

        if config == 0:

            # REF/REF + hapA=ALT
            const = 0.0
            coefA = 0.5
            coefB = 0.0

        elif config == 1:

            # ALT/ALT + hapA=ALT
            const = 0.5
            coefA = 0.5
            coefB = 0.0

        elif config == 2:

            # REF/REF + hapA=REF
            const = 0.0
            coefA = 0.0
            coefB = 0.5

        else:

            # ALT/ALT + hapA=REF
            const = 0.5
            coefA = 0.0
            coefB = 0.5

        # --------------------------------------------------------
        # True t = 0.5
        # --------------------------------------------------------

        t = 0.5

        mu = (
            const
            + coefA * t
            + coefB * (1.0 - t)
        )

        # Random depth.

        n = int(
            rng.integers(
                MIN_DEPTH,
                MAX_DEPTH + 1
            )
        )

        k = beta_binomial_count(
            mu,
            n,
            rho
        )

        snps.append(
            (
                const,
                coefA,
                coefB,
                k,
                n,
            )
        )

    return snps


# ================================================================
# BINOMIAL CONFIDENCE INTERVAL
# ================================================================

def binomial_approx_ci(p, n):
    """
    Simple approximate 95% CI for a proportion.
    Used only for reporting the simulation frequency.
    """

    se = np.sqrt(
        p * (1.0 - p) / n
    )

    return (
        max(0.0, p - 1.96 * se),
        min(1.0, p + 1.96 * se),
    )


# ================================================================
# MAIN
# ================================================================

def main():

    print()
    print("=" * 78)
    print("NULL CALIBRATION TEST FOR fit_block_t.py")
    print("=" * 78)
    print()

    print(f"Random seed:       {SEED}")
    print(f"Number of blocks:  {N_BLOCKS}")
    print(f"SNPs/block:        {MIN_SNPS}-{MAX_SNPS}")
    print(f"Depth/SNP:         {MIN_DEPTH}-{MAX_DEPTH}")
    print(f"rho range:         {MIN_RHO}-{MAX_RHO}")
    print()
    print("TRUE TRANSMISSION RATIO: t = 0.5")
    print()

    # ------------------------------------------------------------
    # Storage
    # ------------------------------------------------------------

    p_values = []

    t_hats = []

    rho_hats = []

    lrts = []

    n_snps_list = []

    # ------------------------------------------------------------
    # Run simulations
    # ------------------------------------------------------------

    for i in range(N_BLOCKS):

        # Random number of SNPs.
        n_snps = int(
            rng.integers(
                MIN_SNPS,
                MAX_SNPS + 1
            )
        )

        # Random true rho.
        rho_true = rng.uniform(
            MIN_RHO,
            MAX_RHO
        )

        snps = simulate_block(
            n_snps,
            rho_true
        )

        # --------------------------------------------------------
        # This is the ACTUAL fitting function from your production
        # script.
        # --------------------------------------------------------

        (
            t_hat,
            rho_hat,
            ll_full,
            ll_null,
            lrt,
            p,
            boundary,
        ) = fit_block(snps)

        p_values.append(p)
        t_hats.append(t_hat)
        rho_hats.append(rho_hat)
        lrts.append(lrt)
        n_snps_list.append(n_snps)

    p_values = np.asarray(p_values)
    t_hats = np.asarray(t_hats)
    rho_hats = np.asarray(rho_hats)
    lrts = np.asarray(lrts)

    # ============================================================
    # SUMMARY
    # ============================================================

    print("=" * 78)
    print("PARAMETER RECOVERY UNDER THE NULL")
    print("=" * 78)
    print()

    print(
        f"t_hat:"
        f"     mean   = {t_hats.mean():.4f}"
        f"     median = {np.median(t_hats):.4f}"
    )

    print(
        f"          10%-90% = "
        f"{np.quantile(t_hats, 0.10):.4f} - "
        f"{np.quantile(t_hats, 0.90):.4f}"
    )

    print()

    print(
        f"rho_hat:"
        f"     mean   = {rho_hats.mean():.4f}"
        f"     median = {np.median(rho_hats):.4f}"
    )

    print(
        f"          10%-90% = "
        f"{np.quantile(rho_hats, 0.10):.4f} - "
        f"{np.quantile(rho_hats, 0.90):.4f}"
    )

    print()

    print(
        f"LRT:"
        f"        mean   = {lrts.mean():.4f}"
        f"        median = {np.median(lrts):.4f}"
    )

    print()

    # ============================================================
    # P-VALUE CALIBRATION
    # ============================================================

    print("=" * 78)
    print("P-VALUE CALIBRATION")
    print("=" * 78)
    print()

    thresholds = [
        0.10,
        0.05,
        0.025,
        0.01,
        0.005,
        0.001,
    ]

    all_pass = True

    for threshold in thresholds:

        observed_fraction = np.mean(
            p_values < threshold
        )

        count = int(
            np.sum(p_values < threshold)
        )

        expected_fraction = threshold

        ci_lo, ci_hi = binomial_approx_ci(
            expected_fraction,
            N_BLOCKS
        )

        # Use a somewhat broader tolerance than the simple
        # binomial CI because this is a diagnostic rather than
        # a formal calibration test.

        tolerance = max(
            CALIBRATION_TOLERANCE,
            3.0 * np.sqrt(
                expected_fraction
                * (1.0 - expected_fraction)
                / N_BLOCKS
            )
        )

        lower = expected_fraction - tolerance
        upper = expected_fraction + tolerance

        if lower <= observed_fraction <= upper:

            status = "OK"

        else:

            status = "CHECK"
            all_pass = False

        print(
            f"p < {threshold:<6g} "
            f"observed = {count:4d}/{N_BLOCKS}"
            f" = {observed_fraction:7.4f}   "
            f"expected ≈ {expected_fraction:7.4f}   "
            f"{status}"
        )

    print()

    # ============================================================
    # P-VALUE QUANTILES
    # ============================================================

    print("=" * 78)
    print("P-VALUE DISTRIBUTION")
    print("=" * 78)
    print()

    print(
        f"min       = {p_values.min():.6g}"
    )

    print(
        f"1%        = {np.quantile(p_values, 0.01):.6g}"
    )

    print(
        f"5%        = {np.quantile(p_values, 0.05):.6g}"
    )

    print(
        f"10%       = {np.quantile(p_values, 0.10):.6g}"
    )

    print(
        f"25%       = {np.quantile(p_values, 0.25):.6g}"
    )

    print(
        f"50%       = {np.quantile(p_values, 0.50):.6g}"
    )

    print(
        f"75%       = {np.quantile(p_values, 0.75):.6g}"
    )

    print(
        f"90%       = {np.quantile(p_values, 0.90):.6g}"
    )

    print(
        f"95%       = {np.quantile(p_values, 0.95):.6g}"
    )

    print(
        f"99%       = {np.quantile(p_values, 0.99):.6g}"
    )

    print(
        f"max       = {p_values.max():.6g}"
    )

    print()

    # ============================================================
    # INTERPRETATION
    # ============================================================

    print("=" * 78)
    print("INTERPRETATION")
    print("=" * 78)
    print()

    fraction_05 = np.mean(
        p_values < 0.05
    )

    fraction_01 = np.mean(
        p_values < 0.01
    )

    print(
        f"At alpha=0.05: "
        f"{fraction_05:.3%} of true-null blocks rejected H0."
    )

    print(
        f"At alpha=0.01: "
        f"{fraction_01:.3%} of true-null blocks rejected H0."
    )

    print()

    if all_pass:

        print(
            "OVERALL: PASS"
        )

        print()
        print(
            "The simulated null p-values are reasonably "
            "consistent with the nominal chi-square calibration."
        )

    else:

        print(
            "OVERALL: CHECK"
        )

        print()
        print(
            "The simulated null p-values show some deviation "
            "from nominal calibration."
        )

        print(
            "This does NOT automatically mean the model is wrong."
        )

        print(
            "With small blocks, estimated rho, and beta-binomial "
            "sampling, finite-sample deviations are possible."
        )

        print(
            "If this happens, the next step should be empirical "
            "null calibration / parametric bootstrap rather than "
            "using the raw chi-square p-values."
        )

    print()
    print("=" * 78)


if __name__ == "__main__":
    main()