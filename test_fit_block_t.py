#!/usr/bin/env python3

"""
Toy end-to-end test for fit_block_t.py.

This creates synthetic parental-table and ANGSD-like files,
runs the real fit_block_t.py script, and checks whether known
transmission ratios are recovered.

Usage:
    python test_fit_block_t.py

Before running:
    Make sure fit_block_t.py is in the same directory.
"""

import gzip
import os
import subprocess
import tempfile

import numpy as np


# ================================================================
# CONFIGURATION
# ================================================================

SCRIPT = "fit_block_t.py"
POOL = "TEST-7x11-PL"

MINREADS = 10
MAXREADS = 100
MIN_SNPS = 5

RNG_SEED = 12345

rng = np.random.default_rng(RNG_SEED)


# ================================================================
# SIMULATION
# ================================================================

def beta_binomial_count(mu, n, rho):
    """Draw one beta-binomial count."""

    if rho <= 0:
        return rng.binomial(n, mu)

    s = (1.0 - rho) / rho

    alpha = mu * s
    beta = (1.0 - mu) * s

    p = rng.beta(alpha, beta)

    return rng.binomial(n, p)


def mu_from_config(t, config):
    """
    Calculate expected ALT frequency using the same parameterization
    as fit_block_t.py.

    Four configurations are tested:

      homref_hapAalt
      homalt_hapAalt
      homref_hapAref
      homalt_hapAref
    """

    if config == "homref_hapAalt":

        # Homozygous parent = REF/REF
        # Segregating parent hapA = ALT

        const = 0.0
        coefA = 0.5
        coefB = 0.0

    elif config == "homalt_hapAalt":

        # Homozygous parent = ALT/ALT
        # Segregating parent hapA = ALT

        const = 0.5
        coefA = 0.5
        coefB = 0.0

    elif config == "homref_hapAref":

        # Homozygous parent = REF/REF
        # Segregating parent hapA = REF

        const = 0.0
        coefA = 0.0
        coefB = 0.5

    elif config == "homalt_hapAref":

        # Homozygous parent = ALT/ALT
        # Segregating parent hapA = REF

        const = 0.5
        coefA = 0.0
        coefB = 0.5

    else:
        raise ValueError(config)

    return const + coefA * t + coefB * (1.0 - t)


# ================================================================
# MAKE FAKE parental_table.tsv
# ================================================================

def make_parental_table(path, blocks):

    header = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "klass",
        "block_id",
        "seg_parent",
        "hapA_allele",
        "exp_alt_freq",
    ]

    with open(path, "w") as out:

        out.write("\t".join(header) + "\n")

        for block in blocks:

            config = block["config"]

            if "hapAalt" in config:
                hapA = "G"
            else:
                hapA = "A"

            if "homalt" in config:
                exp_alt_freq = "0.75"
            else:
                exp_alt_freq = "0.25"

            for i in range(block["m"]):

                pos = block["start"] + i

                out.write(
                    "\t".join([
                        block["chrom"],
                        str(pos),
                        "A",
                        "G",
                        "het_x_hom",
                        block["block_id"],
                        block["seg"],
                        hapA,
                        exp_alt_freq,
                    ])
                    + "\n"
                )


# ================================================================
# MAKE FAKE ANGSD FILES
# ================================================================

def make_pool_files(pooldir, blocks, rho=0.05):

    pos_path = os.path.join(
        pooldir,
        f"freq_{POOL}.pos.gz"
    )

    count_path = os.path.join(
        pooldir,
        f"freq_{POOL}.counts.gz"
    )

    with gzip.open(pos_path, "wt") as pf, \
         gzip.open(count_path, "wt") as cf:

        # Headers.

        pf.write("chromo position\n")
        cf.write("A C G T\n")

        for block in blocks:

            t = block["t"]
            config = block["config"]

            for i in range(block["m"]):

                chrom = block["chrom"]
                pos = block["start"] + i

                # True expected ALT frequency.

                mu = mu_from_config(t, config)

                # Random depth between 35 and 50.

                n = int(rng.integers(35, 51))

                # Simulate ALT reads.

                alt_reads = beta_binomial_count(
                    mu,
                    n,
                    rho
                )

                ref_reads = n - alt_reads

                # REF=A
                # ALT=G
                #
                # ANGSD-like order:
                # A C G T

                pf.write(
                    f"{chrom}\t{pos}\n"
                )

                cf.write(
                    f"{ref_reads} 0 {alt_reads} 0\n"
                )

    return pos_path, count_path


# ================================================================
# MAIN TEST
# ================================================================

def main():

    print()
    print("=" * 72)
    print("TOY TEST OF fit_block_t.py")
    print("=" * 72)
    print()

    # ------------------------------------------------------------
    # Make sure production script exists.
    # ------------------------------------------------------------

    if not os.path.exists(SCRIPT):

        raise SystemExit(
            f"ERROR: Could not find {SCRIPT}\n\n"
            f"Put test_fit_block_t.py in the same directory as "
            f"fit_block_t.py."
        )

    # ------------------------------------------------------------
    # Define synthetic blocks.
    # ------------------------------------------------------------

    blocks = [

        # True t = 0.50
        {
            "chrom": "chr1",
            "start": 1000,
            "block_id": "B1",
            "seg": "PARENT1",
            "config": "homref_hapAalt",
            "t": 0.50,
            "m": 20,
        },

        # True t = 0.75
        {
            "chrom": "chr1",
            "start": 2000,
            "block_id": "B2",
            "seg": "PARENT1",
            "config": "homref_hapAalt",
            "t": 0.75,
            "m": 20,
        },

        # True t = 0.90
        #
        # This tests ALT/ALT homozygous parent.
        {
            "chrom": "chr2",
            "start": 3000,
            "block_id": "B3",
            "seg": "PARENT2",
            "config": "homalt_hapAalt",
            "t": 0.90,
            "m": 20,
        },

        # True t = 0.25
        #
        # This tests ALT/ALT homozygous parent AND
        # hapA = REF.
        {
            "chrom": "chr2",
            "start": 4000,
            "block_id": "B4",
            "seg": "PARENT2",
            "config": "homalt_hapAref",
            "t": 0.25,
            "m": 20,
        },
    ]

    # ------------------------------------------------------------
    # Create temporary fake input data.
    # ------------------------------------------------------------

    with tempfile.TemporaryDirectory(
        prefix="fit_block_t_test_"
    ) as tmp:

        parental = os.path.join(
            tmp,
            "parental_table.tsv"
        )

        pooldir = os.path.join(
            tmp,
            "poolsANGSD"
        )

        output = os.path.join(
            tmp,
            "block_t_perpool.tsv"
        )

        os.makedirs(pooldir)

        print(f"Temporary directory: {tmp}")
        print(f"Random seed:         {RNG_SEED}")
        print()

        # Create fake parental table.

        make_parental_table(
            parental,
            blocks
        )

        # Create fake ANGSD files.

        make_pool_files(
            pooldir,
            blocks,
            rho=0.05
        )

        # --------------------------------------------------------
        # Run the ACTUAL production script.
        # --------------------------------------------------------

        command = [
            "python",
            SCRIPT,

            "--parental",
            parental,

            "--pooldir",
            pooldir,

            "--pools",
            POOL,

            "--minreads",
            str(MINREADS),

            "--maxreads",
            str(MAXREADS),

            "--min-snps",
            str(MIN_SNPS),

            "--out",
            output,
        ]

        print("Running:")
        print(" ".join(command))
        print()

        result = subprocess.run(
            command,
            capture_output=True,
            text=True
        )

        print("Production-script STDERR:")
        print(result.stderr)

        if result.returncode != 0:

            print()
            print("Production-script STDOUT:")
            print(result.stdout)

            raise SystemExit(
                f"\nERROR: {SCRIPT} failed with "
                f"exit code {result.returncode}"
            )

        # --------------------------------------------------------
        # Read production output.
        # --------------------------------------------------------

        with open(output) as fh:

            lines = [
                line.rstrip("\n")
                for line in fh
            ]

        header = lines[0].split("\t")

        rows = []

        for line in lines[1:]:

            values = line.split("\t")

            rows.append(
                dict(zip(header, values))
            )

        # --------------------------------------------------------
        # Display results.
        # --------------------------------------------------------

        print("=" * 72)
        print("RESULTS")
        print("=" * 72)

        print()

        print(
            f"{'Block':<8}"
            f"{'True t':>10}"
            f"{'t_hat':>10}"
            f"{'rho_hat':>10}"
            f"{'LRT':>10}"
            f"{'p-value':>14}"
            f"{'Status':>10}"
        )

        print("-" * 72)

        all_ok = True

        for block, row in zip(blocks, rows):

            true_t = block["t"]

            t_hat = float(row["t_hat"])

            rho_hat = float(row["rho_hat"])

            lrt = float(row["LRT"])

            p = float(row["p_value"])

            error = abs(t_hat - true_t)

            # Toy simulation tolerance.
            #
            # With only 20 SNPs and rho=0.05, estimates will
            # have noticeable sampling noise.

            if error <= 0.15:
                status = "PASS"
            else:
                status = "FAIL"
                all_ok = False

            print(
                f"{block['block_id']:<8}"
                f"{true_t:>10.2f}"
                f"{t_hat:>10.4f}"
                f"{rho_hat:>10.4f}"
                f"{lrt:>10.3f}"
                f"{p:>14.3e}"
                f"{status:>10}"
            )

        # --------------------------------------------------------
        # Basic checks.
        # --------------------------------------------------------

        print()
        print("=" * 72)
        print("SANITY CHECKS")
        print("=" * 72)
        print()

        # Four blocks should be present.

        if len(rows) == 4:

            print("PASS: four blocks were recovered")

        else:

            print(
                f"FAIL: expected 4 blocks, "
                f"found {len(rows)}"
            )

            all_ok = False

        # Every block should contain 20 SNPs.

        if all(
            int(row["n_snps"]) == 20
            for row in rows
        ):

            print(
                "PASS: all blocks contain 20 usable SNPs"
            )

        else:

            print(
                "FAIL: unexpected SNP counts"
            )

            all_ok = False

        # B1 should estimate near 0.50.

        b1 = next(
            row for row in rows
            if row["block_id"] == "B1"
        )

        if abs(float(b1["t_hat"]) - 0.50) <= 0.15:

            print(
                "PASS: B1 recovered t ≈ 0.50"
            )

        else:

            print(
                "FAIL: B1 did not recover t ≈ 0.50"
            )

            all_ok = False

        # B2 should be > 0.50.

        b2 = next(
            row for row in rows
            if row["block_id"] == "B2"
        )

        if float(b2["t_hat"]) > 0.50:

            print(
                "PASS: B2 recovered direction t > 0.50"
            )

        else:

            print(
                "FAIL: B2 estimated t <= 0.50"
            )

            all_ok = False

        # B3 should be > 0.50.

        b3 = next(
            row for row in rows
            if row["block_id"] == "B3"
        )

        if float(b3["t_hat"]) > 0.50:

            print(
                "PASS: B3 recovered direction t > 0.50"
            )

        else:

            print(
                "FAIL: B3 estimated t <= 0.50"
            )

            all_ok = False

        # B4 should be < 0.50.

        b4 = next(
            row for row in rows
            if row["block_id"] == "B4"
        )

        if float(b4["t_hat"]) < 0.50:

            print(
                "PASS: B4 recovered direction t < 0.50"
            )

        else:

            print(
                "FAIL: B4 estimated t >= 0.50"
            )

            all_ok = False

        # --------------------------------------------------------
        # Final result.
        # --------------------------------------------------------

        print()
        print("=" * 72)

        if all_ok:

            print("OVERALL: PASS")

        else:

            print("OVERALL: FAIL")

        print("=" * 72)
        print()


if __name__ == "__main__":
    main()