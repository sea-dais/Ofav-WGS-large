# Understanding `fit_block_t.py` — and learning Python along the way

This guide walks through the script top to bottom. Each section explains **what the code does**, **why it's written that way**, and **which Python concept** you can take away and reuse. No prior stats knowledge is assumed — the biology/statistics is summarized just enough to follow the code.

---

## 0. The 30-second summary of what the whole script does

You have DNA read-count data from several "pools" of organisms. For each genomic **block**, the script asks: *did one parent pass on its two versions (haplotypes) of a chromosome region 50/50, or was there a bias?* That bias is a single number `t` (a fraction between 0 and 1). The script estimates `t` for every block and runs a statistical test on whether `t` differs from `0.5`. Output is one big table.

You don't need to understand the genetics to learn Python from it. Think of it as: *read two files, do math per group, write a results table.*

---

## 1. The very top: the shebang and the docstring

```python
#!/usr/bin/env python3
"""
fit_block_t.py  --  Stage A of the block-level TRD test.
...
"""
```

**What's happening:**
- The first line (`#!/usr/bin/env python3`) is a *shebang*. On Linux/Mac it lets you run the file directly (`./fit_block_t.py`) and have the system pick Python 3. It's ignored when you run `python3 fit_block_t.py`.
- The triple-quoted `"""..."""` block right after is a **module docstring**. It's a string that isn't assigned to anything, so Python just attaches it as documentation to the file. Tools and `help()` can read it.

**Python takeaway:** any string literal placed as the *first statement* of a file, function, or class becomes its docstring. Triple quotes let a string span multiple lines.

---

## 2. Imports

```python
import argparse, gzip, sys, math
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import chi2
from scipy.special import betaln
```

**What each does:**
- `argparse` — builds a command-line interface (the `--parental`, `--out` flags you see later).
- `gzip` — reads `.gz` compressed files without unzipping them first.
- `sys` — access to things like `sys.stderr` (error/log stream) and `sys.exit`.
- `math` — basic math functions (`log`, `lgamma`, ...).
- `numpy as np` — fast numeric arrays; `as np` gives it a short alias.
- The three `from ... import ...` lines pull **specific names** out of SciPy so you can write `chi2` instead of `scipy.stats.chi2`.

**Python takeaway:**
- `import X` gives you `X.thing`.
- `import X as Y` renames it locally (`np` is a universal convention for numpy).
- `from X import thing` pulls one name directly into your namespace.

---

## 3. Module-level constants

```python
BASE = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
BOUNDARY_EPS = 1e-3
```

**What's happening:**
- `BASE` is a **dictionary** mapping each DNA base letter to a column index. Later, `cnt[BASE['A']]` grabs the A-count. Dictionaries are Python's key→value lookup table.
- `BOUNDARY_EPS` is a tiny number (`0.001`). `1e-3` is scientific notation for `0.001`. It's used to decide "is `t` basically at 0 or 1?".

**Python takeaway:** ALL-CAPS names are a *convention* meaning "constant, don't reassign." Python won't stop you, but other programmers read the caps as "leave this alone." `1e-3`, `1e-6` etc. are float literals.

---

## 4. Reading the parental table — `read_parental`

```python
def read_parental(path):
    """(chrom,pos) -> SNP dict, restricted to phased het_x_hom / hom_x_het sites."""
    sites = {}
    with open(path) as fh:
        for i, line in enumerate(fh):
            f = line.rstrip("\n").split("\t")
            if i == 0:
                H = {name: j for j, name in enumerate(f)}
                continue
            ...
```

This is the first real function. Let's unpack the Python patterns, because several of them repeat throughout the file.

### 4a. `def` and the return value
`def read_parental(path):` defines a function taking one argument `path`. It builds a dictionary `sites` and (at the end) `return sites`. A function that doesn't hit a `return` gives back `None`.

### 4b. `with open(path) as fh:`
This is a **context manager**. It opens the file, gives it the name `fh`, and — crucially — *automatically closes it* when the block ends, even if an error happens. Always prefer this over a bare `open()`.

### 4c. `for i, line in enumerate(fh):`
Iterating a file yields one line at a time. `enumerate(...)` wraps that so you also get a counter: `i=0, 1, 2, ...`. So `i == 0` means "this is the header row."

**Python takeaway — `enumerate`:** whenever you catch yourself wanting a manual counter (`i = 0; i += 1`), reach for `enumerate` instead.

### 4d. Parsing a TSV line
```python
f = line.rstrip("\n").split("\t")
```
- `.rstrip("\n")` strips the trailing newline.
- `.split("\t")` splits on tabs into a **list** of column strings.

So `f` is a list like `["chr1", "1200", "A", "G", ...]`.

### 4e. A dict comprehension for the header
```python
H = {name: j for j, name in enumerate(f)}
```
This builds a dictionary in one line: column-name → column-index. So later `f[H["POS"]]` means "the value in the POS column of this row," without hard-coding that POS is column 4. This makes the code robust to column reordering.

**Python takeaway — comprehensions:** `{k: v for ... in ...}` builds a dict, `[x for ... in ...]` builds a list, `(x for ... in ...)` builds a lazy generator. They're compact loops that produce a collection.

### 4f. `continue` to skip rows
```python
if klass not in ("het_x_hom", "hom_x_het"):
    continue
```
`continue` jumps to the next loop iteration, skipping the rest of the body. Here it filters out row types the script doesn't want. `x not in (a, b)` is a clean membership test against a tuple.

### 4g. Type conversion and building the record
```python
chrom = f[H["CHROM"]]; pos = int(f[H["POS"]])
...
exp = float(f[H["exp_alt_freq"]])
```
Everything read from a text file is a **string**. `int(...)` and `float(...)` convert to numbers. Forgetting this is one of the most common beginner bugs (`"12" + 1` errors; `12 + 1` works).

### 4h. f-strings
```python
block=f"{chrom}:{seg}:{block_id}",
```
An **f-string** (`f"..."`) lets you drop variables straight into a string inside `{...}`. This builds a block key like `chr1:P7:3`. The comment explains *why* chrom is included: the same PS block-id number can reappear on different scaffolds, so chrom disambiguates it.

### 4i. Returning a dict of dicts
`sites[(chrom, pos)] = dict(...)` uses a **tuple** `(chrom, pos)` as a dictionary key (tuples are allowed as keys because they're immutable; lists are not). The value is itself a dict built with the `dict(key=value, ...)` constructor. So `sites` is a dict keyed by genomic position, each value a small record.

**Python takeaway:** tuples can be dict keys; lists cannot. `dict(a=1, b=2)` is an alternative to `{"a": 1, "b": 2}`.

---

## 5. A tiny helper — `_count_lines_gz`

```python
def _count_lines_gz(path):
    n = 0
    with gzip.open(path, "rt") as fh:
        for _ in fh:
            n += 1
    return n
```

Counts lines in a gzipped file.

- The leading underscore in `_count_lines_gz` is a **convention** meaning "internal/private helper — not part of the public interface." Python doesn't enforce privacy; it's a signal to readers.
- `gzip.open(path, "rt")` opens a `.gz` file in **text** mode (`"rt"` = read text), so you get strings, not raw bytes.
- `for _ in fh:` — the underscore `_` is the conventional name for "a variable I must name but don't care about." Here we only want the count, not the line.

**Python takeaway:** `_` as a throwaway variable; leading `_name` as "private by convention."

---

## 6. Reading pool counts and zipping two files — `read_pool_counts`

```python
def read_pool_counts(pooldir, pool, sites, minreads, maxreads):
    posf = f"{pooldir}/freq_{pool}.pos.gz"
    cntf = f"{pooldir}/freq_{pool}.counts.gz"

    npos = _count_lines_gz(posf)
    ncnt = _count_lines_gz(cntf)
    if npos != ncnt:
        sys.exit(f"[error] {pool}: pos.gz ({npos} lines) and counts.gz ({ncnt} lines) "
                 f"differ; files are out of sync, refusing to zip.")
```

### 6a. Defensive programming with `sys.exit`
Two parallel files (positions and counts) must have the same number of lines, or zipping them line-by-line would silently pair up the wrong rows. The code counts both and calls `sys.exit(message)` to **abort with an error message** if they differ. Passing a string to `sys.exit` prints it to stderr and exits with a non-zero (failure) status.

**Python takeaway:** validate assumptions early and fail loudly. A silent wrong answer is worse than a crash.

### 6b. Reading two files in lockstep with `zip`
```python
with gzip.open(posf, "rt") as pf, gzip.open(cntf, "rt") as cf:
    pf.readline(); cf.readline()          # skip both headers
    for pline, cline in zip(pf, cf):
        ...
```
- You can open **multiple context managers** in one `with` statement, comma-separated.
- `zip(pf, cf)` walks both files together, giving `(line_from_pos, line_from_counts)` pairs. `zip` stops at the shorter one — which is exactly why the length check above matters.

**Python takeaway — `zip`:** pairs up items from two (or more) iterables position-by-position. Great for "walk these together."

### 6c. `dict.get` with a `None` check
```python
s = sites.get((chrom, pos))
if s is None:
    continue
```
`sites.get(key)` returns the value or `None` if the key is missing — unlike `sites[key]`, which would raise `KeyError`. This is the idiomatic "look it up, skip if absent" pattern. Note `is None`, not `== None`: `is` tests identity and is the correct way to check for `None`.

### 6d. List comprehension + indexing into counts
```python
cnt = [int(x) for x in cc]          # totA totC totG totT
refn = cnt[BASE[s["ref"]]]
altn = cnt[BASE[s["alt"]]]
n = refn + altn
```
`cnt` becomes a list of four integers. `BASE[s["ref"]]` turns the REF base letter into its column index (recall the `BASE` dict), and `cnt[...]` pulls that count. This is why `BASE` existed.

### 6e. Range filtering
```python
if n < minreads or n > maxreads:
    continue
```
Sites with too few reads (noisy) or too many (likely repeat-collapse artifacts) are dropped.

### 6f. `dict.setdefault` to group items
```python
b = blocks.setdefault(s["block"], {"seg": s["seg"], "snps": []})
b["snps"].append((const, coefA, coefB, altn, n))
```
This is a **grouping** idiom. `setdefault(key, default)`:
- if `key` already exists, returns its current value;
- if not, inserts `default` and returns that.

So the first SNP for a block creates the block's record `{"seg": ..., "snps": []}`, and every SNP (including the first) then appends its numbers into that block's `snps` list. The result: `blocks` maps each block → all its SNP tuples.

**Python takeaway — `setdefault`:** a one-liner for "group things into buckets." (`collections.defaultdict(list)` is the other common tool for the same job.)

---

## 7. The statistics functions — reading them as *code*, not math

You don't need the statistics to learn the Python. Here's the code-level story.

### 7a. `_bb_logpmf` — one formula, guarded
```python
def _bb_logpmf(k, n, mu, rho):
    mu = min(max(mu, 1e-6), 1 - 1e-6)
    if rho <= 1e-9:
        return (... binomial formula ...)
    s = (1 - rho) / rho
    a = mu * s
    b = (1 - mu) * s
    return (... beta-binomial formula ...)
```
This computes the log-probability of seeing `k` "alt" reads out of `n`, given a predicted fraction `mu` and an overdispersion `rho`.

Python things worth noticing:
- **Clamping:** `mu = min(max(mu, 1e-6), 1 - 1e-6)` forces `mu` into the open interval `(0, 1)`. `max(mu, 1e-6)` pushes it up off zero; the outer `min(..., 1 - 1e-6)` caps it below one. This avoids `log(0)` blowing up. This min/max clamp idiom appears several times in the file.
- **An early `return` as a branch:** if `rho` is essentially zero, it returns the simpler binomial formula and never reaches the rest. Returning early is a clean alternative to a big `if/else`.
- Parentheses let a single expression span multiple lines — no backslashes needed.

### 7b. `_neg_ll` — sum over the block
```python
def _neg_ll(t, rho, snps):
    ll = 0.0
    for const, coefA, coefB, k, n in snps:
        mu = const + coefA * t + coefB * (1 - t)
        ll += _bb_logpmf(k, n, mu, rho)
    return -ll
```
- **Tuple unpacking in the loop:** `for const, coefA, coefB, k, n in snps:` — each item of `snps` is a 5-tuple, and Python unpacks it into five named variables in one step. Much clearer than `item[0]`, `item[1]`, ...
- It accumulates a total log-likelihood `ll`, then returns `-ll`. Optimizers *minimize*, so returning the **negative** log-likelihood turns "find the maximum likelihood" into "find the minimum." That's why the name starts with `_neg`.

### 7c. `_best_rho` — a lambda and a bounded optimizer
```python
def _best_rho(t, snps, rho0=0.05):
    res = minimize_scalar(lambda r: _neg_ll(t, min(max(r, 1e-6), 0.5), snps),
                          bounds=(1e-6, 0.5), method="bounded",
                          options={"xatol": 1e-5})
    return min(max(res.x, 1e-6), 0.5), res.fun
```
- `rho0=0.05` is a **default argument** — callers may omit it. (Here it's actually unused inside, a harmless leftover — a good example that reading real code means spotting such things.)
- `lambda r: ...` is an **anonymous, one-line function** of `r`. `minimize_scalar` needs "a function of one variable to minimize," and this lambda fixes `t` and `snps` while varying only `r` (rho). Building a small function on the fly like this is extremely common.
- `minimize_scalar(..., bounds=(1e-6, 0.5), method="bounded")` searches for the rho in that range giving the lowest value.
- **Returning two things:** `return best_rho, res.fun` returns a *tuple*. The caller unpacks it with `rho, negll = _best_rho(...)`.

**Python takeaway:** functions return one object, but a comma makes that object a tuple, so you can effectively return several values and unpack them at the call site.

### 7d. `fit_block` — grid search then refine
```python
def fit_block(snps, grid=101):
    ts = np.linspace(0.0, 1.0, grid)
    best_t, best_rho, best_negll = 0.5, 0.05, float("inf")
    for t in ts:
        rho, negll = _best_rho(t, snps)
        if negll < best_negll:
            best_negll, best_t, best_rho = negll, t, rho
    ...
```
- `np.linspace(0, 1, grid)` makes an array of 101 evenly spaced values from 0 to 1 — the candidate `t` values to try.
- `best_negll = float("inf")` initializes the "best score so far" to **positive infinity**, so the first real result always beats it. A standard trick for "find the minimum" loops.
- **Multiple assignment:** `best_negll, best_t, best_rho = negll, t, rho` updates three variables at once. The right side is packed into a tuple, then unpacked into the three names.

The rest of the function does a local refinement around the best grid point, computes a null model (fixing `t = 0.5`), forms a likelihood-ratio test statistic, and returns a **7-tuple** of results. Note:
```python
lrt = max(2.0 * (ll_full - ll_null), 0.0)
```
`max(..., 0.0)` clamps the statistic to be non-negative (numerical noise could otherwise make it slightly negative).
```python
boundary = 1 if (t_hat < BOUNDARY_EPS or t_hat > 1 - BOUNDARY_EPS) else 0
```
This is a **conditional (ternary) expression**: `value_if_true if condition else value_if_false`. It flags whether the estimate landed at the edge (0 or 1), where the usual test math is unreliable.

### 7e. `switch_flag` — counting sign changes
```python
def switch_flag(snps, t_hat):
    signs = []
    for const, coefA, coefB, k, n in snps:
        mu = const + coefA * t_hat + coefB * (1 - t_hat)
        obs = k / n if n else mu
        signs.append(1 if obs > mu else -1)
    return sum(1 for a, b in zip(signs, signs[1:]) if a != b)
```
- `obs = k / n if n else mu` — another ternary. `if n` is truthy when `n != 0`; if `n` were 0 it falls back to `mu` to avoid dividing by zero. (In Python, `0` is "falsy"; nonzero numbers are "truthy.")
- The final line is dense but elegant. `zip(signs, signs[1:])` pairs each element with its **next** neighbor: `(signs[0], signs[1]), (signs[1], signs[2]), ...`. `signs[1:]` is a **slice** — the list from index 1 onward. Then `sum(1 for a, b in ... if a != b)` counts how many adjacent pairs differ, i.e. how many times the sign flips.

**Python takeaway:** `zip(seq, seq[1:])` is the standard idiom for "iterate over consecutive pairs." Slicing `seq[1:]` copies from index 1 to the end.

---

## 8. `main` and the CLI — `argparse`

```python
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parental", required=True)
    ap.add_argument("--pooldir", required=True)
    ap.add_argument("--pools", required=True, help="comma list, e.g. ...")
    ap.add_argument("--minreads", type=int, default=10)
    ap.add_argument("--maxreads", type=int, default=100)
    ap.add_argument("--min-snps", type=int, default=5)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
```

**What's happening:**
- `argparse` turns command-line text into a tidy `args` object. You run the script like:
  ```
  python3 fit_block_t.py --parental parental_table.tsv --pooldir poolsANGSD \
      --pools P1-7x11-PL,P2-7x11-PL --out block_t_perpool.tsv
  ```
- `required=True` means the flag must be supplied.
- `type=int` auto-converts that argument to an integer.
- `default=10` supplies a value when the flag is omitted.
- After `parse_args()`, you access values as `args.parental`, `args.minreads`, etc. Note `--min-snps` becomes `args.min_snps` (dashes become underscores).

**Python takeaway:** `argparse` is the standard-library way to make a real command-line tool. Define arguments, call `parse_args()`, read fields off the result.

### 8a. Writing the output table
```python
with open(args.out, "w") as out:
    out.write("\t".join([...header columns...]) + "\n")
    for pool in args.pools.split(","):
        direction = "7x11" if "7x11" in pool else ("11x7" if "11x7" in pool else "NA")
        blocks = read_pool_counts(...)
        for block_key, d in blocks.items():
            snps = d["snps"]
            if len(snps) < args.min_snps:
                continue
            ...
            out.write("\t".join(str(x) for x in [...]) + "\n")
```
- `open(args.out, "w")` opens for **writing** (`"w"` truncates/creates the file).
- `"\t".join(list_of_strings)` joins a list into one tab-separated string — the inverse of `.split("\t")`. It's how you build a TSV row.
- `str(x) for x in [...]` converts each value to a string first, because `join` only accepts strings. This is a **generator expression** passed straight into `join`.
- Nested ternary picks the direction label from the pool name.
- `for block_key, d in blocks.items():` iterates a dict's **key–value pairs**. `.items()` yields `(key, value)` tuples, unpacked into `block_key` and `d`.
- Formatting: `f"{t_hat:.4f}"` formats a float to 4 decimals; `f"{p:.3e}"` uses scientific notation with 3 digits. Format specs after the colon control the display.

**Python takeaway:**
- `.join()` / `.split()` are the pair you use constantly for delimited text.
- `.items()`, `.keys()`, `.values()` are how you loop over dictionaries.
- `f"{value:.4f}"`-style format specs are worth memorizing: `.4f` = fixed 4 decimals, `.3e` = 3-digit scientific, `,d` = thousands separators, etc.

---

## 9. The `if __name__ == "__main__":` guard

```python
if __name__ == "__main__":
    main()
```

**What's happening:** when you *run* a file directly, Python sets the special variable `__name__` to the string `"__main__"`. When you *import* that file as a module from another script, `__name__` is instead the module's name. So this line means: *"only call `main()` if this file was run directly, not when it's imported."*

This lets you reuse the functions (`import fit_block_t; fit_block_t.fit_block(...)`) without the whole script executing.

**Python takeaway:** this is one of the most important idioms in Python. Put runnable behavior behind this guard, keep reusable functions above it.

---

## 10. Quick-reference: the Python concepts this file teaches

| Concept | Where it appears | One-line reminder |
|---|---|---|
| Docstrings | top of file, each function | first string in a file/function = its documentation |
| `import` styles | imports block | `import x`, `import x as y`, `from x import z` |
| Dictionaries | `BASE`, `sites`, `blocks` | key → value lookups |
| Tuples as keys | `sites[(chrom, pos)]` | immutable, so allowed as dict keys |
| `with` (context managers) | every file open | auto-closes the file |
| `enumerate` | `read_parental` loop | loop with a counter |
| `.split()` / `.join()` | TSV parse & write | text ↔ list of fields |
| Comprehensions | `H = {...}`, `cnt = [...]` | build a collection in one expression |
| `continue` | filters | skip to next iteration |
| Type conversion | `int(...)`, `float(...)` | file text is always strings |
| f-strings + format specs | keys, output | `f"{x:.4f}"` |
| `dict.get` / `is None` | pool reading | safe lookup, correct None check |
| `dict.setdefault` | grouping SNPs | bucket items by key |
| `zip` | reading 2 files, sign changes | walk iterables in parallel |
| Tuple unpacking | loops & returns | `a, b = func()` |
| Returning multiple values | `_best_rho`, `fit_block` | comma makes a tuple |
| `lambda` | `minimize_scalar` calls | tiny inline function |
| Ternary expression | `boundary`, `direction` | `A if cond else B` |
| min/max clamping | `mu`, `rho`, `lrt` | keep a value inside bounds |
| Slicing | `signs[1:]` | sub-sequence from an index |
| `sys.exit` / `sys.stderr` | error handling, logging | fail loudly, log to stderr |
| `argparse` | `main` | build a command-line interface |
| `__name__ == "__main__"` | bottom | run-vs-import guard |

---

## 11. Suggested way to learn from this file

1. **Run it in your head, then for real.** Copy a few lines of a function into a Python REPL (`python3` in a terminal) with fake inputs and watch what they do. E.g. try `"a\tb\tc".split("\t")` and `"\t".join(["x", "y"])`.
2. **Poke the small pieces.** Paste `BASE = {'A':0,'C':1,'G':2,'T':3}` and evaluate `BASE['G']`. Then `cnt = [10, 2, 5, 3]; cnt[BASE['G']]`.
3. **Rewrite one idiom the long way.** Take the dict comprehension `H = {name: j for j, name in enumerate(f)}` and rewrite it as an explicit `for` loop with `H = {}` and `H[name] = j`. Confirm they're equal. Doing this for `setdefault` and `zip(signs, signs[1:])` too will lock in those idioms.
4. **Trace the data shapes.** Sketch what `sites`, `blocks`, and each SNP tuple look like. Most confusion in real code is "what shape is this variable right now?" — answering that repeatedly is the core skill.

---

*This guide focuses on the Python. The statistics (beta-binomial likelihood, likelihood-ratio test, overdispersion) are described only enough to follow the code; if you want, that layer can be explained separately.*