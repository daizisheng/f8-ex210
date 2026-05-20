# TAOCP Pre-Fascicle 8A, Exercise 210 — empirical investigation for m=5

> "Prove or disprove that Q⁺ₘ(z) is a multiple of Qₘ(z)³ when m ≥ 5."
> — Donald Knuth, *TAOCP* Pre-Fascicle 8A, Exercise 210 [HM46]

## What this repository finds

Running the tools in this repository on Knuth's own DYNAHAM as the
source of truth produces, mod p = 10⁹ + 7:

```
deg Q_5(z)         = 8,212        (closed knight tours on 5×n)
deg Q⁺_5(z)        = 39,630       (open knight paths on 5×n)
Q_5(z)^1 divides Q⁺_5(z):   DIVISIBLE
Q_5(z)^2 divides Q⁺_5(z):   DIVISIBLE
Q_5(z)^3 divides Q⁺_5(z):   NOT DIVISIBLE   (remainder degree 24,635)
```

A NOT DIVISIBLE verdict over F_p cannot be a false positive: if
Q^3 divided Q+ over the integers it would also divide it modulo
any prime that does not divide leading coefficients, so a nonzero
remainder mod p implies a nonzero remainder over ℤ as well. The
only remaining concern is that Berlekamp–Massey might have
returned a degenerate divisor of the true minimum polynomial at a
"bad" prime (probability at most deg/p, roughly 10⁻⁵ for our
inputs).

The repository contains a small tool, `verify_bm`, that detects
this: it applies the BM polynomial as a linear recurrence to the
exact integer values produced by Knuth's bignum `dynaham` (which
the `./tools/build_truth.sh` script obtains by downloading and
running Knuth's source unmodified) and checks that every
testable position evaluates to zero modulo p. If even one
position fails, BM lost a factor. For our run, all 40 365
testable positions evaluate to zero, confirming the BM polynomial
is exactly the integer minimum polynomial reduced mod p — not a
degenerate factor.

The 24,635-degree remainder is not a near-miss but a fully
populated nonzero polynomial.

For m = 5, Exercise 210's conjectured cube relation does not hold;
the correct relation appears to be a *square*.

## Quick verification (no external libraries)

If you only want to confirm the divisibility verdict from the
shipped polynomials:

```
make divtest
./divtest data/Q5_closed.txt data/Q5_open.txt 1   # DIVISIBLE
./divtest data/Q5_closed.txt data/Q5_open.txt 2   # DIVISIBLE
./divtest data/Q5_closed.txt data/Q5_open.txt 3   # NOT DIVISIBLE
```

`divtest` is a self-contained 180-line C program; no third-party
libraries are required. Total time: under 1 second.

## Full pipeline from scratch (~3 minutes on 16 cores)

```
./tools/install_sgb.sh        # download + build Stanford GraphBase locally
./tools/build_patched.sh      # download Knuth's dynaham.w + apply patches + build
make                          # build our tools
./run.sh                      # full pipeline
```

No system packages are required; `install_sgb.sh` builds Stanford
GraphBase into `ext/sgb/` and the Makefile picks it up automatically.

`run.sh` performs:

1. `make_knight 5 30 k5x30.gb` — build the 5×30 knight graph (1 second).
2. `./dynaham_modp_dumpT k5x30.gb` and `./dynahamp_modp_dumpT k5x30.gb` — Knuth's DYNAHAM (closed + open), patched to dump the periodic transfer-matrix data to stderr. About 1 minute combined.
3. `./simulate_count ...` — replay the dumped period mod p to produce 40,000 values of S⁺_{5,n}, then Berlekamp–Massey to recover the minimum polynomial. About 75 seconds on 16 cores.
4. Convert the closed u-polynomial to z-space (one-line `awk`).
5. `./divtest` for k = 1, 2, 3 — report the three verdicts. A few milliseconds.

## Independent ground-truth cross-check

If you want to convince yourself the simulator is producing the
same values DYNAHAM would produce, run

```
./tools/build_truth.sh
```

This downloads Knuth's `dynaham.w` and `dynahamp.w` from
`cs.stanford.edu/~knuth/programs/`, ctangles them *unmodified*,
builds the bignum versions, and runs each on a 5×10 board to print
exact integer values of S_{5,n} and S⁺_{5,n} for n = 2..10.

These values can be reduced mod 10⁹+7 and compared with the first
few iterations of `simulate_count`; they match.

## Repository layout

```
src/
  simulate_count.c           Method-3 simulator (~700 lines)
  divtest.c                  self-contained polynomial divisibility (~180 lines)
  make_knight.c              5-line wrapper around SGB's board()
  dynaham_modp_dumpT.c       included for convenience; identical to
  dynahamp_modp_dumpT.c      ./tools/build_patched.sh's output

patches/
  dynaham.patch              unified diff applied to ctangle(dynaham.w)
  dynahamp.patch             (open variant)

data/
  Q5_closed.txt              Q_5(z) mod 10⁹+7, degree 8,212
  Q5_open.txt                Q⁺_5(z) mod 10⁹+7, degree 39,630

tools/
  fetch_knuth.sh             curl Knuth's dynaham*.w into ext/
  build_truth.sh             ctangle + build Knuth's UNPATCHED dynaham, run on 5×10
  build_patched.sh           ctangle + apply patches + build patched binaries
```

## Method summary (one paragraph)

DYNAHAM, when sweeping its sliding window through the interior of an
infinite m × ∞ knight strip, produces a sequence of per-vertex
transition operators that becomes strictly periodic with period 2m
m-steps once the window is far from both column boundaries. We dump
one period (10 m-steps for m = 5) plus the state vector at one
boundary; `simulate_count` then iterates the resulting period
operator in mod-p arithmetic to generate the same count[m,n] values
DYNAHAM would have produced on a much larger board, at about
1.5 ms per outer iteration. Berlekamp–Massey on the resulting
sequence recovers the minimum polynomial. `divtest` performs the
divisibility test with schoolbook polynomial arithmetic mod p.

## Key timings (16-core Xeon Platinum 8375C, AVX-512)

| Stage                                  | Wall time |
|----------------------------------------|-----------|
| `make_knight 5 30`                     | ~1 s |
| `dynaham_modp_dumpT k5x30.gb` (closed) | ~30 s |
| `dynahamp_modp_dumpT k5x30.gb` (open)  | ~60 s |
| `simulate_count` closed, 12,000 iters  | ~3 s |
| `simulate_count` open, 40,000 iters    | ~75 s |
| `divtest ... 3`                        | <0.1 s |
| **total**                              | **~3 minutes** |

For comparison: a full DYNAHAM run on a 5×80,000 board in mod-p
arithmetic takes about 30 hours on the same machine. The
periodicity-and-iterate approach is approximately 600× faster while
producing the same minimum polynomial.

## Matrix structure (informational)

| | closed | open |
|---|---|---|
| stable state-space dim (per phase) | ~17,884 / 17,268 | ~143,448 / 135,913 |
| nnz per period (10 matrices) | ~0.65 M | ~6.83 M |
| matrix-vector multiply (16-core, mod-p) | 0.2 ms | 1.4 ms |
| all entries 0 or 1 | yes | yes |

The transition matrices are *binary*: every nonzero entry is 1.
This is why `bin_csc_apply` in `simulate_count.c` has no
multiplications, only additions mod p.

## Dependencies

- C99 compiler (gcc or clang) and `make`
- `ctangle` from CWEB
  (Linux: `apt install cweb`; macOS: `brew install cweb`)
- An OpenMP-capable compiler (optional; without OpenMP the simulator
  runs single-threaded, about 8× slower)

Stanford GraphBase is fetched and built by `./tools/install_sgb.sh`;
no system package is needed.

macOS notes:

- Apple's bundled clang does not include OpenMP. Install it with
  `brew install libomp`, then run
  `make OMP="-Xpreprocessor -fopenmp -lomp"`.
- The Makefile auto-detects Linux x86_64 vs everything else and
  applies `-mcmodel=large` only where it is meaningful.

## On the choice of prime

The default modulus is p = 10⁹ + 7, a 30-bit prime in common use
in competitive programming and number theory. It is large enough
to make accidental collisions astronomically unlikely yet small
enough that `uint32_t` storage and `uint64_t` accumulation are
sufficient. Modern C compilers turn `% p` into a multiplication
and shift since p is known at compile time.

`simulate_count` and `divtest` both honour a `-DMODP=` override.
For a defensive cross-check against the bad-prime scenario noted
above, rebuild with `-DMODP=1000000009ULL` (the next prime) and
rerun the pipeline; both verdicts should agree.

We measured a Mersenne candidate, 2³¹ − 1 = 2147483647, on the
same hardware as a curiosity. It is actually *slower* than
10⁹ + 7 in both the compiler-generic case and the hand-coded case
(see below); the choice is not material to the verdicts either way.

| operation                          | p = 10⁹+7 | p = 2³¹−1 hand | slowdown |
|------------------------------------|-----------|----------------|----------|
| simulate_count OPEN 1000 iters (1 thread) | 10.19 ms/iter | 10.45 ms/iter | 2.5%   |
| simulate_count CLOSED 5000 iters         | 1.08 ms/iter  | 1.18 ms/iter | 9%     |
| divtest Q³ ÷ Q⁺                          | 0.66 s        | 1.09 s       | 65%    |

GCC compiles `% 1000000007` into a four-instruction Barrett
reduction (one wide `mul`, one shift, one narrow `mul`, one
`sub`), branch-free. The hand-written Mersenne fast reduction
needs two rounds of `(lo + hi)` plus a final correction, which
adds up to seven instructions and a conditional branch. On this
Ice Lake Xeon the compiler's Barrett wins outright.

## Author

Shisheng Li, with implementation assistance from Claude (Anthropic),
May 2026.

This work rests entirely on Knuth's DYNAHAM. The contribution of
this repository is the observation that DYNAHAM's transition
matrix sequence is exactly periodic in the stable interior, plus
the small toolkit needed to act on that observation.
