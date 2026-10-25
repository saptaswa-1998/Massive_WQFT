# Computational checks for massive WQFT elliptic integrals

This repository contains **Wolfram Mathematica notebooks** and a **GiNaC/C++ numerical validator** accompanying the manuscript

> **Geometry of massive loops in gravitational scattering: Maximal cut, Riemann's bilinear relations and genus drop**  
> Arpan Bhattacharyya and Saptaswa Ghosh

The purpose of these files is to make selected symbolic derivations and numerical checks in the manuscript independently inspectable and reproducible. They cover integration-by-parts (IBP) reduction on the maximal cut, intersection theory, differential equations, elliptic periods, monodromy, twisted Riemann period relations, and the dimensional-regularization expansion.

## Background and scope

The manuscript studies a two-loop worldline quantum field theory (WQFT) integral family with four massive quadratic propagators and two worldline delta-function constraints. In the four-dimensional maximal-cut analysis, the residual geometry can be represented by the elliptic curve

\[
E_\tau:\qquad y^2=(1-z^2)(z^2+\tau).
\]

The accompanying calculations investigate the connection between the Baikov maximal-cut representation, a five-master differential system (with three-dimensional elliptic and two-dimensional complementary blocks), and its period and monodromy structures. The manuscript also studies a Bessel-weighted Fourier cut and quotient-curve geometry; **the files in this repository are not a complete automated verification of every statement in the manuscript**.

## Repository contents

| File | Main purpose | Related manuscript material |
| --- | --- | --- |
| [`01-litered-intersection-ibp-crosscheck.nb`](01-litered-intersection-ibp-crosscheck.nb) | Three-way check of maximal-cut ISP-numerator reduction using LiteRed, local intersection-theory reduction, and a global twisted-exact certificate. The example involves \((\ell_1\!\cdot\!v_2)^5\). | Sections 4–5 |
| [`02-orthogonal-dlog-gauss-manin-checks.nb`](02-orthogonal-dlog-gauss-manin-checks.nb) | Constructs the orthogonal dlog basis and its Gauss–Manin connection; checks basis transformations and twisted-exact certificates against the monomial-basis system. | Section 4 |
| [`03-elliptic-fundamental-matrix-checks.nb`](03-elliptic-fundamental-matrix-checks.nb) | Verifies the three-dimensional \(\epsilon=0\) fundamental matrix, complete/incomplete elliptic-integral columns, determinant normalization, and Legendre-type relations. | Sections 6 and 8 |
| [`04-exact-monodromy-bilinear-checks.nb`](04-exact-monodromy-bilinear-checks.nb) | Exact symbolic checks of orbit-basis monodromy matrices, the global monodromy relation, and an invariant twisted bilinear pairing. | Section 4.6 |
| [`05-libra-monodromy-independent-check.nb`](05-libra-monodromy-independent-check.nb) | Independent reconstruction of monodromy using Libra, local series, balance transformations, and Levelt-type analysis; compares with exact global identities. | Section 4.6 |
| [`06-twisted-riemann-period-relations.nb`](06-twisted-riemann-period-relations.nb) | Symbolic and high-precision checks of cocycle/cycle intersection matrices and twisted Riemann period relations, including beta–Appell identities and the elliptic Legendre limit. | Section 6 |
| [`07-exact-epsilon-factorization-checks.nb`](07-exact-epsilon-factorization-checks.nb) | Checks affine \(\epsilon\)-dependence, the fundamental matrix, transformed kernels in both parity blocks, the finite-\(\epsilon\) cut differential equation, and boundary data. | Section 8 |
| [`ginac_empl_validation.cpp`](ginac_empl_validation.cpp) | Independent GiNaC-based evaluation of iterated-integral words and numerical checks of \(\epsilon\)-expansion coefficients, shuffle identities, and differential-equation residuals. | Sections 7–8 |

The filename `01-litered-intersection-ibp-crosscheck.nb` is kept as distributed; the external Mathematica package is **LiteRed**.

## Requirements

**Mathematica notebooks**

- Wolfram Mathematica / Wolfram Language. 
- **LiteRed** for the package-dependent reduction in notebook 01.
- **Libra** for notebook 05. Its saved outputs identify **Libra 1.23**.
- The other notebooks primarily use Wolfram Language's built-in symbolic and numerical functionality; inspect their initialization cells for any local requirements.

**C++ numerical validator**

- A C++17 compiler.
- [GiNaC](https://www.ginac.de/) and its dependencies (including CLN), with development headers and libraries installed.
- `pkg-config`, for the build command below.

The external Mathematica packages and GiNaC libraries are **not bundled** in this repository.

## Running the Mathematica notebooks

1. Open a notebook in Mathematica.
2. Start a **fresh kernel** and evaluate its cells **in the order presented** (for example, using **Evaluation → Evaluate Notebook**).
3. Read the printed checks and residuals. Some notebooks print `CHECK=PASS`, named `PASS`/`FAIL` assertions, or an overall verification summary.
4. If a check fails, inspect the first failed assertion and package-loading messages. Some notebook checks deliberately abort or terminate the kernel on failure.

**LiteRed setup (notebook 01).** Its initialization code reads the environment variable `LITERED_PATH`, which must point to the LiteRed **loader `.m` file**. For example, in the current Mathematica session, before evaluating the notebook:

```wolfram
SetEnvironment["LITERED_PATH" -> "/absolute/path/to/LiteRed2026.m"];
```

Alternatively, set `LITERED_PATH` in the environment **before launching Mathematica**. The notebook contains an author-specific fallback path; do not rely on that fallback on another machine.

**Libra setup (notebook 05).** Its first section uses an author-specific absolute path to `Libra.m`. **Edit that package-loading path** to point to your local Libra installation before evaluating the notebook. The saved outputs correspond to Libra 1.23.

**Saved outputs are not a fresh run.** The `.nb` files may display results from previous evaluations. To independently reproduce them, use a new kernel and re-evaluate the input cells.

## Building and running the GiNaC validator

From the repository root, on a system with GiNaC available through `pkg-config`:

```bash
c++ -std=c++17 -O2 ginac_empl_validation.cpp \
  -o ginac_empl_validation $(pkg-config --cflags --libs ginac)

./ginac_empl_validation 512
```

The optional integer argument is the Gaussian quadrature order **per subinterval**; the default is `512`, and values below `8` are rejected. Larger values can increase computation time.

The program checks, at specified sample kinematics and numerical points:

- The leading, first, and second \(\epsilon\)-expansion coefficients against reference values embedded in the source.
- GiNaC multiple-polylogarithm evaluations, including ordered depth-two words and shuffle identities.
- A native GiNaC `iterated_integral` example against numerical quadrature.
- Differential-equation residuals at sample values of \(\tau\).

Its final line is `VALIDATION: PASS` if **all implemented tolerance tests** succeed, or `VALIDATION: FAIL` otherwise. Exit status `0` means all tests passed, `2` indicates failed validation, and `1` indicates an exception. These are finite-precision checks at specified points, not proofs for arbitrary kinematics.

## Reproducibility notes

- **Conventions matter:** contour orientation, sheet choices, branch prescriptions, and elliptic-integral parameter conventions are fixed in the manuscript and/or individual checks. Apparent sign differences should first be examined against those conventions.
- **Exact versus numerical evidence:** Mathematica performs symbolic equalities in several notebooks; the GiNaC program supplies complementary finite-precision tests. Neither replaces the assumptions made in the analytic derivations.
- **Fresh kernels:** run each notebook separately to avoid leakage of symbols and definitions across sessions.
- **Package paths:** remove or replace machine-specific paths before running on a different system.
- **Scope:** the focus is on maximal-cut/homogeneous systems and associated period identities. This repository does not implement a full calculation of the uncut scattering observable.

## Suggested verification order

For a structured inspection, run **01 → 02 → 03 → 04 → 05 → 06 → 07**, and then compile and run `ginac_empl_validation.cpp`. Notebook 05 provides an external-package cross-check of the monodromy structure also examined in notebook 04; the C++ program independently tests selected Section 7–8 expressions.

## Citation

If these files are useful, please cite the accompanying manuscript:

**Arpan Bhattacharyya and Saptaswa Ghosh**, *Geometry of massive loops in gravitational scattering: Maximal cut, Riemann's bilinear relations and genus drop*.
(https://arxiv.org/abs/2610.07234)

## Acknowledgements and AI-assisted tools

The authors acknowledge the developers of **Wolfram Mathematica, LiteRed, Libra, and GiNaC**, which were used for symbolic and numerical computations in this project.

**ChatGPT (OpenAI)** and **Claude (Anthropic)** were used as  tools during the preparation of the codes. All the logical flow behind the
AI coding has been directed and verified by the authors.
