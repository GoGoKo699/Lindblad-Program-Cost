# Verification

The proofs establish the uniform copy-complexity statements. The finite checks here target tensor-order, complex-conjugation, channel-normalization, and projector-algebra errors that could invalidate those proofs or the explicit processors. The same entrypoint checks document navigation, mathematical delimiters, tables, and milestone structure.

## Run

From the repository root, using Python 3.11 or later:

```bash
python3 -m pip install -r verification/requirements.txt
python3 verification/run_all.py
```

The script prints one JSON report and raises an exception if a checked identity exceeds its tolerance or a document check fails. The document checks can also run without NumPy or SciPy:

```bash
python3 verification/verify_docs.py
```

To save a fresh combined report:

```bash
python3 verification/run_all.py --output verification/expected/results.json
```

The second command writes the named output file. The committed [reference report](expected/results.json) was generated with Python 3.12.14, NumPy 2.3.5, and SciPy 1.17.0. The requirements file pins the numerical dependencies. Floating-point residuals may vary slightly across BLAS implementations; the mathematical acceptance tolerances, rather than exact serialized residuals, determine success.

## What is checked

| Check | Coverage | Scientific risk addressed |
|---|---|---|
| Vectorized hard-family program | Dimensions 4, 5, and 7, with complex directions | Confusing the supplied conjugate vector with a free unconjugated vector |
| Exact damping channel | 15 dimension/time cases, with times 0, 0.01, 0.2, 1, and 5 | Comparing the full Liouvillian exponential with the proposed Kraus map and selected success branch |
| Haar/Choi operator inequality | 8 dimension/copy cases, including zero copies | Symmetric-projector normalization and the coefficient in the success-weighted conjugation inequality |
| Reference-phase sectors | Copy counts 0 through 4 | Retaining cross-sector coherence that should vanish under phase averaging |
| Collective rank-one processor | 38 admissible parameter cases, including 17 with a negative surviving amplitude | Complete positivity, trace preservation, the entire reduced channel, and the attained trace-distance error |
| Exact channel for arbitrary rank-one jumps | 196 dimension, time, and overlap cases, including parallel and orthogonal directions | Comparing the transfer integral with the full Liouvillian exponential and checking trace preservation |
| Long-time hard-family signal | 101 times from 4 to one million, at each of four removal-rate floors including zero | Checking the output weight in the hidden direction and its analytical lower bound across the linear and saturation regimes |
| Symmetric channel average | 54 cases with complex directions and two geodesic displacement sizes | Checking cancellation of the first-order error when the channels for symmetric perturbations are averaged |
| WML invariant algebra and reductions | 18 identity checks at dimensions 2 and 3, with reference dimension 2 | The five-map composition algebra and the reduced maps used in the universal diamond-bound proof |
| Document structure | Relative file links and heading fragments, mathematical delimiters, table widths, and all 30 milestone records | Broken reading routes, malformed Markdown, missing milestone fields, and dependencies on later milestones |

The collective checks include $`\beta=1`$, which is the limiting complete-decay channel for this family. For every admissible case, the script compares the induced map on **all matrix units**, so the comparison covers the entire finite-dimensional linear channel. It also tests an input attaining the claimed normalized diamond error. The proof that this input is optimal, including reference extensions, is analytical.

The phase-sector check uses $`n+1`$ evenly spaced phases. Their discrete average cancels all nonzero Fourier differences between the $`0,\ldots,n`$ sectors exactly in algebra; numerical residuals arise from floating-point representation.

The symmetric-average check applies each channel to one half of a maximally entangled state. For each displacement it compares the largest individual half trace-distance error with the error after equally averaging the four independent sign choices for the two geodesics. Halving the displacement reduces the averaged error by approximately a factor of four. The uniform bound for the full covariant estimator and arbitrary reference inputs is analytical.

The document checker also guards two observed GitHub rendering failures: the rejected operator-name macro and raw angle brackets in mathematical expressions. The supported forms are upright names with `\mathrm{...}` and inequalities with `\lt` or `\gt`. These source checks support a separate inspection of the hosted pages; they do not certify visual layout or the mathematical content of a milestone. External source links are outside the local navigation check.

## Recorded outcome

The committed run passed. The largest core identity residual was approximately $`5.72\times10^{-16}`$, below the core tolerance $`10^{-10}`$. The smallest eigenvalue of the Haar/Choi slack matrix was approximately $`-1.26\times10^{-16}`$, within the same floating-point tolerance. The WML identities had maximum relative residual approximately $`5.16\times10^{-16}`$, below their tolerance $`10^{-12}`$.

The general rank-one matrix-exponential and trace-preservation residuals were at most $`7.16\times10^{-14}`$ and $`4.45\times10^{-16}`$, respectively. The minimum recorded hard-family signal was 0.29968. The ratio of averaged errors before and after halving the displacement ranged from 3.9905 to 4.0000.

These numbers describe finite matrix checks. They are not an estimate of theorem correctness or a certification of all dimensions.

## Implementation map

| File | Purpose |
|---|---|
| [verify_core.py](verify_core.py) | Hard damping family, symmetric tensors, phase averaging, and collective processor |
| [verify_rank_one_time.py](verify_rank_one_time.py) | Arbitrary rank-one transfer integral, long-time hard-family signal, and symmetric channel averaging |
| [verify_wml.py](verify_wml.py) | Invariant-map and data-reference identities for WML |
| [verify_docs.py](verify_docs.py) | Relative links, heading fragments, mathematical delimiters, tables, and the milestone map |
| [run_all.py](run_all.py) | One entrypoint and optional JSON report output |
| [requirements.txt](requirements.txt) | Pinned numerical dependencies |
| [expected/results.json](expected/results.json) | Actual recorded output, including seeds and environment versions |

The scientific formulas are derived in [conjugation](../science/conjugation.md), [dimension cost](../science/dimension_cost.md), [the universal upper bound](../science/upper_bound.md), [the collective processor](../science/rank_one_processor.md), and [the rank-one time law](../science/rank_one_time_law.md).

## Limits of verification

Finite checks do not prove Haar integration for every dimension, the minimax quantifier over all processors, or an asymptotic scaling law. They also do not establish an efficient gate implementation or compare physical devices. Each of those questions is separate from the identities tested here.
