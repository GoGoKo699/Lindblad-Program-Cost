# Verification

The proofs establish the uniform copy-complexity statements. The finite checks here target tensor-order, complex-conjugation, channel-normalization, and projector-algebra errors that could invalidate those proofs or the explicit processors.

## Run

From the repository root, using Python 3.11 or later:

```bash
python3 -m pip install -r verification/requirements.txt
python3 verification/run_all.py
```

The script prints one JSON report and raises an exception if a checked identity exceeds its tolerance. To save a fresh report:

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
| WML invariant algebra and reductions | 18 identity checks at dimensions 2 and 3, with reference dimension 2 | The five-map composition algebra and the reduced maps used in the universal diamond-bound proof |

The collective checks include $`\beta=1`$, which is the limiting complete-decay channel for this family. For every admissible case, the script compares the induced map on **all matrix units**, so the comparison covers the entire finite-dimensional linear channel. It also tests an input attaining the claimed normalized diamond error. The proof that this input is optimal, including reference extensions, is analytical.

The phase-sector check uses $`n+1`$ evenly spaced phases. Their discrete average cancels all nonzero Fourier differences between the $`0,\ldots,n`$ sectors exactly in algebra; numerical residuals arise from floating-point representation.

## Recorded outcome

The committed run passed. The largest core identity residual was approximately $`5.73\times10^{-16}`$, below the core tolerance $`10^{-10}`$. The smallest eigenvalue of the Haar/Choi slack matrix was approximately $`-1.26\times10^{-16}`$, within the same floating-point tolerance. The WML identities had maximum relative residual approximately $`5.17\times10^{-16}`$, below their tolerance $`10^{-12}`$.

These numbers describe finite matrix checks. They are not an estimate of theorem correctness or a certification of all dimensions.

## Implementation map

| File | Purpose |
|---|---|
| [verify_core.py](verify_core.py) | Hard damping family, symmetric tensors, phase averaging, and collective processor |
| [verify_wml.py](verify_wml.py) | Invariant-map and data-reference identities for WML |
| [run_all.py](run_all.py) | One entrypoint and optional JSON report output |
| [requirements.txt](requirements.txt) | Pinned numerical dependencies |
| [expected/results.json](expected/results.json) | Actual recorded output, including seeds and environment versions |

The scientific formulas are derived in [conjugation](../science/conjugation.md), [dimension cost](../science/dimension_cost.md), [the universal upper bound](../science/upper_bound.md), and [the collective processor](../science/rank_one_processor.md).

## Limits of verification

Finite checks do not prove Haar integration for every dimension, the minimax quantifier over all processors, or an asymptotic scaling law. They also do not establish an efficient gate implementation or compare physical devices. Each of those questions is separate from the identities tested here.
