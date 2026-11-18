# What a quantum dissipative program costs

The program $`|L\rangle`$ encodes every coefficient of an unknown jump operator in a normalized quantum state. That compact representation does not imply that a small number of copies can reproduce the operator's action on an independent input. The [universal converse](dimension_cost.md) quantifies this distinction.

## An unknown direction is already hard

The lower bound uses a rank-one jump $`L_\psi=|s\rangle\langle\phi_\psi|`$ with a known sink. The program is a product across its two vectorization registers. The construction therefore needs neither a complicated singular-value spectrum nor entanglement between those registers.

Its difficulty is the orientation of the decay direction. The program carries $`|\bar\phi_\psi\rangle`$, whereas the dissipative action contains $`|\phi_\psi\rangle\langle\phi_\psi|`$. Complex conjugation is basis dependent, and the access model fixes that basis. A simulator must generate the correct action using physical operations on the supplied conjugate programs.

The known reference component in $`\phi_\psi`$ makes this restriction observable on one known input. After evolution, a fixed output projection selects $`|\psi\rangle`$ with probability

```math
p(t)=\frac{(1-e^{-t/2})^2}{4}.
```

Universal conjugation cannot reproduce this branch arbitrarily accurately from a finite number of programs. The [success-weighted inequality](conjugation.md) bounds the average correct overlap relative to the total success probability. It applies to coherent, adaptive, and postselected processing.

## Why the short-time law is quadratic

For small $`t`$,

```math
1-e^{-t/2}=\frac t2+O(t^2),
\qquad
p(t)=\frac{t^2}{16}+O(t^3).
```

The no-jump amplitude in the selected unknown subspace is first order in time, so its probability is second order. The witness in the converse compares this probability-weighted output directly with the simulated output. It neither normalizes a rare branch nor pays a square-root conversion between error measures.

The result is a lower bound proportional to $`dt^2/\varepsilon`$ once the allowed error is sufficiently smaller than this branch weight. The [WML construction](upper_bound.md) supplies the matching upper order for all single jumps. Its independent short steps each have error of order $`d(t/n)^2`$, and channel contraction lets those errors accumulate at most linearly over $`n`$ steps.

This matching law has an explicit accuracy domain. For example, a coarse error tolerance can make the identity channel adequate at very short times, so a positive copy lower bound of the same order cannot be inferred for every possible relationship between $`t`$ and $`\varepsilon`$.

## When long-time copy cost saturates

For a known-sink rank-one jump, the decay has a simple saturation structure. The exact channel depends on time only through $`\beta(t)=1-e^{-t/2}`$. The [collective processor](rank_one_processor.md) reproduces its coherence terms exactly and leaves an error

```math
\frac{\beta(t)^2(d-2)}{n}
```

in normalized diamond distance, whenever its stated CPTP condition holds. This matches the converse in order and gives

```math
n^\star_{d,s}(t,\varepsilon)
=
\Theta\!\left(
\frac{d(1-e^{-t/2})^2}{\varepsilon}
\right)
```

for $`d\ge4`$, $`t\gt 0`$, and $`0\lt \varepsilon\le\beta(t)^2/16`$.

The optimal cost in this family therefore saturates at order $`d/\varepsilon`$ for long times. The construction allows a negative surviving amplitude on one known joint projector. That operation is a valid CPTP processor even when it is not itself a positive-time damping semigroup.

For general rank-one jumps the two directions may overlap. At each jump the same bright component can be populated again. Writing $`q=1-|\langle\phi|s\rangle|^2`$, only a fraction $`q`$ of this component leaves in each jump. The bright population therefore decays as $`e^{-qt}`$. The slow transfer changes the minimax cost.

In the long-time hard family the program has weight $`q=1/t`$ on the unknown conjugated direction, so $`n`$ copies contain $`n/t`$ such factors on average after resolving their occupation sectors. The target output has constant population in the unconjugated direction. The [rank-one proof](rank_one_time_law.md) combines both quantities in a success-weighted conjugation converse to give order $`dt/\varepsilon`$ for all processors.

The matching processor applies the channel for a pair of covariant estimates and averages the result. Unitary symmetry makes the first-order perturbation vanish in the average. The second-order channel perturbation costs at most a factor proportional to $`t`$, uniformly even for arbitrarily slow removal. A positive lower bound $`q_0`$ on the removal rate makes the worst-case time factor saturate at order $`1/q_0`$, within the [stated promise range](rank_one_time_law.md#a-positive-removal-rate-floor). The matching result for arbitrary ranks concerns short times.

## What the quantifiers imply

The converse is a worst-case statement over unknown jumps. It does not say that every jump requires a dimension-dependent number of programs. A known jump requires no informational copies, and additional promises can reduce the task.

Frobenius normalization also matters. It fixes the total squared matrix amplitude, while the operator norm controls the largest local decay rate. The rank-one hard family has operator norm one. The refined WML error estimate retains a factor $`\|L\|_\infty^2`$, so its error can be smaller on other jumps. The associated finite-step condition must still be enforced when turning that estimate into a sample budget.

For every positive time and $`d\ge4`$, the finite-copy lower bound is strictly positive. Thus exact universal simulation of even the hard known-sink family is impossible with a finite deterministic supply of these programs. Approximation is essential.

## Counted and additional resources

| Resource or operation | Role in the results |
|---|---|
| Independent copies of the vectorized jump state | The charged resource |
| Arbitrary fixed CPTP processing and ancillary states | Allowed, independent of the unknown jump |
| Joint processing of supplied copies | Allowed by the converse, used by the orthogonal-sink contraction and covariant measurements |
| Adaptive measurements and classical feed-forward | Allowed within a deterministic copy cap |
| Entanglement with an external reference | Included in the channel error guarantee |
| Program preparation and data loading | Separate costs; programs are supplied explicitly |
| Gate count and coherent workspace | Unrestricted here; no efficiency claim follows from copy sufficiency |
| Classical jump description or preparation-unitary oracle | Different access models |
| A stopping rule charged only by its expected sample count | Outside the deterministic-budget theorem |

A finite-dimensional Stinespring dilation makes each upper-bound processor physically admissible in the model. It does not imply that compiling that dilation is efficient. Conversely, unlimited computational effort cannot remove the converse's copy requirement.

## Contribution and evidence

The universal conjugation fidelity, WML upper scaling, and covariant pure-state estimation are established ingredients. The [related-work account](../references/RELATED_WORK.md) credits them and identifies the previous algorithm-specific dimension lower bound. The dissipative reductions give processor-independent lower bounds on dimension and time. The uniform channel estimate explains why the averaged processor can attain the long-time rank-one law.

The central guarantee is analytical. The [verification suite](../verification/README.md) compares rank-one channel formulas with matrix exponentials, tests the symmetric-projector inequality and reference-phase decomposition, and reconstructs the collective processor and WML algebra in finite dimensions. Separate finite checks cover the long-time signal and symmetric averaging.
