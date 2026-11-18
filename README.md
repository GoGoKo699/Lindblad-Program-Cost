# Lindblad Program Cost

A simulator receives quantum copies encoding an unknown Lindblad jump operator and must reproduce its dissipative action on an independent input. The optimal copy cost grows linearly with dimension. For rank-one jumps, its time dependence crosses from quadratic at short times to linear at long times.

The jump acts on a $`d`$-dimensional system, obeys $`\left\|L\right\|_{\mathrm F}=1`$, and is supplied through copies of

```math
|L\rangle=(L\otimes I)\sum_{j=1}^{d}|j,j\rangle.
```

The target is $`e^{t\mathcal D_L}`$, with $`\mathcal D_L(X)=LXL^\dagger-\{L^\dagger L,X\}/2`$. Error is **one half of the diamond norm**, including arbitrary reference entanglement.

## Rank-one jumps at all times

For $`L=|s\rangle\langle\phi|`$ with arbitrary unknown unit vectors $`s,\phi`$, the optimal copy cost is

```math
\boxed{
n^\star_{\mathrm{rank}\,1}(d,t,\varepsilon)
=\Theta\!\left(\frac{d}{\varepsilon}\min\{t^2,t\}\right)
}
```

for $`d\ge4`$, $`t\gt0`$, and $`0\lt\varepsilon\le\min(t^2,1)/1024`$. At long times, nearly parallel input and output directions make each jump recycle most of the population. When the removal rate is of order $`1/t`$, an input evolves into a state with a constant signal in a direction that each program encodes only weakly. This forces a copy cost linear in time. A processor matches it by estimating the two program factors and applying the channel for those estimates. The first-order errors cancel in the averaged output channel.

The [proof](science/rank_one_time_law.md) includes a rate-floor corollary. If the program is promised to have $`1-|\langle\phi|s\rangle|^2\ge q_0`$ for a supplied $`0\lt q_0\le1/4`$, the time factor saturates at $`1/q_0`$.

## Arbitrary jumps at short times

The minimum worst-case error over all fixed CPTP processors using $`n`$ supplied programs satisfies

```math
\boxed{
\delta_d(n,t)
\ge
\frac{(1-e^{-t/2})^2}{4}
\frac{d-3}{n+d-2}
}
\qquad
(d\ge4,\ n\ge0,\ t\ge0).
```

Together with the universal wave matrix Lindbladization upper bound, this gives the optimal copy law

```math
\boxed{
n_d^\star(t,\varepsilon)
=
\Theta\!\left(\frac{dt^2}{\varepsilon}\right)
}
```

for $`d\ge4`$, $`0\lt t\le1`$, and $`0\lt \varepsilon\le t^2/256`$. Explicit conservative constants are

```math
\frac{dt^2}{512\varepsilon}
\le n_d^\star(t,\varepsilon)
\le \frac{3dt^2}{\varepsilon}.
```

The converse permits arbitrary collective operations, ancillary states, measurements, and adaptive processing within a deterministic copy cap. Its hard instances have rank one and a known sink.

The mechanism is a reduction to quantum-state conjugation. A known reference component makes the desired dissipative channel produce an unknown unconjugated vector from one fixed input, while the supplied programs contain its conjugate. A success-weighted conjugation inequality limits this conversion without dividing by a small branch probability.

The universal conjugation bound and WML upper scaling are [established ingredients](references/RELATED_WORK.md). The dissipative reduction supplies the dimension lower bound for every processor and addresses the gap identified by Park et al. in arXiv:2605.30301v1.

## An orthogonal jump with known sink

For the family $`L=|s\rangle\langle\phi|`$ with known sink $`s`$ and unknown unit vector $`\phi⟂ s`$, an explicit collective processor attains

```math
n_{d,s}^\star(t,\varepsilon)
=
\Theta\!\left(
\frac{d(1-e^{-t/2})^2}{\varepsilon}
\right)
```

for every $`t\gt 0`$, $`d\ge4`$, and $`0\lt \varepsilon\le(1-e^{-t/2})^2/16`$. This family's cost saturates at order $`d/\varepsilon`$ for long times. The construction matches the desired coherences exactly and leaves a quantified population error.

The two rank-one results describe different promises. With a known orthogonal sink, population leaves the decay direction on a fixed timescale and the required copies stay bounded as time grows.

## Read the science

Start with the [model and exact claims](science/model.md), then follow the proof route below.

| To understand or inspect | Read |
|---|---|
| Supplied programs, allowed operations, and error convention | [Model](science/model.md) |
| Universal conjugation with arbitrary success probability | [Conjugation bound](science/conjugation.md) |
| Reference-bearing hard family and the dimension converse | [Dimension cost](science/dimension_cost.md) |
| Universal simulation algorithm and a complete error bound | [WML upper bound](science/upper_bound.md) |
| Arbitrary rank-one jumps, the long-time converse, and the averaged processor | [Rank-one time law](science/rank_one_time_law.md) |
| Known orthogonal sink, collective processor, and exact diamond error | [Orthogonal-sink construction](science/rank_one_processor.md) |
| Time dependence, physical meaning, and resource boundaries | [Interpretation](science/interpretation.md) |
| Precise prior-work comparison and bibliography | [Related work](references/RELATED_WORK.md) · [Bibliography](references/references.bib) |
| Finite computational checks and recorded outputs | [Verification](verification/README.md) |

[The reading guide](llms.txt) maps scientific questions to their authoritative files.

## Learn the project

The single educational anchor is **John Watrous, *The Theory of Quantum Information* (2018)**. The [Watrous bridge](learning/WATROUS_BRIDGE.md) identifies the relevant sections of the author's available manuscript, translates its notation, and develops the missing steps from channel theory and symmetric tensors to this project's dissipative-programming result.

The [30-milestone map](learning/MILESTONES.md) starts with **14 background milestones** on states, channels, norms, and symmetry. The **16 scientific milestones** trace the dimension converse and the two rank-one time regimes.

The map states understanding targets, earlier dependencies, conceptual connections, and source pointers. The bridge includes a complete four-dimensional example and explains which claims come from the anchor, which are established research ingredients, and which are derived in this project.

## Reproduce the finite checks

From the repository root:

```bash
python3 -m pip install -r verification/requirements.txt
python3 verification/run_all.py
```

The deterministic suite checks the rank-one channel against matrix exponentials, the Haar/Choi inequality, reference-phase sectors, the complete orthogonal-sink processor, and the WML algebra including a reference system. It also checks the long-time hard-family signal and first-order cancellation under symmetric estimates. The committed [run report](verification/expected/results.json) records the residuals and environment.

The written proofs establish the uniform results. Finite checks address the formulas and conventions used in those proofs.

## Scope

The cost is the number of supplied program copies. Program preparation, gate complexity, and coherent workspace are separate resources. The result uses no classical description of the unknown jump, preparation-unitary oracle, inverse oracle, or QRAM.

The lower bound is worst case over complex jumps. It does not assign the same cost to every jump or to stronger access models. Adaptive protocols are included when their copy count has a deterministic cap; a cost defined only through an unbounded stopping-time expectation is a different task.

## Purpose and contact

This repository serves as a record of the work and a guide for the author’s self-directed learning. For discussion or potential collaboration, please contact Ruge Lin at [gogoko699@gmail.com](mailto:gogoko699@gmail.com).
