# Lindblad Program Cost

A simulator receives quantum copies encoding an unknown Lindblad jump operator and must reproduce its dissipative action on an independent input. This project establishes an unavoidable copy cost proportional to the system dimension, even when the simulator can process all programs coherently.

The jump acts on a $`d`$-dimensional system, obeys $`\|L\|_{\mathrm F}=1`$, and is supplied through copies of

```math
|L\rangle=(L\otimes I)\sum_{j=1}^{d}|j,j\rangle.
```

The target is $`e^{t\mathcal D_L}`$, with $`\mathcal D_L(X)=LXL^\dagger-\{L^\dagger L,X\}/2`$. Error is **one half of the diamond norm**, including arbitrary reference entanglement.

## Principal result

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

for $`d\ge4`$, $`0<t\le1`$, and $`0<\varepsilon\le t^2/256`$. Explicit conservative constants are

```math
\frac{dt^2}{512\varepsilon}
\le n_d^\star(t,\varepsilon)
\le \frac{3dt^2}{\varepsilon}.
```

The converse permits arbitrary collective operations, ancillary states, measurements, and adaptive processing within a deterministic copy cap. Its hard instances have rank one and a known sink.

The mechanism is a reduction to quantum-state conjugation. A known reference component makes the desired dissipative channel produce an unknown unconjugated vector from one fixed input, while the supplied programs contain its conjugate. A success-weighted conjugation inequality limits this conversion without dividing by a small branch probability.

The universal conjugation bound and WML upper scaling are [established ingredients](references/RELATED_WORK.md). The contribution is the dissipative reduction and the resulting dimension lower bound for every processor. This addresses the unrestricted-processor gap distinguished from the WML-specific lower bound in Park et al., arXiv:2605.30301v1.

## An explicit processor for rank-one decay

For the family $`L=|s\rangle\langle\phi|`$ with known sink $`s`$ and unknown unit vector $`\phi\perp s`$, an explicit collective processor attains

```math
n_{d,s}^\star(t,\varepsilon)
=
\Theta\!\left(
\frac{d(1-e^{-t/2})^2}{\varepsilon}
\right)
```

for every $`t>0`$, $`d\ge4`$, and $`0<\varepsilon\le(1-e^{-t/2})^2/16`$. This family's cost saturates at order $`d/\varepsilon`$ for long times. The construction matches the desired coherences exactly and leaves a quantified population error.

The all-time statement concerns this specified rank-one family. The matching result for arbitrary jumps has the short-time regime stated above.

## Read the science

Start with the [model and exact claims](science/model.md), then follow the proof route below.

| To understand or inspect | Read |
|---|---|
| Supplied programs, allowed operations, and error convention | [Model](science/model.md) |
| Universal conjugation with arbitrary success probability | [Conjugation bound](science/conjugation.md) |
| Reference-bearing hard family and the dimension converse | [Dimension cost](science/dimension_cost.md) |
| Universal simulation algorithm and a complete error bound | [WML upper bound](science/upper_bound.md) |
| Collective rank-one processor and exact diamond error | [Rank-one construction](science/rank_one_processor.md) |
| Time dependence, physical meaning, and resource boundaries | [Interpretation](science/interpretation.md) |
| Precise prior-work comparison and bibliography | [Related work](references/RELATED_WORK.md) · [Bibliography](references/references.bib) |
| Finite computational checks and recorded outputs | [Verification](verification/README.md) |

[The reading guide](llms.txt) maps scientific questions to their authoritative files.

## Learn the project

The single educational anchor is **John Watrous, *The Theory of Quantum Information* (2018)**. The [Watrous bridge](learning/WATROUS_BRIDGE.md) identifies the relevant sections of the author's available manuscript, translates its notation, and develops the missing steps from channel theory and symmetric tensors to this project's dissipative-programming result.

The [30-milestone map](learning/MILESTONES.md) is a common route through the scientific story. Its **14 background/preliminaries milestones** assemble the state, channel, norm, and symmetry tools into the supplied-copy model. Its **16 narrative/innovation milestones** connect the reference-bearing hard family to the universal converse, the two upper constructions, and their resource interpretation.

The map states understanding targets, earlier dependencies, conceptual connections, and source pointers. The bridge includes a complete four-dimensional example and explains which claims come from the anchor, which are established research ingredients, and which are derived in this project.

## Reproduce the finite checks

From the repository root:

```bash
python3 -m pip install -r verification/requirements.txt
python3 verification/run_all.py
```

The deterministic suite checks the hard-family channel against a matrix exponential, the Haar/Choi operator inequality, reference-phase sectors, the complete collective processor, and the WML algebra including a reference system. The committed [run report](verification/expected/results.json) records the actual residuals and environment.

The written proofs establish the uniform results. Finite checks address the formulas and conventions used in those proofs.

## Scope

The cost is the number of supplied program copies. Program preparation, gate complexity, and coherent workspace are separate resources. The result uses no classical description of the unknown jump, preparation-unitary oracle, inverse oracle, or QRAM.

The lower bound is worst case over complex jumps. It does not assign the same cost to every jump or to stronger access models. Adaptive protocols are included when their copy count has a deterministic cap; a cost defined only through an unbounded stopping-time expectation is a different task.

## Purpose and contact

This repository serves as a record of the work and a guide for the author’s self-directed learning. For discussion or potential collaboration, please contact Ruge Lin at [gogoko699@gmail.com](mailto:gogoko699@gmail.com).
