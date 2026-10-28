# Model and scientific route

The resource considered here is a collection of quantum programs. Each program is one copy of a pure state encoding an unknown Lindblad jump operator. A fixed processor must use those programs to simulate the associated dissipative channel on an arbitrary input, including an input entangled with a reference.

The central question is how the required number of programs depends on the system dimension. A lower bound restricted to one simulation algorithm does not determine this cost: a different processor may process all programs coherently. The converse below permits every completely positive, trace-preserving processor.

## 1. Supplied-copy access

Let the system Hilbert space be $`\mathcal H=\mathbb C^d`$, with a known computational basis. An unknown operator $`L\in\mathsf L(\mathcal H)`$ obeys

```math
\|L\|_{\mathrm F}^2=\mathrm{Tr}(L^\dagger L)=1.
```

There is one jump and no additional Hamiltonian. The generator and target channel are

```math
\mathcal D_L(X)
=
LXL^\dagger-\frac12\bigl(L^\dagger LX+XL^\dagger L\bigr),
\qquad
\mathcal E_{L,t}=e^{t\mathcal D_L}.
```

Time $`t\ge0`$ is measured in this normalization. Rescaling a physical jump as $`J=\sqrt{\gamma}L`$ replaces $`t`$ by $`\gamma t_{\mathrm{phys}}`$.

Write

```math
|\Gamma_d\rangle=\sum_{j=1}^d |j\rangle|j\rangle,
\qquad
|L\rangle=(L\otimes I)|\Gamma_d\rangle,
\qquad
\pi_L=|L\rangle\langle L|.
```

The maximally entangled vector in this definition is **unnormalized**. Consequently, $`|L\rangle`$ has norm one. The supplied resource consists of $`n`$ independent copies $`\pi_L^{\otimes n}`$. This is the vectorization convention used by [Park et al. and Go et al.](../references/RELATED_WORK.md).

A processor $`\mathcal P_{d,n,t}`$ is any CPTP map from the system and its program registers to the system. It may depend on $`d,n,t`$, but its operations and ancillary states are independent of the unknown $`L`$. It induces

```math
\widetilde{\mathcal E}_{L,t}^{(n)}(X)
=
\mathcal P_{d,n,t}\bigl(X\otimes\pi_L^{\otimes n}\bigr).
```

Arbitrary coherent operations, measurements, classical processing, and feed-forward are allowed. Any adaptive protocol using at most $`n`$ copies on every branch has this form after adding unused copies and discarding them. Protocols charged only by an unbounded expected stopping time require a different cost definition.

The model supplies copies rather than a preparation circuit. A unitary preparing $`|L\rangle`$, its inverse, a classical description of $`L`$, or copies of an unconjugated hidden vector are additional resources. Their absence is part of the theorem's access assumptions.

## 2. Error and minimax copy cost

For two channels, the error convention is half the diamond norm:

```math
d_\diamond(\mathcal A,\mathcal B)
=
\frac12\|\mathcal A-\mathcal B\|_\diamond.
```

Equivalently, it is the maximum half trace distance between their outputs over joint system-reference density operators. A reference of dimension $`d`$ suffices. This convention gives an error in $`[0,1]`$.

Define

```math
\delta_d(n,t)
=
\inf_{\mathcal P_{d,n,t}}
\sup_{\|L\|_{\mathrm F}=1}
d_\diamond\!\left(
\widetilde{\mathcal E}_{L,t}^{(n)},\mathcal E_{L,t}
\right),
```

and

```math
n^\star_d(t,\varepsilon)
=
\min\{n\in\mathbb N_0:\delta_d(n,t)\le\varepsilon\}.
```

The minimum of an empty feasible set is defined to be $`+\infty`$. The supremum is over all single-jump operators with the stated normalization. Gate count, coherent workspace, and the cost of preparing the supplied programs are separate resources.

## 3. Main results

Set $`\beta(t)=1-e^{-t/2}`$.

### Universal finite-copy converse

For $`d\ge4`$, every $`n\ge0`$ and every $`t\ge0`$,

```math
\boxed{
\delta_d(n,t)
\ge
\frac{\beta(t)^2}{4}
\frac{d-3}{n+d-2}.
}
```

The [dimension-cost proof](dimension_cost.md) uses only rank-one jumps and one fixed, known system input. It therefore constrains the full channel task even though the proof needs no entangled input.

The hard family has a known sink $`|s\rangle`$, a known reference vector $`|r\rangle`$, and an unknown unit vector $`|\psi\rangle`$ in a $`D=d-2`$ dimensional subspace:

```math
|\phi_\psi\rangle
=
\frac{|r\rangle+|\psi\rangle}{\sqrt2},
\qquad
L_\psi=|s\rangle\langle\phi_\psi|.
```

Its program is $`|s\rangle|\overline{\phi_\psi}\rangle`$. Applying the target channel to $`|r\rangle`$ and selecting the unknown subspace produces the unnormalized state

```math
\frac{\beta(t)^2}{4}\,|\psi\rangle\langle\psi|.
```

This turns a simulator into a probabilistic conjugation device. The [conjugation inequality](conjugation.md) bounds its correct overlap without conditioning on success. The factor $`\beta(t)^2`$ is retained throughout the reduction.

### Optimal short-time scaling for arbitrary jumps

For

```math
d\ge4,\qquad
0\lt t\le1,\qquad
0\lt \varepsilon\le\frac{t^2}{256},
```

the minimum copy cost satisfies

```math
\boxed{
n^\star_d(t,\varepsilon)
=
\Theta\!\left(\frac{dt^2}{\varepsilon}\right).
}
```

The constants are independent of $`d,t,\varepsilon`$. A conservative explicit lower bound is

```math
n^\star_d(t,\varepsilon)
\ge
\frac{dt^2}{512\varepsilon}.
```

The [universal upper bound](upper_bound.md) constructs a fixed processor using fresh programs and proves a sufficient integer copy count. Its proof retains the finite-step condition rather than extrapolating a small-step estimate beyond its domain. The universal construction is inherited from wave matrix Lindbladization; the unrestricted-processor converse is the central contribution here.

The exact finite-copy converse holds at all times. The displayed $`t^2`$ matching law is asserted only in the stated regime. In particular, its lower-bound mechanism saturates as $`t\to\infty`$.

### A matching all-time law for known-sink rank-one jumps

If the jump is promised to be $`L=|s\rangle\langle\phi|`$, with known sink $`s`$ and unknown unit vector $`\phi\perp s`$, let $`n^\star_{d,s}(t,\varepsilon)`$ be the corresponding minimax copy cost. The [collective processor](rank_one_processor.md) and the same converse give

```math
\frac{d\,\beta(t)^2}{32\varepsilon}
\le n^\star_{d,s}(t,\varepsilon)
\le
\left\lceil\frac{(d-2)\beta(t)^2}{\varepsilon}\right\rceil,
```

for $`d\ge4`$, $`t\gt 0`$, and $`0\lt \varepsilon\le\beta(t)^2/16`$. Thus this family's optimal cost is $`\Theta(d\beta(t)^2/\varepsilon)`$ at all times. The construction is an arbitrary allowed CPTP processor on the data and symmetric program space. Its coherent-memory and gate costs are separate from the copy count.

## 4. What the reduction must preserve

The phase of $`|\psi\rangle`$ relative to the known $`|r\rangle`$ is part of the hard family. Averaging that phase removes coherence between different occupation sectors of the flagged program. Treating $`|\psi\rangle`$ only as a projective ray with an arbitrarily fixed representative would leave this step unjustified.

The program is conjugated because vectorization maps $`|s\rangle\langle\phi|`$ to $`|s\rangle|\bar\phi\rangle`$. Supplying $`|\phi\rangle`$ instead would change the physical resource. The system input in the reduction is the fixed vector $`|r\rangle`$ and contains no free unknown state.

Finally, the converse uses the full normalized output of the simulator. A witness with spectral diameter one compares the probability-weighted success branch directly. There is no division by its success probability and no conversion from trace distance to a square-root fidelity error.

## 5. Read the full argument

| Scientific question | Authoritative account |
|---|---|
| What does a program contain, and what operations are allowed? | This page |
| How well can a CP branch conjugate unknown state copies? | [Conjugation and reference flags](conjugation.md) |
| Why does dissipative simulation inherit the dimension cost? | [Universal dimension converse](dimension_cost.md) |
| Which algorithm supplies the matching universal upper bound? | [Wave matrix Lindbladization](upper_bound.md) |
| Can the rank-one family attain its time-dependent lower scaling? | [Collective processor](rank_one_processor.md) |
| What is the physical meaning of the resource law? | [Interpretation and scope](interpretation.md) |
| Which ingredients and nearby results are established work? | [Related work](../references/RELATED_WORK.md) |
| Which finite identities were checked computationally? | [Verification](../verification/README.md) |

The written proofs establish the uniform statements. Finite numerical checks test the vectorization, channel formulas, and operator inequalities used in those proofs; they do not establish the minimax theorem by sampling.
