# Conjugation from quantum programs

The lower bound in [Dimension cost](dimension_cost.md) uses a physical restriction on converting copies of a state into its complex conjugate. The restriction must cover a successful branch of a larger protocol, with its success probability retained. This page establishes that form of the bound and then applies it to a program containing a coherent reference component.

All transposes and complex conjugates refer to the fixed basis used to [vectorize the jump operator](model.md). Write $`P_\psi=|\psi\rangle\langle\psi|`$ and $`\overline P_\psi=|\overline\psi\rangle\langle\overline\psi|`$. The unknown vector belongs to a space $`\mathcal H\cong\mathbb C^D`$, with $`D\ge2`$. Integration over $`\psi`$ means normalized Haar measure on the unit sphere, including its phase.

## The inherited transposition bound

[Brzić, Grinko, Studziński and Quintino, Theorem 1 and Appendix B](https://arxiv.org/abs/2603.23628v3), characterize optimal universal pure-state transposition. Their $`k\to1`$ fidelity is

```math
c_{k,D}=\frac{k+1}{k+D}.
```

Their theorem also identifies an optimal estimation channel. The symmetric-subspace and Choi-matrix argument is an inherited ingredient of this project. We give it below in an unnormalized form so it can be applied to a completely positive branch without discarding the probability of that branch. The reference-sector argument in the following section explains how this ingredient enters the particular dissipative program family.

### Lemma 1: conjugation with arbitrary success probability

Let $`k\ge0`$ be an integer, and let $`\mathcal T`$ be a completely positive map from $`k`$ copies of $`\mathcal H`$ to $`\mathcal H`$. Define

```math
\begin{aligned}
s_\psi&=\operatorname{Tr}\mathcal T(\overline P_\psi^{\otimes k}),\\
f_\psi&=\operatorname{Tr}\!\left[P_\psi\,
\mathcal T(\overline P_\psi^{\otimes k})\right].
\end{aligned}
```

Then

```math
\int f_\psi\,d\psi
\le c_{k,D}\int s_\psi\,d\psi.
```

When $`\mathcal T`$ is trace nonincreasing, $`s_\psi`$ is its success probability and $`f_\psi`$ is the probability weight of the correct output. The lemma allows the success probability to depend on $`\psi`$ and to be arbitrarily small. It imposes no covariance or measure-and-prepare restriction on $`\mathcal T`$.

#### Proof

Let $`P_k`$ be the projector onto the symmetric subspace $`\operatorname{Sym}^k(\mathcal H)`$, and let

```math
h_k=\dim\operatorname{Sym}^k(\mathcal H)
=\binom{D+k-1}{k}.
```

For $`k=0`$, set $`P_0=1`$ and $`h_0=1`$. Haar invariance and normalization give the standard moment identity

```math
\int P_\psi^{\otimes k}\,d\psi=\frac{P_k}{h_k}.
```

Use the unnormalized Choi convention

```math
\mathcal T(X)=\operatorname{Tr}_{\mathrm{in}}
\left[J_{\mathcal T}(I_{\mathrm{out}}\otimes X^{\mathsf T})\right],
\qquad J_{\mathcal T}\ge0.
```

The transpose in this identity removes the conjugation on the supplied program:
$`(\overline P_\psi^{\otimes k})^{\mathsf T}=P_\psi^{\otimes k}`$.
Consequently,

```math
\begin{aligned}
\int f_\psi\,d\psi
&=\frac{1}{h_{k+1}}\operatorname{Tr}(J_{\mathcal T}P_{k+1}),\\
\int s_\psi\,d\psi
&=\frac{1}{h_k}\operatorname{Tr}
\left[J_{\mathcal T}(I_{\mathrm{out}}\otimes P_k)\right].
\end{aligned}
```

Here $`P_{k+1}`$ acts on the output factor followed by the $`k`$ input factors, identified in the fixed basis. Its range lies inside $`\mathcal H\otimes\operatorname{Sym}^k(\mathcal H)`$, so

```math
P_{k+1}\le I_{\mathrm{out}}\otimes P_k.
```

Taking the trace against the positive matrix $`J_{\mathcal T}`$ proves the result, because

```math
\frac{h_k}{h_{k+1}}=\frac{k+1}{k+D}.
```

Trace preservation was not used. This is why the same inequality applies to postselected branches. For $`k=0`$, the argument reduces to the identity $`\int P_\psi\,d\psi=I/D`$.

## A coherent reference in the program

The hard family in [Dimension cost](dimension_cost.md#the-hard-family) supplies copies of

```math
|\chi_\psi\rangle
=\frac{|r\rangle+|\overline\psi\rangle}{\sqrt2},
\qquad \langle r|\overline\psi\rangle=0,
```

where $`|r\rangle`$ is a fixed basis vector outside $`\mathcal H`$. The relative phase between the reference component and $`|\overline\psi\rangle`$ is part of the input. Lemma 1 cannot be applied merely by treating this state as a copy of $`|\overline\psi\rangle`$. The following occupation-sector argument accounts for the reference explicitly.

### Lemma 2: conjugation from a reference-bearing program

Let $`n\ge0`$ be an integer, and let $`\mathcal T`$ be a completely positive map from $`n`$ copies of $`\mathbb C|r\rangle\oplus\mathcal H`$ to $`\mathcal H`$. Define

```math
\begin{aligned}
s_\psi&=\operatorname{Tr}\mathcal T
\left((|\chi_\psi\rangle\langle\chi_\psi|)^{\otimes n}\right),\\
f_\psi&=\operatorname{Tr}\!\left[P_\psi\,\mathcal T
\left((|\chi_\psi\rangle\langle\chi_\psi|)^{\otimes n}\right)\right].
\end{aligned}
```

Then

```math
\boxed{\displaystyle
\int f_\psi\,d\psi
\le \frac{n+1}{n+D}\int s_\psi\,d\psi.}
```

As in Lemma 1, the statement includes trace-nonincreasing branches and state-dependent success probabilities.

#### Proof

There is an orthogonal occupation decomposition

```math
\operatorname{Sym}^n(\mathbb C|r\rangle\oplus\mathcal H)
\cong\bigoplus_{k=0}^{n}\operatorname{Sym}^k(\mathcal H).
```

Let $`V_k`$ be the normalized isometry embedding the sector with exactly $`k`$ factors in $`\mathcal H`$ and $`n-k`$ reference factors. The program expands as

```math
|\chi_\psi\rangle^{\otimes n}
=\sum_{k=0}^{n}\sqrt{w_k}\,
V_k|\overline\psi\rangle^{\otimes k},
\qquad w_k=2^{-n}\binom nk.
```

Replacing $`\psi`$ by $`e^{i\theta}\psi`$ multiplies the $`k`$th program component by $`e^{-ik\theta}`$. The desired output projector $`P_\psi`$ is unchanged. Therefore averaging over $`\theta`$ eliminates every term between distinct sectors in both the success integral and the correct-overlap integral. This step concerns the Haar average and does not assume that the processor itself is covariant.

Define the CP maps

```math
\mathcal T_k(X)=\mathcal T(V_kXV_k^\dagger).
```

Here $`\mathcal T_k`$ is initially defined on $`\operatorname{Sym}^k(\mathcal H)`$. Composing with the projection onto that subspace extends it to a CP map on all $`k`$ input factors, so Lemma 1 applies directly.

After phase averaging,

```math
\begin{aligned}
\int f_\psi\,d\psi
&=\sum_{k=0}^{n}w_k
\int\operatorname{Tr}
\left[P_\psi\,\mathcal T_k(\overline P_\psi^{\otimes k})\right]d\psi,\\
\int s_\psi\,d\psi
&=\sum_{k=0}^{n}w_k
\int\operatorname{Tr}\mathcal T_k(\overline P_\psi^{\otimes k})\,d\psi.
\end{aligned}
```

Apply Lemma 1 in each sector. Since

```math
\frac{k+1}{k+D}\le\frac{n+1}{n+D}
\qquad(0\le k\le n),
```

and every success contribution is nonnegative, the claimed inequality follows.

## Why the unnormalized form matters

A simulator may use intermediate measurements, retain classical outcomes, or direct only some outcomes into the subspace examined by the proof. Fixing the simulator's data input and then selecting that output subspace defines a CP trace-nonincreasing map from its program to the selected output. Lemma 2 therefore applies directly to the quantities that can be compared with a trace-distance guarantee.

No successful branch is normalized in the argument. In particular, a very small success probability does not make the inequality inapplicable. [The spectral-diameter witness](dimension_cost.md#converting-the-conjugation-bound-into-channel-error) combines the success weight and correct overlap into one observable and produces a bound on the simulator's unconditional error.

The transposition fidelity itself is prior work. The contribution developed in this repository is the use of the reference-bearing program and a fixed-input dissipative evolution to make this restriction a dimension-dependent lower bound for every allowed sample-based Lindblad processor. See [Related work](../references/RELATED_WORK.md) for the full comparison.

---

[Repository overview](../README.md) · [Access model](model.md) · [Dimension cost](dimension_cost.md) · [Upper bound](upper_bound.md)
