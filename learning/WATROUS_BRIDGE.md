# From Watrous to the cost of a Lindblad program

The proof uses three connections: an operator becomes a supplied quantum state through vectorization; a desired dissipative channel becomes a probabilistic state-conjugation task; and a restriction on that task becomes a channel-error lower bound. This page develops those connections from the finite-dimensional language of John Watrous, *The Theory of Quantum Information* (Cambridge University Press, 2018).

The [author's book page](https://jhwatrous.github.io/books-and-courses.html) provides a manuscript and individual chapters. The locations below refer to that author-hosted text. The subsequent derivations use this repository's notation. The [model](../science/model.md) fixes the resource assumptions and theorem statements; the [milestones](MILESTONES.md) organize the background and proof narrative.

## Verified anchor locations

| Location in Watrous | Role in this bridge |
|---|---|
| [§1.1.2, Eqs. (1.127)–(1.134)](https://jhwatrous.github.io/TQI.pdf) | Operator-vector correspondence, its Hilbert–Schmidt isometry, and partial-trace identities |
| [§1.1.3](https://jhwatrous.github.io/TQI.pdf) | Operator decompositions and Schatten norms |
| [§§2.1.2–2.1.3](https://jhwatrous.github.io/TQI.2.pdf) | Density operators, pure states and their phase ambiguity, tensor products, reductions, and purifications |
| [§§2.2.1–2.2.2](https://jhwatrous.github.io/TQI.2.pdf) | Channels, Choi operators, Kraus representations, and Stinespring representations; Eq. (2.66) is the Choi reconstruction formula, and Theorem 2.22 characterizes complete positivity |
| [§3.1](https://jhwatrous.github.io/TQI.3.pdf) | State discrimination and the operational meaning of trace norm |
| [§§3.3.2–3.3.3](https://jhwatrous.github.io/TQI.3.pdf) | Completely bounded trace norm and channel discrimination |
| [§7.1.1, Proposition 7.1](https://jhwatrous.github.io/TQI.7.pdf) | Symmetric tensors and their projector |
| [§7.2.1, Lemma 7.24, Eq. (7.159)](https://jhwatrous.github.io/TQI.7.pdf) | Uniform spherical measure, Haar measure, and the normalized symmetric-projector moment identity |

The Lindblad evolution, the supplied-copy model, the reference-flag reduction, and the two simulation constructions are developed locally below. The [related-work account](../references/RELATED_WORK.md) identifies the primary research sources for the inherited conjugation and wave matrix Lindbladization results.

## Notation and resource model

A quantum state is a positive semidefinite matrix of trace one. A unit vector $`|\psi\rangle`$ represents the pure state $`P_\psi=|\psi\rangle\langle\psi|`$. Multiplying the vector by a phase leaves this projector unchanged. In a coherent sum with a fixed reference vector, however, changing only the unknown component's phase generally changes the state. That distinction will be used in the hard family.

For a joint state on registers $`A,B`$, discarding $`B`$ means applying the partial trace $`\operatorname{Tr}_B`$. On matrix units,

```math
\operatorname{Tr}_B
\bigl(|i\rangle\langle j|\otimes|k\rangle\langle\ell|\bigr)
=\delta_{k\ell}|i\rangle\langle j|.
```

This operation preserves the statistics of all measurements on $`A`$. Appending a state, applying a joint channel, and taking a partial trace are therefore the basic operations used to specify a processor.

Watrous writes the adjoint of a matrix as $`A^*`$; this repository writes $`A^\dagger`$. An overline denotes entrywise complex conjugation and $`A^{\mathsf T}`$ denotes transpose. These operations refer to the fixed computational basis. In particular, $`A^\dagger=\overline A^{\mathsf T}`$.

For singular values $`s_j`$, the three norms needed here are

```math
\|A\|_1=\sum_j s_j,
\qquad
\|A\|_{\mathrm F}
=\sqrt{\operatorname{Tr}(A^\dagger A)}
=\sqrt{\sum_j s_j^2},
\qquad
\|A\|_\infty=\max_j s_j.
```

The Frobenius norm is the Schatten 2-norm. We write its subscript as $`\mathrm F`$ so that the normalization of a program cannot be confused with the operator norm controlling a channel estimate. The inequality $`\|A\|_\infty\le\|A\|_{\mathrm F}`$ follows immediately from these formulas.

The unknown jump $`L`$ acts on $`\mathbb C^d`$ and satisfies $`\|L\|_{\mathrm F}=1`$. Its vectorized program is

```math
|L\rangle
=\operatorname{vec}(L)
=\sum_{i,j}L_{ij}|i\rangle|j\rangle
=(L\otimes I)|\Gamma_d\rangle,
\qquad
|\Gamma_d\rangle=\sum_{j=1}^{d}|j,j\rangle.
```

The vector $`\Gamma_d`$ is unnormalized. Directly summing matrix entries gives

```math
\langle L|L\rangle
=\sum_{i,j}|L_{ij}|^2
=\|L\|_{\mathrm F}^2
=1.
```

Using a normalized maximally entangled vector here would instead produce a vector of squared norm $`1/d`$. The unnormalized convention removes that extra factor from both the program and the Choi formulas.

A processor receives $`n`$ independent copies of $`\pi_L=|L\rangle\langle L|`$ and acts on an independent data input:

```math
\widetilde{\mathcal E}_{L,t}^{(n)}(X)
=\mathcal P_{d,n,t}(X\otimes\pi_L^{\otimes n}).
```

The channel $`\mathcal P_{d,n,t}`$ may depend on the known dimension, copy budget, and time. Its gates, ancillas, and measurements are independent of $`L`$. Joint operations on all programs are allowed. A protocol using measurements and feed-forward with a deterministic cap of $`n`$ copies is included by supplying those copies initially and keeping its classical outcomes inside the processor. The [resource model](../science/model.md#1-supplied-copy-access) explains why a preparation unitary, its inverse, and a classical description of $`L`$ are different interfaces.

## Channels and the Choi transpose

A channel is completely positive and trace preserving. Complete positivity ensures that the map remains positive when an untouched reference register is present. A selected outcome of a physical instrument is instead described by a completely positive, trace-nonincreasing map. Its output trace is the probability of that outcome.

For a linear map $`\mathcal T`$ from an input space to an output space, use the unnormalized, output-first Choi matrix

```math
J_{\mathcal T}
=\sum_{i,j}\mathcal T(|i\rangle\langle j|)
\otimes|i\rangle\langle j|.
```

The reconstruction formula corresponding to Watrous Eq. (2.66) is

```math
\mathcal T(X)
=\operatorname{Tr}_{\mathrm{in}}
\left[J_{\mathcal T}
\bigl(I_{\mathrm{out}}\otimes X^{\mathsf T}\bigr)\right].
```

The transpose is required: $`\operatorname{Tr}(|i\rangle\langle j|X^{\mathsf T})=X_{ij}`$, so the right-hand side is $`\sum_{ij}X_{ij}\mathcal T(|i\rangle\langle j|)`$. Omitting the transpose would reconstruct the action on $`X^{\mathsf T}`$.

The matrix $`J_{\mathcal T}`$ is positive exactly when $`\mathcal T`$ is CP. Indeed, a Kraus representation gives

```math
\mathcal T(X)=\sum_a A_aXA_a^\dagger,
\qquad
J_{\mathcal T}
=\sum_a|\operatorname{vec}(A_a)\rangle
\langle\operatorname{vec}(A_a)|.
```

Conversely, a positive-matrix decomposition of $`J_{\mathcal T}`$ can be reshaped into those Kraus operators. Trace preservation is $`\sum_aA_a^\dagger A_a=I`$, equivalently $`\operatorname{Tr}_{\mathrm{out}}J_{\mathcal T}=I`$. Replacing equality by an operator inequality gives the trace-nonincreasing condition.

The supplied program has a useful Choi interpretation:

```math
J_{X\mapsto LXL^\dagger}=|L\rangle\langle L|.
```

This is the Choi matrix of the CP jump branch. It is a normalized state because of Frobenius normalization, but it is generally not the Choi matrix of a trace-preserving channel. The dissipative target includes additional terms that restore trace preservation.

For a rank-one operator, vectorization gives

```math
\operatorname{vec}(|s\rangle\langle\phi|)
=|s\rangle\otimes|\bar\phi\rangle.
```

Thus the supplied hidden direction is complex conjugated. Later, when $`\overline P_\psi`$ enters the Choi reconstruction formula, its transpose becomes $`P_\psi`$. That cancellation places the output and input factors in the same symmetric tensor space.

## Distinguishability and the factor of two

For density operators, define trace distance by

```math
d_1(\rho,\sigma)=\frac12\|\rho-\sigma\|_1.
```

If $`0\le F\le I`$ is a measurement effect, then

```math
|\operatorname{Tr}[F(\rho-\sigma)]|
\le d_1(\rho,\sigma).
```

To derive the bound, write the traceless Hermitian difference as $`\rho-\sigma=\Delta_+-\Delta_-`$, where the positive and negative parts have orthogonal support. Each has trace $`\|\rho-\sigma\|_1/2`$. An effect cannot collect more weight than the positive part or remove more weight than the negative part. The projector onto the positive support attains the upper bound. With equal prior probabilities, optimal binary discrimination therefore succeeds with probability $`(1+d_1)/2`$.

The completely bounded trace norm in Watrous is the diamond norm used here. For the difference of two channels,

```math
d_\diamond(\mathcal A,\mathcal B)
=\frac12\|\mathcal A-\mathcal B\|_\diamond
=\max_{\rho_{SR}}
\frac12\left\|
\bigl[(\mathcal A-\mathcal B)\otimes\operatorname{id}_R\bigr]
(\rho_{SR})\right\|_1.
```

A reference of dimension equal to the input dimension suffices. Watrous's norm is the full norm; this repository's error is half that norm. Every conversion from a norm estimate to a stated error must retain this factor.

An upper bound must control every input and its reference correlations. For a lower bound, one fixed input and one measurement test already constrain a uniform channel guarantee. If a Hermitian observable $`X`$ has spectral diameter at most one, subtracting its smallest eigenvalue times the identity makes it an effect. Hence

```math
|\operatorname{Tr}[X(\rho-\sigma)]|
\le\frac12\|\rho-\sigma\|_1.
```

This form permits the witness used below to contain a negative eigenvalue without paying an additional factor of two.

## Symmetric tensors and Haar moments

Identical pure-state copies lie in the symmetric subspace: permuting their tensor factors leaves the vector unchanged. Write

```math
\mathsf S_k=\operatorname{Sym}^k(\mathbb C^D),
\qquad
\Pi_k=\frac1{k!}\sum_{\pi\in S_k}U_\pi,
\qquad
h_k=\dim\mathsf S_k=\binom{D+k-1}{k}.
```

The dimension counts occupation vectors $`(n_1,\ldots,n_D)`$ with nonnegative entries summing to $`k`$. Set $`\mathsf S_0=\mathbb C`$, $`\Pi_0=1`$, and $`h_0=1`$. The projector formula is Watrous Proposition 7.1. His Lemma 7.24, Eq. (7.159), gives

```math
\int P_\psi^{\otimes k}\,d\psi
=\frac{\Pi_k}{h_k},
```

where the integral uses normalized uniform measure on the complex unit sphere. Haar-random unitaries applied to any fixed unit vector produce this measure.

The normalization can also be checked directly with Gaussian coordinates. Let the entries of $`g\in\mathbb C^D`$ be independent standard complex Gaussians. Their moments satisfy $`\mathbb E[z^p\bar z^q]=0`$ for $`p\ne q`$ and $`\mathbb E|z|^{2p}=p!`$. Expanding matrix entries therefore gives

```math
\mathbb E\bigl[(gg^\dagger)^{\otimes k}\bigr]
=k!\Pi_k.
```

Writing $`g=R\psi`$, the radius and uniform direction are independent, and

```math
\mathbb E R^{2k}=D(D+1)\cdots(D+k-1)=k!h_k.
```

Dividing the two identities recovers the normalized moment. This also fixes the convention at $`k=0`$.

A vector symmetric in all $`k+1`$ registers is symmetric in its last $`k`$ registers. Consequently,

```math
\Pi_{k+1}\le I\otimes\Pi_k.
```

The two facts needed for conjugation are now visible: this inclusion of subspaces, and the dimension ratio

```math
\frac{h_k}{h_{k+1}}=\frac{k+1}{k+D}.
```

## Single-jump Lindblad dynamics

The target has one jump and no additional Hamiltonian:

```math
\mathcal D_L(X)
=LXL^\dagger-\frac12\{L^\dagger L,X\},
\qquad
\mathcal E_{L,t}=e^{t\mathcal D_L}.
```

Its local meaning follows from a short physical step. For sufficiently small $`h\ge0`$, take Kraus operators

```math
V_0=\sqrt{I-hL^\dagger L},
\qquad V_1=\sqrt h\,L.
```

They obey $`V_0^\dagger V_0+V_1^\dagger V_1=I`$. Expanding the square root shows that their channel is $`\mathcal I+h\mathcal D_L+O(h^2)`$. Repeating steps of duration $`t/m`$ and taking $`m\to\infty`$ gives the exponential in finite dimensions. The jump term transfers population; the anticommutator supplies the corresponding loss and coherence decay. Its trace cancels the jump term's trace.

For $`L=|s\rangle\langle\phi|`$ with unit vectors $`s\perp\phi`$, let $`P=|\phi\rangle\langle\phi|`$ and $`P_s=|s\rangle\langle s|`$. Then $`L^\dagger L=P`$, $`L^2=0`$, and the exact channel is

```math
\begin{aligned}
K_0&=I-\beta(t)P,\\
K_1&=\sqrt{1-e^{-t}}\,|s\rangle\langle\phi|,
\qquad \beta(t)=1-e^{-t/2},\\
\mathcal E_{\phi,t}(X)&=K_0XK_0^\dagger+K_1XK_1^\dagger.
\end{aligned}
```

The identity $`2\beta-\beta^2=1-e^{-t}`$ verifies Kraus normalization. To verify the dynamics itself, decompose operators into blocks along $`\phi`$ and its orthogonal complement. The $`P`$ population decays as $`e^{-t}`$ and is transferred to $`s`$. The off-diagonal blocks decay as $`e^{-t/2}`$. The pre-existing block fully orthogonal to $`\phi`$ is retained and additionally receives $`(1-e^{-t})\operatorname{Tr}(PX)P_s`$. These are precisely the solutions of the generator equation.

Expanding the same channel gives a useful form for the later construction:

```math
\mathcal E_{\phi,t}(X)
=X-\beta(PX+XP)
+\beta^2\operatorname{Tr}(PX)P
+(2\beta-\beta^2)\operatorname{Tr}(PX)P_s.
```

These formulas are derived here for the project. They are not being attributed to a Lindblad chapter of the selected book.

## The reference-flag program

Choose a known sink $`|s\rangle`$ and a known reference $`|r\rangle`$ to be distinct computational-basis vectors. Let the remaining basis vectors span a hidden subspace $`\mathcal H`$ of dimension $`D=d-2`$:

```math
\mathbb C^d
=\mathbb C|s\rangle\oplus\mathbb C|r\rangle\oplus\mathcal H.
```

For an arbitrary unit vector $`\psi\in\mathcal H`$, define

```math
|\phi_\psi\rangle
=\frac{|r\rangle+|\psi\rangle}{\sqrt2},
\qquad
L_\psi=|s\rangle\langle\phi_\psi|.
```

This family already satisfies the required Frobenius normalization. Its supplied state is

```math
|L_\psi\rangle
=|s\rangle\otimes|\chi_\psi\rangle,
\qquad
|\chi_\psi\rangle
=\frac{|r\rangle+|\bar\psi\rangle}{\sqrt2}.
```

The fixed sink factors can be absorbed into the processor. The informative part is a coherent reference plus the conjugated unknown vector.

Apply the desired channel to the fixed, known data input $`|r\rangle`$. The two Kraus vectors are

```math
K_0|r\rangle
=\left(1-\frac\beta2\right)|r\rangle
-\frac\beta2|\psi\rangle,
\qquad
K_1|r\rangle
=\sqrt{\frac{1-e^{-t}}2}\,|s\rangle.
```

Let $`Q`$ project onto $`\mathcal H`$. Selecting that fixed subspace gives

```math
Q\mathcal E_{L_\psi,t}(P_r)Q
=p(t)P_\psi,
\qquad p(t)=\frac{\beta(t)^2}{4}.
```

The input and the selecting projection are independent of $`\psi`$. The desired dynamics therefore produces an unconjugated unknown state in a known output subspace, with a known probability. The projection is a useful probe of the full channel; the simulation task remains unconditional.

## Phase averaging and occupation sectors

Replacing $`\psi`$ by $`e^{i\theta}\psi`$ does not change $`P_\psi`$, but it changes the relative phase in $`\phi_\psi`$ and $`\chi_\psi`$. All such vectors belong to the hard family. This permits a phase average inside a worst-case lower-bound proof.

The symmetric program space splits according to the number of factors in the hidden subspace:

```math
\operatorname{Sym}^n(\mathbb C|r\rangle\oplus\mathcal H)
\cong\bigoplus_{k=0}^n\operatorname{Sym}^k(\mathcal H).
```

Let $`V_k`$ insert $`n-k`$ reference factors and symmetrize with the normalization making it an isometry. Then

```math
|\chi_\psi\rangle^{\otimes n}
=\sum_{k=0}^n\sqrt{w_k}\,
V_k|\bar\psi\rangle^{\otimes k},
\qquad
w_k=2^{-n}\binom nk.
```

For example, with two copies the sectors have weights $`1/4,1/2,1/4`$: both factors are references, one is hidden, or both are hidden. In the middle sector the normalized vector is $`(|r\rangle|\bar\psi\rangle+|\bar\psi\rangle|r\rangle)/\sqrt2`$.

Under the phase change, the $`k`$th amplitude acquires $`e^{-ik\theta}`$. A density-matrix term between sectors $`k`$ and $`\ell`$ acquires $`e^{-i(k-\ell)\theta}`$. Its average vanishes unless $`k=\ell`$. Thus

```math
\int_0^{2\pi}
\bigl(P_{\chi_{e^{i\theta}\psi}}\bigr)^{\otimes n}
\frac{d\theta}{2\pi}
=\sum_{k=0}^n w_kV_k\overline P_\psi^{\otimes k}V_k^\dagger.
```

This is an identity for the averaged input ensemble. It places no symmetry assumption on the processor. Because the desired projector is phase invariant, it applies to both its average correct-output weight and its average success probability.

## Success-weighted conjugation

Consider a CP map $`\mathcal T`$ from $`k`$ copies of $`\mathcal H`$ to $`\mathcal H`$. Define

```math
f_\psi
=\operatorname{Tr}\bigl[P_\psi\mathcal T(\overline P_\psi^{\otimes k})\bigr],
\qquad
s_\psi
=\operatorname{Tr}\mathcal T(\overline P_\psi^{\otimes k}).
```

For a trace-nonincreasing branch, $`s_\psi`$ is success probability. If it is positive, the conditional probability of passing the $`P_\psi`$ test is $`f_\psi/s_\psi`$. The quantities $`f_\psi`$ and $`s_\psi`$ themselves are linear in the CP map, which is what a channel-error comparison needs.

Choi reconstruction, the transpose cancellation, and the Haar moment give

```math
\begin{aligned}
\int f_\psi\,d\psi
&=\frac1{h_{k+1}}\operatorname{Tr}(J_{\mathcal T}\Pi_{k+1}),\\
\int s_\psi\,d\psi
&=\frac1{h_k}\operatorname{Tr}
\bigl[J_{\mathcal T}(I\otimes\Pi_k)\bigr].
\end{aligned}
```

Taking the trace of $`\Pi_{k+1}\le I\otimes\Pi_k`$ against the positive matrix $`J_{\mathcal T}`$ yields

```math
\boxed{
\int f_\psi\,d\psi
\le\frac{k+1}{k+D}\int s_\psi\,d\psi.
}
```

The argument uses complete positivity and does not require trace preservation, covariance, or state-independent success probability. It includes $`k=0`$, where the ratio is $`1/D`$. It does not assert the same bound pointwise for every $`\psi`$.

For the actual flagged programs, compose the processor's selected branch with each $`V_k`$. The resulting map is CP on that sector; projecting onto the symmetric input space extends it to all $`k`$ input factors. Apply the displayed inequality sector by sector. Since

```math
\frac{k+1}{k+D}\le\frac{n+1}{n+D}
\quad\text{for }0\le k\le n,
```

and the sector success contributions are nonnegative, the flagged-program branch obeys

```math
\int f_\psi\,d\psi
\le c\int s_\psi\,d\psi,
\qquad c=\frac{n+1}{n+D}.
```

This is the [conjugation lemma](../science/conjugation.md) in the form used by the simulator converse. The universal transposition bound and symmetric-subspace argument are established ingredients, as described by [Brzić et al., Theorem 1 and Appendix B](https://arxiv.org/abs/2603.23628v3). The role of the present reduction is to identify a dissipative output branch to which that argument applies.

## The witness and the minimax bound

Fix any processor and supply its known data input $`P_r`$. Let its full output be $`\widetilde\rho_\psi`$. Selecting $`Q`$ defines a CP trace-nonincreasing map from its programs, so set

```math
s_\psi=\operatorname{Tr}(Q\widetilde\rho_\psi),
\qquad
f_\psi=\operatorname{Tr}(P_\psi\widetilde\rho_\psi).
```

The selected target has both success and correct-output weight equal to $`p(t)`$. The witness

```math
X_\psi=P_\psi-cQ
```

has eigenvalues $`1-c`$ on $`\psi`$, $`-c`$ on the other hidden directions, and zero on the reference and sink. For $`D\ge2`$, its spectral diameter is exactly one. Equivalently, $`X_\psi+cI=P_\psi+c(I-Q)`$ is an effect.

Let $`\varepsilon_{\mathcal P}`$ be this processor's worst-case half-diamond error. On the fixed input, the observable bound gives

```math
\begin{aligned}
\varepsilon_{\mathcal P}
&\ge\operatorname{Tr}
\left[X_\psi\left(\mathcal E_{L_\psi,t}(P_r)
-\widetilde\rho_\psi\right)\right]\\
&=p(t)(1-c)-(f_\psi-cs_\psi).
\end{aligned}
```

Average over $`\psi`$. The conjugation inequality makes the average of the parenthesized expression nonpositive. Therefore every processor satisfies

```math
\varepsilon_{\mathcal P}
\ge p(t)(1-c)
=\frac{\beta(t)^2}{4}\frac{D-1}{n+D}.
```

Taking the infimum over processors and substituting $`D=d-2`$ proves the [finite-copy converse](../science/dimension_cost.md#theorem-1-finite-copy-dimension-bound):

```math
\boxed{
\delta_d(n,t)
\ge\frac{(1-e^{-t/2})^2}{4}
\frac{d-3}{n+d-2}
}\qquad(d\ge4,\ n\ge0,\ t\ge0).
```

The unknown-dependent witness is a mathematical test of a uniform channel guarantee. The processor is not supplied that test or an unknown state. No branch has been divided by its success probability, and no square-root conversion between fidelity and trace distance has entered the argument.

## A four-dimensional worked example

Use the computational basis $`|0\rangle,|1\rangle,|2\rangle,|3\rangle`$ and set

```math
|s\rangle=|0\rangle,
\qquad |r\rangle=|1\rangle,
\qquad |\psi\rangle=a|2\rangle+b|3\rangle,
\qquad |a|^2+|b|^2=1.
```

The hidden subspace has dimension $`D=2`$. The decay direction, jump matrix, and supplied program are

```math
\begin{aligned}
|\phi_\psi\rangle
&=\frac{|1\rangle+a|2\rangle+b|3\rangle}{\sqrt2},\\
L_\psi
&=\frac1{\sqrt2}
\begin{pmatrix}
0&1&\bar a&\bar b\\
0&0&0&0\\
0&0&0&0\\
0&0&0&0
\end{pmatrix},\\
|L_\psi\rangle
&=|0\rangle\otimes
\frac{|1\rangle+\bar a|2\rangle+\bar b|3\rangle}{\sqrt2}.
\end{aligned}
```

The Frobenius norm is one because the sum of squared entries is $`(1+|a|^2+|b|^2)/2=1`$. For example, if $`a=1/\sqrt2`$ and $`b=i/\sqrt2`$, the hidden direction in the program contains $`-i`$, while the desired output direction contains $`+i`$.

Choose

```math
t=2\log2,
\qquad e^{-t/2}=\frac12,
\qquad e^{-t}=\frac14,
\qquad \beta=\frac12.
```

The exact Kraus operators are $`K_0=I-P_{\phi_\psi}/2`$ and $`K_1=(\sqrt3/2)L_\psi`$. On the known input,

```math
K_0|1\rangle=\frac34|1\rangle-\frac14|\psi\rangle,
\qquad
K_1|1\rangle=\sqrt{\frac38}|0\rangle.
```

The full output is consequently

```math
\rho_\psi
=\frac1{16}
\begin{pmatrix}
6&0&0&0\\
0&9&-3\bar a&-3\bar b\\
0&-3a&|a|^2&a\bar b\\
0&-3b&b\bar a&|b|^2
\end{pmatrix}.
```

Its trace is $`(6+9+1)/16=1`$. Its positivity also follows directly from the two Kraus vectors. With $`Q=|2\rangle\langle2|+|3\rangle\langle3|`$,

```math
Q\rho_\psi Q=\frac1{16}P_\psi.
```

For any processor using $`n`$ supplied programs, the conjugation coefficient is $`c=(n+1)/(n+2)`$. The witness has eigenvalues $`1/(n+2)`$ and $`-(n+1)/(n+2)`$ on the hidden plane, and zero on the remaining two directions. Its target expectation is

```math
\operatorname{Tr}(X_\psi\rho_\psi)
=\frac1{16(n+2)}.
```

Its average simulator expectation is nonpositive. Thus

```math
\boxed{
\delta_4(n,2\log2)\ge\frac1{16(n+2)}.
}
```

For a requested error $`\varepsilon`$, this implies $`n\ge1/(16\varepsilon)-2`$ whenever the right-hand side is positive. For example, error at most $`1/1600`$ requires at least $`98`$ copies. This is the all-time finite-copy bound. Since $`2\log2>1`$, this example does not invoke the separate $`0<t\le1`$ matching corollary.

## From WML to the matching short-time scale

The converse constrains every allowed processor. To determine the optimal order of the copy cost, a construction must also simulate every Frobenius-normalized jump. Wave matrix Lindbladization (WML) does so by consuming one fresh program in each short step. The method and its dimension-dependent analysis are inherited from the [primary WML literature](../references/RELATED_WORK.md#sample-access-and-the-dimension-gap). The repository's [upper-bound proof](../science/upper_bound.md) reconstructs a conservative sufficient coefficient.

### A known interaction reads one program

Let register 1 contain the data and registers 2,3 contain one program. Define the known joint operator

```math
M=\frac1{\sqrt d}
\bigl(I_1\otimes|\Gamma_d\rangle\langle\Gamma_d|_{23}\bigr)
\bigl(\operatorname{SWAP}_{12}\otimes I_3\bigr),
```

and the known generator

```math
\mathcal G(Y)=MYM^\dagger-\frac12\{M^\dagger M,Y\}.
```

One step is the channel

```math
\mathcal W_{L,\Delta}(X)
=\operatorname{Tr}_{23}
e^{\Delta\mathcal G}(X\otimes\pi_L).
```

Only the appended state depends on the unknown jump. Watrous's channel and dilation framework explains why this is an admissible fixed processor: append a supplied state, apply a known channel, and discard its environment. The copy model does not require the known channel's gate decomposition to be efficient.

The first-order identity can be understood from

```math
M(X\otimes\pi_L)M^\dagger
=LXL^\dagger\otimes
\frac{|\Gamma_d\rangle\langle\Gamma_d|}{d}.
```

The program contraction produces precisely the jump term. The two reduced anticommutator terms produce $`L^\dagger LX/2`$ and $`XL^\dagger L/2`$. Hence

```math
\operatorname{Tr}_{23}\mathcal G(X\otimes\pi_L)
=\mathcal D_L(X),
```

so the WML step and the desired semigroup have the same derivative at zero duration. These identities hold with an untouched reference because they are linear operator identities, not statements restricted to pure data states.

### Why the accumulated error is of order one over the copy count

Matching the first derivative makes the local error second order. Its dimension coefficient requires controlling higher terms after tracing out the programs. The explicit invariant-span calculation in [the upper-bound proof](../science/upper_bound.md#2-algebra-that-controls-the-higher-order-terms) gives, for $`a=\|L\|_\infty^2\le1`$,

```math
\left\|
\operatorname{Tr}_{23}\circ\mathcal G^k
\circ(X\mapsto X\otimes\pi_L)
\right\|_\diamond
\le2a(2d)^{k-1},
\qquad k\ge1.
```

The first application reduces to data maps with norm at most $`a`$. Subsequent applications remain in a finite operator span whose coefficient growth is at most $`2d`$ per application. This is the step that keeps the leading reduced remainder linear in dimension.

Summing the exponential series yields

```math
\left\|\mathcal W_{L,\Delta}
-\mathcal I-\Delta\mathcal D_L\right\|_\diamond
\le\frac ad\bigl(e^{2d\Delta}-1-2d\Delta\bigr)
\le4ad\Delta^2
```

when $`0\le\Delta<1/(2d)`$. The last inequality uses $`e^x-1-x\le x^2`$ for $`0\le x\le1`$. For the exact channel, $`\|\mathcal D_L\|_\diamond\le2a`$ and the integral Taylor remainder give

```math
\left\|e^{\Delta\mathcal D_L}
-\mathcal I-\Delta\mathcal D_L\right\|_\diamond
\le2a^2\Delta^2.
```

Combining the two and dividing by two to obtain the repository's error convention gives

```math
d_\diamond(\mathcal W_{L,\Delta},e^{\Delta\mathcal D_L})
\le(2d+1)\Delta^2.
```

Now set $`\Delta=t/n`$ and use a fresh copy at every step. A telescoping expansion of the difference between two $`n`$-fold channel compositions has $`n`$ terms. Every surrounding channel has diamond norm one, so

```math
d_\diamond\left(\mathcal W_{L,t/n}^{\,n},e^{t\mathcal D_L}\right)
\le n(2d+1)(t/n)^2
=\frac{(2d+1)t^2}{n},
\qquad n>2dt.
```

An explicit sufficient integer count is therefore

```math
N_d(t,\varepsilon)
=\max\left\{
\lfloor2dt\rfloor+1,
\left\lceil\frac{(2d+1)t^2}{\varepsilon}\right\rceil
\right\}.
```

The first term enforces the step-size condition. The second enforces the requested accuracy. A small-step expansion without that first condition would not establish a uniform sample guarantee.

### Matching the converse in its stated regime

For $`0<t\le1`$, the elementary bound $`\beta(t)\ge t/4`$ holds. If also $`0<\varepsilon\le t^2/256`$, then $`\varepsilon\le\beta(t)^2/16`$. Set $`A=\beta(t)^2/(4\varepsilon)\ge4`$. The finite-copy converse gives

```math
n\ge A(D-1)-D
\ge\frac{A(D-1)}2
=\frac{\beta(t)^2(D-1)}{8\varepsilon}.
```

Here $`D\ge2`$ implies $`D\le2(D-1)`$. Using $`\beta(t)^2\ge t^2/16`$ and $`D-1=d-3\ge d/4`$ now gives

```math
n_d^\star(t,\varepsilon)
\ge\frac{dt^2}{512\varepsilon}.
```

In this regime the ceiling term in $`N_d`$ already enforces the strict step condition, and $`N_d\le3dt^2/\varepsilon`$. Consequently,

```math
\boxed{
\frac{dt^2}{512\varepsilon}
\le n_d^\star(t,\varepsilon)
\le\frac{3dt^2}{\varepsilon}
}
\quad
\left(d\ge4,\ 0<t\le1,\ 0<\varepsilon\le\frac{t^2}{256}\right).
```

The upper construction is sequential, but its existence is enough for the unrestricted minimax upper bound. The converse already permits collective processing. Together they show that allowing more general processors cannot remove the dimension factor in this regime.

## A collective processor for a known sink

The hard family also suggests a separate constructive question. Suppose the jump is promised to be $`L=|s\rangle\langle\phi|`$, where the sink is known and $`\phi`$ is any unknown unit vector orthogonal to it. This is a larger family than the flagged hard instances, but smaller than the set of all single jumps. Write

```math
\mathcal A=s^\perp,
\qquad m=\dim\mathcal A=d-1,
\qquad P=|\phi\rangle\langle\phi|.
```

The symbol $`m`$ refers to this entire active space. It differs from the converse's $`D=d-2`$, which reserves one active direction as a known reference.

For the construction, take $`m\ge2`$ and an integer $`n\ge1`$. After discarding the fixed sink factors, the programs are $`|\bar\phi\rangle^{\otimes n}`$. They lie in $`\operatorname{Sym}^n(\overline{\mathcal A})`$, where $`\overline{\mathcal A}`$ is the coordinatewise conjugate of the known active space. Choose a known orthonormal data basis $`\{|i\rangle\}_{i=1}^m`$ for $`\mathcal A`$ and the conjugate basis $`\{|\bar i\rangle\}_{i=1}^m`$ for the program modes. Thus $`\phi=\sum_i\phi_i|i\rangle`$ corresponds to $`\bar\phi=\sum_i\bar\phi_i|\bar i\rangle`$. The following operations are known linear maps on that symmetric space. An occupation basis makes them explicit; no physical bosonic encoding is assumed.

### A known contraction creates the unknown overlap

For an occupation vector $`\mathbf k=(k_1,\ldots,k_m)`$ with total occupation $`n`$, define

```math
a_i|\mathbf k\rangle
=\sqrt{k_i}|\mathbf k-\mathbf e_i\rangle.
```

The map $`a_i`$ lowers the degree by one. On a coherent tensor power it obeys

```math
a_i|\bar\phi\rangle^{\otimes n}
=\sqrt n\,\bar\phi_i
|\bar\phi\rangle^{\otimes(n-1)}.
```

Define a contraction from data and programs to degree $`n-1`$ memory:

```math
C=\sum_{i=1}^m\langle i|\otimes a_i,
\qquad N=n+m-1.
```

It acts as zero on sink inputs. For an active data vector $`v`$,

```math
C\bigl(|v\rangle\otimes|\bar\phi\rangle^{\otimes n}\bigr)
=\sqrt n\,\langle\phi|v\rangle
|\bar\phi\rangle^{\otimes(n-1)}.
```

The unknown overlap appears through the supplied program, while the matrix $`C`$ is independent of $`\phi`$. On degree $`n-1`$ memory, the occupation rules imply

```math
CC^\dagger
=\sum_i a_i a_i^\dagger
=(n-1+m)I=NI.
```

Thus

```math
\Pi=\frac{C^\dagger C}{N},
\qquad
J=\frac{|s\rangle\otimes C}{\sqrt N}
```

satisfy $`\Pi^2=\Pi`$ and $`J^\dagger J=\Pi`$. The projector $`\Pi`$ acts jointly on data and programs; it is known even though the desired data projector $`P`$ is unknown.

### A valid joint channel and its reduced action

For $`0\le\alpha\le2`$, define Kraus operators

```math
A_0=I-\alpha\Pi,
\qquad
A_1=\sqrt{\alpha(2-\alpha)}\,J.
```

Their output memories have degrees $`n`$ and $`n-1`$, respectively; embed those spaces in an orthogonal direct sum. Since $`\Pi`$ is a projector,

```math
A_0^\dagger A_0+A_1^\dagger A_1
=I-2\alpha\Pi+\alpha^2\Pi
+\alpha(2-\alpha)\Pi
=I.
```

Discarding the output memory therefore defines a channel. A CPTP extension to invalid, nonsymmetric program inputs is obtained by first projecting onto the known symmetric subspace and using any fixed trace-and-prepare channel on its complement. Valid supplied programs enter the symmetric outcome with certainty.

The interval $`1<\alpha\le2`$ is physically allowed: the surviving amplitude on the range of $`\Pi`$ is then negative. Complete positivity depends on the Kraus representation and normalization, not on that amplitude being positive. The allowed processor class does not require this joint operation to be a positive-time damping semigroup.

Let $`\sigma_\phi=\overline P_\phi^{\otimes n}`$ and $`q=n/N`$. Contraction with the coherent program gives

```math
\begin{aligned}
\operatorname{Tr}_{\mathrm{prog}}
\bigl[\Pi(X\otimes\sigma_\phi)\bigr]&=qPX,\\
\operatorname{Tr}_{\mathrm{prog}}
\bigl[(X\otimes\sigma_\phi)\Pi\bigr]&=qXP,\\
\operatorname{Tr}_{\mathrm{prog}}
\bigl[J(X\otimes\sigma_\phi)J^\dagger\bigr]
&=q\operatorname{Tr}(PX)P_s,\\
\operatorname{Tr}_{\mathrm{prog}}
\bigl[\Pi(X\otimes\sigma_\phi)\Pi\bigr]
&=\operatorname{Tr}(PX)
\left[q^2P+\frac{q^2}{n}(I_{\mathcal A}-P)\right].
\end{aligned}
```

The quadratic expression is the source of the residual error. Its occupation-operator moment is

```math
\left\langle a_k^\dagger a_\ell a_i^\dagger a_j\right\rangle
=n(n-1)\phi_k\phi_i\bar\phi_\ell\bar\phi_j
+n\delta_{\ell i}\phi_k\bar\phi_j.
```

The first term comes from lowering two occupied factors; the second comes from commuting $`a_\ell`$ past $`a_i^\dagger`$. Substitution into the matrix entries of the reduced projector expression gives the stated formula. These representation contractions are local derivations, not results assigned to a Watrous section.

Choose

```math
\alpha=\frac\beta q
=\beta\frac{n+m-1}{n}.
```

The admissibility condition is

```math
n\ge\frac{\beta(m-1)}{2-\beta}.
```

In particular, $`n\ge m-1`$ suffices at every time. Substituting the contractions into the two-Kraus channel makes the terms linear in $`P`$ agree exactly with the target. The remaining difference is

```math
\widetilde{\mathcal E}^{(n)}_{\phi,t}(X)
-\mathcal E_{\phi,t}(X)
=\frac{\beta^2}{n}\operatorname{Tr}(PX)
\bigl[I_{\mathcal A}-P-(m-1)P_s\bigr].
```

The error redirects part of the intended sink population into the active subspace perpendicular to $`\phi`$.

### Exact error with an arbitrary reference

Put

```math
\eta=\frac{\beta^2(m-1)}n,
\qquad
\tau_\perp=\frac{I_{\mathcal A}-P}{m-1}.
```

For any joint data-reference density operator, the difference channel produces

```math
\eta(\tau_\perp-P_s)\otimes\omega_R,
\qquad
\omega_R
=(\langle\phi|\otimes I_R)\rho_{SR}
(|\phi\rangle\otimes I_R).
```

The operator $`\omega_R`$ is positive and has trace at most one. Since $`\tau_\perp`$ and $`P_s`$ have orthogonal supports, their difference has trace norm two. The output difference therefore has trace norm at most $`2\eta`$, and the data input $`P`$ attains equality. Restoring the half-diamond convention gives the exact error

```math
\boxed{
d_\diamond\bigl(
\widetilde{\mathcal E}^{(n)}_{\phi,t},\mathcal E_{\phi,t}
\bigr)
=\frac{\beta(t)^2(d-2)}n.
}
```

This construction uses $`n`$ supplied copies coherently. Its symmetric memory dimension is $`\binom{n+m-1}{n}`$. Copy sufficiency does not by itself give a low gate count or a small coherent memory implementation.

### The all-time law and its domain

For $`d\ge4`$, $`t>0`$, and $`0<\varepsilon\le\beta(t)^2/16`$, take

```math
n=\left\lceil\frac{(d-2)\beta(t)^2}{\varepsilon}\right\rceil.
```

This automatically meets the CPTP condition. The flagged hard family is contained in this known-sink family, so its converse also applies. The resulting [all-time copy law](../science/rank_one_processor.md#5-matching-copy-law-for-this-family-at-all-times) is

```math
\frac{d\beta(t)^2}{32\varepsilon}
\le n_{d,s}^\star(t,\varepsilon)
\le
\left\lceil\frac{(d-2)\beta(t)^2}{\varepsilon}\right\rceil.
```

As time grows, $`\beta(t)`$ approaches one, and this family's cost saturates at order $`d/\varepsilon`$. At $`t=0`$ the identity channel needs no programs. The construction assumes a known sink and a rank-one jump; the matching result for arbitrary jumps remains the separate short-time theorem.

## Scientific dependencies and contribution

Watrous supplies the general framework of states, CP maps, Choi operators, operational norms, and symmetric tensors. The universal transposition limit and WML construction are established research ingredients. The reference-flag reduction identifies a fixed-input dissipative output whose success-weighted correct overlap is restricted by that limit. The witness then forces the dimension factor for every processor in the supplied-copy model.

The known-sink construction provides explicit collective sufficiency for the promised family and identifies its all-time dependence. Its formulas are proved in [the construction page](../science/rank_one_processor.md); the [related-work page](../references/RELATED_WORK.md) bounds the attribution claims. The [verification suite](../verification/README.md) checks finite-dimensional identities and normalizations. The general minimax conclusions follow from the analytical arguments and their quantifiers.

---

[Milestones](MILESTONES.md) · [Model](../science/model.md) · [Dimension converse](../science/dimension_cost.md) · [WML upper bound](../science/upper_bound.md) · [Collective processor](../science/rank_one_processor.md)
