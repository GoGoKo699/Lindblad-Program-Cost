# Related work

## Sample access and the dimension gap

The resource model supplies independent copies of the normalized vectorized jump state

```math
|L\rangle=(L\otimes I)|\Gamma\rangle,
\qquad |\Gamma\rangle=\sum_{j=1}^{d}|j,j\rangle,
\qquad \|L\|_F=1.
```

A processor independent of the unknown jump operator uses at most $`n`$ copies to approximate $`e^{t\mathcal D_L}`$, where $`\mathcal D_L(\rho)=L\rho L^\dagger-\{L^\dagger L,\rho\}/2`$. Error is half the diamond norm. The quantities $`\delta_d(n,t)`$ and $`n_d^*(t,\varepsilon)`$ refer to the minimax error and corresponding sample complexity in this model. A classical description of $`L`$, a state-preparation unitary and its inverse, and additional complex-conjugated versions of the supplied programs are different resources.

[Patel and Wilde, *Wave Matrix Lindbladization I: Quantum Programs for Simulating Markovian Dynamics*](https://doi.org/10.1142/S1230161223500105) introduced this pure-state encoding and the WML simulation method. The subsequent dimension-dependent analyses below provide the comparison for the present result. See bibliography key `PatelWilde2023`.

### Go et al.: unrestricted lower bound without a growing dimension factor

[Go et al., arXiv:2412.02134v3](https://arxiv.org/abs/2412.02134v3), published in *Quantum Science and Technology* **10**, 045058 (2025), define sample complexity over arbitrary fixed quantum processors. Their Theorem 12 proves

```math
n_d^*(t,\varepsilon)\ge 10^{-4}t^2/\varepsilon
```

for $`0<\varepsilon\le\min\{0.039,0.013t\}`$. The proof uses two qubit jumps $`L_\varphi=\operatorname{diag}(1,e^{i\varphi})/\sqrt2`$, together with a reduction from simulated-channel discrimination to discrimination of the supplied programs. It therefore already permits collective processing. Its hard instances have fixed dimension. Their WML analysis gives an $`O(d^2t^2/\varepsilon)`$ upper bound. Theorem numbers and conditions here refer to v3, dated 29 September 2025. Bibliography key: `Go2025`.

### Park et al.: linear-dimensional WML performance

[Park et al., arXiv:2605.30301v1](https://arxiv.org/abs/2605.30301v1), dated 28 May 2026, refine WML to the scaling

```math
O\!\left(d\|L\|_\infty^2t^2/\varepsilon\right).
```

Theorem 3 and its supporting Lemma 2 supply the upper bound; the lemma assumes the step condition $`n>2dt`$. Frobenius normalization implies $`\|L\|_\infty^2\le1`$, giving worst-case order $`O(dt^2/\varepsilon)`$ in the relevant accuracy regime. Theorem 10 supplies a matching dimension-dependent lower bound for WML on an explicit rank-one jump. It analyzes that algorithm's output. Section VII separately conjectures that the remaining factor of $`d`$ is also necessary for arbitrary processors. Their typical-case analysis concerns normalized random jumps with small operator norm; a worst-case converse does not turn that statement into a typical-case lower bound. Bibliography key: `Park2026`.

## Transposition and conjugation as an inherited ingredient

[Brzić, Grinko, Studziński, and Quintino, arXiv:2603.23628v3](https://arxiv.org/abs/2603.23628v3), dated 26 September 2026, prove in Theorem 1 that optimal deterministic $`N`$-to-$`K`$ universal pure-state transposition has worst-case fidelity

```math
\frac{\dim\operatorname{Sym}^N(\mathbb C^D)}
{\dim\operatorname{Sym}^{N+K}(\mathbb C^D)}.
```

Their Appendix B derives the Haar-average formulation and its Choi-operator optimization. For $`K=1`$, the fidelity is $`(N+1)/(N+D)`$. This formula and the underlying symmetric-subspace argument are inherited. Bibliography key: `Brzic2026`.

Earlier results include the optimal qubit universal-NOT fidelity $`(N+1)/(N+2)`$ of [Bužek, Hillery, and Werner (1999)](https://doi.org/10.1103/PhysRevA.60.R2626), equivalent up to a fixed unitary to qubit conjugation, and the optimal one-copy qudit transposition fidelity $`2/(D+1)`$ of [Buscemi et al. (2003)](https://doi.org/10.1016/S0375-9601(03)00954-X), Section 2, Eqs. (14)–(16). Bibliography keys: `Buzek1999`, `Buscemi2003`.

The proof in this repository uses the same positive-operator comparison for a completely positive, trace-nonincreasing map $`T`$:

```math
\int\langle\psi|T(\overline{P_\psi}^{\otimes k})|\psi\rangle\,d\psi
\le \frac{k+1}{k+D}
\int\operatorname{Tr}T(\overline{P_\psi}^{\otimes k})\,d\psi.
```

This is the success-weighted formulation needed for a postselected branch, including state-dependent success probabilities. It follows directly from positivity of the Choi operator and $`\Pi_{\mathrm{sym}}^{(k+1)}\le\Pi_{\mathrm{sym}}^{(k)}\otimes I`$. It is not claimed as a separate new conjugation limit.

## Other programming and sample-complexity models

### Jing et al.: programming a fixed semigroup through time

[Jing et al., *Programmable Open Quantum Systems*](https://doi.org/10.1103/yqlr-2dhr), *Physical Review Letters* **137**, 040403 (23 July 2026), study a fixed generator with a time-varying family of programs $`\pi_t`$ and one retrieval map acting throughout a time interval. The detailed model in [arXiv:2512.08279v1](https://arxiv.org/abs/2512.08279v1), dated 9 December 2025, includes finite-dimensional exact programmability and Hermitian-preserving retrieval maps simulated through quasiprobability sampling. Equation (1) and supplementary Definition S6 quantify programming cost through the logarithm of the minimum retrieval-map diamond norm. This resource differs from the number of independent vectorized copies encoding an unknown jump. Their amplitude-damping programmability obstruction does not supply a lower bound on $`n_d^*(t,\varepsilon)`$. Bibliography key: `Jing2026`; the equation and definition labels refer to the stated preprint version.

### Wei et al.: exponentiating the output of a known linear map

[Wei et al., *Simulating non-completely positive actions via exponentiation of Hermitian-preserving maps*](https://doi.org/10.1038/s41534-024-00949-z), *npj Quantum Information* **10**, 134 (30 December 2024), consume copies of an unknown state $`\rho`$ to implement the unitary generated by $`\mathcal N(\rho)`$, for a specified Hermitian-preserving map $`\mathcal N`$. Their published Theorem 5 gives an arbitrary-processor lower bound $`\Omega(R_*^2t^2/\epsilon)`$, where $`R_*`$ is an optimized output spectral diameter over a stated feasible set. The proof uses state and channel discrimination. This is an antecedent for sample-based simulation converses; its target is Hamiltonian evolution, rather than the dissipative channel generated by an unknown vectorized jump. Bibliography key: `Wei2024`.

The original [arXiv:2308.07956v1](https://arxiv.org/abs/2308.07956v1), dated 15 August 2023, has the title *Realizing Non-Physical Actions through Hermitian-Preserving Map Exponentiation*, an additional author, and labels the corresponding result Theorem 4. It has a separate bibliography entry, `Wei2023Preprint`, so that theorem numbering and author lists remain version-specific.

## Contribution and scope

The additional step is a reduction from simulated dissipation to the success-weighted conjugation task. For $`D=d-2`$, a fixed sink $`|s\rangle`$, a fixed reference $`|r\rangle`$, and an unknown unit vector $`|\psi\rangle\in\mathbb C^D`$, use

```math
L_\psi=|s\rangle\langle\phi_\psi|,
\qquad
|\phi_\psi\rangle=(|r\rangle+|\psi\rangle)/\sqrt2.
```

The supplied program contains $`|\overline{\phi_\psi}\rangle`$. Evolving the fixed input $`|r\rangle`$ and projecting onto the unknown-vector subspace produces the unnormalized state $`(1-e^{-t/2})^2P_\psi/4`$. Averaging the program's occupation sectors over the phase of $`\psi`$, then using the conjugation inequality and a spectral-diameter-one observable, gives

```math
\delta_d(n,t)\ge
\frac{(1-e^{-t/2})^2}{4}\,
\frac{d-3}{n+d-2},
\qquad d\ge4.
```

This applies to every fixed processor with a deterministic budget of $`n`$ program copies, including collective and adaptive processing. It supplies the missing dimension factor in the stated short-time, sufficiently accurate regime. It does not assert matching quadratic time dependence at arbitrarily long times, a typical-case lower bound, or a result for stronger oracle-access models. Relative to the precise results above, the contribution is the reference-flag dissipative reduction and the resulting processor-independent dimension lower bound.
