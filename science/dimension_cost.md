# The dimension cost of a dissipative program

A quantum program for a jump operator contains the coefficients of that operator as a state. To reproduce the jump's dynamics, a processor must convert those coefficients into the correct action on an independent input. For a rank-one jump, the vectorized program encodes the decay direction in its complex-conjugate form. This page shows that this encoding imposes a sample cost proportional to the system dimension, even when the processor can use arbitrary collective operations.

The proof has three parts. A rank-one damping family converts a known input into a selected unknown output with a calculable probability. The [conjugation inequality](conjugation.md) limits how accurately a processor can obtain that output from its programs. A trace-distance witness retains the success probability and transfers the restriction to the full channel approximation.

## Statement in the sample-access model

Use the [access model](model.md): the system has dimension $`d`$, the jump satisfies $`\|L\|_{\mathrm F}=1`$, and its program is

```math
\pi_L=|L\rangle\langle L|,
\qquad
|L\rangle=(L\otimes I)|\Gamma\rangle,
\qquad
|\Gamma\rangle=\sum_{j=1}^{d}|j\rangle|j\rangle.
```

Here $`\|\cdot\|_{\mathrm F}`$ is the Frobenius (Hilbert-Schmidt) norm. The target generator and channel are

```math
\mathcal D_L(\rho)=L\rho L^\dagger
-\frac12\{L^\dagger L,\rho\},
\qquad \mathcal E_{L,t}=e^{t\mathcal D_L}.
```

A processor $`\mathcal P`$ may depend on $`d,n,t`$, but its dependence on the unknown $`L`$ is through the supplied $`n`$ program copies. Write

```math
\widetilde{\mathcal E}_{L,t}(\rho)
=\mathcal P(\rho\otimes\pi_L^{\otimes n}).
```

The error is the normalized diamond distance. Define

```math
\delta_d(n,t)
=\inf_{\mathcal P\,\mathrm{CPTP}}
\sup_{\|L\|_{\mathrm F}=1}
\frac12\left\|\widetilde{\mathcal E}_{L,t}
-\mathcal E_{L,t}\right\|_\diamond.
```

The minimum copy cost at accuracy $`\epsilon`$ is denoted by $`n_d^\star(t,\epsilon)`$. This definition is the unrestricted processor model used for sample-based simulation in [Park et al., Section II](https://arxiv.org/abs/2605.30301v1); it does not restrict the processor to wave matrix Lindbladization.

### Theorem 1: finite-copy dimension bound

For every integer $`d\ge4`$, every integer $`n\ge0`$, and every $`t\ge0`$,

```math
\boxed{\displaystyle
\delta_d(n,t)
\ge
\frac{(1-e^{-t/2})^2}{4}\,
\frac{d-3}{n+d-2}.}
```

The same bound already holds when the unknown jump is promised to have rank one and a fixed, known sink state. If a processor achieves error at most $`\epsilon>0`$, then necessarily

```math
n\ge
\frac{(1-e^{-t/2})^2(d-3)}{4\epsilon}-(d-2).
```

The integer sample count is also nonnegative, so a negative right-hand side imposes no additional restriction.

## The hard family

Decompose the system as

```math
\mathcal H_S
=\mathbb C|s\rangle\oplus\mathbb C|r\rangle\oplus\mathcal H,
\qquad \dim\mathcal H=D=d-2.
```

Choose the known sink $`|s\rangle`$ and reference $`|r\rangle`$ to be basis vectors. For every unit vector $`|\psi\rangle\in\mathcal H`$, set

```math
|\phi_\psi\rangle=\frac{|r\rangle+|\psi\rangle}{\sqrt2},
\qquad
L_\psi=|s\rangle\langle\phi_\psi|.
```

Every member has $`\|L_\psi\|_{\mathrm F}=1`$, $`L_\psi^\dagger L_\psi=P_{\phi_\psi}`$, and $`L_\psi^2=0`$. Its vectorized program is

```math
|L_\psi\rangle
=|s\rangle\otimes|\overline{\phi_\psi}\rangle
=|s\rangle\otimes
\frac{|r\rangle+|\overline\psi\rangle}{\sqrt2}.
```

The first program factor is fixed. The second is exactly the reference-bearing program treated in [Lemma 2 of Conjugation](conjugation.md#lemma-2-conjugation-from-a-reference-bearing-program).

The phase of $`\psi`$ relative to $`r`$ is allowed to vary. It changes the decay direction $`\phi_\psi`$ and is therefore a physical parameter of this family. This is why the Haar average includes the phase even though the selected output projector $`P_\psi`$ is phase invariant.

### Exact evolution of the known input

Let

```math
\beta(t)=1-e^{-t/2}.
```

The damping channel has the exact Kraus representation

```math
\begin{aligned}
K_{0,\psi}&=I-\beta(t)P_{\phi_\psi},\\
K_{1,\psi}&=\sqrt{1-e^{-t}}\,
|s\rangle\langle\phi_\psi|.
\end{aligned}
```

Indeed, the population in $`\phi_\psi`$ decays at rate one, coherences between that direction and its orthogonal complement decay at rate one half, and the lost population enters the sink. The Kraus normalization follows from $`2\beta-\beta^2=1-e^{-t}`$.

Apply this channel to the fixed state $`|r\rangle\langle r|`$. Its two Kraus vectors are

```math
\begin{aligned}
K_{0,\psi}|r\rangle
&=\left(1-\frac{\beta}{2}\right)|r\rangle
-\frac{\beta}{2}|\psi\rangle,\\
K_{1,\psi}|r\rangle
&=\sqrt{\frac{1-e^{-t}}2}\,|s\rangle.
\end{aligned}
```

Let $`Q`$ project onto $`\mathcal H`$. The selected output is exactly

```math
Q\mathcal E_{L_\psi,t}(|r\rangle\langle r|)Q
=p(t)P_\psi,
\qquad
p(t)=\frac{\beta(t)^2}{4}.
```

Thus the target evolution converts a known input into the unknown state $`\psi`$ in a fixed output subspace with probability $`p(t)`$. The supplied programs contain the conjugate direction. Both this probability and the program's coherent reference are retained in the lower bound.

## Converting the conjugation bound into channel error

Fix any allowed processor using $`n`$ programs and denote its output on the known input by

```math
\widetilde\rho_\psi
=\mathcal P\left(|r\rangle\langle r|\otimes
\pi_{L_\psi}^{\otimes n}\right).
```

Absorbing the known data input and the fixed sink factors into $`\mathcal P`$, the map from its remaining programs to $`Q\widetilde\rho_\psi Q`$ is CP and trace nonincreasing. Define

```math
s_\psi=\operatorname{Tr}(Q\widetilde\rho_\psi),
\qquad
f_\psi=\operatorname{Tr}(P_\psi\widetilde\rho_\psi),
\qquad
c=\frac{n+1}{n+D}.
```

Lemma 2 gives

```math
\int(f_\psi-cs_\psi)\,d\psi\le0.
```

Now choose the Hermitian witness

```math
X_\psi=P_\psi-cQ.
```

It has eigenvalue $`1-c`$ on $`\psi`$, eigenvalue $`-c`$ on the orthogonal part of $`\mathcal H`$, and eigenvalue zero on the sink and reference. Since $`D\ge2`$, its spectral diameter is one.

For any two density matrices $`\rho,\sigma`$ and any Hermitian $`X`$ with spectral diameter one,

```math
|\operatorname{Tr}[X(\rho-\sigma)]|
\le\frac12\|\rho-\sigma\|_1.
```

To see this, subtract the smallest eigenvalue of $`X`$ times the identity. The result is an effect between zero and identity, and $`\operatorname{Tr}(\rho-\sigma)=0`$.

Suppose the processor's worst-case normalized diamond error is $`\epsilon_{\mathcal P}`$. Diamond distance bounds the trace distance on the particular input $`|r\rangle`$, so for every $`\psi`$,

```math
\begin{aligned}
\epsilon_{\mathcal P}
&\ge
\operatorname{Tr}\!\left[
X_\psi\left(
\mathcal E_{L_\psi,t}(|r\rangle\langle r|)
-\widetilde\rho_\psi\right)\right]\\
&=p(t)(1-c)-(f_\psi-cs_\psi).
\end{aligned}
```

Averaging over the unit sphere and applying the conjugation inequality yields

```math
\epsilon_{\mathcal P}
\ge p(t)(1-c)
=\frac{\beta(t)^2}{4}\frac{D-1}{n+D}.
```

Substitute $`D=d-2`$ and take the infimum over all processors. This proves Theorem 1.

The witness depends on the unknown state only as a mathematical test used to lower-bound trace distance. The protocol is not granted this state-dependent measurement. A diamond-distance guarantee must control every such test, so its use in the converse does not add an access assumption.

## A matching short-time sample law

### Corollary 2

For $`d\ge4`$, $`0<t\le1`$, and

```math
0<\epsilon\le\frac{t^2}{256},
```

the minimum copy cost obeys

```math
\boxed{\displaystyle
n_d^\star(t,\epsilon)\ge\frac{dt^2}{512\epsilon}.}
```

Together with the universal construction in [Upper bound](upper_bound.md), this gives

```math
\boxed{\displaystyle
n_d^\star(t,\epsilon)=\Theta\!\left(\frac{dt^2}{\epsilon}\right)}
```

throughout the stated regime, with constants independent of $`d,t,\epsilon`$. The upper construction and its sufficient sample count are given on that page; the converse here applies to every processor in the model.

#### Proof of the explicit lower constant

For $`0<t\le1`$,

```math
\beta(t)=1-e^{-t/2}\ge\frac{t}{4}.
```

The assumed error therefore satisfies $`\epsilon\le\beta(t)^2/16`$. Set $`A=\beta(t)^2/(4\epsilon)\ge4`$. Theorem 1 implies

```math
n\ge A(D-1)-D
\ge\frac{A(D-1)}2
=\frac{\beta(t)^2(D-1)}{8\epsilon},
```

where the second inequality uses $`D\ge2`$ and $`D\le2(D-1)`$. Finally, $`\beta(t)^2\ge t^2/16`$ and $`D-1=d-3\ge d/4`$ give the displayed constant.

## Scope and limiting cases

The finite-copy theorem applies for every $`t\ge0`$. At $`t=0`$, the target is the identity channel and zero copies suffice; the bound is exactly zero. For positive time and $`d\ge4`$, its right-hand side is positive at every finite $`n`$, excluding an exact finite-copy processor for the entire promised family. The proof includes $`n=0`$ through the zero-copy case of the conjugation lemma.

At long times, $`\beta(t)^2`$ approaches one, and the lower bound approaches $`(d-3)/[4(n+d-2)]`$. The matching universal $`dt^2/\epsilon`$ result above is stated for $`0<t\le1`$. The finite-copy theorem does not establish quadratic growth in arbitrarily long evolution time. The construction uses an unknown subspace of dimension at least two, which is why the theorem begins at $`d=4`$.

Adaptive protocols are included when they have a deterministic cap of $`n`$ source copies. Such a protocol can be supplied all $`n`$ copies at the outset, leave unused copies untouched, and discard them at the end. Its measurement outcomes and feedforward are part of one CPTP processor. Internal heralding is also allowed, with all outcomes included in the channel whose error is assessed. The theorem does not replace this deterministic cap with an expected-copy budget.

The resource measured is the number of supplied program copies. Gate counts, preparation time and a classical description of the jump are separate access questions addressed in the [model](model.md). The result establishes a cost for this specific sample interface; it does not assign the same cost to a jump given through a different interface.

## What the dimension factor establishes

The hard family has a known sink and a rank-one jump, and its program is a product across the two vectorization registers. Its difficulty therefore persists without entanglement between those registers or a complicated jump spectrum. The restriction comes from converting the conjugate decay direction in the program into an action on an independent input.

The fixed reference component makes that restriction observable: the no-jump part of the evolution places the unconjugated unknown vector in a specified output subspace. Complete positivity limits the accuracy of producing that vector from conjugate programs. The witness turns this limit into an error bound for the entire simulator.

[Park et al.](https://arxiv.org/abs/2605.30301v1) establish a dimension-dependent lower bound for the WML algorithm. The theorem here places a matching dimension dependence on the unrestricted sample-access task in the stated regime. The conjugation fidelity used in its proof is credited to [Brzić et al.](https://arxiv.org/abs/2603.23628v3); the new implication concerns the dissipative programming task. [Related work](../references/RELATED_WORK.md) separates these contributions and the earlier sample-complexity results.

---

[Repository overview](../README.md) · [Access model](model.md) · [Conjugation bound](conjugation.md) · [Upper bound](upper_bound.md)
