# A collective processor for a known sink

For a rank-one jump with a known sink, the dimension cost can also be attained by an explicit collective processor. This construction gives a matching time dependence for this family at arbitrary evolution times. The general-jump upper bound remains the separate [wave matrix Lindbladization result](upper_bound.md).

## 1. Family and target channel

Let

```math
\mathcal H=\mathbb C|s\rangle\oplus\mathcal A,
\qquad
\dim\mathcal A=m=d-1,
```

where the sink $`|s\rangle`$ and the subspace $`\mathcal A`$ are known. The unknown jump is

```math
L_\phi=|s\rangle\langle\phi|,
\qquad
|\phi\rangle\in\mathcal A,\quad
\langle\phi|\phi\rangle=1.
```

A program is $`|s\rangle|\bar\phi\rangle`$. Its fixed sink factor can be discarded. Write $`P=|\phi\rangle\langle\phi|`$, $`P_s=|s\rangle\langle s|`$, and $`\beta=1-e^{-t/2}`$. The exact channel is

```math
\mathcal E_{\phi,t}(X)
=
X-\beta(PX+XP)+\beta^2PXP
+(2\beta-\beta^2)\operatorname{Tr}(PX)P_s.
```

This follows from the Kraus operators $`I-\beta P`$ and $`\sqrt{2\beta-\beta^2}\,|s\rangle\langle\phi|`$. Since $`P`$ has rank one, $`PXP=\operatorname{Tr}(PX)P`$.

The basis used to conjugate $`|\phi\rangle`$ is the same known basis used in the program vectorization.

## 2. A parameter-independent contraction

Take an integer $`n\ge1`$. The $`n`$ program copies lie in $`\operatorname{Sym}^n(\overline{\mathcal A})`$. This symmetric subspace is known independently of $`\phi`$. It can be used as a mathematical description of the original registers; no efficient compression circuit is assumed.

To define the processor on the full tensor-product program space, first project onto this known symmetric subspace and its orthogonal complement. Use the construction below on the symmetric outcome and any fixed trace-and-prepare channel on the other outcome. This is a CPTP extension, and every valid supplied program enters the symmetric outcome with certainty.

Let $`a_i`$ be the annihilation maps from the degree-$`n`$ symmetric space to degree $`n-1`$. On data in $`\mathcal A`$, define

```math
C=\sum_{i=1}^{m}\langle i|\otimes a_i,
\qquad
N=n+m-1.
```

Extend $`C`$ by zero on sink inputs. The bosonic commutation relations give, on the degree-$`n-1`$ space,

```math
CC^\dagger
=
\sum_{i=1}^{m}a_i a_i^\dagger
=
N I.
```

Consequently,

```math
\Pi=\frac{C^\dagger C}{N}
```

is an orthogonal projector. The partial isometry

```math
J=\frac{|s\rangle\otimes C}{\sqrt N}
```

satisfies $`J^\dagger J=\Pi`$. Its output program register has degree $`n-1`$.

For any real $`0\le\alpha\le2`$, the two operators

```math
A_0=I-\alpha\Pi,
\qquad
A_1=\sqrt{\alpha(2-\alpha)}\,J
```

define a trace-preserving operation because

```math
A_0^\dagger A_0+A_1^\dagger A_1
=
I.
```

The degree-$`n`$ and degree-$`n-1`$ output memory spaces can be placed in an orthogonal direct sum. Discarding that memory gives the system output. Thus this is a finite, fixed CPTP processor using exactly the supplied $`n`$ copies.

The interval $`\alpha\in[0,2]`$ is material. For $`\alpha>1`$, the surviving amplitude on $`\Pi`$ is negative, but the operation remains CPTP. Requiring the joint processor itself to be a positive-time amplitude-damping semigroup would impose an unnecessary restriction on the allowed processor class.

## 3. Program contractions

Let

```math
\sigma_\phi
=
\bigl(|\bar\phi\rangle\langle\bar\phi|\bigr)^{\otimes n},
\qquad
q=\frac{n}{N}.
```

The annihilation rule is

```math
a_i|\bar\phi\rangle^{\otimes n}
=
\sqrt n\,\bar\phi_i
|\bar\phi\rangle^{\otimes(n-1)}.
```

For every system operator $`X`$, this gives

```math
\operatorname{Tr}_{\mathrm{prog}}
\bigl[\Pi(X\otimes\sigma_\phi)\bigr]
=
qPX,
\qquad
\operatorname{Tr}_{\mathrm{prog}}
\bigl[(X\otimes\sigma_\phi)\Pi\bigr]
=
qXP.
```

The quadratic contraction is

```math
\operatorname{Tr}_{\mathrm{prog}}
\bigl[\Pi(X\otimes\sigma_\phi)\Pi\bigr]
=
\operatorname{Tr}(PX)
\left[
q^2P+\frac{q^2}{n}(I_{\mathcal A}-P)
\right].
```

To verify the last identity, its $`i,\ell`$ matrix element is

```math
\frac1{N^2}
\sum_{j,k}X_{jk}
\langle a_k^\dagger a_\ell a_i^\dagger a_j\rangle_{\sigma_\phi}.
```

Commuting $`a_\ell`$ past $`a_i^\dagger`$ gives

```math
\langle a_k^\dagger a_\ell a_i^\dagger a_j\rangle_{\sigma_\phi}
=
n(n-1)\phi_k\phi_i\bar\phi_\ell\bar\phi_j
+
n\delta_{\ell i}\phi_k\bar\phi_j.
```

Summing over $`j,k`$ proves the formula. These equalities include sink coherences: the projector acts as zero on the sink, while the linear terms retain the appropriate coherences with the active space.

## 4. Exact error channel

Choose

```math
\alpha=\frac{\beta}{q}
=
\beta\,\frac{n+m-1}{n}.
```

The CPTP condition $`\alpha\le2`$ is equivalent to

```math
n\ge\frac{\beta(m-1)}{2-\beta}.
```

In particular, $`n\ge m-1`$ is sufficient at every time.

The preceding contractions show that the reduced processor channel $`\widetilde{\mathcal E}^{(n)}_{\phi,t}`$ obeys the exact identity

```math
\widetilde{\mathcal E}^{(n)}_{\phi,t}(X)
-
\mathcal E_{\phi,t}(X)
=
\frac{\beta^2}{n}\operatorname{Tr}(PX)
\bigl[I_{\mathcal A}-P-(m-1)P_s\bigr].
```

Every linear coherence term matches the target. The remaining error moves a small amount of the intended sink population into the active subspace orthogonal to $`\phi`$.

For $`m\ge2`$, define

```math
\delta_{\mathrm{coll}}
=
\frac{\beta^2(m-1)}{n},
\qquad
\tau_\perp=\frac{I_{\mathcal A}-P}{m-1}.
```

The error map is $`X\mapsto\delta_{\mathrm{coll}}\operatorname{Tr}(PX)(\tau_\perp-P_s)`$. The two output states have orthogonal supports. Therefore,

```math
\boxed{
\frac12
\left\|
\widetilde{\mathcal E}^{(n)}_{\phi,t}
-
\mathcal E_{\phi,t}
\right\|_\diamond
=
\frac{\beta^2(m-1)}{n}.
}
```

For completeness, on a joint system-reference density operator the error output is $`\delta_{\mathrm{coll}}(\tau_\perp-P_s)\otimes\omega_R`$, where

```math
\omega_R
=
\operatorname{Tr}_{\mathcal H}
\bigl[(P\otimes I_R)\rho_{\mathcal HR}\bigr]
\ge0,\qquad
\operatorname{Tr}\omega_R\le1.
```

Its trace norm is at most $`2\delta_{\mathrm{coll}}`$. A system input $`P`$ attains equality. Thus the error calculation includes arbitrary reference entanglement.

## 5. Matching copy law for this family at all times

Let $`n^\star_{d,s}(t,\varepsilon)`$ denote the minimax copy count over this known-sink rank-one family. The processor may know $`|s\rangle`$ and its orthogonal complement, but not $`|\phi\rangle`$.

For $`d\ge4`$, $`t>0`$, and

```math
0<\varepsilon\le\frac{\beta(t)^2}{16},
```

the construction uses

```math
n
=
\left\lceil
\frac{(d-2)\beta(t)^2}{\varepsilon}
\right\rceil
```

copies. This automatically satisfies $`n\ge d-2=m-1`$, so its CPTP condition holds.

The [hard family in the universal converse](dimension_cost.md) is a subset of this known-sink family: choose a known reference in $`\mathcal A`$ and let the unknown component range through its orthogonal complement. Hence

```math
n^\star_{d,s}(t,\varepsilon)
\ge
\frac{\beta(t)^2(d-3)}{4\varepsilon}-(d-2)
\ge
\frac{d\,\beta(t)^2}{32\varepsilon}.
```

The last inequality uses $`\beta(t)^2/\varepsilon\ge16`$ and $`d\ge4`$. Combining both directions gives

```math
\boxed{
n^\star_{d,s}(t,\varepsilon)
=
\Theta\!\left(
\frac{d(1-e^{-t/2})^2}{\varepsilon}
\right).
}
```

This is an all-time statement for the specified rank-one family. The unrestricted single-jump problem has the separate short-time matching theorem in [the model page](model.md). At $`t=0`$, the target is the identity and needs no program copies.

## 6. Resource meaning and verification

The construction is collective: it uses a projector on the system and the symmetric program space. It establishes copy sufficiency under arbitrary CPTP processing. It does not provide a polynomial gate decomposition or a small coherent-memory implementation.

The symmetric input memory has dimension $`\binom{n+m-1}{n}`$, and the output includes a degree-$`n-1`$ sector of dimension $`\binom{n+m-2}{n-1}`$. These are dimensions of an explicit mathematical implementation, not additional supplied copies.

The [finite verification](../verification/README.md) builds the annihilation operators, constructs this CPTP map directly, and compares its entire superoperator against the derived target-plus-error identity. It includes $`\alpha>1`$ cases so that the full admissible interval is tested.
