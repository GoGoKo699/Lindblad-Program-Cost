# Universal upper bound from wave matrix Lindbladization

Wave matrix Lindbladization (WML) gives a universal processor using one supplied jump-program copy per step. This page reconstructs a conservative explicit bound from the invariant operator space used by Park et al. The upper bound is an inherited result; the constants below are chosen for a short, self-contained proof.

For $`d\ge 2`$, $`t\gt 0`$, and $`\varepsilon\gt 0`$, define

```math
N_d(t,\varepsilon)
=\max\!\left\{
 \lfloor 2dt\rfloor+1,\,
 \left\lceil\frac{(2d+1)t^2}{\varepsilon}\right\rceil
\right\}.
```

There is a fixed processor using exactly $`N_d(t,\varepsilon)`$ program copies such that

```math
\sup_{\|L\|_{\mathrm F}=1}
\frac12
\left\|
\mathcal W_{L,t/N_d}^{\,N_d}
-e^{t\mathcal D_L}
\right\|_\diamond
\le \varepsilon.
```

Here $`\mathcal D_L(X)=LXL^\dagger-\tfrac12\{L^\dagger L,X\}`$, and each $`\mathcal W_{L,\Delta}`$ is implemented by the same channel acting on the data and one supplied program. The channel's dependence on $`L`$ comes entirely from that program.

In the matching regime $`d\ge4`$, $`0\lt t\le1`$, and $`0\lt \varepsilon\le t^2/256`$, the ceiling term already exceeds $`2dt`$. Consequently,

```math
n_d^\star(t,\varepsilon)
\le
\left\lceil\frac{(2d+1)t^2}{\varepsilon}\right\rceil
\le \frac{3dt^2}{\varepsilon}.
```

The proof below is in the full diamond norm, with arbitrary reference systems. It does not assume a product input, knowledge of $`L`$, a state-preparation oracle, or access to an inverse program-preparation circuit.

## 1. One program copy implements one WML step

Let registers $`1,2,3`$ each have dimension $`d`$. Register $`1`$ is the data, and registers $`2,3`$ hold

```math
|L\rangle_{23}=(L\otimes I)|\Gamma\rangle_{23},
\qquad
|\Gamma\rangle=\sum_{j=1}^{d}|j,j\rangle,
\qquad
\pi_L=|L\rangle\langle L|.
```

The condition $`\|L\|_{\mathrm F}=1`$ makes $`\pi_L`$ a normalized state. Define the known operator

```math
M=\frac1{\sqrt d}
\left(I_1\otimes|\Gamma\rangle\langle\Gamma|_{23}\right)
\left(\mathrm{SWAP}_{12}\otimes I_3\right)
```

and its Lindblad generator

```math
\mathcal G(Y)=MYM^\dagger-\frac12\{M^\dagger M,Y\}.
```

For a duration $`\Delta\ge0`$, set

```math
\mathcal W_{L,\Delta}(X)
=\mathrm{Tr}_{23}
 \left[e^{\Delta\mathcal G}(X\otimes\pi_L)\right].
```

This is a CPTP map. Its implementation appends one program, applies the fixed channel $`e^{\Delta\mathcal G}`$, and discards registers $`2,3`$. The interaction depends on $`d`$ and $`\Delta`$, and is independent of the unknown jump operator. Composing $`n`$ steps consumes $`n`$ fresh copies.

These are the WML definitions in Park et al., Eqs. (10)-(15). The construction also appears in the earlier WML work cited there. The resource counted here is supplied program copies. A finite-dimensional Stinespring implementation of the known channel is permitted by the processor model; this argument does not assert a gate-complexity bound or free preparation of the unknown programs.

## 2. Algebra that controls the higher-order terms

Write $`A=M^\dagger M`$. Directly from the definition of $`M`$,

```math
M|a,b,c\rangle
=\frac{\delta_{a,c}}{\sqrt d}
 \sum_{j=1}^{d}|b,j,j\rangle.
```

It follows that

```math
M^2=\frac{M}{\sqrt d},
\qquad
MM^\dagger M=dM,
\qquad
A^2=dA.
```

For an arbitrary operator $`Y`$ on the three registers, introduce five linear maps:

```math
\begin{aligned}
\mathcal S_0(Y)&=Y,\\
\mathcal S_1(Y)&=MYM^\dagger,\\
\mathcal S_2(Y)&=\tfrac12(AY+YA),\\
\mathcal S_3(Y)&=\tfrac{\sqrt d}{2}
                  (AYM^\dagger+MYA),\\
\mathcal S_4(Y)&=\tfrac1d AYA.
\end{aligned}
```

The preceding operator identities give the closed composition rules

```math
\begin{aligned}
\mathcal G\mathcal S_0
 &=\mathcal S_1-\mathcal S_2,\\
\mathcal G\mathcal S_1
 &=d^{-1}\mathcal S_1-d^{-1}\mathcal S_3,\\
\mathcal G\mathcal S_2
 &=d\mathcal S_1-\tfrac d2\mathcal S_2-\tfrac d2\mathcal S_4,\\
\mathcal G\mathcal S_3
 &=d\mathcal S_1-\tfrac d2\mathcal S_3-\tfrac d2\mathcal S_4,\\
\mathcal G\mathcal S_4
 &=d\mathcal S_1-d\mathcal S_4.
\end{aligned}
```

For example, the jump part of $`\mathcal G\mathcal S_2`$ is $`d\mathcal S_1`$, using $`MA=dM`$ and $`AM^\dagger=dM^\dagger`$. Its anticommutator part is

```math
-\frac14(A^2Y+2AYA+YA^2)
=-\frac d2\mathcal S_2(Y)-\frac d2\mathcal S_4(Y).
```

The remaining rows follow by the same substitutions. This is the invariant span in Park et al., Appendix B, Eqs. (B1)-(B19). The present proof uses its composition rules directly, without needing the closed-form exponential coefficients or the norm estimate asserted in Lemma 6.

## 3. Reduction to data maps, including a reference

Let $`\mathcal A_{\pi_L}(X)=X\otimes\pi_L`$. Reducing each basis map gives

```math
\begin{aligned}
\mathcal R_1(X)
&:=\mathrm{Tr}_{23}\mathcal S_1
       \mathcal A_{\pi_L}(X)
 =LXL^\dagger,\\
\mathcal R_2(X)
&:=\mathrm{Tr}_{23}\mathcal S_2
       \mathcal A_{\pi_L}(X)
 =\tfrac12\{L^\dagger L,X\},\\
\mathcal R_3(X)
&:=\mathrm{Tr}_{23}\mathcal S_3
       \mathcal A_{\pi_L}(X)
 =LXL^\dagger,\\
\mathcal R_4(X)
&:=\mathrm{Tr}_{23}\mathcal S_4
       \mathcal A_{\pi_L}(X)
 =\frac{I}{d}\mathrm{Tr}(LXL^\dagger).
\end{aligned}
```

These identities can be checked on $`X=|a\rangle\langle a'|`$ using the displayed action of $`M`$. For example,

```math
M(X\otimes\pi_L)M^\dagger
=LXL^\dagger\otimes
  \frac{|\Gamma\rangle\langle\Gamma|}{d},
```

while

```math
\mathrm{Tr}_{23}
\!\left[A(X\otimes\pi_L)M^\dagger\right]
=\frac1{\sqrt d}LXL^\dagger.
```

The adjoint cross term has the same reduction, which verifies $`\mathcal R_3`$. Applying $`A`$ on both sides and tracing registers $`2,3`$ gives $`I\mathrm{Tr}(LXL^\dagger)`$, which verifies $`\mathcal R_4`$ after its factor $`1/d`$. The corresponding source identities are in Appendix B, Eqs. (B10)-(B12) and the text following Eq. (B28).

All these identities remain valid after tensoring with the identity on any reference $`R`$: expand an operator on $`R1`$ into reference matrix units and apply the identities to each data block. In the last line, the reference version replaces only the data system by $`I/d`$ and retains the reference's unnormalized marginal.

Set $`a=\|L\|_\infty^2\le1`$. Each reduced map has diamond norm at most $`a`$:

```math
\|\mathcal R_j\|_\diamond\le a,
\qquad j=1,2,3,4.
```

For $`\mathcal R_1=\mathcal R_3`$, the norm equals $`\|L^\dagger L\|_\infty=a`$. For $`\mathcal R_2`$, left and right multiplication each have completely bounded trace norm at most $`a`$. Finally, $`\mathcal R_4`$ is $`\mathcal R_1`$ followed by the completely depolarizing channel, whose diamond norm is one. These arguments bound arbitrary ancilla extensions, not just ordinary output trace distance.

In particular,

```math
\mathrm{Tr}_{23}\mathcal G
 \mathcal A_{\pi_L}
=\mathcal R_1-\mathcal R_2
=\mathcal D_L.
```

Thus WML has the required first-order generator.

## 4. A conservative WML remainder bound

For $`k\ge1`$, define

```math
\mathcal B_k
=\mathrm{Tr}_{23}\circ\mathcal G^k
 \circ\mathcal A_{\pi_L}.
```

After the first application of $`\mathcal G`$, the coefficient vector in the ordered basis $`(\mathcal S_1,\mathcal S_2,\mathcal S_3,\mathcal S_4)`$ is $`(1,-1,0,0)`$, with sum of absolute values two. Every subsequent row in the composition rules has sum of absolute coefficients at most $`2d`$. There is no return to $`\mathcal S_0`$. Induction therefore gives

```math
\|\mathcal B_k\|_\diamond
\le 2a(2d)^{k-1},
\qquad k\ge1.
```

This is the source of the linear dimension cost. The first generator application reduces to data maps of norm at most $`a`$; each additional application costs at most $`2d`$ in this coefficient bound.

Expanding the finite-dimensional exponential and using $`\mathcal B_1=\mathcal D_L`$,

```math
\begin{aligned}
\left\|
\mathcal W_{L,\Delta}-\mathcal I-\Delta\mathcal D_L
\right\|_\diamond
&\le
\sum_{k=2}^{\infty}\frac{\Delta^k}{k!}
       2a(2d)^{k-1}\\
&=\frac ad\left(e^{2d\Delta}-1-2d\Delta\right).
\end{aligned}
```

For $`0\le x\le1`$,

```math
e^x-1-x
=\sum_{k=2}^\infty\frac{x^k}{k!}
\le x^2\sum_{k=2}^\infty\frac1{k!}
=(e-2)x^2
\le x^2.
```

Consequently, whenever $`0\le\Delta\lt 1/(2d)`$,

```math
\left\|
\mathcal W_{L,\Delta}-\mathcal I-\Delta\mathcal D_L
\right\|_\diamond
\le 4ad\Delta^2.
```

No sign cancellation between different basis maps is used in this estimate.

## 5. Exact-channel remainder and one-step error

The target generator satisfies

```math
\|\mathcal D_L\|_\diamond
\le \|\mathcal R_1\|_\diamond+\|\mathcal R_2\|_\diamond
\le 2a.
```

Its exponential is a channel for every nonnegative time, so its diamond norm is one. The integral Taylor remainder gives

```math
e^{\Delta\mathcal D_L}-\mathcal I-\Delta\mathcal D_L
=\int_0^\Delta(\Delta-s)
 e^{s\mathcal D_L}\mathcal D_L^2\,ds.
```

By submultiplicativity,

```math
\left\|
e^{\Delta\mathcal D_L}-\mathcal I-\Delta\mathcal D_L
\right\|_\diamond
\le \frac{\Delta^2}{2}(2a)^2
=2a^2\Delta^2.
```

Combining the two remainders and restoring the normalization factor $`1/2`$ yields

```math
\begin{aligned}
\frac12\left\|
\mathcal W_{L,\Delta}-e^{\Delta\mathcal D_L}
\right\|_\diamond
&\le (2ad+a^2)\Delta^2\\
&\le (2d+1)a\Delta^2\\
&\le (2d+1)\Delta^2.
\end{aligned}
```

The first inequality retains the dependence on $`\|L\|_\infty`$. The final inequality supplies a single budget valid for all Frobenius-normalized jumps.

## 6. Composition and the integer sample budget

Let $`\Delta=t/n`$. The difference between two $`n`$-fold channel compositions telescopes into $`n`$ terms. Since every surrounding channel has diamond norm one,

```math
\frac12\left\|
\mathcal W_{L,t/n}^{\,n}-e^{t\mathcal D_L}
\right\|_\diamond
\le
n\,\frac12\left\|
\mathcal W_{L,t/n}-e^{(t/n)\mathcal D_L}
\right\|_\diamond.
```

Provided $`n\gt 2dt`$, the one-step condition holds and therefore

```math
\frac12\left\|
\mathcal W_{L,t/n}^{\,n}-e^{t\mathcal D_L}
\right\|_\diamond
\le \frac{(2d+1)a t^2}{n}
\le \frac{(2d+1)t^2}{n}.
```

The announced $`N_d(t,\varepsilon)`$ enforces both this strict step-size condition and the desired error. Park et al.'s Lemma 2 also imposes $`n\gt 2dt`$; the integer maximum makes that condition explicit here.

For $`d\ge4`$, $`0\lt t\le1`$, and $`\varepsilon\le t^2/256`$, let $`x=t^2/\varepsilon\ge256`$. Then $`(2d+1)x\gt 2dt`$, so

```math
N_d(t,\varepsilon)
=\lceil(2d+1)x\rceil
\le(2d+1)x+1
\le3dx.
```

Together with the converse in this repository, this establishes the matching worst-case order $`\Theta(dt^2/\varepsilon)`$ in the stated regime. No matching dependence on unbounded evolution time is asserted.

## Source and attribution

Siheon Park, Youngjin Seo, Byeongseon Go, Dhrumil Patel, Mark M. Wilde, and Hyukjoon Kwon, **Improved sample complexity bound for sample-based Lindbladian simulation**, [arXiv:2605.30301v1](https://arxiv.org/abs/2605.30301v1), 28 May 2026; [primary PDF](https://arxiv.org/pdf/2605.30301v1).

The relevant dependencies are the supplied-program model in Definition 1; the WML construction in Eqs. (10)-(15); the first-order identities in Eqs. (16)-(17); the step-size condition in Lemma 2; and the invariant operator space and reduction identities in Appendix B, Eqs. (B1)-(B19) and (B26)-(B28). Their Theorem 3 establishes the linear-dimension WML upper scaling. The proof on this page derives a sufficient larger coefficient directly from the displayed algebra and standard channel norm properties, and does not claim a new upper-scaling result.
