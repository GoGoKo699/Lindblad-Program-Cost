# The copy cost of rank-one decay

A rank-one jump $`L=|s\rangle\langle\phi|`$ can return population to its own decay direction when $`s`$ and $`\phi`$ overlap. Each jump then removes only part of the remaining population from that direction. When the removal rate is of order $`1/t`$, a time-$`t`$ evolution can move constant population into a direction that was only weakly represented in the supplied programs.

Use the [supplied-copy model](model.md) with $`d\ge4`$, unit vectors $`s,\phi`$, and target $`\mathcal E_t=e^{t\mathcal D_L}`$ generated only by the indicated jump. Each program copy is the product $`|s\rangle|\bar\phi\rangle`$. The error is half the diamond norm for the output channel on arbitrary data and its reference. The processor may know $`d,t`$ and the specified accuracy, and its deterministic copy budget covers all branches.

Write $`q=1-|\langle\phi|s\rangle|^2`$. For all $`t\gt0`$ and $`0\lt\varepsilon\le\min(t^2,1)/1024`$, the minimax copy count over all rank-one jumps is

```math
n^\star_{\mathrm{rank}\,1}(d,t,\varepsilon)
=\Theta\!\left(\frac{d}{\varepsilon}\min\{t^2,t\}\right),
```

with universal constants. The copy cost crosses from quadratic to linear time dependence. The [orthogonal known-sink construction](rank_one_processor.md) gives its own sharper finite-copy guarantee for that special case. At $`t=0`$ the channel is the identity and uses no programs.

## Exact channel

Let $`P=|\phi\rangle\langle\phi|`$, $`P_s=|s\rangle\langle s|`$, and $`K_u=I-(1-e^{-u/2})P`$. For any initial operator $`X`$, write $`X_u=\mathcal E_u(X)`$. The generator gives

```math
\frac{dX_u}{du}
= -\frac{1}{2}(PX_u+X_uP)+\mathrm{Tr}(PX_u)P_s,
\qquad
\frac{d}{du}\mathrm{Tr}(PX_u)=-q\,\mathrm{Tr}(PX_u).
```

Solving the second equation and varying constants in the first gives

```math
\mathcal E_t(X)
=K_tXK_t+\mathrm{Tr}(PX)R_t,
\qquad
R_t=\int_0^t e^{-q(t-u)}g_ug_u^\dagger\,du,
\qquad
g_u=K_us.
```

The operator $`R_t`$ is positive and has trace $`1-e^{-t}`$, since $`K_t^2=I-(1-e^{-t})P`$ and $`\mathcal E_t`$ preserves trace. For $`q=1`$ the transfer is $`R_t=(1-e^{-t})P_s`$. For $`q=0`$ it gives the population term of pure dephasing along $`P`$. The formula is continuous in between, including at vanishing $`q`$.

## Processor-independent lower bound

For $`t\ge4`$, choose $`q=1/t`$. Then $`q\le1/4`$ and $`qt\ge1`$. Choose a known basis vector $`r`$ and write $`\mathcal H=\mathbb C r\oplus\mathcal U`$ with $`D=\dim\mathcal U=d-1`$. Take a Haar unit vector $`\psi\in\mathcal U`$ including its relative phase and set

```math
s=r,\qquad
\phi_\psi=\sqrt{1-q}\,r+\sqrt q\,\psi.
```

The input is the fixed state $`|r\rangle\langle r|`$. The output remains in the span of $`r,\psi`$. Its population in the direction orthogonal to $`\phi_\psi`$ within that span is $`f=1-(1-q)e^{-qt}`$. By positivity, the magnitude of the bright-dark coherence is at most $`\sqrt{f(1-f)}`$. The $`\psi`$ output population therefore obeys

```math
p\ge
\left(\sqrt{(1-q)f}-\sqrt{q(1-f)}\right)^2
\ge\frac{(\sqrt{15}-\sqrt3)^2}{32}\gt\frac18.
```

Here $`f\gt5/8`$ follows from $`qt\ge1`$. If $`Q`$ is the projection onto $`\mathcal U`$ and $`P_\psi=|\psi\rangle\langle\psi|`$, the target output has $`Q`$ mass $`p`$ and zero $`Q-P_\psi`$ mass.

### How many informative factors reach a successful output?

Discard the fixed first program factors. Expand the remaining $`n`$ copies of $`\sqrt{1-q}r+\sqrt q\,\bar\psi`$ into sectors with $`k`$ factors in $`\overline{\mathcal U}`$. Their squared weights are $`p_k=\binom nk q^k(1-q)^{n-k}`$. Each normalized sector is the image of $`|\bar\psi\rangle^{\otimes k}`$ under the known isometric embedding that appends $`n-k`$ copies of $`r`$ and symmetrizes with the sector normalization. Rotating the phase of $`\psi`$ multiplies sector $`k`$ by $`e^{-ik\vartheta}`$. Both output effects are phase invariant, so their joint Haar averages are unchanged if all cross-sector terms are deleted.

The processor with its data input fixed and the output projected onto $`Q`$ induces a CP trace-nonincreasing map $`\mathcal T_k`$ on each sector. Let $`s_k`$ denote its Haar-average success weight and $`h_k`$ its Haar-average overlap weight with $`P_\psi`$. The [success-weighted conjugation inequality](conjugation.md) gives

```math
h_k\le \frac{k+1}{k+D}s_k.
```

Indeed, with $`J_k\ge0`$ the Choi operator of $`\mathcal T_k`$, the integrals defining $`h_k`$ and $`s_k`$ become the traces against $`\Pi_{k+1}/\dim\mathrm{Sym}^{k+1}(\mathcal U)`$ and $`(I\otimes\Pi_k)/\dim\mathrm{Sym}^k(\mathcal U)`$. The transpose in the Choi identity turns the input conjugated projectors into unconjugated ones. The inclusion $`\Pi_{k+1}\le I\otimes\Pi_k`$ and the ratio of the two symmetric dimensions prove the inequality, also for $`k=0`$.

Write $`S=\sum_k p_k s_k`$ for the averaged $`Q`$ mass. A worst-case half-diamond error of at most $`\varepsilon`$ requires $`S\ge1/8-\varepsilon`$ and bounds the erroneous $`Q-P_\psi`$ output mass by $`\varepsilon`$. Cauchy–Schwarz retains the actual sector weights, including state-dependent success and possible concentration on rare sectors,

```math
\begin{aligned}
\varepsilon
&\ge\sum_k p_k(s_k-h_k)
\ge(D-1)\sum_k\frac{p_ks_k}{k+D}\\
&\ge\frac{(D-1)S^2}{\sum_k p_ks_k(k+D)}
\ge\frac{(D-1)S^2}{nq+D}.
\end{aligned}
```

The last step uses $`s_k\le1`$ and $`\sum_k k p_k=nq`$. If $`\varepsilon\le1/1024`$, then $`S\ge1/16`$ and $`D\le(D-1)/(512\varepsilon)`$. Thus

```math
n\ge\frac{d-2}{512q\varepsilon}
=\frac{(d-2)t}{512\varepsilon}.
```

The lower bound allows the processor to know both $`r`$ and $`q`$. It only hides the direction and relative phase of $`\psi`$. At $`0\lt t\le1`$, the [existing dimension lower bound](dimension_cost.md) applies to its orthogonal-sink subclass. At $`1\le t\le4`$, the same subclass gives the lower order because its $`1-e^{-t/2}`$ factor is bounded below by $`1/4`$.

## Processor for arbitrary rank-one jumps

At $`t\le1`$ apply the [existing WML upper bound](upper_bound.md). For $`t\ge1`$, measure the two product factors from each program separately across copies. For $`n`$ copies of a pure direction $`v`$, use the covariant state-estimation POVM on the symmetric subspace. Its effect density on rays $`x`$ is $`m_{n,d}|x\rangle\langle x|^{\otimes n}`$ with respect to Haar measure, where $`m_{n,d}=\dim\mathrm{Sym}^{n}(\mathbb C^d)`$. On the orthogonal subspace use an arbitrary fixed outcome. Conjugate the classical ray estimate for $`\bar\phi`$ to obtain the estimate of $`\phi`$. Run the exact finite-dimensional channel for the estimated $`s,\phi`$ using the channel formula above, and discard the two internal outcomes.

### Why averaging cancels the first-order error

Relative to the true direction choose the estimate representative $`\widehat v=\cos\theta\,v+\sin\theta\,z`$ with $`z\perp v`$ unit and $`0\le\theta\le\pi/2`$. Given $`\theta`$, the direction $`z`$ is uniform on that sphere, including phase. The POVM and the $`(n+1)`$-copy Haar integral give

```math
\mathbb E[\sin^2\theta]=1-\frac{m_{n,d}}{m_{n+1,d}}
=\frac{d-1}{n+d},
\qquad
\mathbb E[z\mid\theta]=0,
\qquad
\theta^2\le\frac{\pi^2}{4}\sin^2\theta.
```

Apply this to both estimates and interpolate each along $`v_\lambda=\cos(\lambda\theta)v+\sin(\lambda\theta)z`$ with $`\lambda\in[0,1]`$. Let $`a=\theta_s+\theta_\phi`$. Then $`\mathbb E[a^2]\le\pi^2(d-1)/(n+d)`$. Differentiation of the channel at $`\lambda=0`$ is linear in the two initial tangents, so its mean vanishes, also conditional on the two angles. The quantity to bound is the diamond distance of the averaged channel, where the measurement is internal, and the entire second-order remainder.

### Uniform second-order bound

Primes denote differentiation in $`\lambda`$. The projector and the two vectors in the exact channel satisfy

```math
\begin{gathered}
\|P'\|_\infty,\ \|g_u'\|\le2a,
\qquad
\|P''\|_\infty,\ \|g_u''\|\le4a^2,\\
|q'|\le4a\sqrt q,
\qquad
|q''|\le16a^2,
\qquad
\|g_u\|^2=q+(1-q)e^{-u}.
\end{gathered}
```

The $`q`$ bounds follow by differentiating $`q=\|(I-P)s\|^2`$, since the first and second derivatives of $`(I-P)s`$ have norms at most $`2a`$ and $`4a^2`$. With $`w=t-u`$, the trace norm of the second derivative of the $`R_t`$ integrand is at most

```math
a^2e^{-qw}
\left[(16w^2q+16w)(q+e^{-u})
+32w\sqrt q(\sqrt q+e^{-u/2})+16\right].
```

The terms containing $`(wq)^2`$ and $`wq`$ times $`e^{-qw}`$ are bounded pointwise by constants. For the terms containing $`e^{-u}`$, use $`w^2q e^{-qw}\le w/e`$ and $`w\le t`$. The cross term containing $`e^{-u/2}`$ also integrates to order $`t`$ since $`w\sqrt q\le t`$. This gives $`\|R_t''\|_1\le144ta^2`$. Applying the same inequalities to the first derivative gives $`\|R_t'\|_1\le16ta`$.

The map $`X\mapsto A X B`$ has diamond norm at most $`\|A\|_\infty\|B\|_\infty`$. The map $`X\mapsto\mathrm{Tr}(AX)B`$ has diamond norm at most $`\|A\|_\infty\|B\|_1`$. Since $`\|K_t\|_\infty\le1`$ and $`\|R_t\|_1\le1`$, the total diamond norm of the channel's second derivative is at most $`256 t a^2`$ for $`t\ge1`$. Taylor's integral remainder and the cancellation of its linear term yield

```math
\frac12\left\|\mathbb E[\mathcal E_{t,\widehat s,\widehat\phi}]
-\mathcal E_{t,s,\phi}\right\|_\diamond
\le \frac{64\pi^2 t(d-1)}{n+d}.
```

The bound is on the complete output channel with internal estimates discarded. A finite partition of the compact ray space with one representative per outcome approximates the continuous outcome rule to arbitrary additional diamond precision with the same number of supplied copies. Gate count, measurement implementation, and classical synthesis precision are separate from copy count.

## A positive removal rate floor

Let $`n^\star_{\mathrm{rank}\,1}(d,t,\varepsilon;q_0)`$ be the minimax count under the promise $`q\ge q_0`$. With the same time and error ranges and $`0\le q_0\le1/4`$, the proof extends to

```math
n^\star_{\mathrm{rank}\,1}(d,t,\varepsilon;q_0)
=\Theta\!\left(\frac{d}{\varepsilon}\min\{t^2,t,q_0^{-1}\}\right),
```

with $`q_0^{-1}=+\infty`$ when $`q_0=0`$. To obtain the lower bound for $`t\ge4`$, choose $`q=\max(q_0,1/t)`$ in the same hard family. Then $`q\le1/4`$ and $`qt\ge1`$, so the output signal remains above $`1/8`$. The sector argument gives $`n\ge(d-2)/(512q\varepsilon)`$ as above. The orthogonal-sink hard set used for $`t\le4`$ has $`q=1`$ and satisfies the rate floor.
For a supplied $`q_0\gt0`$, consider the event $`a\le\sqrt{q_0}/4`$. The derivative bound for $`(I-P)s`$ implies $`q_\lambda\ge q_0/4`$ along the complete interpolation on this event. There the same integrand estimates give $`\|R_t''\|_1=O(a^2/q_0)`$ and $`\|R_t'\|_1=O(a/q_0)`$. For terms containing only $`wq`$ and $`e^{-qw}`$, substitute $`x=wq`$ in the integral. For terms also containing $`e^{-u}`$ or $`e^{-u/2}`$, use $`w e^{-qw}=O(1/q)`$, $`w^2q e^{-qw}=O(1/q)`$, and $`w\sqrt q\,e^{-qw}=O(1/\sqrt q)`$ before integrating the remaining exponential.

The event depends only on the two angles, so the first-order term still averages to zero there. The complementary event has probability at most $`16\mathbb E[a^2]/q_0`$ by Markov's inequality and half-diamond error at most one. The same processor therefore also gives error $`O(d/[(n+d)q_0])`$ for $`t\ge1`$. Taking the two valid bounds together gives the claimed $`\min\{t,q_0^{-1}\}`$ long-time factor.

## What changes with time

The long-time hard family has a fixed, known output state $`r`$ and a removal rate of order $`1/t`$ in the absence of a floor. The effective number of informative conjugated vectors in the programs is of order $`n/t`$, even though the target transfers constant population to the hidden direction. A uniform positive removal rate makes that rare-program mechanism saturate. The exactly parallel case has dephasing dynamics and is not the hard family.

The upper construction uses an inherited covariant state estimator and an exact rank-one channel formula. The contribution is the uniform second-order control of the averaged channel and the matching conjugation converse that counts informative factors. The result concerns rank-one jumps. The general-jump time lower bound of [Go et al.](https://arxiv.org/abs/2412.02134v3) uses rank-two jumps. The rank-one lower bound of [Park et al.](https://arxiv.org/abs/2605.30301v1) analyzes WML, which can overspend copies on rank-one inputs. The covariant estimation ingredient is treated in [Hayashi, Hashimoto and Horibe](https://arxiv.org/abs/quant-ph/0410207v3).
