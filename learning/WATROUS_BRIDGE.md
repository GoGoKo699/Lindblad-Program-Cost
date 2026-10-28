# Watrous bridge

**Single anchor:** John Watrous, *The Theory of Quantum Information* (2018), https://cs.uwaterloo.ca/~watrous/TQI/TQI.pdf.

Read the operator and tensor-product foundations (Chapter 1), states, channels, Kraus and Choi representations (Chapter 2), and distinguishability and completely bounded trace norms (Chapter 3). These provide the foundations, not the Lindblad-specific research theorem.

## Translation into the project

The unnormalized maximally entangled vector `|Γ_d⟩=Σ_j|j,j⟩` encodes a Frobenius-normalized jump as `|L⟩=(L⊗I)|Γ_d⟩`, with `⟨L|L⟩=Tr(L†L)=1`. In particular, `L=|s⟩⟨φ|` gives the program `|s⟩|bar φ⟩`: the unknown decay direction is conjugated in the program.

A CPTP processor can jointly process the data and n copies of that state. The error is `½||E−F||◇`, including reference-entangled inputs. Watrous's channel formalism explains why this is stronger than checking selected ordinary states.

For `L=|s⟩⟨φ|` with orthogonal sink, the Kraus operators are `K₀=I−βP_φ` and `K₁=√(2β−β²)|s⟩⟨φ|`, where `β=1−e^(−t/2)`. Check `K₀†K₀+K₁†K₁=I`. This particular Lindblad exponential is derived in `science/dimension_cost.md`, not taken from the anchor.

**Worked example.** Take `s=|0⟩`, `r=|1⟩`, `ψ=|2⟩`, and `φ=(r+ψ)/√2`. On input r, the no-jump vector is `(1−β/2)r−(β/2)ψ`; projecting onto ψ gives probability `β²/4`. The actual hard family varies an arbitrary *complex* ψ in dimension `D=d−2`; the processor cannot use the known rotation of this toy example.

For a CP map T, Choi positivity and the symmetric Haar moment `∫P_ψ^⊗k dψ=P_sym,k/binom(D+k−1,k)` give the success-weighted conjugation inequality

`∫Tr[P_ψ T(bar P_ψ^⊗k)]dψ ≤ (k+1)/(k+D) ∫Tr[T(bar P_ψ^⊗k)]dψ`.

The project-specific reference-phase sector argument extends it to the flagged program. A spectral-diameter-one witness then yields the finite-copy lower bound `δ_d(n,t)≥β²(d−3)/[4(n+d−2)]`. See `science/conjugation.md` and `science/dimension_cost.md`.

The general matching upper bound uses wave matrix Lindbladization with the step condition `n>2dt`; see `science/upper_bound.md`. The known-sink rank-one collective processor has error `β²(d−2)/n` under its CPTP condition; see `science/rank_one_processor.md`. These are research constructions, not Watrous textbook results. The all-time family statement must not be transferred to arbitrary jumps.

The selected anchor is the sole educational source; primary research citations remain in `references/RELATED_WORK.md` for attribution.