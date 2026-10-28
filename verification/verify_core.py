"""Finite matrix checks for the converse and collective rank-one processor.

The tests check identities, tensor conventions, and specific parameter values.
The general minimax statements are established by the written proofs.
"""
from itertools import product
from math import comb, factorial
import json

import numpy as np
from scipy.linalg import expm

TOLERANCE = 1e-10


def projector(v):
    return np.outer(v, v.conj())


def unit_vector(rng, dimension):
    v = rng.normal(size=dimension) + 1j * rng.normal(size=dimension)
    return v / np.linalg.norm(v)


def residual(actual, expected):
    value = float(np.linalg.norm(actual - expected) / max(1.0, np.linalg.norm(expected)))
    if not np.isfinite(value) or value > TOLERANCE:
        raise AssertionError(f"matrix residual {value} exceeds {TOLERANCE}")
    return value


def occupation_tuples(dimension, degree):
    if dimension == 1:
        return [(degree,)]
    return [
        (first,) + tail
        for first in range(degree + 1)
        for tail in occupation_tuples(dimension - 1, degree - first)
    ]


def symmetric_isometry(dimension, degree):
    occupations = occupation_tuples(dimension, degree)
    columns = {occ: j for j, occ in enumerate(occupations)}
    v = np.zeros((dimension**degree, len(occupations)), dtype=complex)
    for flat, word in enumerate(product(range(dimension), repeat=degree)):
        occ = tuple(word.count(i) for i in range(dimension))
        multiplicity = factorial(degree)
        for count in occ:
            multiplicity //= factorial(count)
        v[flat, columns[occ]] = 1.0 / np.sqrt(multiplicity)
    return v


def annihilation_maps(dimension, degree):
    source = occupation_tuples(dimension, degree)
    target = occupation_tuples(dimension, degree - 1)
    target_index = {occ: j for j, occ in enumerate(target)}
    result = []
    for mode in range(dimension):
        a = np.zeros((len(target), len(source)), dtype=complex)
        for j, occ in enumerate(source):
            if occ[mode]:
                lower = list(occ)
                lower[mode] -= 1
                a[target_index[tuple(lower)], j] = np.sqrt(occ[mode])
        result.append(a)
    return source, target, result


def symmetric_power(v, degree):
    answer = np.array([1.0 + 0.0j])
    for _ in range(degree):
        answer = np.kron(answer, v)
    return answer


def symmetric_program(v, degree, occupations):
    coefficients = []
    for occ in occupations:
        multiplicity = factorial(degree)
        for count in occ:
            multiplicity //= factorial(count)
        amplitude = np.sqrt(multiplicity)
        for component, count in zip(v, occ):
            amplitude *= component**count
        coefficients.append(amplitude)
    return np.array(coefficients, dtype=complex)


def trace_memory(x, system_dimension, memory_dimension):
    return np.trace(
        x.reshape(system_dimension, memory_dimension,
                  system_dimension, memory_dimension),
        axis1=1,
        axis2=3,
    )


def hard_family_checks(rng):
    worst_program = worst_channel = worst_branch = 0.0
    cases = 0
    for d in (4, 5, 7):
        s = np.eye(d, dtype=complex)[:, 0]
        r = np.eye(d, dtype=complex)[:, 1]
        psi = np.zeros(d, dtype=complex)
        psi[2:] = unit_vector(rng, d - 2)
        phi = (r + psi) / np.sqrt(2)
        jump = np.outer(s, phi.conj())
        gamma = np.eye(d).reshape(-1)
        encoded = np.kron(jump, np.eye(d)) @ gamma
        worst_program = max(
            worst_program, residual(encoded, np.kron(s, phi.conj()))
        )
        a = jump.conj().T @ jump
        generator = (
            np.kron(jump.conj(), jump)
            - 0.5 * (np.kron(np.eye(d), a) + np.kron(a.T, np.eye(d)))
        )
        q = np.diag([0, 0] + [1] * (d - 2))
        for t in (0.0, 0.01, 0.2, 1.0, 5.0):
            beta = -np.expm1(-t / 2)
            k0 = np.eye(d) - beta * projector(phi)
            k1 = np.sqrt(-np.expm1(-t)) * jump
            kraus_super = np.kron(k0.conj(), k0) + np.kron(k1.conj(), k1)
            exact_super = expm(t * generator)
            worst_channel = max(worst_channel, residual(exact_super, kraus_super))
            output = (exact_super @ projector(r).reshape(-1, order="F")).reshape(
                (d, d), order="F"
            )
            worst_branch = max(
                worst_branch,
                residual(q @ output @ q, beta**2 * projector(psi) / 4),
            )
            cases += 1
    return {
        "cases": cases,
        "max_vectorization_residual": worst_program,
        "max_kraus_vs_exponential_residual": worst_channel,
        "max_success_branch_residual": worst_branch,
    }


def conjugation_operator_checks():
    minimum_eigenvalue = 0.0
    worst_projector = 0.0
    cases = 0
    for dimension in (2, 3):
        for degree in range(4):
            v = symmetric_isometry(dimension, degree)
            w = symmetric_isometry(dimension, degree + 1)
            p = v @ v.conj().T
            pnext = w @ w.conj().T
            h = comb(dimension + degree - 1, degree)
            hnext = comb(dimension + degree, degree + 1)
            worst_projector = max(
                worst_projector,
                residual(p @ p, p),
                residual(pnext @ pnext, pnext),
                residual(np.array([np.trace(p)]), np.array([h])),
            )
            average_overlap = pnext / hnext
            average_success = np.kron(p, np.eye(dimension)) / h
            coefficient = (degree + 1) / (degree + dimension)
            eigenvalues = np.linalg.eigvalsh(
                coefficient * average_success - average_overlap
            )
            minimum_eigenvalue = min(minimum_eigenvalue, float(eigenvalues.min()))
            if eigenvalues.min() < -TOLERANCE:
                raise AssertionError("Haar/Choi operator inequality violated")
            cases += 1
    return {
        "cases": cases,
        "max_symmetric_projector_residual": worst_projector,
        "minimum_choi_slack_eigenvalue": minimum_eigenvalue,
    }


def phase_sector_checks(rng):
    worst = 0.0
    for n in range(5):
        psi = np.r_[0.0, unit_vector(rng, 2)]
        reference = np.array([1.0, 0.0, 0.0], dtype=complex)
        average = np.zeros((3**n, 3**n), dtype=complex)
        for j in range(n + 1):
            theta = 2 * np.pi * j / (n + 1)
            flag = (reference + np.exp(-1j * theta) * psi.conj()) / np.sqrt(2)
            average += projector(symmetric_power(flag, n)) / (n + 1)
        sectors = np.zeros_like(average)
        for k in range(n + 1):
            sector = np.zeros(3**n, dtype=complex)
            for word in product((0, 1), repeat=n):
                if sum(word) != k:
                    continue
                term = np.array([1.0 + 0.0j])
                for selected in word:
                    term = np.kron(term, psi.conj() if selected else reference)
                sector += term
            sector /= np.sqrt(comb(n, k))
            sectors += comb(n, k) / 2**n * projector(sector)
        worst = max(worst, residual(average, sectors))
    return {"cases": 5, "max_phase_dephasing_residual": worst}


def collective_processor_checks(rng):
    worst_projection = worst_tp = worst_channel = worst_attained = 0.0
    cases = 0
    negative_survival_cases = 0
    for m in (2, 3):
        d = m + 1
        sink = np.eye(d, dtype=complex)[:, 0]
        for n in (1, 2, 3, 4):
            source, target, maps = annihilation_maps(m, n)
            h, hprev = len(source), len(target)
            c = np.zeros((hprev, d * h), dtype=complex)
            for i, a in enumerate(maps):
                c[:, (i + 1) * h:(i + 2) * h] = a
            normalization = n + m - 1
            pi = c.conj().T @ c / normalization
            jmap = np.kron(sink[:, None], np.eye(hprev)) @ c / np.sqrt(normalization)
            worst_projection = max(
                worst_projection,
                residual(c @ c.conj().T, normalization * np.eye(hprev)),
                residual(pi @ pi, pi),
                residual(jmap.conj().T @ jmap, pi),
            )
            phi_active = unit_vector(rng, m)
            phi = np.r_[0.0, phi_active]
            p = projector(phi)
            ps = projector(sink)
            ia = np.diag([0] + [1] * m)
            program = symmetric_program(phi_active.conj(), n, source)
            sigma = projector(program)
            residual(np.array([np.vdot(program, program)]), np.array([1.0]))
            q = n / normalization
            for beta in (0.0, 0.2, 0.6, 0.9, 1.0):
                alpha = beta / q
                if alpha > 2 + 1e-14:
                    continue
                a0 = np.eye(d * h) - alpha * pi
                a1 = np.sqrt(max(0.0, alpha * (2 - alpha))) * jmap
                worst_tp = max(
                    worst_tp,
                    residual(a0.conj().T @ a0 + a1.conj().T @ a1, np.eye(d * h)),
                )

                def apply(x):
                    joint = np.kron(x, sigma)
                    return (
                        trace_memory(a0 @ joint @ a0.conj().T, d, h)
                        + trace_memory(a1 @ joint @ a1.conj().T, d, hprev)
                    )

                # Matrix-unit comparison checks the entire linear map.
                for row in range(d):
                    for col in range(d):
                        x = np.zeros((d, d), dtype=complex)
                        x[row, col] = 1.0
                        u = np.trace(p @ x)
                        target_map = (
                            x - beta * (p @ x + x @ p)
                            + beta**2 * p @ x @ p
                            + (2 * beta - beta**2) * u * ps
                        )
                        error_map = beta**2 / n * u * (ia - p - (m - 1) * ps)
                        worst_channel = max(
                            worst_channel, residual(apply(x), target_map + error_map)
                        )
                delta = beta**2 * (m - 1) / n
                target_p = (1 - beta)**2 * p + (2 * beta - beta**2) * ps
                actual_error = (apply(p) - target_p)
                trace_error = np.linalg.svd(actual_error, compute_uv=False).sum() / 2
                worst_attained = max(
                    worst_attained,
                    residual(np.array([trace_error]), np.array([delta])),
                )
                cases += 1
                negative_survival_cases += int(alpha > 1 + 1e-14)
    return {
        "cases": cases,
        "alpha_greater_than_one_cases": negative_survival_cases,
        "max_contraction_projector_residual": worst_projection,
        "max_trace_preservation_residual": worst_tp,
        "max_full_channel_identity_residual": worst_channel,
        "max_attained_error_residual": worst_attained,
    }


def run_checks():
    rng = np.random.default_rng(20261028)
    return {
        "seed": 20261028,
        "tolerance": TOLERANCE,
        "hard_family": hard_family_checks(rng),
        "conjugation_operator": conjugation_operator_checks(),
        "reference_phase_sectors": phase_sector_checks(rng),
        "collective_rank_one_processor": collective_processor_checks(rng),
        "status": "passed",
    }


if __name__ == "__main__":
    print(json.dumps(run_checks(), indent=2, sort_keys=True))
