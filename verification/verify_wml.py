"""Finite algebra checks for science/upper_bound.md.

These deterministic checks test the invariant-map identities and their reduction
with an arbitrary reference. They do not prove the uniform diamond-norm bound.
Only NumPy is needed.
"""

from __future__ import annotations

import json

import numpy as np


SEED = 207
TOLERANCE = 1.0e-12
SYSTEM_DIMENSIONS = (2, 3)
REFERENCE_DIMENSION = 2


def _wml_operator(dimension: int) -> np.ndarray:
    """Build M from its action on the ordered registers (data, program, program)."""
    size = dimension**3
    operator = np.zeros((size, size), dtype=np.complex128)

    def index(a: int, b: int, c: int) -> int:
        return (a * dimension + b) * dimension + c

    for a in range(dimension):
        for b in range(dimension):
            for j in range(dimension):
                operator[index(b, j, j), index(a, b, a)] += (
                    1.0 / np.sqrt(dimension)
                )
    return operator


def _invariant_maps(
    argument: np.ndarray,
    operator: np.ndarray,
    effect: np.ndarray,
    dimension: int,
) -> list[np.ndarray]:
    adjoint = operator.conj().T
    return [
        argument,
        operator @ argument @ adjoint,
        (effect @ argument + argument @ effect) / 2,
        np.sqrt(dimension)
        * (effect @ argument @ adjoint + operator @ argument @ effect)
        / 2,
        effect @ argument @ effect / dimension,
    ]


def _generator(
    argument: np.ndarray, operator: np.ndarray, effect: np.ndarray
) -> np.ndarray:
    return (
        operator @ argument @ operator.conj().T
        - (effect @ argument + argument @ effect) / 2
    )


def _trace_program(
    argument: np.ndarray, dimension: int, reference_dimension: int
) -> np.ndarray:
    """Trace the final two registers in the order (reference, data, 2, 3)."""
    tensor = argument.reshape(
        reference_dimension,
        dimension,
        dimension,
        dimension,
        reference_dimension,
        dimension,
        dimension,
        dimension,
    )
    return np.einsum("iabcjxbc->iajx", tensor).reshape(
        reference_dimension * dimension,
        reference_dimension * dimension,
    )


def run_checks() -> dict[str, object]:
    """Run the finite identities and return a JSON-serializable summary.

    Failure or a nonfinite residual raises AssertionError. Complex, generally
    non-Hermitian inputs are intentional: the tested equalities are identities
    of linear maps, including reference matrix blocks.
    """
    rng = np.random.default_rng(SEED)
    records: list[dict[str, object]] = []

    def check(
        name: str, dimension: int, actual: np.ndarray, expected: np.ndarray
    ) -> None:
        residual = float(
            np.linalg.norm(actual - expected)
            / max(1.0, float(np.linalg.norm(expected)))
        )
        if not np.isfinite(residual) or residual > TOLERANCE:
            raise AssertionError(
                f"{name}, d={dimension}: relative residual {residual:.6g} "
                f"exceeds {TOLERANCE:.6g}"
            )
        records.append(
            {
                "name": name,
                "dimension": dimension,
                "relative_residual": residual,
            }
        )

    for dimension in SYSTEM_DIMENSIONS:
        size = dimension**3
        operator = _wml_operator(dimension)
        effect = operator.conj().T @ operator

        # The five rows of G S_j, expressed in S_0,...,S_4.
        composition = np.array(
            [
                [0, 1, -1, 0, 0],
                [0, 1 / dimension, 0, -1 / dimension, 0],
                [0, dimension, -dimension / 2, 0, -dimension / 2],
                [0, dimension, 0, -dimension / 2, -dimension / 2],
                [0, dimension, 0, 0, -dimension],
            ],
            dtype=float,
        )
        argument = rng.normal(size=(size, size)) + 1j * rng.normal(
            size=(size, size)
        )
        basis = _invariant_maps(argument, operator, effect, dimension)
        for row in range(5):
            expected = sum(
                composition[row, column] * basis[column]
                for column in range(5)
            )
            check(
                f"invariant_map_row_{row}",
                dimension,
                _generator(basis[row], operator, effect),
                expected,
            )

        # Row-major vectorization is exactly (L tensor I)|Gamma>.
        jump = rng.normal(size=(dimension, dimension)) + 1j * rng.normal(
            size=(dimension, dimension)
        )
        jump /= np.linalg.norm(jump)
        program_vector = jump.reshape(-1, order="C")
        program = np.outer(program_vector, program_vector.conj())

        reference_dimension = REFERENCE_DIMENSION
        input_size = reference_dimension * dimension
        data_reference = rng.normal(
            size=(input_size, input_size)
        ) + 1j * rng.normal(size=(input_size, input_size))
        full_input = np.kron(data_reference, program)
        reference_identity = np.eye(reference_dimension)
        extended_operator = np.kron(reference_identity, operator)
        extended_effect = np.kron(reference_identity, effect)
        full_maps = _invariant_maps(
            full_input, extended_operator, extended_effect, dimension
        )[1:]

        extended_jump = np.kron(reference_identity, jump)
        extended_jump_effect = np.kron(
            reference_identity, jump.conj().T @ jump
        )
        jump_output = (
            extended_jump @ data_reference @ extended_jump.conj().T
        )
        reference_marginal = np.einsum(
            "iaja->ij",
            jump_output.reshape(
                reference_dimension,
                dimension,
                reference_dimension,
                dimension,
            ),
        )
        expected_maps = [
            jump_output,
            (
                extended_jump_effect @ data_reference
                + data_reference @ extended_jump_effect
            )
            / 2,
            jump_output,
            np.kron(reference_marginal, np.eye(dimension) / dimension),
        ]
        for number, (full_map, expected) in enumerate(
            zip(full_maps, expected_maps), start=1
        ):
            check(
                f"reference_reduction_R{number}",
                dimension,
                _trace_program(
                    full_map, dimension, reference_dimension
                ),
                expected,
            )

    return {
        "passed": True,
        "check_count": len(records),
        "seed": SEED,
        "tolerance": TOLERANCE,
        "system_dimensions": list(SYSTEM_DIMENSIONS),
        "reference_dimension": REFERENCE_DIMENSION,
        "maximum_relative_residual": max(
            float(record["relative_residual"]) for record in records
        ),
        "checks": records,
        "scope": (
            "Finite invariant-map and data-reference reduction identities; "
            "the analytic proof establishes the uniform diamond-norm bound."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(run_checks(), indent=2, sort_keys=True))
