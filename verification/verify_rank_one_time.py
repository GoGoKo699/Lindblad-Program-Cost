"""Finite checks for the long-time rank-one copy law."""

from itertools import product
import json

import numpy as np
from scipy.linalg import expm

from verify_core import projector, unit_vector


SEED = 73191


def vectorize(matrix):
    return matrix.reshape(-1, order="F")


def orthogonal_unit(rng, vector):
    direction = unit_vector(rng, vector.size)
    direction -= np.vdot(vector, direction) * vector
    return direction / np.linalg.norm(direction)


def decay_integral(time, q, rate):
    """Integrate exp(-q*(time-u)-rate*u) without subtracting close values."""
    gap = abs(q - rate)
    if gap == 0:
        return time * np.exp(-q * time)
    return np.exp(-min(q, rate) * time) * (-np.expm1(-gap * time)) / gap


def exact_channel(sink, bright, time):
    dimension = sink.size
    bright_projector = projector(bright)
    bright_part = bright_projector @ sink
    dark_part = sink - bright_part
    q = np.vdot(dark_part, dark_part).real
    cross = np.outer(dark_part, bright_part.conj())
    transfer = (
        decay_integral(time, q, 0) * projector(dark_part)
        + decay_integral(time, q, 0.5) * (cross + cross.conj().T)
        + decay_integral(time, q, 1) * projector(bright_part)
    )
    no_jump = np.eye(dimension) + np.expm1(-time / 2) * bright_projector
    return np.kron(no_jump.conj(), no_jump) + np.outer(
        vectorize(transfer), vectorize(bright_projector).conj()
    )


def generator(sink, bright):
    identity = np.eye(sink.size)
    jump = np.outer(sink, bright.conj())
    bright_projector = projector(bright)
    return np.kron(jump.conj(), jump) - (
        np.kron(identity, bright_projector) + np.kron(bright_projector.T, identity)
    ) / 2


def check_exact_channel():
    rng = np.random.default_rng(SEED)
    channel_tolerance = 1e-11
    trace_tolerance = 1e-12
    max_residual = 0.0
    max_trace_residual = 0.0
    cases = 0
    for dimension in (2, 3, 4, 6):
        identity_row = vectorize(np.eye(dimension))
        for q in (0, 1e-8, 0.01, 0.25, 0.5, 0.75, 1):
            sink = unit_vector(rng, dimension)
            bright = np.sqrt(1 - q) * sink + np.sqrt(q) * orthogonal_unit(rng, sink)
            lindbladian = generator(sink, bright)
            for time in (0, 1e-4, 0.1, 1, 4, 32, 512):
                actual = exact_channel(sink, bright, time)
                expected = expm(time * lindbladian)
                residual = float(np.max(np.abs(actual - expected)))
                trace_residual = float(np.max(np.abs(identity_row @ actual - identity_row)))
                if not (residual < channel_tolerance and trace_residual < trace_tolerance):
                    raise AssertionError(
                        f"rank-one channel d={dimension}, q={q}, t={time}: "
                        f"matrix residual {residual}, trace residual {trace_residual}"
                    )
                max_residual = max(residual, max_residual)
                max_trace_residual = max(trace_residual, max_trace_residual)
                cases += 1
    return {
        "cases": cases,
        "channel_tolerance": channel_tolerance,
        "trace_tolerance": trace_tolerance,
        "max_matrix_exponential_residual": max_residual,
        "max_trace_preservation_residual": max_trace_residual,
    }


def check_hard_family():
    times = np.geomspace(4, 1e6, 101)
    rate_floors = (0, 1e-4, 0.01, 0.25)
    tolerance = 1e-12
    min_signal = 1.0
    min_margin = 1.0
    sink = np.array([1.0, 0, 0, 0])
    initial = vectorize(projector(sink))
    for floor, time in product(rate_floors, times):
        q = max(floor, 1 / time)
        bright = np.array([np.sqrt(1 - q), np.sqrt(q), 0, 0])
        output = exact_channel(sink, bright, time) @ initial
        signal = float(output.reshape((4, 4), order="F")[1, 1].real)
        dark_population = 1 - (1 - q) * np.exp(-q * time)
        bound = (
            np.sqrt((1 - q) * dark_population)
            - np.sqrt(q * (1 - dark_population))
        ) ** 2
        margin = float(signal - bound)
        if not (signal > 1 / 8 and margin > -tolerance):
            raise AssertionError(
                f"hard-family signal {signal}, lower bound {bound}, q0={floor}, t={time}"
            )
        min_signal = min(signal, min_signal)
        min_margin = min(margin, min_margin)
    return {
        "time_points": len(times),
        "rate_floors": rate_floors,
        "minimum_signal": min_signal,
        "minimum_signal_minus_analytic_bound": min_margin,
    }


def normalized_choi_half_trace_norm(channel, dimension):
    tensor = channel.reshape((dimension,) * 4, order="F")
    choi = tensor.transpose(2, 0, 3, 1).reshape((dimension**2,) * 2) / dimension
    return float(np.linalg.svd(choi, compute_uv=False).sum() / 2)


def check_symmetric_channel_average():
    """Check Taylor cancellation on a maximally entangled input."""
    rng = np.random.default_rng(SEED)
    ratios = []
    max_scaled_error = 0.0
    for dimension, time in product((3, 4, 6), (1, 4, 32, 512)):
        for q in sorted({0, 1e-8, 1 / time, 0.25, 1}):
            sink = unit_vector(rng, dimension)
            bright = np.sqrt(1 - q) * sink + np.sqrt(q) * orthogonal_unit(rng, sink)
            sink_tangent = orthogonal_unit(rng, sink)
            bright_tangent = orthogonal_unit(rng, bright)
            target = exact_channel(sink, bright, time)
            errors = []
            for angle in (2**-9, 2**-10):
                cosine, sine = np.cos(angle), np.sin(angle)
                differences = [
                    exact_channel(
                        cosine * sink + sink_sign * sine * sink_tangent,
                        cosine * bright + bright_sign * sine * bright_tangent,
                        time,
                    ) - target
                    for sink_sign, bright_sign in product((-1, 1), repeat=2)
                ]
                largest_error = max(
                    normalized_choi_half_trace_norm(difference, dimension)
                    for difference in differences
                )
                average_error = normalized_choi_half_trace_norm(
                    np.mean(differences, axis=0), dimension
                )
                errors.append((largest_error, average_error))
                max_scaled_error = max(max_scaled_error, average_error / (time * angle**2))
            # Halving the displacement gives a factor two at first order and
            # four at second order. The finite windows allow higher-order terms.
            individual_ratio = errors[0][0] / errors[1][0]
            average_ratio = errors[0][1] / errors[1][1]
            if not (1.7 < individual_ratio < 2.5 and 3.8 < average_ratio < 4.2):
                raise AssertionError(
                    f"Taylor check d={dimension}, q={q}, t={time}: "
                    f"individual ratio {individual_ratio}, average ratio {average_ratio}"
                )
            ratios.append((individual_ratio, average_ratio))
    return {
        "cases": len(ratios),
        "individual_error_halving_ratio_range": [min(r[0] for r in ratios), max(r[0] for r in ratios)],
        "averaged_error_halving_ratio_range": [min(r[1] for r in ratios), max(r[1] for r in ratios)],
        "max_averaged_error_over_time_times_angle_squared": max_scaled_error,
    }


def run_checks():
    return {
        "seed": SEED,
        "exact_channel": check_exact_channel(),
        "hard_family": check_hard_family(),
        "symmetric_channel_average": check_symmetric_channel_average(),
        "status": "passed",
    }


if __name__ == "__main__":
    print(json.dumps(run_checks(), indent=2, sort_keys=True))
