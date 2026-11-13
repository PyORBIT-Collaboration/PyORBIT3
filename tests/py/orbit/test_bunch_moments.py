"""Regression tests for normalized transverse bunch moments."""

import math
import random

import pytest

from orbit.core.bunch import Bunch, BunchTwissAnalysis


@pytest.mark.parametrize("weighted", [False, True])
@pytest.mark.parametrize("dispersion", [0, 1])
@pytest.mark.parametrize("emitnorm", [0, 1])
def test_bunch_moments(weighted, dispersion, emitnorm):
    rng = random.Random(1729)
    particles = [
        (rng.uniform(-0.01, 0.02), rng.uniform(-0.003, 0.002),
         rng.uniform(-0.02, 0.01), rng.uniform(-0.002, 0.003),
         rng.uniform(-0.1, 0.1), rng.uniform(-0.01, 0.02))
        for _ in range(64)
    ]
    weights = [1.0 + k % 5 for k in range(len(particles))] if weighted else [1.0] * len(particles)
    bunch = Bunch()
    bunch.mass(0.938272)
    bunch.getSyncParticle().kinEnergy(1.0)
    if weighted:
        bunch.addPartAttr("macrosize")
    for k, particle in enumerate(particles):
        bunch.addParticle(*particle)
        if weighted:
            bunch.partAttrValue("macrosize", k, 0, weights[k])

    total = sum(weights)
    means = [sum(w * p[k] for w, p in zip(weights, particles)) / total for k in range(6)]
    energy_variance = sum(w * (p[5] - means[5]) ** 2 for w, p in zip(weights, particles)) / total
    covariance = sum(w * (p[0] - means[0]) * (p[5] - means[5]) for w, p in zip(weights, particles)) / total
    coefficient = covariance / energy_variance if dispersion else 0.0
    x_mean = means[0] - coefficient * means[5]
    y_mean = means[2]

    analysis = BunchTwissAnalysis()
    analysis.computeBunchMoments(bunch, 4, dispersion, emitnorm)
    x_scale = math.sqrt(analysis.getBeta(0) * (analysis.getEmittance(0) if emitnorm else 1.0))
    y_scale = math.sqrt(analysis.getBeta(1) * (analysis.getEmittance(1) if emitnorm else 1.0))
    expected_moments = {}
    for j in range(5):
        for i in range(5 - j):
            if (i, j) == (1, 0):
                expected = x_mean
            elif (i, j) == (0, 1):
                expected = y_mean
            else:
                expected = sum(
                    w * ((p[0] - coefficient * p[5] - x_mean) / x_scale) ** i
                    * ((p[2] - y_mean) / y_scale) ** j
                    for w, p in zip(weights, particles)
                ) / total
            expected_moments[i, j] = expected
            assert analysis.getBunchMoment(i, j) == pytest.approx(expected, rel=1e-11, abs=1e-12)

    # Reuse the same analysis while shrinking and growing its contiguous storage.
    for order in (2, 1, 0, 4):
        analysis.computeBunchMoments(bunch, order, dispersion, emitnorm)
        assert analysis.getGlobalCount() == len(particles)
        assert analysis.getGlobalMacrosize() == pytest.approx(total)
        for j in range(order + 1):
            for i in range(order + 1):
                expected = expected_moments[i, j] if i + j <= order else 0.0
                assert analysis.getBunchMoment(i, j) == pytest.approx(expected, rel=1e-11, abs=1e-12)


@pytest.mark.parametrize("order", [0, 1, 2])
def test_empty_bunch_low_order(order):
    bunch = Bunch()
    bunch.getSyncParticle().kinEnergy(1.0)
    analysis = BunchTwissAnalysis()
    analysis.computeBunchMoments(bunch, order, 0, 0)
    assert analysis.getGlobalCount() == 0
    assert analysis.getGlobalMacrosize() == 0.0
    for j in range(order + 1):
        for i in range(order + 1):
            assert analysis.getBunchMoment(i, j) == (1.0 if i == j == 0 else 0.0)
