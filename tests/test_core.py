import numpy as np
import pytest
from nonlinear_oscillator_lab import (Model, simulate, stationary_points,
    turning_points, classify_energy, period, energy_drift)

def test_harmonic_period():
    m=Model(a=4)
    assert abs(period(m,0.4,0)-np.pi)<1e-7

def test_harmonic_solution_and_energy():
    m=Model(a=4)
    s=simulate(m,0.4,0,10)
    assert np.max(np.abs(s.y[0]-0.4*np.cos(2*s.t)))<2e-8
    assert energy_drift(m,s)<2e-9

def test_double_well_equilibria():
    x=stationary_points(Model(a=-1,b=1))
    assert np.allclose(x,[-1,0,1])

def test_turning_points():
    x=turning_points(Model(a=1),0.5)
    assert np.allclose(x,[-1,1])

def test_separatrix():
    m=Model(a=-1,b=1)
    assert classify_energy(m,1,2**-0.5)['on_separatrix_energy']
    with pytest.raises(ValueError): period(m,1,2**-0.5)

def test_cubic_quintic():
    m=Model(a=1,b=0.4,c=0.2)
    s=simulate(m,0.8,0.1,t_end=20)
    assert energy_drift(m,s)<5e-8
    assert period(m,0.8,0.1)>0

def test_equilibrium_period():
    assert np.isinf(period(Model(),0,0))


# --- closed-form Jacobi-elliptic references (independent of the quadrature) ---
from scipy.special import ellipk


def _duffing_exact_period(a, b, A):
    E = a * A * A / 2 + b * A ** 4 / 4
    if a < 0 and b > 0 and E < 0:          # intra-well orbit: x = A dn(w t | m)
        return 2 * ellipk(2 + 2 * a / (b * A * A)) / (np.sqrt(b / 2) * A)
    w2 = a + b * A * A                       # x = A cn(w t | m)
    return 4 * ellipk(b * A * A / (2 * w2)) / np.sqrt(w2)


@pytest.mark.parametrize("a,b,A", [(1, 1, 1.0), (1, -0.5, 0.9), (-1, 1, 1.2), (-1, 1, 1.6), (1, 0.4, 2.0)])
def test_period_matches_elliptic_closed_form(a, b, A):
    T = period(Model(a, b), A, 0.0)
    assert abs(T / _duffing_exact_period(a, b, A) - 1) < 1e-12


@pytest.mark.parametrize("dE", [1e-4, 1e-6])
def test_period_near_separatrix(dE):
    # double well: E = 0 - dE (intra-well) and E = 0 + dE (inter-well)
    m = Model(-1, 1)
    xin = np.sqrt(1 + np.sqrt(1 - 4 * dE))
    assert abs(period(m, xin, 0.0) / _duffing_exact_period(-1, 1, xin) - 1) < 1e-10
    xout = np.sqrt(1 + np.sqrt(1 + 4 * dE))
    assert abs(period(m, xout, 0.0) / _duffing_exact_period(-1, 1, xout) - 1) < 1e-10


def test_triple_well_equilibria():
    x = stationary_points(Model(1, -4, 3))
    assert np.allclose(x, [-1, -3 ** -0.5, 0, 3 ** -0.5, 1])
