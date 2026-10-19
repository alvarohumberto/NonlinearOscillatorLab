"""Period-energy grids computed with the Python package.

Output: benchmarks/python_period_energy.csv (read by mathematica/MakeFigures.wls,
which recomputes every point in the Wolfram Language and against closed-form or
40-digit references).
"""
from pathlib import Path
import csv
import numpy as np
from nonlinear_oscillator_lab import Model, period

def approach(es, n=40, lo=0.02, hi=6.0):
    """Energies approaching es from below, geometrically closer to the barrier."""
    return es * (1.0 - 10.0 ** -np.linspace(lo, hi, n))

FAMILIES = {
    # name: (a, b, c, x0, energies)
    "hardening": (1.0, 1.0, 0.0, 0.0, np.linspace(0.05, 5.0, 40)),
    "softening": (1.0, -0.5, 0.0, 0.0, approach(0.5)),
    "dw_inner": (-1.0, 1.0, 0.0, 1.0, -0.25 * 10.0 ** -np.linspace(0.005, 6.0, 40)),
    "dw_outer": (-1.0, 1.0, 0.0, 0.0, 10.0 ** np.linspace(-6.0, 0.7, 40)),
    "cq_demo": (1.0, 0.4, 0.2, 0.0, np.linspace(0.05, 5.0, 40)),
    "tw_center": (1.0, -4.0, 3.0, 0.0, approach(2.0 / 27.0)),
    "tw_side": (1.0, -4.0, 3.0, 1.0, approach(2.0 / 27.0)),
    "tw_outer": (1.0, -4.0, 3.0, 0.0, 2.0 / 27.0 + 10.0 ** np.linspace(-6.0, 0.0, 40)),
}

out = Path(__file__).resolve().parent / "python_period_energy.csv"
with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["family", "a", "b", "c", "x0", "E", "period"])
    for name, (a, b, c, x0, Es) in FAMILIES.items():
        m = Model(a, b, c)
        for E in Es:
            v0 = float(np.sqrt(2.0 * (E - float(m.potential(x0)))))
            w.writerow([name, a, b, c, x0, repr(float(E)), repr(period(m, x0, v0))])
print("wrote", out)
