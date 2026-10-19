"""Write Python-side benchmark periods for cross-implementation checks.

Output: benchmarks/python_periods.csv, read by mathematica/CrossValidation.wls.
"""
from pathlib import Path
import csv
from nonlinear_oscillator_lab import Model, period, simulate, energy_drift

CASES = [
    # id, a, b, c, x0, v0, description
    ("harmonic",        4.0,  0.0, 0.0, 0.4, 0.0, "harmonic reference, T = pi"),
    ("hardening",       1.0,  1.0, 0.0, 1.0, 0.0, "hardening Duffing"),
    ("softening",       1.0, -0.5, 0.0, 0.9, 0.0, "softening Duffing"),
    ("dw_inner",       -1.0,  1.0, 0.0, 1.2, 0.0, "double well, intra-well orbit"),
    ("dw_outer",       -1.0,  1.0, 0.0, 1.6, 0.0, "double well, inter-well orbit"),
    ("cq_demo",         1.0,  0.4, 0.2, 0.8, 0.1, "cubic-quintic demo case"),
    ("tw_center",       1.0, -4.0, 3.0, 0.3, 0.0, "triple well, central well"),
    ("tw_side",         1.0, -4.0, 3.0, 1.1, 0.0, "triple well, side well"),
    ("tw_outer",        1.0, -4.0, 3.0, 1.2, 0.3, "triple well, above barriers"),
]

out = Path(__file__).resolve().parent / "python_periods.csv"
with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["id", "a", "b", "c", "x0", "v0", "period", "energy_drift_T40", "description"])
    for cid, a, b, c, x0, v0, d in CASES:
        m = Model(a, b, c)
        T = period(m, x0, v0)
        s = simulate(m, x0, v0, t_end=40.0, n=4001)
        w.writerow([cid, a, b, c, x0, v0, repr(T), repr(energy_drift(m, s)), d])
print("wrote", out)
