# NonlinearOscillatorLab (v0.1.0)

A small, research-oriented, MIT-licensed toolkit -- a Python package and a
Wolfram Language (Mathematica) package -- for the conservative Duffing and
cubic--quintic oscillator

    x'' + a*x + b*x**3 + c*x**5 = 0.

It computes equilibria, energies, polynomial turning points, periods of bounded
non-separatrix trajectories (singularity-free factored quadrature), and IVP
simulations. The Wolfram Language package adds exact/arbitrary-precision
arithmetic, closed-form Jacobi-elliptic periods and the closed-form solution
x = A y / sqrt(1 + lambda y^2 + mu y^4) of the cubic-quintic equation.

## Installation

From the extracted project folder:

```bash
python -m pip install -e ".[test]"
python -m pytest -q
python examples/demo.py
```

## API

```python
from nonlinear_oscillator_lab import Model, simulate, period, energy_drift
m = Model(a=1, b=0.4, c=0.2)
s = simulate(m, x0=0.8, v0=0.1)
print(period(m, 0.8, 0.1), energy_drift(m, s))
```

## Wolfram Language (Mathematica 14+, tested on 15.0)

```
wolframscript -file mathematica/RunTests.wls      # 31 verification tests
wolframscript -file mathematica/MakeFigures.wls   # all paper figures and tables
wolframscript -file mathematica/RunAll.wls        # both
```

Interactive use: open `mathematica/NonlinearOscillatorLab_Demo.nb`, or

```mathematica
Get["mathematica/NonlinearOscillatorLab.wl"]
NOLPeriod[NOLModel[1, 2/5, 1/5], 4/5, 1/10, WorkingPrecision -> 50]
sol = NOLExactSolution[NOLModel[1, 0.4, 0.2], 0.8, 0.1]; sol["x"][3.]
```

Outputs: `figures/` (PDF + 300 dpi PNG), `tables/` (CSV and LaTeX).
Python-side reference data used by the figures is produced by
`benchmarks/python_benchmarks.py` and `benchmarks/python_period_energy.py`.

## Computational definition

Potential V(x) = a*x^2/2 + b*x^4/4 + c*x^6/6.
The trajectory is confined only if the admissible energy interval containing
x0 is bounded by two consecutive real roots of V(x)=E. The period uses the exact factorisation
E - V(x) = (x - x_-)(x_+ - x) Q(x); with x = midpoint + halfwidth*sin(theta) the
integrand becomes 1/sqrt(2 Q(x)), smooth on the whole interval. The quadrature is
split where the orbit crosses a barrier.
The code detects stationary separatrix energies but does not assert universal
symbolic closed forms or stability certification.

## Known limitations

* Conservative, autonomous single-degree-of-freedom systems only.
* No forcing, damping or events/impact dynamics.
* Near a separatrix, quadrature is ill-conditioned; an explicit error is raised
  when energy matches a local-maximum barrier at the chosen tolerance.
* The closed-form solver exists only in the Wolfram Language package; its
  coverage is empirical (see the randomised campaign in tables/).
* Scientific verification here is limited to supplied automated tests and
  benchmarks; independently reviewed benchmarks are still needed.
* Period computation does not cover all degenerate multiple-root cases.

## Repository and citation

Source: https://github.com/alvarohumberto/NonlinearOscillatorLab
Paper: A. H. Salas, NonlinearOscillatorLab (SoftwareX, submitted 2026).

Contact: ahsalass@unal.edu.co
