"""NonlinearOscillatorLab: conservative polynomial-oscillator diagnostics."""
from .core import Model, simulate, stationary_points, classify_energy, turning_points, period, energy_drift
__version__ = "0.1.0"
__all__ = ["Model","simulate","stationary_points","classify_energy","turning_points","period","energy_drift"]
