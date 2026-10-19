"""Produce a sample orbit and its energy diagnostic."""
from pathlib import Path
import csv
from nonlinear_oscillator_lab import Model, simulate, energy_drift, period
m=Model(a=1,b=0.4,c=0.2)
r=simulate(m,0.8,0.1,t_end=30.,n=3001)
out=Path(__file__).resolve().parent/'orbit.csv'
with out.open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['t','x','velocity','energy'])
    for t,x,v in zip(r.t,r.y[0],r.y[1]): w.writerow([t,x,v,float(m.energy(x,v))])
print('period=',period(m,0.8,0.1),'relative_energy_drift=',energy_drift(m,r))
print('wrote',out)
