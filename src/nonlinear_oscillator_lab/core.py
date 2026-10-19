"""Energy methods for x'' + a*x + b*x**3 + c*x**5 = 0.

No assertion of analytic closed-form trajectories is made. Period is computed
by numerical quadrature of the factored, singularity-free integrand
obtained from E - V(x) = (x - x_-)(x_+ - x) Q(x) and x = mid + half*sin(theta).
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp, quad

@dataclass(frozen=True)
class Model:
    a: float = 1.0
    b: float = 0.0
    c: float = 0.0

    def potential(self, x):
        x = np.asarray(x)
        return 0.5*self.a*x*x + 0.25*self.b*x**4 + self.c*x**6/6

    def force(self, x):
        return -(self.a*x+self.b*x**3+self.c*x**5)

    def energy(self, x, v):
        return 0.5*np.asarray(v)**2+self.potential(x)


def stationary_points(model, tol=1e-10):
    """Return sorted real roots of dV/dx, including zero."""
    roots=[0.0]
    # c y^2 + b y + a = 0, with y=x^2>0
    coefficients=[model.c,model.b,model.a]
    while len(coefficients)>1 and abs(coefficients[0])<tol:
        coefficients.pop(0)
    if len(coefficients)>1:
        for y in np.roots(coefficients):
            if abs(y.imag)<tol and y.real>tol:
                z=float(np.sqrt(y.real));roots.extend([-z,z])
    return np.array(sorted(set(round(x,12) for x in roots)))


def classify_energy(model, x0, v0, tol=1e-9):
    """Label energy surface relative to equilibrium barrier energies.

    Labels are descriptive; `bounded_component` is not inferred merely from E.
    """
    E=float(model.energy(x0,v0))
    eq=stationary_points(model)
    barriers=[(float(x),float(model.potential(x))) for x in eq if model.a+3*model.b*x*x+5*model.c*x**4 < -tol]
    at_barrier=any(abs(E-b)<tol*(1+abs(b)) for _,b in barriers)
    return {'energy': E,'equilibria':eq.tolist(),'maxima':barriers,
            'on_separatrix_energy': bool(at_barrier)}


def turning_points(model, E, tol=1e-9):
    """All real turning points V(x)=E, without assuming confinement."""
    E=float(E)
    coef=[model.c/6,model.b/4,model.a/2,-E]
    while len(coef)>1 and abs(coef[0])<1e-14: coef.pop(0)
    yy=np.roots(coef) if len(coef)>1 else np.array([])
    xs=[]
    for y in yy:
        if abs(y.imag) < tol*(1+abs(y.real)) and y.real >= -tol:
            z=np.sqrt(max(0.,float(y.real)))
            xs.extend([-z,z] if z>tol else [0.])
    # Newton polishing of each root of V(x) - E (skipped at near-double roots)
    polished=[]
    for z in xs:
        for _ in range(3):
            dv=-float(model.force(z))
            if abs(dv) < 1e-8*(1+abs(z)): break
            z=z-(float(model.potential(z))-E)/dv
        polished.append(float(z))
    out=[]
    for z in sorted(polished):
        if not out or abs(z-out[-1]) > 1e-10*(1+abs(z)): out.append(z)
    return np.array(out)


def period(model, x0, v0, *, tol=1e-9):
    """Period in the connected bounded admissible component containing x0.

    Returns np.inf for equilibrium and raises ValueError for non-closed
    or separatrix trajectories. Period=2*integral dx/sqrt(2*(E-V(x))).
    """
    E=float(model.energy(x0,v0))
    if abs(v0)<tol and abs(float(model.force(x0)))<tol:
        return float('inf')
    if classify_energy(model,x0,v0,tol)['on_separatrix_energy']:
        raise ValueError('Separatrix energy: infinite period')
    roots=turning_points(model,E)
    pairs=[]
    for l,r in zip(roots[:-1],roots[1:]):
        mid=(l+r)/2
        if E-float(model.potential(mid))>tol and l-tol<=x0<=r+tol:
            pairs.append((float(l),float(r)))
    if len(pairs)!=1:
        raise ValueError('Initial state is not in a verified bounded component')
    left,right=pairs[0]; center=(right+left)/2; half=(right-left)/2
    # Factor E - V(x) = (x - left)(right - x) Q(x).  With x = center + half*sin(theta),
    # dx/sqrt(2(E-V)) = dtheta/sqrt(2 Q(x)): the integrand is smooth and nonsingular.
    p=np.array([-model.c/6,0.,-model.b/4,0.,-model.a/2,0.,E])
    while len(p)>1 and p[0]==0: p=p[1:]
    q,_=np.polydiv(p,np.array([-1.,left+right,-left*right]))
    def integrand(theta):
        qx=float(np.polyval(q,center+half*np.sin(theta)))
        if qx<=0:
            raise ValueError('Degenerate turning point: period quadrature not regular')
        return 1.0/np.sqrt(2*qx)
    # Near a separatrix the integrand peaks where the orbit crosses a barrier
    # (a local maximum of V); split the quadrature there.
    cuts=sorted(float(np.arcsin((xe-center)/half)) for xe in stationary_points(model)
                if left<xe<right)
    nodes=[-np.pi/2]+cuts+[np.pi/2]
    val=sum(quad(integrand,lo,hi,epsabs=1e-14,epsrel=1e-13,limit=400)[0]
            for lo,hi in zip(nodes[:-1],nodes[1:]) if hi>lo)
    if not np.isfinite(val) or val<=0:
        raise ValueError('Period quadrature failed')
    return float(2*val)


def simulate(model, x0, v0, t_end=40., n=2001, rtol=1e-10, atol=1e-12):
    """Integrate IVP using adaptive DOP853; returns scipy OdeResult."""
    if t_end<=0 or n<2: raise ValueError('t_end must be positive and n>=2')
    def f(t,z): return [z[1],float(model.force(z[0]))]
    t=np.linspace(0,t_end,n)
    result=solve_ivp(f,(0.,t_end),[x0,v0],method='DOP853',t_eval=t,rtol=rtol,atol=atol)
    if not result.success: raise RuntimeError(result.message)
    return result


def energy_drift(model, solution):
    energies=model.energy(solution.y[0],solution.y[1])
    E0=float(energies[0]); return float(np.max(np.abs(energies-E0))/(1+abs(E0)))
