(* ::Package:: *)

(* ::Title:: *)
(*NonlinearOscillatorLab - Wolfram Language companion implementation*)


(* ::Text:: *)
(*Version 0.1.0. MIT License. Author: Alvaro H. Salas (ahsalass@unal.edu.co).*)
(*Conservative oscillator  x'' + a x + b x^3 + c x^5 == 0.*)
(**)
(*This package mirrors the Python API (Model, stationary_points, classify_energy,*)
(*turning_points, period, simulate, energy_drift) and adds:*)
(*  - exact / arbitrary-precision arithmetic for every routine;*)
(*  - a singularity-free period integrand obtained by polynomial factorisation*)
(*    of E - V(x) = (x - x_-)(x_+ - x) Q(x);*)
(*  - closed-form Jacobi-elliptic periods for the cubic (Duffing) case, used as an*)
(*    independent analytic reference;*)
(*  - the logarithmic separatrix asymptotics of the period.*)


BeginPackage["NonlinearOscillatorLab`"];


NOLModel::usage =
"NOLModel[a, b, c] represents the conservative oscillator x'' + a x + b x^3 + c x^5 == 0. \
NOLModel[a, b] sets c = 0 and NOLModel[a] sets b = c = 0.";
NOLPotential::usage = "NOLPotential[model, x] gives V(x) = a x^2/2 + b x^4/4 + c x^6/6.";
NOLForce::usage = "NOLForce[model, x] gives -(a x + b x^3 + c x^5).";
NOLEnergy::usage = "NOLEnergy[model, x, v] gives v^2/2 + V(x).";
NOLCurvature::usage = "NOLCurvature[model, x] gives V''(x) = a + 3 b x^2 + 5 c x^4.";
NOLStationaryPoints::usage =
"NOLStationaryPoints[model] gives the sorted real roots of V'(x) = 0 (exact when the coefficients are exact).";
NOLClassifyEnergy::usage =
"NOLClassifyEnergy[model, x0, v0] returns an Association with the energy, the equilibria, the \
minima and maxima {x, V(x)}, and whether E coincides with a barrier (separatrix) energy.";
NOLTurningPoints::usage = "NOLTurningPoints[model, E] gives all real solutions of V(x) == E, sorted.";
NOLBoundedComponent::usage =
"NOLBoundedComponent[model, x0, v0] gives {x-, x+}, the consecutive turning points enclosing x0 \
inside which E - V(x) > 0, or a Failure object.";
NOLPeriod::usage =
"NOLPeriod[model, x0, v0] gives the period of the closed orbit through (x0, v0) by quadrature of the \
regularised integrand 2 Integrate[1/Sqrt[2 Q(c + h Sin[t])], {t, -Pi/2, Pi/2}]. Returns Infinity at an \
equilibrium and a Failure object for separatrix or unbounded motion.";
NOLPeriodFromEnergy::usage =
"NOLPeriodFromEnergy[model, E, x0] gives the period of the energy-E orbit through the well that contains x0.";
NOLDuffingExactPeriod::usage =
"NOLDuffingExactPeriod[a, b, A] gives the closed-form Jacobi-elliptic period of x'' + a x + b x^3 == 0 \
for the orbit with turning point A (cn branch, or dn branch for intra-well orbits of the double well).";
NOLDuffingExactSolution::usage =
"NOLDuffingExactSolution[a, b, A, t] gives the closed-form Jacobi-elliptic solution with x(0) = A, x'(0) = 0.";
NOLSeparatrixPeriodAsymptotics::usage =
"NOLSeparatrixPeriodAsymptotics[model, xs, dE, k] gives the leading term (k/lambda) Log[1/|dE|] of the period \
near a barrier at xs, lambda = Sqrt[-V''(xs)], where k is the number of passages near the saddle per period.";
NOLSimulate::usage =
"NOLSimulate[model, x0, v0, tEnd] integrates the IVP with an 8th-order explicit Runge-Kutta method \
(the analogue of SciPy's DOP853) and returns an Association with keys \"x\", \"v\", \"tEnd\".";
NOLEnergyDrift::usage =
"NOLEnergyDrift[model, sol] gives max_j |E_j - E_0|/(1 + |E_0|) on a uniform grid of n points (option \"Samples\").";
NOLPhasePortrait::usage = "NOLPhasePortrait[model, {xmin, xmax}, {vmin, vmax}] draws energy level sets and separatrices.";
NOLDuffingSolution::usage =
"NOLDuffingSolution[p, q, y0, y1] gives an Association with the closed-form Jacobi-elliptic (or hyperbolic) \
solution \"y\" of y'' + p y + q y^3 == 0, y(0) = y0, y'(0) = y1, its case label, modulus, frequency and \
period or blow-up time (A. H. Salas).";
NOLExactSolution::usage =
"NOLExactSolution[model, x0, v0] or NOLExactSolution[a, b, c, x0, v0] gives the closed-form solution of \
x'' + a x + b x^3 + c x^5 == 0 in the form x = A y/Sqrt[1 + lambda y^2 + mu y^4], where y solves a cubic \
Duffing equation (A. H. Salas). The result is an Association with key \"x\" (a pure function of t) and the \
representation data, or a Failure object when no real representation is found.";
NOLExactSolutionError::usage =
"NOLExactSolutionError[model, x0, v0, sol, tmax] gives max |x_exact(t) - x_ref(t)| on a uniform grid of [0, tmax], \
where x_ref is a 32-digit NDSolve reference solution.";
NOLVersion::usage = "NOLVersion gives the package version string.";

NOLPeriod::sep = "Initial energy `1` coincides with a separatrix (barrier) energy: the period is infinite.";
NOLPeriod::unb = "Initial state is not in a verified bounded component.";
NOLDuffingExactPeriod::dom = "No bounded periodic orbit with turning point `1` for a = `2`, b = `3`.";
NOLSimulate::args = "tEnd must be positive.";


Begin["`Private`"];


NOLVersion = "0.1.0";

NOLModel[a_] := NOLModel[a, 0, 0];
NOLModel[a_, b_] := NOLModel[a, b, 0];

NOLPotential[NOLModel[a_, b_, c_], x_] := a x^2/2 + b x^4/4 + c x^6/6;
NOLForce[NOLModel[a_, b_, c_], x_] := -(a x + b x^3 + c x^5);
NOLEnergy[m_NOLModel, x_, v_] := v^2/2 + NOLPotential[m, x];
NOLCurvature[NOLModel[a_, b_, c_], x_] := a + 3 b x^2 + 5 c x^4;


(* ---------- helpers ---------- *)

exactQ[m_NOLModel] := And @@ (Precision[#] === Infinity & /@ List @@ m);

(* real roots of a univariate polynomial, exact when possible;
   inexact roots are Newton-polished (skipped at near-double roots) *)
newtonPolish[poly_, x_, r_] := Module[{d = D[poly, x], z = r, dz},
  Do[dz = d /. x -> z; If[Abs[dz] < 10^-8 (1 + Abs[z]), Break[]]; z = z - (poly /. x -> z)/dz, {3}];
  z];
realRoots[poly_, x_, wp_] := Module[{sol, rts},
  sol = If[Precision[poly] === Infinity,
    Quiet@Solve[poly == 0, x, Reals],
    NSolve[poly == 0, x, Reals, WorkingPrecision -> wp]];
  If[sol === {} || ! ListQ[sol], Return[{}]];
  rts = x /. sol;
  If[Precision[poly] =!= Infinity, rts = newtonPolish[poly, x, #] & /@ rts];
  SortBy[DeleteDuplicates[rts, Abs[N[#1 - #2]] < 10^-10 (1 + Abs[N[#1]]) &], N]];


(* ---------- energy geometry ---------- *)

NOLStationaryPoints[m : NOLModel[a_, b_, c_], opts : OptionsPattern[{WorkingPrecision -> MachinePrecision}]] :=
 Module[{x}, realRoots[a x + b x^3 + c x^5, x, OptionValue[WorkingPrecision]]];

Options[NOLClassifyEnergy] = {"Tolerance" -> 10^-9};
NOLClassifyEnergy[m_NOLModel, x0_, v0_, OptionsPattern[]] :=
 Module[{tol = OptionValue["Tolerance"], en, eq, mx, mn},
  en = NOLEnergy[m, x0, v0];
  eq = NOLStationaryPoints[m];
  mx = {#, NOLPotential[m, #]} & /@ Select[eq, N[NOLCurvature[m, #]] < -tol &];
  mn = {#, NOLPotential[m, #]} & /@ Select[eq, N[NOLCurvature[m, #]] > tol &];
  <|"Energy" -> en, "Equilibria" -> eq, "Minima" -> mn, "Maxima" -> mx,
    "OnSeparatrixEnergy" -> AnyTrue[mx, Abs[N[en - #[[2]]]] <= tol (1 + Abs[N[#[[2]]]]) &]|>];

NOLTurningPoints[m_NOLModel, en_, opts : OptionsPattern[{WorkingPrecision -> MachinePrecision}]] :=
 Module[{x}, realRoots[NOLPotential[m, x] - en, x, OptionValue[WorkingPrecision]]];

Options[NOLBoundedComponent] = {"Tolerance" -> 10^-9, WorkingPrecision -> MachinePrecision};
NOLBoundedComponent[m_NOLModel, x0_, v0_, OptionsPattern[]] :=
 Module[{tol = OptionValue["Tolerance"], en, r, pairs},
  en = NOLEnergy[m, x0, v0];
  r = NOLTurningPoints[m, en, WorkingPrecision -> OptionValue[WorkingPrecision]];
  pairs = Select[Partition[r, 2, 1],
    N[en - NOLPotential[m, Mean[#]]] > tol && N[#[[1]]] - tol <= N[x0] <= N[#[[2]]] + tol &];
  If[Length[pairs] == 1, First[pairs],
   Failure["Unbounded", <|"MessageTemplate" :> NOLPeriod::unb|>]]];


(* ---------- period ---------- *)

(* E - V(x) = (x - xl)(xr - x) Q(x); with x = c + h Sin[t], dx/Sqrt[2(E-V)] = dt/Sqrt[2 Q(x)] *)
regularisedPeriod[m_NOLModel, en_, {xl_, xr_}, wp_, pg_] :=
 Module[{x, t, q, ctr = (xl + xr)/2, h = (xr - xl)/2, f, nodes},
  q = PolynomialQuotient[Expand[en - NOLPotential[m, x]], Expand[(x - xl) (xr - x)], x];
  f = Function[{tt}, Evaluate[1/Sqrt[2 q /. x -> ctr + h Sin[tt]]]];
  (* split where the orbit crosses a stationary point (integrand peaks near separatrices) *)
  nodes = Join[{-Pi/2},
    Sort[ArcSin[(# - ctr)/h] & /@ Select[NOLStationaryPoints[m, WorkingPrecision -> wp], xl < # < xr &]], {Pi/2}];
  nodes = If[wp === MachinePrecision, N[nodes], SetPrecision[nodes, wp]];
  2 NIntegrate[f[t], Evaluate@Prepend[nodes, t], WorkingPrecision -> wp, PrecisionGoal -> pg,
    Method -> "GaussKronrodRule", MaxRecursion -> 30]];

Options[NOLPeriod] = {"Tolerance" -> 10^-9, WorkingPrecision -> MachinePrecision, PrecisionGoal -> Automatic};
NOLPeriod[m_NOLModel, x0_, v0_, OptionsPattern[]] :=
 Module[{tol = OptionValue["Tolerance"], wp = OptionValue[WorkingPrecision], pg, en, comp},
  pg = Replace[OptionValue[PrecisionGoal], Automatic :> If[wp === MachinePrecision, 13, wp - 6]];
  If[Abs[N[v0]] < tol && Abs[N[NOLForce[m, x0]]] < tol, Return[Infinity]];
  en = NOLEnergy[m, x0, v0];
  If[NOLClassifyEnergy[m, x0, v0, "Tolerance" -> tol]["OnSeparatrixEnergy"],
   Message[NOLPeriod::sep, N[en]];
   Return[Failure["Separatrix", <|"MessageTemplate" -> "Separatrix energy: infinite period", "Energy" -> en|>]]];
  comp = NOLBoundedComponent[m, x0, v0, "Tolerance" -> tol,
    WorkingPrecision -> If[wp === MachinePrecision, wp, wp + 20]];
  If[FailureQ[comp], Message[NOLPeriod::unb]; Return[comp]];
  (* 20 guard digits absorb the cancellation in E - V(x) near separatrices *)
  {en, comp} = If[wp === MachinePrecision, N[{en, comp}], N[{en, comp}, wp + 20]];
  regularisedPeriod[If[wp === MachinePrecision, N[m], N[m, wp + 20]], en, comp, wp, pg]];

NOLPeriodFromEnergy[m_NOLModel, en_, x0_, opts : OptionsPattern[NOLPeriod]] :=
 Module[{v2 = 2 (en - NOLPotential[m, x0])},
  If[N[v2] < 0, Failure["Inaccessible", <|"MessageTemplate" -> "x0 is not accessible at this energy"|>],
   NOLPeriod[m, x0, Sqrt[v2], opts]]];


(* ---------- closed-form Duffing references (c = 0) ---------- *)

NOLDuffingExactPeriod[a_, b_, amp_] :=
 Module[{en = a amp^2/2 + b amp^4/4, w2, mm},
  Which[
   (* intra-well orbit of the double well: x = A dn(w t | m), period 2K/w *)
   N[a] < 0 && N[b] > 0 && N[en] < 0,
    w2 = b amp^2/2; mm = 2 + 2 a/(b amp^2); 2 EllipticK[mm]/Sqrt[w2],
   (* cn branch: x = A cn(w t | m), period 4K/w *)
   N[a + b amp^2] > 0 && N[b amp^2/(2 (a + b amp^2))] < 1 && ! (N[a] < 0 && N[b] > 0 && N[en] == 0),
    w2 = a + b amp^2; mm = b amp^2/(2 w2); 4 EllipticK[mm]/Sqrt[w2],
   True, Message[NOLDuffingExactPeriod::dom, amp, a, b]; $Failed]];

NOLDuffingExactSolution[a_, b_, amp_, t_] :=
 Module[{en = a amp^2/2 + b amp^4/4, w2, mm},
  If[N[a] < 0 && N[b] > 0 && N[en] < 0,
   w2 = b amp^2/2; mm = 2 + 2 a/(b amp^2); amp JacobiDN[Sqrt[w2] t, mm],
   w2 = a + b amp^2; mm = b amp^2/(2 w2); amp JacobiCN[Sqrt[w2] t, mm]]];

NOLSeparatrixPeriodAsymptotics[m_NOLModel, xs_, dE_, k_ : 1] :=
 k/Sqrt[-NOLCurvature[m, xs]] Log[1/Abs[dE]];


(* ---------- closed-form solutions (A. H. Salas) ---------- *)
(* The three routines below are the author's original closed-form solver, kept verbatim. *)

(* =================Cubic Duffing solver for y=================*)
duffingSolution[p_?NumericQ, q_?NumericQ, y0_?NumericQ, 
   y1_?NumericQ] := 
  Catch[Module[{h, X, dd, roots, pos, A2, A, B2, b2, B, m, w, c1, k, 
     R, sg, ua, below, r, s, sn0, a}, h = y1^2 + p y0^2 + q y0^4/2; 
    X = y0^2;
    If[q == 0, 
     Throw[<|"Case" -> "Linear", 
       "y" -> Function[t, 
         Evaluate[
          Re[y0 Cos[Sqrt[p + 0. I] t] + 
            y1 Sin[Sqrt[p + 0. I] t]/Sqrt[p + 0. I]]]]|>]];
    dd = p^2 + 2 q h;
    roots = If[dd < 0, {}, Sort[{(-p + Sqrt[dd])/q, (-p - Sqrt[dd])/q}]];
    If[y1 == 0 && roots =!= {}, 
     roots = ReplacePart[roots, 
       First@Ordering[Abs[roots - X], 1] -> X]];
    pos = Select[roots, # > 0 &];
    (*orbit around 0:A cn*)
    If[h > 0 && pos =!= {} && X <= Min[pos], A2 = Min[pos]; 
     A = Sqrt[A2]; w = p + q A2; m = q A2/(2 w);
     c1 = If[y1 >= 0, -1, 1] InverseJacobiCN[Clip[y0/A, {-1, 1}], m];
     Throw[<|"Case" -> "cn: orbit around 0", "A" -> A, "m" -> m, 
       "\[Omega]" -> w, "c1" -> c1, "Period" -> 4 EllipticK[m]/Sqrt[w],
        "y" -> Function[t, Evaluate[A JacobiCN[c1 + Sqrt[w] t, m]]]|>]];
    If[q > 0,(*separatrix*)
     If[h == 0, 
      If[y0 == 0, 
       Throw[<|"Case" -> "Equilibrium", "y" -> Function[t, 0]|>]];
      a = Sqrt[-2 p/q]; c1 = If[y0 y1 > 0, -1, 1] ArcSech[Abs[y0]/a];
      Throw[<|"Case" -> "Separatrix (sech)", 
        "y" -> Function[t, Evaluate[Sign[y0] a Sech[Sqrt[-p] t + c1]]]|>]];
     (*one well:B dn*)B2 = Max[roots]; b2 = Min[roots]; B = Sqrt[B2]; 
     m = 1 - b2/B2; w = q B2/2;
     c1 = 
      If[y0 y1 > 0, -1, 1] InverseJacobiDN[
        Clip[Abs[y0]/B, {Sqrt[1 - m], 1}], m];
     Throw[<|"Case" -> "dn: oscillation inside one well", "B" -> B, 
       "m" -> m, "\[Omega]" -> w, "c1" -> c1, 
       "Period" -> 2 EllipticK[m]/Sqrt[w], 
       "y" -> Function[t, 
         Evaluate[Sign[y0] B JacobiDN[c1 + Sqrt[w] t, m]]]|>]];
    (*q<0*)k = -q;
    If[h > 0 && pos === {}, R = Sqrt[2 h/k]; m = 1/2 + p/(2 k R); 
     w = 2 k R; sg = If[y1 >= 0, 1, -1];
     c1 = Sign[y0] InverseJacobiCN[(R - X)/(R + X), m];
     Throw[<|"Case" -> "Crosses 0 and blows up", "m" -> m, 
       "\[Omega]" -> w, "c1" -> c1, 
       "BlowUpTime" -> (2 EllipticK[m] - sg c1)/Sqrt[w], 
       "y" -> Function[t, 
         Evaluate[
          Sqrt[R] JacobiSN[c1 + sg Sqrt[w] t, 
             m]/(1 + JacobiCN[c1 + sg Sqrt[w] t, m])]]|>]];
    ua = Max[Select[pos, # <= X &]];
    below = Sort[Select[Join[roots, {0}], # < ua &]];
    If[Length[below] < 2, 
     Throw[Failure["NotCovered", <|"Message" -> "Degenerate case"|>]]];
    r = below[[-1]]; s = below[[-2]];
    m = (r - s)/(ua - s); w = k/2 (ua - s);
    sn0 = Sqrt[Clip[(X - ua)/(X - r), {0, 1}]];
    c1 = If[y0 y1 >= 0, 1, -1] InverseJacobiSN[sn0, m];
    <|"Case" -> "Outside barrier: blows up", "m" -> m, 
     "\[Omega]" -> w, "c1" -> c1, 
     "BlowUpTime" -> (EllipticK[m] - c1)/Sqrt[w], 
     "y" -> Function[t, 
       Evaluate[
        Sign[y0] Sqrt[ua - r JacobiSN[c1 + Sqrt[w] t, m]^2]/
          JacobiCN[c1 + Sqrt[w] t, m]]]|>]];

(* ==========Cubic-quintic:x=A y/Sqrt[1+\[Lambda] y^2+\[Mu] y^4]==========*)

cqCandidates[\[Alpha]_, \[Beta]_, \[Gamma]_, x0_, v0_] := 
  Module[{h, X = x0^2, P, roots, out = {}, s, \[Lambda], \[Mu], p, q, 
    y1, \[Sigma], pp, j, kl, tol = 10^-9}, 
   h = v0^2 + \[Alpha] x0^2 + \[Beta] x0^4/2 + \[Gamma] x0^6/3;
   roots = 
    P /. NSolve[h P^3 - \[Alpha] P^2 - \[Beta]/2 P - \[Gamma]/3 == 0, 
      P];
   If[\[Gamma] == 0, 
    roots = Append[P /. NSolve[h P^2 - \[Alpha] P - \[Beta]/2 == 0, P],
       0.]];
   Do[If[Abs[Im[roots[[j]]]] < tol, 
     With[{pj = Re[roots[[j]]]},(*Family I:\[Mu]=0*)
      If[Abs[1 - pj X] > tol, s = 1/(1 - pj X); \[Lambda] = s pj;
       AppendTo[
        out, <|"Family" -> "I", "s" -> s, "\[Lambda]" -> \[Lambda], 
         "\[Mu]" -> 0, "p" -> \[Alpha] - 3 h pj, 
         "q" -> s (\[Beta] + 4 \[Alpha] pj - 6 h pj^2), "y1" -> s v0|>]];
      (*Family II:\[Mu]\[NotEqual]0*)kl = Delete[roots, j];
      \[Sigma] = Total[kl]; pp = Times @@ kl;
      If[
       Abs[Im[\[Sigma]]] < tol && Abs[Im[pp]] < tol && 
        1 - Re[\[Sigma]] X + Re[pp] X^2 >= 0, \[Sigma] = Re[\[Sigma]];
        pp = Re[pp];
       s = 4/(2 - \[Sigma] X + 2 Sqrt[1 - \[Sigma] X + pp X^2]);
       \[Lambda] = s \[Sigma]/2; \[Mu] = s^2 (\[Sigma]^2 - 4 pp)/16;
       If[Abs[1 - \[Mu] X^2] > tol, 
        AppendTo[
         out, <|"Family" -> "II", "s" -> s, "\[Lambda]" -> \[Lambda], 
          "\[Mu]" -> \[Mu], "p" -> h (pj - \[Sigma]/2), 
          "q" -> -h s (\[Sigma]^2 - 4 pp)/8, 
          "y1" -> s v0/(1 - \[Mu] X^2)|>]]]]], {j, Length[roots]}];
   out];

cqBuild[\[Alpha]_, \[Beta]_, \[Gamma]_, x0_, v0_, c_] := 
  Module[{ys, yf, xf, signed, h, Pf, ts, yv, Dv, dl = 10^-3, acc, tay,
     dt = 10^-6, xp, xv, res, tw = 6.}, 
   ys = duffingSolution[c["p"], c["q"], x0, c["y1"]];
   If[FailureQ[ys], Return[$Failed]];
   yf = ys["y"];
   signed = 
    StringStartsQ[ys["Case"], "cn"] || 
     StringStartsQ[ys["Case"], 
      "Linear"] || (c["\[Lambda]"] == 0 && c["\[Mu]"] == 0);
   xf = If[signed, 
     Function[t, 
      Evaluate[
       Sqrt[c["s"]] yf[t]/
         Sqrt[1 + c["\[Lambda]"] yf[t]^2 + c["\[Mu]"] yf[t]^4]]], 
     Function[t, 
      Evaluate[
       Sign[x0] Sqrt[c["s"]] Abs[yf[t]]/
         Sqrt[1 + c["\[Lambda]"] yf[t]^2 + c["\[Mu]"] yf[t]^4]]]];
   If[c["\[Lambda]"] == 0 && c["\[Mu]"] == 0 && 
     KeyExistsQ[ys, "BlowUpTime"], tw = Min[tw, 0.9 ys["BlowUpTime"]]];
   h = v0^2 + \[Alpha] x0^2 + \[Beta] x0^4/2 + \[Gamma] x0^6/3;
   Pf[u_] := h - \[Alpha] u - \[Beta] u^2/2 - \[Gamma] u^3/3;
   ts = Range[tw/400, tw, tw/400];
   yv = Select[Quiet[yf /@ ts], NumericQ[#] && Abs[#] < 10^6 &];
   Dv = 1 + c["\[Lambda]"] yv^2 + c["\[Mu]"] yv^4;
   If[Min[Dv] <= 0, Return[$Failed]];
   If[Abs[xf[0.] - x0] > 10^-8, Return[$Failed]];
   acc = -\[Alpha] x0 - \[Beta] x0^3 - \[Gamma] x0^5;
   tay = x0 + v0 dl + acc dl^2/2;
   If[Abs[tay - x0] > 0 && Abs[xf[dl] - tay] > 0.05 Abs[tay - x0], 
    Return[$Failed]];
   res = Quiet@Table[xv = xf[t]; xp = (xf[t + dt] - xf[t - dt])/(2 dt);
      If[NumericQ[xp] && Abs[yf[t]] < 10^6, 
       Abs[xp^2 - Pf[xv^2]]/(1 + Abs[Pf[xv^2]]), 0], {t, ts}];
   If[Max[res] > 10^-3, Return[$Failed]];
   <|"Representation" -> c["Family"] <> " / y: " <> ys["Case"], 
    "A" -> Sqrt[c["s"]], "\[Lambda]" -> c["\[Lambda]"], 
    "\[Mu]" -> c["\[Mu]"], "Duffing p" -> c["p"], "Duffing q" -> c["q"],
     "y(0)" -> x0, "y'(0)" -> c["y1"], "ySolution" -> KeyDrop[ys, "y"],
     "y" -> yf, "x" -> xf, 
    "Bounded y" -> ! KeyExistsQ[ys, "BlowUpTime"]|>];

cqSolution[\[Alpha]_?NumericQ, \[Beta]_?NumericQ, \[Gamma]_?NumericQ, 
   x0_?NumericQ, v0_?NumericQ] := 
  Catch[Module[{good, h, roots, u, Pf, others, xo, so, tb}, 
    good = DeleteCases[
      cqBuild[\[Alpha], \[Beta], \[Gamma], x0, v0, #] & /@ 
       cqCandidates[\[Alpha], \[Beta], \[Gamma], x0, v0], $Failed];
    If[good =!= {}, 
     Throw[First[SortBy[good, If[#["Bounded y"], 0, 1] &]]]];
    (*start at rest at a turning point that is a pole of every y:
    start from the other end of the motion and shift time*)
    If[v0 == 0, h = \[Alpha] x0^2 + \[Beta] x0^4/2 + \[Gamma] x0^6/3;
     Pf[w_] := h - \[Alpha] w - \[Beta] w^2/2 - \[Gamma] w^3/3;
     roots = 
      Cases[u /. 
        NSolve[\[Gamma]/3 u^3 + \[Beta]/2 u^2 + \[Alpha] u - h == 0, 
         u, Reals], _?NumericQ];
     others = 
      SortBy[Select[
        roots, # > 0 && Abs[# - x0^2] > 10^-9 && 
          Pf[(# + x0^2)/2] > 0 &], Abs[# - x0^2] &];
     Do[xo = Sign[x0] Sqrt[uo];
      so = 
       DeleteCases[
        cqBuild[\[Alpha], \[Beta], \[Gamma], xo, 0, #] & /@ 
         cqCandidates[\[Alpha], \[Beta], \[Gamma], xo, 0], $Failed];
      so = Select[so, ! #["Bounded y"] &];
      so = If[so === {}, $Failed, First[so]];
      If[so =!= $Failed, tb = so["ySolution"]["BlowUpTime"];
       Throw[
        Join[so, <|
          "Representation" -> so["Representation"] <> " (time-shifted)",
           "x" -> Function[t, Evaluate[so["x"][t + tb]]]|>]]], {uo, 
       others}]];
    Failure[
     "NotCovered", <|
      "Message" -> "No real representation (unbounded motion)"|>]]];


NOLDuffingSolution[p_?NumericQ, q_?NumericQ, y0_?NumericQ, y1_?NumericQ] := duffingSolution[p, q, y0, y1];
NOLExactSolution[NOLModel[a_, b_, c_], x0_, v0_] := NOLExactSolution[a, b, c, x0, v0];
NOLExactSolution[a_?NumericQ, b_?NumericQ, c_?NumericQ, x0_?NumericQ, v0_?NumericQ] :=
  cqSolution[N[a], N[b], N[c], N[x0], N[v0]];

Options[NOLExactSolutionError] = {"Samples" -> 400, WorkingPrecision -> 32};
NOLExactSolutionError[m : NOLModel[a_, b_, c_], x0_, v0_, sol_Association, tmax_, OptionsPattern[]] :=
 Module[{z, s, ref, wp = OptionValue[WorkingPrecision], ts, ap, bp, cp},
  {ap, bp, cp} = SetPrecision[{a, b, c}, wp];
  ref = NDSolveValue[{z''[s] + ap z[s] + bp z[s]^3 + cp z[s]^5 == 0,
     z[0] == SetPrecision[x0, wp], z'[0] == SetPrecision[v0, wp]}, z, {s, 0, tmax},
    WorkingPrecision -> wp, PrecisionGoal -> wp - 8, AccuracyGoal -> wp - 8, MaxSteps -> Infinity];
  ts = Subdivide[0, tmax, OptionValue["Samples"]];
  Max[Select[Quiet[Abs[(sol["x"] /@ N[ts]) - N[ref /@ ts]]], NumericQ]]];

(* ---------- simulation ---------- *)

Options[NOLSimulate] = {WorkingPrecision -> MachinePrecision, PrecisionGoal -> 10, AccuracyGoal -> 12,
   Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 8}, MaxSteps -> Infinity};
NOLSimulate[m_NOLModel, x0_, v0_, tEnd_ : 40, OptionsPattern[]] :=
 Module[{x, v, t, sol},
  If[! (N[tEnd] > 0), Message[NOLSimulate::args]; Return[$Failed]];
  sol = First@NDSolve[{x'[t] == v[t], v'[t] == NOLForce[m, x[t]], x[0] == x0, v[0] == v0},
     {x, v}, {t, 0, tEnd},
     WorkingPrecision -> OptionValue[WorkingPrecision], PrecisionGoal -> OptionValue[PrecisionGoal],
     AccuracyGoal -> OptionValue[AccuracyGoal], Method -> OptionValue[Method],
     MaxSteps -> OptionValue[MaxSteps], InterpolationOrder -> All];
  <|"x" -> (x /. sol), "v" -> (v /. sol), "tEnd" -> tEnd, "Model" -> m, "Initial" -> {x0, v0}|>];

Options[NOLEnergyDrift] = {"Samples" -> 2001};
NOLEnergyDrift[m_NOLModel, sol_Association, OptionsPattern[]] :=
 Module[{ts = Subdivide[0, sol["tEnd"], OptionValue["Samples"] - 1], es, e0},
  es = NOLEnergy[m, sol["x"][#], sol["v"][#]] & /@ ts;
  e0 = First[es];
  Max[Abs[es - e0]]/(1 + Abs[e0])];


(* ---------- visualisation ---------- *)

Options[NOLPhasePortrait] = Join[{"Levels" -> 18, "SeparatrixStyle" -> Directive[Red, Thick]}, Options[ContourPlot]];
NOLPhasePortrait[m_NOLModel, {x1_, x2_}, {v1_, v2_}, opts : OptionsPattern[]] :=
 Module[{x, v, cls, sepE, base},
  cls = NOLClassifyEnergy[m, 0, 0];
  sepE = N[cls["Maxima"][[All, 2]]];
  base = ContourPlot[NOLEnergy[m, x, v], {x, x1, x2}, {v, v1, v2},
    Evaluate@FilterRules[{opts}, Options[ContourPlot]],
    Contours -> OptionValue["Levels"], ContourShading -> None,
    ContourStyle -> GrayLevel[0.45], FrameLabel -> {"x", "v"}];
  If[sepE === {}, base,
   Show[base, ContourPlot[NOLEnergy[m, x, v], {x, x1, x2}, {v, v1, v2},
     Contours -> sepE, ContourShading -> None, ContourStyle -> OptionValue["SeparatrixStyle"],
     PlotPoints -> 60]]]];


End[];
EndPackage[];
