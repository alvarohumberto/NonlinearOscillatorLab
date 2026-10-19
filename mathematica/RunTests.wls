#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Automated verification suite for the Wolfram Language companion of NonlinearOscillatorLab.
   Run:  wolframscript -file RunTests.wls
   Writes test_report.txt and test_report.csv next to this script. *)

$dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{$dir, "NonlinearOscillatorLab.wl"}]];
$root = ParentDirectory[$dir];

relErr[x_, y_] := Abs[x - y]/Abs[y];

tests = {
  (* ---- the seven Python regression tests, re-implemented ---- *)
  VerificationTest[Abs[NOLPeriod[NOLModel[4], 0.4, 0] - Pi] < 10^-12, True, TestID -> "harmonic_period"],
  VerificationTest[
   Module[{s = NOLSimulate[NOLModel[4], 0.4, 0, 10]},
    {Max[Abs[s["x"][#] - 0.4 Cos[2 #]] & /@ Subdivide[0., 10., 2000]] < 2*^-8,
     NOLEnergyDrift[NOLModel[4], s] < 2*^-9}], {True, True}, TestID -> "harmonic_solution_and_energy"],
  VerificationTest[NOLStationaryPoints[NOLModel[-1, 1]], {-1, 0, 1}, TestID -> "double_well_equilibria_exact"],
  VerificationTest[NOLTurningPoints[NOLModel[1], 1/2], {-1, 1}, TestID -> "turning_points_exact"],
  VerificationTest[
   {NOLClassifyEnergy[NOLModel[-1, 1], 1, 2^(-1/2)]["OnSeparatrixEnergy"],
    FailureQ[Quiet@NOLPeriod[NOLModel[-1, 1], 1, 2^(-1/2)]]}, {True, True}, TestID -> "separatrix_detection"],
  VerificationTest[
   Module[{m = NOLModel[1, 0.4, 0.2], s},
    s = NOLSimulate[m, 0.8, 0.1, 20];
    {NOLEnergyDrift[m, s] < 5*^-8, NOLPeriod[m, 0.8, 0.1] > 0}], {True, True}, TestID -> "cubic_quintic"],
  VerificationTest[NOLPeriod[NOLModel[1], 0, 0], Infinity, TestID -> "equilibrium_period"],

  (* ---- new: arbitrary precision ---- *)
  VerificationTest[Abs[NOLPeriod[NOLModel[4], 2/5, 0, WorkingPrecision -> 40] - Pi] < 10^-30, True,
   TestID -> "harmonic_period_40_digits"],
  VerificationTest[
   Module[{hp = NOLPeriod[NOLModel[1, 2/5, 1/5], 4/5, 1/10, WorkingPrecision -> 40]},
    {Precision[hp] >= 30, relErr[NOLPeriod[NOLModel[1, 0.4, 0.2], 0.8, 0.1], hp] < 10^-13}], {True, True},
   TestID -> "cubic_quintic_machine_vs_40_digits"],

  (* ---- new: independent closed-form (Jacobi elliptic) references ---- *)
  VerificationTest[relErr[NOLPeriod[NOLModel[1, 1], 1., 0], NOLDuffingExactPeriod[1, 1, 1.]] < 10^-13, True,
   TestID -> "elliptic_hardening"],
  VerificationTest[relErr[NOLPeriod[NOLModel[1, -1/2], 0.9, 0], NOLDuffingExactPeriod[1, -1/2, 0.9]] < 10^-13, True,
   TestID -> "elliptic_softening"],
  VerificationTest[relErr[NOLPeriod[NOLModel[-1, 1], 1.2, 0], NOLDuffingExactPeriod[-1, 1, 1.2]] < 10^-13, True,
   TestID -> "elliptic_double_well_inner_dn"],
  VerificationTest[relErr[NOLPeriod[NOLModel[-1, 1], 1.6, 0], NOLDuffingExactPeriod[-1, 1, 1.6]] < 10^-13, True,
   TestID -> "elliptic_double_well_outer_cn"],
  VerificationTest[
   Module[{p = NOLPeriod[NOLModel[-1, 1], 12/10, 0, WorkingPrecision -> 50],
     q = N[NOLDuffingExactPeriod[-1, 1, 12/10], 50]}, Abs[p - q] < 10^-40], True,
   TestID -> "elliptic_double_well_50_digits"],
  VerificationTest[
   Module[{t, res},
    res = Table[Module[{xs = NOLDuffingExactSolution[a, b, amp, t]},
       Quiet[N[(D[xs, {t, 2}] + a xs + b xs^3) /. t -> 7/10, 30], {N::meprec, General::stop}]],
      {p, {{1, 1, 1}, {1, -1/2, 9/10}, {-1, 1, 6/5}, {-1, 1, 8/5}}}, {a, {p[[1]]}}, {b, {p[[2]]}}, {amp, {p[[3]]}}];
    Max[Abs[Flatten[res]]] < 10^-25], True, TestID -> "elliptic_solutions_satisfy_ODE"],
  VerificationTest[
   Module[{s = NOLSimulate[NOLModel[-1, 1], 1.2, 0, 20]},
    Max[Abs[s["x"][#] - NOLDuffingExactSolution[-1, 1, 1.2, #]] & /@ Subdivide[0., 20., 400]] < 10^-8], True,
   TestID -> "simulation_vs_elliptic_trajectory"],

  (* ---- new: geometry of the cubic-quintic triple well ---- *)
  VerificationTest[NOLStationaryPoints[NOLModel[1, -4, 3]], {-1, -1/Sqrt[3], 0, 1/Sqrt[3], 1},
   TestID -> "triple_well_equilibria_exact"],
  VerificationTest[Simplify[NOLPotential[NOLModel[1, -4, 3], 1/Sqrt[3]]], 2/27, TestID -> "triple_well_barrier"],

  (* ---- new: period/trajectory consistency ---- *)
  VerificationTest[
   Module[{m = NOLModel[1, 0.4, 0.2], T, s},
    T = NOLPeriod[m, 0.8, 0.1]; s = NOLSimulate[m, 0.8, 0.1, 1.5 T];
    Norm[{s["x"][T] - 0.8, s["v"][T] - 0.1}] < 10^-8], True, TestID -> "period_closes_orbit"],
  VerificationTest[FailureQ[Quiet@NOLPeriod[NOLModel[1, -1], 2., 0]], True, TestID -> "unbounded_detected"],

  (* ---- new: separatrix asymptotics T ~ (1/lambda) Log[1/dE] + C ---- *)
  VerificationTest[
   Module[{m = NOLModel[-1, 1], dE = 10^-10, ti, to},
    ti = NOLPeriodFromEnergy[m, -dE, 1, WorkingPrecision -> 40, "Tolerance" -> 10^-30];
    to = NOLPeriodFromEnergy[m, dE, 0, WorkingPrecision -> 40, "Tolerance" -> 10^-30];
    {Abs[ti - NOLSeparatrixPeriodAsymptotics[m, 0, dE, 1] - Log[16]] < 10^-6,
     Abs[to - NOLSeparatrixPeriodAsymptotics[m, 0, dE, 2] - 2 Log[16]] < 10^-6}], {True, True},
   TestID -> "separatrix_log_asymptotics"],

  (* ---- new: closed-form solutions x = A y/Sqrt[1 + lambda y^2 + mu y^4] (A. H. Salas) ---- *)
  Sequence @@ Table[With[{cs = cs},
     VerificationTest[
      Module[{m = NOLModel @@ cs[[2 ;; 4]], sol},
       sol = NOLExactSolution[m, cs[[5]], cs[[6]]];
       If[FailureQ[sol], sol, NOLExactSolutionError[m, cs[[5]], cs[[6]], sol, 20] < 10^-8]], True,
      TestID -> "closed_form_" <> cs[[1]]]],
    {cs, {{"hardening", 1, 1, 0, 1., 0.}, {"softening", 1, -1/2, 0, 0.9, 0.},
      {"dw_inner", -1, 1, 0, 1.2, 0.}, {"dw_outer", -1, 1, 0, 1.6, 0.},
      {"cq_demo", 1, 0.4, 0.2, 0.8, 0.1}, {"tw_center", 1, -4, 3, 0.3, 0.},
      {"tw_side", 1, -4, 3, 1.1, 0.}, {"tw_outer", 1, -4, 3, 1.2, 0.3}}}],
  VerificationTest[
   Module[{d = NOLDuffingSolution[1., 1., 1., 0.]}, relErr[d["Period"], NOLDuffingExactPeriod[1, 1, 1.]] < 10^-14],
   True, TestID -> "closed_form_duffing_period"],

  (* ---- new: cross-implementation agreement with the Python package ---- *)
  VerificationTest[
   Module[{rows = Rest@Import[FileNameJoin[{$root, "benchmarks", "python_periods.csv"}], "CSV"]},
    Max[Table[relErr[NOLPeriod[NOLModel @@ r[[2 ;; 4]], r[[5]], r[[6]]], r[[7]]], {r, rows}]] < 10^-9], True,
   TestID -> "python_cross_validation"]
};

report = TestReport[tests];
lines = Table[{t["TestID"], ToString[t["Outcome"]], Round[1000 QuantityMagnitude[t["AbsoluteTimeUsed"], "Seconds"], 0.1]},
   {t, Values[report["TestResults"]]}];
Export[FileNameJoin[{$dir, "test_report.csv"}], Prepend[lines, {"test", "outcome", "ms"}], "CSV"];
summary = StringJoin[
   "NonlinearOscillatorLab (Wolfram Language) ", NOLVersion, " - ", $Version, "\n",
   "Passed ", ToString[report["TestsSucceededCount"]], " / ", ToString[Length[tests]], "\n\n",
   StringRiffle[StringRiffle[ToString /@ #, "  "] & /@ lines, "\n"], "\n"];
Export[FileNameJoin[{$dir, "test_report.txt"}], summary, "Text"];
Print[summary];
Do[If[t["Outcome"] =!= "Success", Print[t["TestID"], ": ", t["ActualOutput"]]], {t, Values[report["TestResults"]]}];
