#!/usr/bin/env wolframscript
(* Writes NonlinearOscillatorLab_Demo.nb, an unevaluated tutorial notebook. *)
$dir = DirectoryName[$InputFileName];
in[s_String] := Cell[s, "Input"];
tx[s_String] := Cell[s, "Text"];
nb = Notebook[{
   Cell["NonlinearOscillatorLab \[LongDash] Wolfram Language companion", "Title"],
   tx["Conservative oscillator x'' + a x + b x^3 + c x^5 = 0. Evaluate the cells in order (Shift+Enter)."],
   in["SetDirectory[NotebookDirectory[]]; Get[\"NonlinearOscillatorLab.wl\"]; NOLVersion"],
   Cell["Model, energy geometry", "Section"],
   in["m = NOLModel[1, -4, 3];   (* triple well *)\nNOLStationaryPoints[m]"],
   in["NOLClassifyEnergy[m, 1/2, 0]"],
   in["NOLTurningPoints[m, 1/10]"],
   in["Plot[NOLPotential[m, x], {x, -1.4, 1.4}, Frame -> True, FrameLabel -> {\"x\", \"V(x)\"}]"],
   in["NOLPhasePortrait[m, {-1.35, 1.35}, {-0.75, 0.75}, \"Levels\" -> {0.01, 0.04, 0.07, 0.1, 0.15, 0.2}]"],
   Cell["Periods: machine and arbitrary precision", "Section"],
   in["NOLPeriod[NOLModel[1, 0.4, 0.2], 0.8, 0.1]"],
   in["NOLPeriod[NOLModel[1, 2/5, 1/5], 4/5, 1/10, WorkingPrecision -> 50]"],
   in["{NOLPeriod[NOLModel[-1, 1], 1.2, 0], NOLDuffingExactPeriod[-1, 1, 1.2]}   (* quadrature vs elliptic K *)"],
   in["NOLPeriod[NOLModel[-1, 1], 1, 2^(-1/2)]   (* separatrix: Failure *)"],
   Cell["Closed-form solution x = A y / Sqrt[1 + \[Lambda] y^2 + \[Mu] y^4] (A. H. Salas)", "Section"],
   in["sol = NOLExactSolution[NOLModel[1, 0.4, 0.2], 0.8, 0.1];\nKeyDrop[sol, {\"x\", \"y\"}]"],
   in["sim = NOLSimulate[NOLModel[1, 0.4, 0.2], 0.8, 0.1, 30];\nPlot[{sol[\"x\"][t], sim[\"x\"][t]}, {t, 0, 30}, PlotStyle -> {Thick, Dashed}, Frame -> True,\n PlotLegends -> {\"closed form\", \"RK8\"}]"],
   in["NOLExactSolutionError[NOLModel[1, 0.4, 0.2], 0.8, 0.1, sol, 30]"],
   Cell["Random example (as in the author's original script)", "Section"],
   in["{a0, b0, c0, x00, v00} = {RandomReal[{-12, 12}], RandomReal[{-12, 12}], RandomReal[], RandomReal[{-2, 0}], RandomReal[{-2, 2}]};\nsol = NOLExactSolution[a0, b0, c0, x00, v00];\nIf[FailureQ[sol], sol, {sol[\"Representation\"], NOLExactSolutionError[NOLModel[a0, b0, c0], x00, v00, sol, 20]}]"],
   Cell["Simulation and energy drift", "Section"],
   in["s = NOLSimulate[NOLModel[-1, 1], 1.6, 0, 100];\nNOLEnergyDrift[NOLModel[-1, 1], s]"]
   }, WindowTitle -> "NonlinearOscillatorLab demo"];
Export[FileNameJoin[{$dir, "NonlinearOscillatorLab_Demo.nb"}], nb];
Print["wrote demo notebook"];
