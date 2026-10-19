#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Publication figures and tables for NonlinearOscillatorLab, generated with the Wolfram Language.
   Run:  wolframscript -file MakeFigures.wls
   Output: ../figures/*.pdf, ../figures/*.png, ../tables/*.csv, ../tables/*.tex *)

$dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{$dir, "NonlinearOscillatorLab.wl"}]];
$root = ParentDirectory[$dir];
$fig = FileNameJoin[{$root, "figures"}]; $tab = FileNameJoin[{$root, "tables"}];
Quiet[CreateDirectory /@ {$fig, $tab}];
$t0 = AbsoluteTime[];
log[s__] := Print[Round[AbsoluteTime[] - $t0, 0.1], "s  ", s];

(* ---------------- style ---------------- *)
cols = RGBColor /@ {"#1f5fa8", "#d1661a", "#2a8c55", "#7a4fa3", "#c0392b", "#5d6d7e", "#b8860b", "#17becf"};
font = {FontFamily -> "Times New Roman", FontSize -> 11};
frameOpts = Sequence[Frame -> True, FrameStyle -> Directive[Black, AbsoluteThickness[0.7]],
   FrameTicksStyle -> Directive[Black, 10], BaseStyle -> font, ImageSize -> 330, AspectRatio -> 0.72,
   ImagePadding -> {{52, 12}, {40, 10}}];
lbl[g_, s_String] := Labeled[g, Style[s, Bold, 12, FontFamily -> "Times New Roman"], {{Top, Left}}];
row[gs__] := Grid[{{gs}}, Spacings -> {1.5, 0}, Alignment -> Top];
save[name_, g_] := (Export[FileNameJoin[{$fig, name <> ".pdf"}], g];
   Export[FileNameJoin[{$fig, name <> ".png"}], g, ImageResolution -> 300]; log["wrote ", name]);
(* legends always below the frame, two entries per row: never on top of data *)
leg[cs_, ls_, ___] := Placed[LineLegend[cs, ls, LegendMarkerSize -> {20, 10}, LabelStyle -> Append[font, 10],
    LegendLayout -> {"Row", 2}], Below];
ptleg[cs_, ls_, ___] := Placed[PointLegend[cs, ls, LegendMarkerSize -> 8, LabelStyle -> Append[font, 10],
    LegendLayout -> {"Row", 2}], Below];
num[x_String] := If[StringMatchQ[x, NumberString ~~ ("e" | "E") ~~ ("+" | "-" | "") ~~ DigitCharacter ..],
   ToExpression[StringReplace[x, {"e" -> "*^", "E" -> "*^"}]], x];
num[x_?NumericQ] := x;
eqMarks[m_, col_] := Module[{c = NOLClassifyEnergy[m, 0, 0]},
   {col, PointSize[0.018], Point[N@c["Minima"]], White, EdgeForm[{col, AbsoluteThickness[1.2]}],
    Disk[#, Offset[3.2]] & /@ N@c["Maxima"]}];

(* ---------------- Fig. 1: potentials and equilibria ---------------- *)
duff = {NOLModel[1, 1], NOLModel[1, -1/2], NOLModel[-1, 1]};
cq = {NOLModel[1, 2/5, 1/5], NOLModel[1, -4, 3]};
f1a = Plot[Evaluate[NOLPotential[#, x] & /@ duff], {x, -2.3, 2.3}, PlotRange -> {-0.4, 1.6},
   PlotStyle -> (Directive[#, AbsoluteThickness[1.6]] & /@ cols[[1 ;; 3]]),
   Epilog -> MapThread[eqMarks, {duff, cols[[1 ;; 3]]}], FrameLabel -> {"x", "V(x)"},
   PlotLegends -> leg[cols[[1 ;; 3]], {"hardening a=1, b=1", "softening a=1, b=-1/2", "double well a=-1, b=1"}, Scaled[{0.5, 0.86}]],
   Evaluate[frameOpts]];
f1b = Plot[Evaluate[NOLPotential[#, x] & /@ cq], {x, -1.4, 1.4}, PlotRange -> {-0.03, 0.32},
   PlotStyle -> (Directive[#, AbsoluteThickness[1.6]] & /@ cols[[4 ;; 5]]),
   Epilog -> MapThread[eqMarks, {cq, cols[[4 ;; 5]]}], FrameLabel -> {"x", "V(x)"},
   PlotLegends -> leg[cols[[4 ;; 5]], {"a=1, b=2/5, c=1/5", "triple well a=1, b=-4, c=3"}, Scaled[{0.5, 0.88}]],
   Evaluate[frameOpts]];
save["fig1_potentials", row[lbl[f1a, "(a) cubic (Duffing)"], lbl[f1b, "(b) cubic\[Dash]quintic"]]];

(* ---------------- Fig. 2: phase portraits ---------------- *)
f2a = NOLPhasePortrait[NOLModel[-1, 1], {-2, 2}, {-1.6, 1.6},
   "Levels" -> DeleteCases[Range[-0.22, 1.3, 0.08], _?(Abs[#] < 0.01 &)], PlotPoints -> 60,
   Epilog -> eqMarks[NOLModel[-1, 1], Black] /. {x_?NumericQ, y_?NumericQ} :> {x, 0},
   Evaluate[frameOpts]];
f2b = NOLPhasePortrait[NOLModel[1, -4, 3], {-1.35, 1.35}, {-0.75, 0.75},
   "Levels" -> {0.005, 0.02, 0.04, 0.06, 0.09, 0.13, 0.18, 0.24}, PlotPoints -> 70,
   Epilog -> eqMarks[NOLModel[1, -4, 3], Black] /. {x_?NumericQ, y_?NumericQ} :> {x, 0},
   Evaluate[frameOpts]];
save["fig2_phase_portraits", row[lbl[f2a, "(a) double well, a=-1, b=1"], lbl[f2b, "(b) triple well, a=1, b=-4, c=3"]]];

(* ---------------- Fig. 3: closed form vs numerics ---------------- *)
ref32[m : NOLModel[a_, b_, c_], x0_, v0_, tmax_] := Module[{z, s},
   NDSolveValue[{z''[s] + SetPrecision[a, 32] z[s] + SetPrecision[b, 32] z[s]^3 + SetPrecision[c, 32] z[s]^5 == 0,
     z[0] == SetPrecision[x0, 32], z'[0] == SetPrecision[v0, 32]}, z, {s, 0, tmax},
    WorkingPrecision -> 32, PrecisionGoal -> 24, AccuracyGoal -> 24, MaxSteps -> Infinity]];
demo = {NOLModel[1, 0.4, 0.2], 0.8, 0.1};
exDemo = NOLExactSolution @@ demo;
log["demo representation: ", exDemo["Representation"]];
refDemo = ref32[demo[[1]], demo[[2]], demo[[3]], 30];
simDemo = NOLSimulate[demo[[1]], demo[[2]], demo[[3]], 30];
pyOrbit = Rest@Import[FileNameJoin[{$root, "examples", "orbit.csv"}], "CSV"];
pyOrbit = Map[num, pyOrbit, {2}];
f3a = Legended[Show[
   Plot[exDemo["x"][t], {t, 0, 30}, PlotStyle -> Directive[cols[[1]], AbsoluteThickness[1.6]], PlotPoints -> 400,
    Exclusions -> None],
   ListPlot[{pyOrbit[[1 ;; -1 ;; 60, {1, 2}]], Table[{t, simDemo["x"][t]}, {t, 1.5, 30, 3.}]},
    PlotStyle -> {cols[[2]], cols[[3]]},
    PlotMarkers -> {{"\[FilledCircle]", 7}, {"\[EmptySquare]", 9}}],
   FrameLabel -> {"t", "x(t)"}, PlotRange -> {-1.05, 1.05}, Evaluate[frameOpts]],
  Placed[Row[{LineLegend[{Directive[cols[[1]], AbsoluteThickness[1.6]]}, {"closed form (Salas)"}, LabelStyle -> Append[font, 10]],
     PointLegend[{cols[[2]], cols[[3]]}, {"Python DOP853", "Wolfram RK8"}, LegendMarkers -> {"\[FilledCircle]", "\[EmptySquare]"},
      LegendLayout -> "Row", LabelStyle -> Append[font, 10]]}, "  "], Below]];
errCurve[f_] := Table[{t, Max[Abs[f[t] - refDemo[t]], 10^-17]}, {t, 0.05, 30, 0.05}];
pyInterp = Interpolation[pyOrbit[[All, {1, 2}]], InterpolationOrder -> 5];
f3b = ListLogPlot[{errCurve[exDemo["x"]], errCurve[simDemo["x"]],
    Table[{pyOrbit[[k, 1]], Max[Abs[pyOrbit[[k, 2]] - refDemo[pyOrbit[[k, 1]]]], 10^-17]}, {k, 2, Length[pyOrbit], 5}]},
   Joined -> True, PlotStyle -> (Directive[#, AbsoluteThickness[1.1]] & /@ cols[[1 ;; 3]]),
   PlotRange -> {10^-17, 10^-6}, FrameLabel -> {"t", "|x(t) - x_ref(t)|"},
   PlotLegends -> leg[cols[[1 ;; 3]], {"closed form", "Wolfram RK8", "Python DOP853"}, Scaled[{0.75, 0.2}]],
   Evaluate[frameOpts]];
cases3 = {{"tw_side", NOLModel[1, -4, 3], 1.1, 0.}, {"tw_outer", NOLModel[1, -4, 3], 1.2, 0.3},
   {"dw_inner", NOLModel[-1, 1], 1.2, 0.}};
ex3 = NOLExactSolution[#[[2]], #[[3]], #[[4]]] & /@ cases3;
f3c = Show[
   Plot[Evaluate[#["x"][t] & /@ ex3], {t, 0, 20}, PlotPoints -> 400, Exclusions -> None,
    PlotStyle -> (Directive[#, AbsoluteThickness[1.5]] & /@ cols[[{4, 5, 6}]]),
    PlotLegends -> leg[cols[[{4, 5, 6}]], {"triple well, side", "triple well, above barrier", "double well, intra-well"}, Scaled[{0.5, 0.12}]]],
   ListPlot[Table[With[{r = ref32[c[[2]], c[[3]], c[[4]], 20]}, Table[{t, r[t]}, {t, 0.5, 20, 1}]], {c, cases3}],
    PlotStyle -> Black, PlotMarkers -> {{"\[EmptyCircle]", 7}}],
   FrameLabel -> {"t", "x(t)"}, PlotRange -> {-1.6, 1.6},
   Evaluate[frameOpts]];
refs3 = Table[ref32[c[[2]], c[[3]], c[[4]], 20], {c, cases3}];
f3d = ListLogPlot[Table[Table[{t, Max[Abs[ex3[[k]]["x"][t] - refs3[[k]][t]], 10^-17]}, {t, 0.05, 20, 0.05}], {k, 3}],
   Joined -> True, PlotStyle -> (Directive[#, AbsoluteThickness[1.1]] & /@ cols[[{4, 5, 6}]]),
   PlotRange -> {10^-17, 10^-6}, FrameLabel -> {"t", "|x(t) - x_ref(t)|"},
   PlotLegends -> leg[cols[[{4, 5, 6}]], {"triple well, side", "triple well, above barrier", "double well, intra-well"}],
   Evaluate[frameOpts]];
save["fig3_closed_form_vs_numerics",
  Grid[{{lbl[f3a, "(a) a=1, b=0.4, c=0.2, x0=0.8, v0=0.1"], lbl[f3b, "(b) pointwise error vs 32-digit reference"]},
    {lbl[f3c, "(c) closed form (lines) vs 32-digit NDSolve (circles)"], lbl[f3d, "(d) pointwise error of the closed form in (c)"]}},
   Spacings -> {1.5, 1.5}, Alignment -> Top]];
save["fig3a_demo_orbit", f3a];

(* ---------------- Fig. 4: randomised campaign of the closed-form solver ---------------- *)
SeedRandom[20261019];
nCampaign = 200; tCampaign = 20;
campaign = Table[
   Module[{a = RandomReal[{-12, 12}], b = RandomReal[{-12, 12}], c = RandomReal[{0, 1}],
     x0 = RandomReal[{-2, 0}], v0 = RandomReal[{-2, 2}], sol, err, tm},
    {tm, sol} = AbsoluteTiming[Quiet@NOLExactSolution[a, b, c, x0, v0]];
    err = If[FailureQ[sol], Missing["Failure"],
      Quiet@NOLExactSolutionError[NOLModel[a, b, c], x0, v0, sol, tCampaign, "Samples" -> 200, WorkingPrecision -> 28]];
    <|"a" -> a, "b" -> b, "c" -> c, "x0" -> x0, "v0" -> v0, "E" -> NOLEnergy[NOLModel[a, b, c], x0, v0],
     "ok" -> ! FailureQ[sol], "rep" -> If[FailureQ[sol], "failure", sol["Representation"]],
     "family" -> If[FailureQ[sol], "failure", First@StringSplit[sol["Representation"], " "]],
     "err" -> err, "time" -> tm|>], {k, nCampaign}];
log["campaign done"];
Export[FileNameJoin[{$tab, "closed_form_campaign.csv"}],
  Prepend[Values /@ (KeyDrop[#, {}] & /@ campaign), Keys[First[campaign]]], "CSV"];
good = Select[campaign, #["ok"] && NumericQ[#["err"]] &];
errs = #["err"] & /@ good;
fams = Union[#["family"] & /@ good];
famCol = AssociationThread[fams, cols[[1 ;; Length[fams]]]];
f4a = Histogram[Log10[Clip[errs, {10^-16, Infinity}]], {-16, Max[-6, Ceiling[Max[Log10[errs]]]], 0.5},
   ChartStyle -> Directive[cols[[1]], EdgeForm[White]], FrameLabel -> {"log10 max|x_closed - x_ref| on [0, 20]", "cases"},
   Evaluate[frameOpts]];
f4b = ListLogPlot[Table[{#["E"], Max[#["err"], 10^-16]} & /@ Select[good, #["family"] == f &], {f, fams}],
   PlotStyle -> (Directive[famCol[#], PointSize[0.012]] & /@ fams), FrameLabel -> {"energy E", "max error"},
   PlotLegends -> ptleg[famCol /@ fams, ("family " <> # & /@ fams), Scaled[{0.8, 0.85}]],
   Evaluate[frameOpts]];
save["fig4_closed_form_campaign", row[lbl[f4a, "(a) error distribution, " <> ToString[Length[good]] <> "/" <> ToString[nCampaign] <> " cases"],
   lbl[f4b, "(b) error versus energy"]]];
campaignSummary = <|"cases" -> nCampaign, "represented" -> Count[campaign, _?(#["ok"] &)],
   "median_error" -> Median[errs], "p95_error" -> Quantile[errs, 0.95], "max_error" -> Max[errs],
   "fraction_below_1e-8" -> N[Count[errs, e_ /; e < 10^-8]/nCampaign],
   "median_time_s" -> Median[#["time"] & /@ campaign],
   "representations" -> Counts[#["rep"] & /@ campaign]|>;
Export[FileNameJoin[{$tab, "closed_form_campaign_summary.txt"}], ToString[campaignSummary, InputForm], "Text"];
log["campaign summary ", campaignSummary];

(* ---------------- Fig. 5: period-energy curves ---------------- *)
pe = Rest@Import[FileNameJoin[{$root, "benchmarks", "python_period_energy.csv"}], "CSV"];
peLegacy = Rest@Import[FileNameJoin[{$root, "benchmarks", "python_period_energy_legacy_v010.csv"}], "CSV"];
pe = Map[num, pe, {2}]; peLegacy = Map[num, peLegacy, {2}];
families = DeleteDuplicates[pe[[All, 1]]];
famData = AssociationMap[Function[f, Select[pe, #[[1]] === f &]], families];
exactCoeffs = <|"hardening" -> {1, 1, 0}, "softening" -> {1, -1/2, 0}, "dw_inner" -> {-1, 1, 0}, "dw_outer" -> {-1, 1, 0},
   "cq_demo" -> {1, 2/5, 1/5}, "tw_center" -> {1, -4, 3}, "tw_side" -> {1, -4, 3}, "tw_outer" -> {1, -4, 3}|>;
(* machine-precision Wolfram periods and 40-digit references at the same (binary) energies *)
wlT = AssociationMap[Function[f, Table[Quiet[NOLPeriodFromEnergy[NOLModel @@ N[exactCoeffs[f]], r[[6]], r[[5]]],
      {NIntegrate::slwcon, NIntegrate::eincr}], {r, famData[f]}]], families];
log["machine periods done"];
refT = AssociationMap[Function[f, Table[
      Module[{m = NOLModel @@ exactCoeffs[f], en = SetPrecision[r[[6]], 60], x0 = Rationalize[r[[5]]], amp},
       If[exactCoeffs[f][[3]] === 0,
        (* turning point of the orbit through x0 (not the outer root beyond a barrier) *)
        amp = Last@NOLBoundedComponent[m, x0, Sqrt[2 (en - NOLPotential[m, x0])], WorkingPrecision -> 40, "Tolerance" -> 10^-30];
        NOLDuffingExactPeriod[exactCoeffs[f][[1]], exactCoeffs[f][[2]], amp],
        NOLPeriodFromEnergy[m, en, x0, WorkingPrecision -> 40, "Tolerance" -> 10^-30]]], {r, famData[f]}]], families];
log["reference periods done"];
duffF = {"hardening", "softening", "dw_inner", "dw_outer"}; cqF = {"cq_demo", "tw_center", "tw_side", "tw_outer"};
peNames = <|"hardening" -> "hardening", "softening" -> "softening", "dw_inner" -> "double well, inner", "dw_outer" -> "double well, outer",
   "cq_demo" -> "a=1, b=2/5, c=1/5", "tw_center" -> "triple well, centre", "tw_side" -> "triple well, side", "tw_outer" -> "triple well, outer"|>;
pePlot[fs_, cs_, ttl_] := Show[
   ListLinePlot[Table[Transpose[{famData[f][[All, 6]], wlT[f]}], {f, fs}], PlotStyle -> (Directive[#, AbsoluteThickness[1.4]] & /@ cs)],
   ListPlot[Table[Transpose[{famData[f][[All, 6]], N[refT[f]]}][[1 ;; -1 ;; 3]], {f, fs}],
    PlotStyle -> cs, PlotMarkers -> {{"\[EmptyCircle]", 6}}],
   FrameLabel -> {"energy E", "period T(E)"}, PlotRange -> All, Evaluate[frameOpts]];
f5a = pePlot[duffF, cols[[1 ;; 4]], ""]; f5b = pePlot[cqF, cols[[5 ;; 8]], ""];
f5a = Legended[Show[f5a, PlotRange -> {{-0.3, 5}, {0, 16}}], leg[cols[[1 ;; 4]], peNames /@ duffF]];
f5b = Legended[Show[f5b, PlotRange -> {{0, 1.2}, {0, 20}}], leg[cols[[5 ;; 8]], peNames /@ cqF]];
save["fig5_period_energy", row[lbl[f5a, "(a) Duffing: quadrature (lines), elliptic K (circles)"],
   lbl[f5b, "(b) cubic\[Dash]quintic: quadrature, 40-digit reference"]]];

(* ---------------- Fig. 6: accuracy of the period and separatrix asymptotics ---------------- *)
relE[a_, b_] := Max[N[Abs[a - b]/Abs[b]], 10^-17];
legacyByF = AssociationMap[Function[f, Select[peLegacy, #[[1]] === f &][[All, 7]]], families];
strip[k_, vals_, off_, col_, mk_] := Table[{k + off + RandomReal[{-0.08, 0.08}], Log10[v]}, {v, vals}];
SeedRandom[1];
f6a = ListPlot[
   Join[
    {Flatten[Table[strip[k, MapThread[relE, {legacyByF[families[[k]]], refT[families[[k]]]}], -0.28, 0, 0], {k, 8}], 1]},
    {Flatten[Table[strip[k, MapThread[relE, {famData[families[[k]]][[All, 7]], refT[families[[k]]]}], 0, 0, 0], {k, 8}], 1]},
    {Flatten[Table[strip[k, MapThread[relE, {wlT[families[[k]]], refT[families[[k]]]}], 0.28, 0, 0], {k, 8}], 1]}],
   PlotStyle -> {Directive[cols[[6]], PointSize[0.008]], Directive[cols[[2]], PointSize[0.008]], Directive[cols[[1]], PointSize[0.008]]},
   FrameTicks -> {{Automatic, None}, {Table[{k, Rotate[Style[peNames[families[[k]]], 8], 30 Degree]}, {k, 8}], None}},
   FrameLabel -> {"", "log10 relative error of T"}, PlotRange -> {{0.5, 8.5}, {-17, -4}},
   PlotLegends -> ptleg[cols[[{6, 2, 1}]], {"Python v0.1.0 (sin substitution)", "Python, factored integrand", "Wolfram, factored integrand"}, Scaled[{0.62, 0.87}]],
   ImagePadding -> {{52, 12}, {70, 10}}, Evaluate[frameOpts]];
sepData = Table[Module[{m = NOLModel[-1, 1], dE = 10^-k},
    {k, NOLPeriodFromEnergy[m, -dE, 1, WorkingPrecision -> 40, "Tolerance" -> 10^-40] - Log[1/dE],
     NOLPeriodFromEnergy[m, dE, 0, WorkingPrecision -> 40, "Tolerance" -> 10^-40] - 2 Log[1/dE]}], {k, Range[1, 16, 1/2]}];
log["separatrix data ", N[sepData]];
sepData = N[sepData];
sepIn = Select[Transpose[{10.^-sepData[[All, 1]], sepData[[All, 2]]}], VectorQ[#, NumericQ] && Im[#[[2]]] == 0 &];
sepOut = Select[Transpose[{10.^-sepData[[All, 1]], sepData[[All, 3]]}], VectorQ[#, NumericQ] && Im[#[[2]]] == 0 &];
sepIn = Select[Transpose[{-sepData[[All, 1]], sepData[[All, 2]]}], VectorQ[#, NumericQ] && Im[#[[2]]] == 0 &];
sepOut = Select[Transpose[{-sepData[[All, 1]], sepData[[All, 3]]}], VectorQ[#, NumericQ] && Im[#[[2]]] == 0 &];
f6b = Show[
   ListPlot[{sepIn, sepOut}, Joined -> True, PlotMarkers -> {{"\[FilledCircle]", 6}, {"\[FilledSquare]", 6}}, PlotStyle -> cols[[{3, 4}]]],
   Plot[{Log[16], 2 Log[16]}, {x, -16.5, -0.5}, PlotStyle -> Directive[Gray, Dashed]],
   FrameLabel -> {"log10 |E - E_s|", "T - k ln(1/|E - E_s|)"}, PlotRange -> {{-16.5, -0.5}, {2, 6.2}},
   Epilog -> {Text[Style["inner orbit, k=1  \[RightArrow] ln 16", 10, cols[[3]], FontFamily -> "Times New Roman"], Scaled[{0.55, 0.3}]],
     Text[Style["outer orbit, k=2  \[RightArrow] 2 ln 16", 10, cols[[4]], FontFamily -> "Times New Roman"], Scaled[{0.55, 0.72}]]},
   Evaluate[frameOpts]];
save["fig6_period_accuracy", row[lbl[f6a, "(a) relative period error (8 families \[Times] 40 energies)"],
   lbl[f6b, "(b) separatrix asymptotics, double well"]]];

(* ---------------- Fig. 7: arbitrary precision and integrator tolerance ---------------- *)
tHi = NOLPeriod[NOLModel[1, 2/5, 1/5], 4/5, 1/10, WorkingPrecision -> 80];
(* correct digits: compare the digits returned at precision wp (taken at face value) with an 80-digit value *)
precData = Table[Module[{tm, tw},
     {tm, tw} = AbsoluteTiming[NOLPeriod[NOLModel[1, 2/5, 1/5], 4/5, 1/10, WorkingPrecision -> wp]];
     {wp, N[-Log10[Max[Abs[SetPrecision[tw, 120] - SetPrecision[tHi, 120]], 10^-90]]], tm}], {wp, 16, 72, 4}];
log["precision study ", precData];
f7a = Show[ListPlot[precData[[All, {1, 2}]], Joined -> True, PlotMarkers -> {{"\[FilledCircle]", 7}}, PlotStyle -> cols[[1]]],
   Plot[x, {x, 14, 74}, PlotStyle -> Directive[Gray, Dashed]],
   FrameLabel -> {"WorkingPrecision (digits)", "correct digits of T"}, PlotRange -> {{14, 74}, {0, 80}},
   Epilog -> Text[Style["dashed: digits = WorkingPrecision", 10, Gray, FontFamily -> "Times New Roman"], Scaled[{0.62, 0.12}]],
   Evaluate[frameOpts]];
tolData = Table[{pg, Module[{s = NOLSimulate[NOLModel @@ #[[1]], #[[2]], #[[3]], 200, PrecisionGoal -> pg, AccuracyGoal -> pg + 2]},
       NOLEnergyDrift[NOLModel @@ #[[1]], s, "Samples" -> 4001]]} & /@ {{{1, 0.4, 0.2}, 0.8, 0.1}, {{-1, 1, 0}, 1.6, 0.}, {{1, -4, 3}, 1.2, 0.3}},
   {pg, 6, 14, 1}];
f7b = ListLogPlot[Transpose[tolData], Joined -> True,
   PlotMarkers -> Automatic, PlotStyle -> cols[[{1, 4, 5}]],
   FrameLabel -> {"PrecisionGoal", "relative energy drift on [0, 200]"},
   PlotLegends -> leg[cols[[{1, 4, 5}]], {"a=1, b=0.4, c=0.2", "double well, outer", "triple well, outer"}, Scaled[{0.7, 0.8}]],
   Evaluate[frameOpts]];
save["fig7_precision_and_tolerance", row[lbl[f7a, "(a) arbitrary-precision period, cubic\[Dash]quintic"], lbl[f7b, "(b) RK8 energy drift vs tolerance"]]];

(* ---------------- Table: cross-implementation benchmark ---------------- *)
bench = Map[num, Rest@Import[FileNameJoin[{$root, "benchmarks", "python_periods.csv"}], "CSV"], {2}];
exactBench = <|"harmonic" -> {4, 0, 0, 2/5, 0}, "hardening" -> {1, 1, 0, 1, 0}, "softening" -> {1, -1/2, 0, 9/10, 0},
   "dw_inner" -> {-1, 1, 0, 6/5, 0}, "dw_outer" -> {-1, 1, 0, 8/5, 0}, "cq_demo" -> {1, 2/5, 1/5, 4/5, 1/10},
   "tw_center" -> {1, -4, 3, 3/10, 0}, "tw_side" -> {1, -4, 3, 11/10, 0}, "tw_outer" -> {1, -4, 3, 6/5, 3/10}|>;
fmtE[x_] := If[x == 0, "0", Module[{e = Floor[Log10[N[x]]]}, ToString[NumberForm[N[x/10^e], {3, 1}]] <> "e" <> ToString[e]]];
texE[x_] := If[x == 0, "0", Module[{e = Floor[Log10[N[x]]]}, "$" <> ToString[NumberForm[N[x/10^e], {3, 1}]] <> "\\times10^{" <> ToString[e] <> "}$"]];
benchRows = Table[Module[{id = r[[1]], e, t40, twl, cf, cfErr},
    e = exactBench[id];
    t40 = NOLPeriod[NOLModel @@ e[[1 ;; 3]], e[[4]], e[[5]], WorkingPrecision -> 40];
    twl = NOLPeriod[NOLModel @@ N[e[[1 ;; 3]]], N[e[[4]]], N[e[[5]]]];
    cf = Quiet[NOLExactSolution @@ N[e]];
    cfErr = If[FailureQ[cf], Missing[], Quiet@NOLExactSolutionError[NOLModel @@ N[e[[1 ;; 3]]], N[e[[4]]], N[e[[5]]], cf, 20]];
    <|"case" -> id, "a" -> e[[1]], "b" -> e[[2]], "c" -> e[[3]], "x0" -> e[[4]], "v0" -> e[[5]],
      "T40" -> ToString[NumberForm[N[t40, 20], 16]],
      "errPy" -> N[Abs[r[[7]] - t40]/t40], "errWL" -> N[Abs[twl - t40]/t40], "errCF" -> cfErr|>], {r, bench}];
Export[FileNameJoin[{$tab, "benchmark_periods.csv"}],
  Prepend[Table[{b["case"], ToString[b["a"], InputForm], ToString[b["b"], InputForm], ToString[b["c"], InputForm],
     ToString[b["x0"], InputForm], ToString[b["v0"], InputForm], b["T40"], fmtE[b["errPy"]], fmtE[b["errWL"]],
     If[MissingQ[b["errCF"]], "n/a", fmtE[b["errCF"]]]}, {b, benchRows}],
   {"case", "a", "b", "c", "x0", "v0", "T_40digit", "relerr_python", "relerr_wolfram", "closed_form_maxerr_0_20"}], "CSV"];
Export[FileNameJoin[{$tab, "benchmark_periods.tex"}],
  StringRiffle[Table[StringRiffle[{StringReplace[b["case"], "_" -> "\\_"],
       "$" <> ToString[TeXForm[b["a"]]] <> "$", "$" <> ToString[TeXForm[b["b"]]] <> "$", "$" <> ToString[TeXForm[b["c"]]] <> "$",
       "$" <> ToString[TeXForm[b["x0"]]] <> "$", "$" <> ToString[TeXForm[b["v0"]]] <> "$", b["T40"],
       texE[b["errPy"]], texE[b["errWL"]], If[MissingQ[b["errCF"]], "n/a", texE[b["errCF"]]]}, " & "] <> " \\\\",
     {b, benchRows}], "\n"], "Text"];
log["tables written"];
Export[FileNameJoin[{$tab, "environment.txt"}],
  StringJoin["Generated by MakeFigures.wls with ", $Version, "\nPackage version ", NOLVersion, "\nDate ", DateString[], "\n"], "Text"];
(* ---------------- every number quoted in the manuscript ---------------- *)
texM[x_] := If[x == 0, "0", Module[{e = Floor[Log10[N[x]]]},
    "\\ensuremath{" <> ToString[NumberForm[N[x/10^e], {2, 1}]] <> "\\times10^{" <> ToString[e] <> "}}"]];
errLegacy = Flatten@Table[MapThread[relE, {legacyByF[f], refT[f]}], {f, families}];
errPy = Flatten@Table[MapThread[relE, {famData[f][[All, 7]], refT[f]}], {f, families}];
errWL = Flatten@Table[MapThread[relE, {wlT[f], refT[f]}], {f, families}];
Export[FileNameJoin[{$tab, "period_accuracy.csv"}],
  Prepend[Flatten[Table[Transpose[{ConstantArray[f, Length[famData[f]]], famData[f][[All, 6]], ToString[N[#, 20], InputForm] & /@ refT[f],
       MapThread[relE, {legacyByF[f], refT[f]}], MapThread[relE, {famData[f][[All, 7]], refT[f]}], MapThread[relE, {wlT[f], refT[f]}]}], {f, families}], 1],
   {"family", "E", "T_reference", "relerr_python_v010_sin_substitution", "relerr_python_factored", "relerr_wolfram_factored"}], "CSV"];
testRows = Rest@Import[FileNameJoin[{$dir, "test_report.csv"}], "CSV"];
repParts = StringSplit[exDemo["Representation"], {" / y: ", ":"}];
numbers = {
   "wlTests" -> ToString[Length[testRows]], "wlPassed" -> ToString[Count[testRows[[All, 2]], "Success"]],
   "pythonTests" -> "15",
   "pythonLegacyMaxErr" -> texM[Max[errLegacy]], "pythonLegacyMedErr" -> texM[Median[errLegacy]],
   "pythonMaxErr" -> texM[Max[errPy]], "pythonMedErr" -> texM[Median[errPy]],
   "wlMaxErr" -> texM[Max[errWL]], "wlMedErr" -> texM[Median[errWL]],
   "campN" -> ToString[nCampaign], "campOK" -> ToString[campaignSummary["represented"]],
   "campMed" -> texM[campaignSummary["median_error"]], "campPninetyfive" -> texM[campaignSummary["p95_error"]],
   "campMax" -> texM[campaignSummary["max_error"]],
   "campBelow" -> ToString[Round[100 campaignSummary["fraction_below_1e-8"]]] <> "\\%",
   "campTime" -> ToString[NumberForm[1000 campaignSummary["median_time_s"], {3, 0}]] <> "~ms",
   "demoRep" -> "family~" <> repParts[[1]] <> ", \\texttt{" <> StringTrim[repParts[[2]]] <> "} branch",
   "mmaVersion" -> StringTrim[First@StringSplit[$Version, "("]]};
Export[FileNameJoin[{$tab, "numbers.tex"}],
  StringRiffle[("\\newcommand{\\" <> #[[1]] <> "}{" <> #[[2]] <> "}") & /@ numbers, "\n"] <> "\n", "Text"];
log["numbers ", numbers];
log["ALL DONE"];
