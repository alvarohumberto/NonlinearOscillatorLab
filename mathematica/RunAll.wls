#!/usr/bin/env wolframscript
(* Runs the test suite and regenerates every figure and table.
   Usage (from any folder):  wolframscript -file RunAll.wls *)
$dir = DirectoryName[$InputFileName];
Print["== NonlinearOscillatorLab: tests =="];
Get[FileNameJoin[{$dir, "RunTests.wls"}]];
Print["== NonlinearOscillatorLab: figures and tables =="];
Get[FileNameJoin[{$dir, "MakeFigures.wls"}]];
Print["== finished =="];
