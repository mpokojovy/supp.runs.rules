%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on the original article                                                                       %
%               Detecting Selected Non-Random Patterns with Individuals Control Charts                %
% by J. Marcus Jobe & Michael Pokojovy                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Copyright (C): Michael Pokojovy and J. Marcus Jobe (2020)
See license.txt for license conditions

This R code implements Western Electric runs rules, Nelson's runs rules, Shewhart X and CUSUM charts.

Instructions:

A) In-control (IC) simulation:
  Run IC.sim/shell.R

  Warning: Use HPC computing or reduce the value of nrep (to, say, 10000) in IC.sim/run.Nelson.R and 
           IC.sim/run.other.R

B) In-control (OC) simulation:
  Run IC.sim/shell.R

B) Evaluate out-of-control performance of the TS, CUSUM, CUSUM head start and Shewhart X charts
  Run evaluate_chart_performance.m (Note that the program will use precomputed control limits, etc.)

C) Run ds.plots/ds.3plots.init.state.R and ds.plots/ds.3plots.steady.state.R to display sample OC 
   streams for the initial state (steady state size = 0) and steady state (steady state size = 25).



