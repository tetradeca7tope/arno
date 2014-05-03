% A Unit test in 1 dimension for mbpObjective.m, mbpGradient.m, maxBandPoint.m
% and alMaxBandPoint.m

close all;
clear all;
addpath ../UCLibkky/
addpath ../GPLibkky/
addpath ../helper/

% For the experiments
num_pts = 10;
resolution = 100;

% Test 1: One Dimension
% ---------------------
fprintf('Test in 1 Dimension:\n');
f = @(t) t.^2 /4 + sin(7*t);

% TEST for alMaxBandPoint
% =======================
numALIters = 100;
ucParams.numALCandidates = 10;
ucParams.lowestLogliklVal = -1;
ucParams.alBandwidth = 0.35;
ucParams.logLiklRange = 5;
bounds = [0 1];

[ucPts, ucVals] = alGPUncertaintyReduction(f, [], [], bounds, numALIters, ...
  ucParams);

figure;
th = linspace(0, 1, 100)';
plot(th, f(th), 'b'); hold on,
plot(ucPts, ucVals, 'rx');

