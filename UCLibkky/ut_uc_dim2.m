% A Unit test in 2 dimensions for mbpObjective.m, mbpGradient.m, maxBandPoint.m
% and alMaxBandPoint.m

close all;
clear all;
addpath ../LipLibkky/

% For the experiments
resolution = 100;

% Test 2 : Two Dimensions
% -----------------------
fprintf('Test in 2 Dimensions\n');

% specify the function 
f = @(T) diag(T*T') /4 + sin(2*T(:,1) + 2*T(:, 2));
% Plot the function
t1 = linspace(0,1,resolution); [T1 T2] = meshgrid(t1, t1); t = [T1(:), T2(:)];
ft = f(t); fT = reshape(ft, resolution, resolution);
figure;
mesh(T1, T2, fT);

% test for Uncertainty Reduction
% ==============================
numALIters = 300;
ucParams.numALCandidates = 40;
ucParams.lowestLogliklVal = -1;
ucParams.alBandwidth = 0.01;
ucParams.logLiklRange = 5;
bounds = [0 1; 0 1];

[ucPts, ucVals] = alGPUncertaintyReduction(f, [], [], bounds, numALIters, ...
  ucParams);
figure;
contour(T1, T2, fT); hold on,
plot(ucPts(:,1), ucPts(:,2), 'rx');

