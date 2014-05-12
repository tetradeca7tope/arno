% v15: Galaxy Colouring

clear all;
close all;
rng('shuffle'); % Shuffle the seed for Matlab's RNG
addpath ../helper/
addpath ../GPLibkky/
addpath ../LipLibkky/
addpath ../UCLibkky/
addpath ../MCMCLibkky/

% Load all constants
loadConstants;

% First plot out the function
if numDims ==2
  figure;
  f = @(x) - (2*x(:,1).^2 + 3*x(:,2).^2) - 25;
  plot2DFunction(f, [problemSpaceBounds(1,:), ...
        problemSpaceBounds(2,:)], 'mesh');
end

for experiment_iter = 1:NUM_EXPERIMENTS

  fprintf('MaxBandPointForGalaxyColoring\n');
  MaxBandPointForGalaxyColoring;

  fprintf('Uncertainty Reduction\n');
  UncertaintyReductionGalaxyColoring;

  fprintf('MCMC \n');
  MCMCGalaxyColoring;

end

