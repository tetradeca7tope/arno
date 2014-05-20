% v16: SDSS

clear all;
close all;
rng('shuffle');
addpath ../GPLibkky/
addpath ../UCLibkky/
addpath ../LipLibkky/
addpath ../MCMCLibkky/
addpath ../ABCLibkky/
addpath ../helper/

% Load all constants
loadConstants;

for experimentIter = 1:NUM_EXPERIMENTS

  fprintf('Uncertainty Reduction\n');
  UncertaintyReductionForSDSS;

  fprintf('Max Band Point\n');
  MaxBandPointForSDSS;

  fprintf('MCMC\n');
  MCMCForSDSS;

  fprintf('RAND\n');
  RANDSDSS;

end

