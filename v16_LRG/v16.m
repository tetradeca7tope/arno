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

% Obtain an estimate of the ground truth
obtainGroundTruth = true;
% obtainGroundTruth = false;
if obtainGroundTruth
  tic,
  obtainGroundTruthKL;
  toc,
else 
  load(gtFile);
end

% For saving results


for experimentIter = 1:NUM_EXPERIMENTS

  fprintf('Uncertainty Reduction\n');
  UncertaintyReductionForLRG;

  % Before running MBP reduce numALCandidates to 1000.
  fprintf('Max Band Point\n');
  MaxBandPointForLRG;

  fprintf('MCMC\n');
  MCMCForLRG;

  fprintf('RAND\n');
  RANDForLRG;

end

