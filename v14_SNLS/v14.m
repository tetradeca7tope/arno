% Script for running experiments on the SNLS dataset

% Prelims
clear all;
close all;
rng('shuffle'); % shuffle the seed for the RNG
addpath ../helper/
addpath ../GPLibkky/
addpath ../UCLibkky/
addpath ../LipLibkky/
addpath ../MCMCLibkky/
addpath ../ABCLibkky/

% Constants
numDims = 3;
lowestLogLiklVal = - 2000;
logLiklRange = 2000;
paramSpaceBounds = [0 1; 0 1; 0 1]; % we are normalizing them later on
% Constants for the tests
alphaL = 0.05; % Testing at 95% confidence
partitionRes = 100; % For integrating and estimating the KL etc.

% First set the SNLS Experiment up
% ==============================================================================
% These are the points at which we'll be evaluating the KL for comparison
[g1, g2, g3] = gen3DCentrePointGrid(partitionRes, paramSpaceBounds);
evalPts = [g1(:), g2(:), g3(:)];
% Create function handle for Evaluating Joing Log Likelihood and post prob
fprintf('First obtaining True Posterior\n');
snlsParamSpaceBounds = [60 80; 0 1; 0 1];
snls = SNExperiment('davisdata.txt', snlsParamSpaceBounds);
evalLogJoint = @(arg) snls.normCoordLogJointProbs(arg);
truePostFile = sprintf('truePost_%d_%d_%d.txt', ...
  snlsParamSpaceBounds(1,1), snlsParamSpaceBounds(1,2),  partitionRes);
% Now evaluate the posterior at evalPts
    loadTruePostFromFile = false;
    % loadTruePostFromFile = true;
if ~loadTruePostFromFile
  [truePosterior, trueLogPostValsAtEvalPts, normConst] = ...
    obtainProbHandle(evalLogJoint, evalPts, paramSpaceBounds);
  save(truePostFile, 'trueLogPostValsAtEvalPts', '-ascii');
else
  trueLogPostValsAtEvalPts = load(truePostFile);
end
% ==============================================================================

% params for fitting the GP
noiseLevelGP = 5;
cvCostFunc = @(y1, y2) (exp(y1) - exp(y2)).^2;

DEBUG_MODE = false;
% DEBUG_MODE = true;
if ~DEBUG_MODE
  NUM_AL_ITERS = 600;
  NUM_MCMC_ITERS_FOR_EST = 100000;
  NUM_MCMC_BURNIN_FOR_EST = 100;
  NUM_EXPERIMENTS = 1;
  STORE_RESULTS_EVERY = NUM_AL_ITERS/5;
else
  NUM_AL_ITERS = 6;
  NUM_MCMC_ITERS_FOR_EST = 10;
  NUM_MCMC_BURNIN_FOR_EST = 10;
  NUM_EXPERIMENTS = 2;
  STORE_RESULTS_EVERY = 3;
end
numResultsToBeStored = NUM_AL_ITERS / STORE_RESULTS_EVERY;

% Parameters for Active Learning
numALCandidates = 1000;

% Constants for MCMC
NUM_MCMC_SAMPLES = 3*NUM_AL_ITERS;
numMCMCResultsToBeStored = NUM_MCMC_SAMPLES / STORE_RESULTS_EVERY;

% Constants for the tests

% For Saving results
[~, hostname] = system('hostname');
saveFileName = sprintf('snlsResults/snls_%s.mat', ...
                        datestr(now, 'mm:dd-HH:MM:SS'));
uc_errs = zeros(NUM_EXPERIMENTS, numResultsToBeStored);
mbp_errs = zeros(NUM_EXPERIMENTS, numResultsToBeStored);
mcmcReg_errs = zeros(NUM_EXPERIMENTS, numResultsToBeStored);
% Larger arrays for MCMC and ABC
mcmc_errs = zeros(NUM_EXPERIMENTS, numMCMCResultsToBeStored);
abc_errs = zeros(NUM_EXPERIMENTS, numMCMCResultsToBeStored);

% Now run the experiments
for experimentIter = 1:NUM_EXPERIMENTS

  fprintf('Host: %s, date/time: %s\n', hostname, datestr(now, 'mm:dd-HH:MM'));
  fprintf('EXPERIMENT: %d\n==============================\n\n', experimentIter);
  
  fprintf('UNCERTAINTY REDUCTION\n');
  UncertaintyReductionSNLS;

  fprintf('MAX-BAND-POINT\n');
  MaxBandPointSNLS;

  fprintf('MCMC\n');
  MCMCSNLS;

  % Save the results
  save(saveFileName, 'uc_errs', 'mbp_errs', 'mcmc_errs', 'mcmcReg_errs', ...
                     'abc_errs');

end

