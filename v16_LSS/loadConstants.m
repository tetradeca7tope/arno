% Sets all parameters for the problem

numDims = 6;
lowestLogliklVal = -200;
logLiklRange = 200;

% Set up problem class
paramSpaceBounds = repmat([0 1], numDims, 1);
sdss6 = SDSS6Experiment();
evalLogJoint = @(arg) sdss6.normCoordLogJointProbs(arg);

DEBUG_MODE = false;
% DEBUG_MODE = true;
if ~DEBUG_MODE
  NUM_AL_ITERS = 2500;
  NUM_EXPERIMENTS = 1;
else
  NUM_AL_ITERS = 2;
  NUM_EXPERIMENTS = 1; 
end

% Params for fitting the GP
noiseLevelGP = logLiklRange/100;
cvCostFunc = @(y1, y2) (exp(y1) - exp(y2)).^2;

% Parameters for active learning
numALCandidates = 2000;

% Constants for MaxBandPoint
initLipschitzConstant = logLiklRange/ (sqrt(numDims) * 2);
initLipschitzConstant = 600; %logLiklRange/ (sqrt(numDims) * 2);

% Parameters for Uncertainty Reduction
alBandwidth = 2.5 * NUM_AL_ITERS ^ (-1 /(1.5 + numDims));
alScale = logLiklRange/2;

% Parameters for MCMC
NUM_MCMC_SAMPLES = 10 * NUM_AL_ITERS;
mcmcProposalStd = 5; % after the logit transform
mcmcInitPt = 0.5 * ones(numDims, 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other Ancillary Variables

timestamp = datestr(now, 'mm-dd_HH-MM-SS');
UC_EXP_FILE = sprintf('results/ucexp_%s.mat', timestamp);
MBP_EXP_FILE = sprintf('results/mbpexp_%s.mat', timestamp);
MCMC_EXP_FILE = sprintf('results/mcmcexp_%s.mat', timestamp);
RAND_EXP_FILE = sprintf('results/rand_%s.mat', timestamp);

