% A file to set all parameters for the problem

% Problem Dependent Constants
% numDims = 1; % problemSpaceBounds = [-1 1];
numDims = 2; problemSpaceBounds = [-1 1; -2 1];
% numDims = 3; problemSpaceBounds = [-1 1; -2 1; -1 0];
% Would be best if all bounds are in [0 1]. Lets work on this later.
lowestLogliklVal = -30; % TODO @Ying Need an estimate for this
logLiklRange = 8; % TODO @Ying: need an estimate for this.

% Set up problem class 
paramSpaceBounds = repmat([0 1], numDims, 1);
gcExp = GCExperiment(problemSpaceBounds, lowestLogliklVal);
evalLogJoint = @(arg) gcExp.normCoordLogJointProbs(arg);

DEBUG_MODE = false;
% DEBUG_MODE = true;
if ~DEBUG_MODE
  NUM_AL_ITERS = 50;
  NUM_EXPERIMENTS = 1;
else
  NUM_AL_ITERS = 5;
  NUM_EXPERIMENTS = 1;
end

% Parameters for fitting the GP
noiseLevelGP = logLiklRange / 100;
cvCostFunc = @(y1, y2) (exp(y1) - exp(y2)).^2;

% Parameters for ACtive Learning
numALCandidates = 1000;

% Constants for MaxBandPoint
INIT_LIPSCHITZ_CONST = 6;

% Parameters for Uncertainty Reduction
alBandwidth = 0.4 * NUM_AL_ITERS ^ (-1 /(1.3 + numDims));
alScale = logLiklRange / 2;

% Parameters for MCMC
NUM_MCMC_SAMPLES = 10 * NUM_AL_ITERS;
mcmcProposalStd = 0.3;
mcmcInitPt = 0.5*ones(numDims, 1);

