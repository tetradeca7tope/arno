% Sets all parameters for the problem

numDims = 8;
lowestLogliklVal = -200;
logLiklRange = 400;

% Set up problem class
paramSpaceBounds = repmat([0 1], numDims, 1);
lrgExp = LRGExperiment();
evalLogJoint = @(arg) lrgExp.normCoordLogJointProbs(arg);

DEBUG_MODE = false;
% DEBUG_MODE = true;
if ~DEBUG_MODE
  NUM_AL_ITERS = 2500;
  NUM_EXPERIMENTS = 30;
  STORE_RESULTS_EVERY = NUM_AL_ITERS/10;
else
  NUM_AL_ITERS = 6;
  NUM_EXPERIMENTS = 2; 
  STORE_RESULTS_EVERY = 3;
end
numResultsToBeStored = NUM_AL_ITERS / STORE_RESULTS_EVERY;

% Params for fitting the GP
noiseLevelGP = logLiklRange/100;
cvCostFunc = @(y1, y2) (exp(y1) - exp(y2)).^2;
% Create this structure for convenience
gpFitParams.noiseLevelGP = noiseLevelGP;
gpFitParams.cvCostFunc = cvCostFunc;
gpFitParams.lowestLogliklVal = lowestLogliklVal;
gpFitParams.logLiklRange = logLiklRange;

% Parameters for active learning
numALCandidates = 2000;
alInitPts = lrgExp.getNormCoords([0 0.762 0.1045 0.02233 0.951 0.6845 0 1.908]);
alInitLogProbs = evalLogJoint(alInitPts);
fprintf('Initing with logp = %s at %s\n', ...
  mat2str(alInitLogProbs), mat2str(alInitPts) );

% Constants for MaxBandPoint
initLipschitzConstant = logLiklRange/ (sqrt(numDims) * 2);
initLipschitzConstant = 500; %logLiklRange/ (sqrt(numDims) * 2);

% Parameters for Uncertainty Reduction
alBandwidth = 1.5 * NUM_AL_ITERS ^ (-1 /(1.3 + numDims));
% the constant 2.5 works well in front
alScale = logLiklRange/2;

% Parameters for MCMC
NUM_MCMC_SAMPLES = 6 * NUM_AL_ITERS;
mcmcProposalStd = 0.05; % after the logit transform
mcmcInitPt = logit( ...
  lrgExp.getNormCoords([0 0.762 0.1045 0.02233 0.951 0.6845 0 1.908])  )';
numMCMCResultsToBeStored = NUM_MCMC_SAMPLES / STORE_RESULTS_EVERY;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other Ancillary Variables

% For saving the results
timestamp = datestr(now, 'mm-dd_HH-MM-SS');
saveFileName = sprintf('lrgResults/lrg_%s.mat', timestamp);
uc_errs = zeros(NUM_EXPERIMENTS, numResultsToBeStored);
mbp_errs = zeros(NUM_EXPERIMENTS, numResultsToBeStored);
mcmcReg_errs = zeros(NUM_EXPERIMENTS, numResultsToBeStored);
rand_errs = zeros(NUM_EXPERIMENTS, numResultsToBeStored);
% Larger arrays for MCM
mcmc_errs = zeros(NUM_EXPERIMENTS, numMCMCResultsToBeStored);

% For evaluation purposes
numBurninSampleEstEval = 1000;
numSamplesEstEval = 1e5;
evalMCMCProposalStd = mcmcProposalStd; % No reason to use anything else
evalMCMCInitPt = logit( ...
  lrgExp.getNormCoords([0 0.762 0.1045 0.02233 0.951 0.6845 0 1.908])  )';
evalMCMCLogJoint = @(t) evalLogJoint(logitinv(t));
gtFile = sprintf('groundTruthKL_%d_%.3f.mat', numSamplesEstEval, ...
                 evalMCMCProposalStd);
% Create this structure for convenience
evalMCMCParams.evalMCMCProposalStd = evalMCMCProposalStd;
evalMCMCParams.evalMCMCInitPt = evalMCMCInitPt;
evalMCMCParams.numBurninSampleEstEval = numBurninSampleEstEval;
evalMCMCParams.numSamplesEstEval = numSamplesEstEval;
% Use another file to store all the ground truth details - just in case
allGtFile = sprintf('allGroundTruth_%d.mat', numSamplesEstEval);

