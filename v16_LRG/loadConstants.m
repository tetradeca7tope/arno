% Sets all parameters for the problem

% First determine the tests you are going to perfrm
DOING_MMD = true;
DOING_LOGIT = false;

numDims = 8;
lowestLogliklVal = -2000;
logLiklRange = 2000;

% Set up problem class
paramSpaceBounds = repmat([0 1], numDims, 1);
lrgExp = LRGExperiment(lowestLogliklVal);
evalLogJoint = @(arg) lrgExp.normCoordLogJointProbs(arg);
MLEPoint = lrgExp.getNormCoords([0 0.762 0.1045 0.02233 0.951 0.6845 0 1.908]);
MLELogP = evalLogJoint(MLEPoint);

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
% alInitPts=lrgExp.getNormCoords([0 0.762 0.1045 0.02233 0.951 0.6845 0 1.908]);
% alInitLogProbs = evalLogJoint(alInitPts);
alInitPts = [];
alInitLogProbs = [];
fprintf('Initing with logp = %s at %s MLE-LOGP = %0.5f\n', ...
  mat2str(alInitLogProbs), mat2str(alInitPts), MLELogP );

% Constants for MaxBandPoint
initLipschitzConstant = logLiklRange/ (sqrt(numDims) * 2);
initLipschitzConstant = 500; %logLiklRange/ (sqrt(numDims) * 2);

% Parameters for Uncertainty Reduction
% the KL was stuck at ~ 3-5 here for all iterations.
alBandwidth = 0.5 * NUM_AL_ITERS ^ (-1 /(1.3 + numDims));
% So I am going to try reducing it to see if we can have more exploration
% alBandwidth = 5 * NUM_AL_ITERS ^ (-1 /(1.3 + numDims));
% the constant 2.5 works well in front
alScale = logLiklRange/2;

% Log on alBandwidth's that worked

% Parameters for MCMC
NUM_MCMC_SAMPLES = 4 * NUM_AL_ITERS;
numMCMCResultsToBeStored = NUM_MCMC_SAMPLES / STORE_RESULTS_EVERY;
if DOING_LOGIT
  mcmcProposalStd = 2; % after the logit transform
  mcmcInitPt = zeros(numDims, 1);
  mcmcLogJoint = @(t) evalLogJoint(logitinv(t));
else
  mcmcProposalStd = 0.002;
  mcmcInitPt = 0.5*ones(numDims, 1);
  mcmcLogJoint = @(t) evalLogJoint(t);
end

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
numSamplesEstEval = 1e5;
numBurninSampleEstEval = numSamplesEstEval/10;
evalMCMCProposalStd = mcmcProposalStd;
evalMCMCLogJoint = mcmcLogJoint;
if DOING_LOGIT
  evalMCMCInitPt = logit(MLEPoint )';
else
  evalMCMCInitPt = MLEPoint';
end
gtFile = sprintf('groundTruth_logit%d_%d_%.3f.mat', DOING_LOGIT, ...
  numSamplesEstEval, evalMCMCProposalStd);
% Create this structure for convenience
evalMCMCParams.evalMCMCProposalStd = evalMCMCProposalStd;
evalMCMCParams.evalMCMCInitPt = evalMCMCInitPt;
evalMCMCParams.numBurninSampleEstEval = numBurninSampleEstEval;
evalMCMCParams.numSamplesEstEval = numSamplesEstEval;
evalMCMCParams.doingLogit = DOING_LOGIT;
% Use another file to store all the ground truth details - just in case
allGtFile = sprintf('allGroundTruth_%d.mat', numSamplesEstEval);

