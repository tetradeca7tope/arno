% Sets all parameters for the problem

% First determine the tests you are going to perfrm
DOING_LOGIT = false;

numDims = 8;
lowestLogliklVal = -200;
logLiklRange = 200;

% Set up problem class
paramSpaceBounds = repmat([0 1], numDims, 1);
lrgExp = LRGExperiment(lowestLogliklVal);
evalLogJoint = @(arg) lrgExp.normCoordLogJointProbs(arg);
MLEPoint = lrgExp.getNormCoords([0 0.762 0.1045 0.02233 0.951 0.6845 0 1.908]);
MLELogP = evalLogJoint(MLEPoint);

DEBUG_MODE = false;
% DEBUG_MODE = true;
if ~DEBUG_MODE
  NUM_AL_ITERS = 12000;
  NUM_EXPERIMENTS = 10;
%   STORE_RESULTS_AT = round( logspace(0, log10(NUM_AL_ITERS), 15) )';
  STORE_RESULTS_AT = [10 50 100 500 1000 2000 4000 6000 8000 10000 12000]';
else
  NUM_AL_ITERS = 6;
  NUM_EXPERIMENTS = 2; 
  STORE_RESULTS_AT = [2 4 6]';
end
numResultsToBeStored = numel(STORE_RESULTS_AT);

% Params for fitting the GP
noiseLevelGP = logLiklRange/100;
cvCostFunc = @(y1, y2) (exp(y1) - exp(y2)).^2;
% Create this structure for convenience
gpFitParams.noiseLevelGP = noiseLevelGP;
gpFitParams.cvCostFunc = cvCostFunc;
gpFitParams.lowestLogliklVal = lowestLogliklVal;
gpFitParams.logLiklRange = logLiklRange;

% Parameters for active learning
numALCandidates = 4000;
% alInitPts = lrgExp.getNormCoords([0 0.762 0.1045 0.02233 0.951 0.6845 0 1.908]);
alInitPts = 0.5 * ones(1, numDims); 
alInitLogProbs = evalLogJoint(alInitPts);
fprintf('Initing with logp = %s at %s MLE-LOGP = %0.5f\n', ...
  mat2str(alInitLogProbs), mat2str(alInitPts), MLELogP );

% Constants for MaxBandPoint
initLipschitzConstant = logLiklRange/ (sqrt(numDims) * 2);
initLipschitzConstant = 500; %logLiklRange/ (sqrt(numDims) * 2);

% Parameters for Uncertainty Reduction
% the KL was stuck at ~ 3-5 here for all iterations.
alBandwidth = 0.5 * NUM_AL_ITERS ^ (-1 /(1.3 + numDims));
alScale = 4*logLiklRange;

% Parameters for MCMC
mcmcMultipleFactor = 2;
NUM_MCMC_SAMPLES = mcmcMultipleFactor * NUM_AL_ITERS;
if ~DEBUG_MODE
  STORE_MCMC_RESULTS_AT = mcmcMultipleFactor * STORE_RESULTS_AT;
else
  STORE_MCMC_RESULTS_AT = (2:2:NUM_MCMC_SAMPLES)';
end
numMCMCResultsToBeStored = numel(STORE_MCMC_RESULTS_AT);
if DOING_LOGIT
  mcmcProposalStd = 2; % after the logit transform
  mcmcInitPt = zeros(numDims, 1);
  mcmcLogJoint = @(t) evalLogJoint(logitinv(t));
else
  mcmcProposalStd = 0.05 * sqrt(numDims);
  mcmcInitPt = rand(numDims, 1);
  mcmcLogJoint = @(t) evalLogJoint(t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other Ancillary Variables

% For saving the results
timestamp = datestr(now, 'mmddHHMMSS');
saveFileName = sprintf('lrgResults/lrgnew_%s.mat', timestamp);
savePtsFilePrefix = sprintf('lrgResults/lrgPts_%s', timestamp);
uc_errs = zeros(NUM_EXPERIMENTS, numResultsToBeStored);
mbp_errs = zeros(NUM_EXPERIMENTS, numResultsToBeStored);
mcmcReg_errs = zeros(NUM_EXPERIMENTS, numResultsToBeStored);
rand_errs = zeros(NUM_EXPERIMENTS, numResultsToBeStored);
% Larger arrays for MCM
mcmc_errs = zeros(NUM_EXPERIMENTS, numMCMCResultsToBeStored);

% For evaluation purposes
randTestFiles = {'gt250K.mat'};
gtPts = zeros(0, numDims); gtLogProbs = zeros(0, 1);
for i = 1:numel(randTestFiles)
  gt = load(randTestFiles{i});
  gtPts = [gtPts; gt.randPts]; gtLogProbs = [gtLogProbs; gt.randLogProbs];
end
numGTPts = size(gtPts, 1);
fprintf('Total Num Eval Pts: %d\n', numGTPts);
fprintf('Evaluating at %s\n', mat2str(STORE_RESULTS_AT));

