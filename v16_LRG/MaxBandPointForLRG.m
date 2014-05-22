% Performs Max Band Point for SDSS

% Initialization
initialPts = [];
initialLogProbs = [];

% Set parameters for MBP
phi = @exp; gradPhi = @exp; % Use exponential Transformation
alMbpParams.num_iters = 10;
alMbpParams.init_step_size = 1;
% First obtain the points via MBP
[mbpPts, mbpLogProbs, mbpLipConst] = alMaxBandPoint(evalLogJoint, ...
  initialPts, initialLogProbs, phi, gradPhi, initLipschitzConstant, ...
  paramSpaceBounds, NUM_AL_ITERS, alMbpParams);

for mbpResIter = 1:numResultsToBeStored

  currNumMBPPts = mbpResIter * STORE_RESULTS_EVERY;
  % Obtain regressors and regressands
  Xtr = mbpPts(1:currNumMBPPts, :);
  Ytr = mbpLogProbs(1:currNumMBPPts);

  % Obtain errors and other estimates
  [currKL, mbpLogJointEst, mbpProbEst] = evalRegMethodKLProgress( Xtr, Ytr, ...
    gpFitParams, klEvalPts, truePAtEvalPts, evalMCMCParams, optKDEBandWidth);

  % record & report
  mbp_errs(experimentIter, mbpResIter) = currKL;
  fprintf(' MBP iter: %d, #pts: %d, KL: %.4f\n', mbpResIter, currNumMBPPts, ...
          currKL);

end

