% MaxBanPoint for SNLS

initLipschitzConst = 12000; % I think this is an under estimate

% Set params for alMaxBandPoint
phi = @exp; gradPhi = @exp;
alMbpParams.num_grad_desc_init_pts = 300;
alMbpParams.init_step_size = 0.1;
[mbpPts, mbpLogProbs, mbpLipConst] = alMaxBandPoint(evalLogJoint, [], [], ...
  phi, gradPhi, initLipschitzConst, paramSpaceBounds, NUM_AL_ITERS, ...
  alMbpParams);

for mbpResIter = 1:numResultsToBeStored

  currNumMBPpts = mbpResIter * STORE_RESULTS_EVERY;
  % Obtain Regressors and Regressands
  Xtr = ucPts(1:currNumMBPpts, :);
  Ytr = ucLogProbs(1:currNumMBPpts);

  % Regress on these points
  mbpLogJointEst = regressionWrap(Xtr, Ytr, noiseLevelGP, lowestLogLiklVal, ...
    logLiklRange, cvCostFunc);
  % obtain the estimted posterior
  [postEst, mbpEstLogPostAtEvalPts] = ...
    obtainProbHandle(mbpLogJointEst, evalPts, paramSpaceBounds);

  % Evaluate the KL
  currKL = eval3DKLFromLogProbs(trueLogPostValsAtEvalPts, ...
    mbpEstLogPostAtEvalPts, paramSpaceBounds);
  mbp_errs(experimentIter, mbpResIter) = currKL;
  fprintf(' MBP Iter: %d, num-points: %d, KL: %0.4f\n', mbpResIter, ...
    currNumMBPpts, currKL);
end
