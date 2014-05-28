% Script runs UC for SDSS

% Initialization
ucParams.numALCandidates = numALCandidates;
ucParams.lowestLogliklVal = lowestLogliklVal;
ucParams.alBandwidth = alBandwidth;
ucParams.alScale = alScale;
ucParams.gpNoiseLevel = noiseLevelGP;
[ucPts, ucLogProbs] = alGPUncertaintyReduction(evalLogJoint, alInitPts, ...
  alInitLogProbs, paramSpaceBounds, NUM_AL_ITERS, ucParams);

for ucResIter = 1:numResultsToBeStored

  currNumUCPts = ucResIter * STORE_RESULTS_EVERY;
  % Obtain regressors and regressands
  Xtr = ucPts(1:currNumUCPts, :);
  Ytr = ucLogProbs(1:currNumUCPts);

  % Obtain the errors and estimtes
  [currKL, ucLogJointEst, ucProbEst]= evalRegMethodKLProgress( ...
    Xtr, Ytr, ...
    gpFitParams, klEvalPts, truePAtEvalPts, evalMCMCParams, optKDEBandWidth);
  prs = probRatioStatistic(prsPts, prsLogProbs, ucLogJointEst, MLEPoint, ...
        MLELogP);

  % record and report
  uc_errs(experimentIter, ucResIter) = currKL;
  fprintf(' UC Iter: %d, #pts: %d, KL: %0.4f, prs: %0.4f\n', ucResIter, ...
          currNumUCPts, currKL, prs);

end

