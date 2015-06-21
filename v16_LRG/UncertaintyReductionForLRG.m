% Script runs UC for SDSS

% %Initialization
ucParams.numALCandidates = numALCandidates;
ucParams.lowestLogliklVal = lowestLogliklVal;
ucParams.alBandwidth = alBandwidth;
ucParams.alScale = alScale;
ucParams.gpNoiseLevel = noiseLevelGP;
[ucPts, ucLogProbs] = alGPUncertaintyReduction(evalLogJoint, alInitPts, ...
  alInitLogProbs, paramSpaceBounds, NUM_AL_ITERS, ucParams);

for ucResIter = 1:numResultsToBeStored

  currNumUCPts = STORE_RESULTS_AT(ucResIter);
  % Obtain regressors and regressands
  Xtr = ucPts(1:currNumUCPts, :);
  Ytr = ucLogProbs(1:currNumUCPts);

  % Obtain the errors and estimtes
  [currErr, ucLogJointEst]= evalRegMethodProgress(Xtr, Ytr, ...
    gtPts, gtLogProbs, gpFitParams);

  % record and report
  uc_errs(experimentIter, ucResIter) = currErr;
  fprintf(' UC Iter: %d, #pts: %d, Err: %e\n', ucResIter, ...
          currNumUCPts, currErr);

end

