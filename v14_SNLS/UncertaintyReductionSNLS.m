% Uncertainty Reduction for SNLS

% Params
alBandwidth = 0.5 * NUM_AL_ITERS ^ (-1/(1.3 + numDims) );
alScale = 2000;

% First obtain the points via Active Learning
numALIters = NUM_AL_ITERS;
ucParams.numALCandidates = numALCandidates;
ucParams.lowestLogliklVal = lowestLogLiklVal;
ucParams.alBandwidth = alBandwidth;
ucParams.alScale = alScale;
ucParams.gpNoiseLevel = noiseLevelGP;
[ucPts, ucLogProbs] = alGPUncertaintyReduction(evalLogJoint, [], [], ...
  paramSpaceBounds, numALIters, ucParams);

% Now Evaluate the progress
for ucResIter = 1:numResultsToBeStored

  currNumUCPts = ucResIter * STORE_RESULTS_EVERY;
  % Obtain the regressors and regressands
  Xtr = ucPts(1:currNumUCPts, :);
  Ytr = ucLogProbs(1:currNumUCPts);

  % Regress on these points
  ucLogJointEst = regressionWrap(Xtr, Ytr, noiseLevelGP, lowestLogLiklVal, ...
    logLiklRange, cvCostFunc);
  % and obtain a handle on the estimated posterior
  [postEst, ucEstLogPostValsAtEvalPts] = ...
    obtainProbHandle(ucLogJointEst, evalPts, paramSpaceBounds);
  
  % Evaluate the KL
  currKL = eval3DKLFromLogProbs(trueLogPostValsAtEvalPts, ...
    ucEstLogPostValsAtEvalPts, paramSpaceBounds);
  uc_errs(experimentIter, ucResIter) = currKL;
  fprintf(' UC Iter: %d, num-points: %d, KL: %0.4f\n', ucResIter, ...
    currNumUCPts, currKL);

end

