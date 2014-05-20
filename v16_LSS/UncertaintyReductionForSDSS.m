% Script runs UC for SDSS

% Initialization
ucParams.numALCandidates = numALCandidates;
ucParams.lowestLogliklVal = lowestLogliklVal;
ucParams.alBandwidth = alBandwidth;
ucParams.alScale = alScale;
ucParams.gpNoiseLevel = noiseLevelGP;
[ucPts, ucLogProbs] = alGPUncertaintyReduction(evalLogJoint, [], [], ...
  paramSpaceBounds, NUM_AL_ITERS, ucParams);

% Now do the regression
ucLogJointEst = regressionWrap(ucPts, ucLogProbs, noiseLevelGP, ...
  lowestLogliklVal, logLiklRange, cvCostFunc);

% Save the Results
save(UC_EXP_FILE, 'ucPts', 'ucLogProbs', 'ucParams', 'noiseLevelGP', ...
  'lowestLogliklVal', 'logLiklRange', 'cvCostFunc');

