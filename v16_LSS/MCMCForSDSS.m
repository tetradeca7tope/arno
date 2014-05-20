% Performs MCMC for SDSS

mcmcEvalLogJoint = @(t) evalLogJoint(logitinv(t));
[logitMcmcSamples, logitMcmcQueries, mcmcLogProbs] = CustomMCMC( ...
  NUM_MCMC_SAMPLES, mcmcProposalStd, mcmcInitPt, mcmcEvalLogJoint);
mcmcSamples = logitinv(logitMcmcSamples);
mcmcQueries = logitinv(logitMcmcQueries);

% Perform MCMC Density Estimation
[~, mcmcProbEst] = kde01(mcmcSamples);

% Now perform MCMC Regression
mrLogJointEst = regressionWrap(mcmcQueries, mcmcLogProbs, noiseLevelGP, ...
  lowestLogliklVal, logLiklRange, cvCostFunc);

% Save results
save(MCMC_EXP_FILE, 'mcmcSamples', 'mcmcQueries', 'mcmcLogProbs', ...
  'noiseLevelGP', 'lowestLogliklVal', 'logLiklRange', 'cvCostFunc');

