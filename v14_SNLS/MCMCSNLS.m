% Implements MCMC & MCMC Regression for the SNLS Experiment

mcmcProposalStd = 0.5;
mcmcInitPt = 0.5*ones(numDims, 1); % init at centre point
mcmcInitPt = [0.45; 0.24; 0.68]; % init at centre point

% Don't do Logit
%%%%%%%%%%%%%%%%
% % modify the log likelihood oracle to take the logit transform into account
% mcmcEvalLogJoint = @(t) evalLogJoint(logitinv(t));
% [logitMcmcSamples, logitMcmcQueries, mcmcLogProbs] = ...
%   CustomMCMC(NUM_MCMC_SAMPLES, mcmcProposalStd, mcmcInitPt, mcmcEvalLogJoint);
% mcmcSamples = logitinv(logitMcmcSamples);
% mcmcQueries = logitinv(logitMcmcQueries);
[mcmcSamples, mcmcQueries, mcmcLogProbs] = CustomMCMC(NUM_MCMC_SAMPLES, ...
  mcmcProposalStd, mcmcInitPt, evalLogJoint);

% Now perform MCMC
for mcmcResIter = 1:numMCMCResultsToBeStored

  currNumMcmcPts = mcmcResIter * STORE_RESULTS_EVERY;
  % Obtain the current set of points
  currMCMCSamples = mcmcSamples(1:currNumMcmcPts, :);
  
  % Perform KDE
  [~, mcmcProbEst] = kde(currMCMCSamples);
  mcmcLogProbEstAtEvalPts = log(mcmcProbEst(evalPts));

  % Evaluate the KL
  currKL = eval3DKLFromLogProbs(trueLogPostValsAtEvalPts, ...
    mcmcLogProbEstAtEvalPts, paramSpaceBounds);
  mcmc_errs(experimentIter, mcmcResIter) = currKL;
  fprintf(' MCMC Iter: %d, num-points: %d, KL: %0.4f\n', mcmcResIter, ...
    currNumMcmcPts, currKL);

end

% Now perform MCMC Regression
fprintf('MCMC Regression\n====================\n');
for mcmcRegResIter = 1:numResultsToBeStored

  currNumMRPts = mcmcRegResIter * STORE_RESULTS_EVERY;
  % Obtain Regressors and Regressands
  Xtr = mcmcQueries(1:currNumMRPts, :);
  Ytr = mcmcLogProbs(1:currNumMRPts);

  % Regress on these points
  mrLogJointEst = regressionWrap(Xtr, Ytr, noiseLevelGP, lowestLogLiklVal, ...
    logLiklRange, cvCostFunc);
  % Now obtain estimates at evalPts
  [mrPostEst, mrEstLogPostValsAtEvalPts] = ...
    obtainProbHandle(mrLogJointEst, evalPts, paramSpaceBounds);

  % Evaluate the KL
  currKL = eval3DKLFromLogProbs(trueLogPostValsAtEvalPts, ...
    mrEstLogPostValsAtEvalPts, paramSpaceBounds);
  mcmcReg_errs(experimentIter, mcmcRegResIter) = currKL;
  fprintf(' MCMC-Reg Iter: %d, num-pts: %d, KL: %0.4f\n', mcmcRegResIter, ...
    currNumMRPts, currKL);

end

