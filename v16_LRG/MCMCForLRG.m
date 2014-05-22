% Performs MCMC for SDSS

mcmcEvalLogJoint = @(t) evalLogJoint(logitinv(t));
[logitMcmcSamples, logitMcmcQueries, mcmcLogProbs] = CustomMCMC( ...
  NUM_MCMC_SAMPLES, mcmcProposalStd, mcmcInitPt, mcmcEvalLogJoint);
mcmcSamples = logitinv(logitMcmcSamples);
mcmcQueries = logitinv(logitMcmcQueries);

% Now perform MCMC
for mcmcResIter = 1:numMCMCResultsToBeStored

  currNumMcmcPts = mcmcResIter * STORE_RESULTS_EVERY;
  % Obtain the current set of points
  currMCMCSamples = mcmcSamples(1:currNumMcmcPts, :);

  % Perform KDE and obtain the results
  [~, mcmcProbEst] = kde01(currMCMCSamples);
  currKL = estimKLForLRG( klEvalPts, truePAtEvalPts, mcmcProbEst);

  % record and report
  mcmc_errs(experimentIter, mcmcResIter) = currKL;  
  fprintf(' MCMC iter: %d, #pts: %d, KL: %.4f\n', mcmcResIter, ...
    currNumMcmcPts, currKL);
  
end


% Now perform MCMC REgression
% ===========================
fprintf('MCMC Regression\n==================================');

for mcmcRegResIter = 1:numResultsToBeStored

  currNumMRPts = mcmcRegResIter * STORE_RESULTS_EVERY;
  % Obtain Regressors and Regressands
  Xtr = mcmcQueries(1:currNumMRPts, :);
  Ytr = mcmcLogProbs(1:currNumMRPts);

  % obtain errors and estimates
  [currKL, mrLogJointEst, mrProbEst] = evalRegMethodKLProgress( Xtr, Ytr, ...
    gpFitParams, klEvalPts, truePAtEvalPts, evalMCMCParams, optKDEBandWidth);

  % record and report
  mcmcReg_errs(experimentIter, mcmcRegResIter) = currKL;
  fprintf(' MR Iter: %d, #pts: %d, KL: %.4f\n', mcmcRegResIter, ...
    currNumMRPts, currKL);  

end

