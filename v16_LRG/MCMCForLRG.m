% Performs MCMC for SDSS

[logitMcmcSamples, logitMcmcQueries, mcmcLogProbs] = CustomMCMC( ...
  NUM_MCMC_SAMPLES, mcmcProposalStd, mcmcInitPt, mcmcLogJoint);
max(max(logitMcmcSamples)), min(min(logitMcmcSamples)),
if DOING_LOGIT
  mcmcSamples = logitinv(logitMcmcSamples);
  mcmcQueries = logitinv(logitMcmcQueries);
else
  mcmcSamples = logitMcmcSamples;
  mcmcQueries = logitMcmcQueries;
end

% Now perform MCMC
for mcmcResIter = 1:numMCMCResultsToBeStored
% for mcmcResIter = 1:0 %numMCMCResultsToBeStored

  currNumMcmcPts = STORE_MCMC_RESULTS_AT(mcmcResIter);
  % Obtain the current set of points
  currMCMCSamples = mcmcSamples(1:currNumMcmcPts, :);

  % Perform KDE and obtain the results
  [~, mcmcProbEst] = kde01(currMCMCSamples);
  mcmcLogProbEst = @(t) mcmcProbEst(t);
  estGTLogProbs = mcmcLogProbEst(gtPts);
  gtProbs = exp(gtLogProbs);
  estGTProbs = exp(estGTLogProbs);
  k = estGTProbs(1)/gtProbs(1);
  estGTProbs = estGTProbs / k;
  currErr = 1/numGTPts * ( sum(gtProbs.^2) - ...
                           (estGTProbs'*gtProbs)^2/sum(estGTProbs.^2) );
  currErr = sqrt( mean ( (estGTProbs - gtProbs).^2 ) );

  % record and report
  mcmc_errs(experimentIter, mcmcResIter) = currErr;  
  fprintf(' MCMC iter: %d, #pts: %d, err: %e\n', mcmcResIter, ...
    currNumMcmcPts, currErr);
  
end

% Now perform MCMC REgression
% ===========================
fprintf('MCMC Regression\n==================================');
mrPts = mcmcQueries(1:NUM_AL_ITERS, :);
mrLogProbs = mcmcLogProbs(1:NUM_AL_ITERS, :);

for mcmcRegResIter = 1:numResultsToBeStored

  currNumMRPts = STORE_RESULTS_AT(mcmcRegResIter);
  % Obtain Regressors and Regressands
  Xtr = mcmcQueries(1:currNumMRPts, :);
  Ytr = mcmcLogProbs(1:currNumMRPts);

  % obtain errors and estimates
  [currErr, mrLogJointEst]= evalRegMethodProgress(Xtr, Ytr, ...
    gtPts, gtLogProbs, gpFitParams);

  % record and report
  mcmcReg_errs(experimentIter, mcmcRegResIter) = currErr;
  fprintf(' MR Iter: %d, #pts: %d, err: %e\n', mcmcRegResIter, ...
    currNumMRPts, currErr);

end

