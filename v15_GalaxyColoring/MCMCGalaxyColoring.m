% Implements MCMC for the problem

% Don't do logit for now
[mcmcSamples, mcmcQueries, mcmcLogProbs] = CustomMCMC(NUM_MCMC_SAMPLES, ...
  mcmcProposalStd, mcmcInitPt, evalLogJoint);

% Now perform KDE on the collected Samples
[~, mcmcProbEst] = kde(mcmcSamples);

% Now do regression on the first NUM_AL_ITERS points
fprintf('MCMC Regression\n===========================\n');
Xtr = mcmcQueries(mcmcLogProbs(1:NUM_AL_ITERS) > -inf, : );
Ytr = mcmcLogProbs(mcmcLogProbs(1:NUM_AL_ITERS) > -inf );
mrLogJointEst = regressionWrap(Xtr, Ytr, noiseLevelGP, lowestLogliklVal, ...
  logLiklRange, cvCostFunc);

  % Now plot the values out
  PLOT_LOCAL_OK = false;
  PLOT_LOCAL_OK = true;
  if PLOT_LOCAL_OK
    % Plot the points picked by MBP
    figure;
    if numDims == 1
      plot(mcmcQueries, mcmcLogProbs, 'rx'); hold on;
      th = linspace(paramSpaceBounds(1), paramSpaceBounds(2), 100)';
      plot(th, logJointEst(th), 'g-');
    elseif numDims == 2
      plot3(mcmcQueries(1:NUM_AL_ITERS,1), mcmcQueries(1:NUM_AL_ITERS,2), ...
        mcmcLogProbs(1:NUM_AL_ITERS, :), 'kx', 'MarkerSize', 10);
      hold on;
      plot2DFunction(mrLogJointEst, [paramSpaceBounds(1,:), ...
        paramSpaceBounds(2,:)], 'mesh');
      axis([0, 1, 0, 1, (min(Ytr)-1), (max(Ytr)+1)]);
      figure;
      plot(mcmcSamples(:,1), mcmcSamples(:,2), 'rx');
    else
      plot(mcmcQueries(:,1), mcmcQueries(:,2), 'rx'); hold on;
      axis([paramSpaceBounds(1,:), paramSpaceBounds(2, :)]);
    end
  end

