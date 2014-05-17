% This script implements ABC for the SNLS experiment

% Randomly sample points from the prior
randPts = rand(NUM_MCMC_SAMPLES, 3);
fprintf('Evaluating the Log Likelihood for Randomly selected Points..\n');
[randLogProbs, randNormSamples] = evalLogJoint(randPts);

for abcIter = 1:numMCMCResultsToBeStored

  currNumABCPts = abcIter * STORE_RESULTS_EVERY;
  currABCPts = randPts(1:currNumABCPts, :);
  currABCNormSamples = randNormSamples(1:currNumABCPts);

  % Perform ABC
  selABCIdxs = abs(currABCNormSamples) < 0.05;
  numSelPts = sum(selABCIdxs);
  selABCPts = currABCPts(selABCIdxs, :);

  % Do a KDE
  if numSelPts > 0
    [~, abcProbEst] = kde01(selABCPts);
    abcLogProbEstAtEvalPts = log(abcProbEst(evalPts));

    % Evaluate the KL
    currKL = eval3DKLFromLogProbs(trueLogPostValsAtEvalPts, ...
      abcLogProbEstAtEvalPts, paramSpaceBounds);
  else
    currKL = Inf;
  end

  abc_errs(experimentIter, abcIter) = currKL;
  fprintf(' ABC Iter: %d, num-pts: %d, num-sel-pts: %d, KL: %0.4f\n', ...
    abcIter, currNumABCPts, numSelPts, currKL);

end

