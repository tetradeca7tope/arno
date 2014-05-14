% Performs a Random sampling for SNLS

% Use the points in ABC
% randPts = rand(NUM_AL_ITERS, 3);
% fprintf('Evaluating the Log Likelihood for Randomly selected Points..\n');
% randLogProbs = evalLogJoint(randPts);

for randIter = 1:numResultsToBeStored

  currNumRandPts = randIter * STORE_RESULTS_EVERY;
  % Obtain regressors and regressands
  Xtr = randPts(1:currNumRandPts, :);
  Ytr = randLogProbs(1:currNumRandPts);
  
  % Regress on these points
  randLogJointEst = regressionWrap(Xtr, Ytr, noiseLevelGP, lowestLogLiklVal, ...
    logLiklRange, cvCostFunc);
  % obtain the estimated posterior at the grind pts
  [~, randEstLogPostValsAtEvalPts] = ...
    obtainProbHandle(randLogJointEst, evalPts, paramSpaceBounds);

  % Evaluate the KL
  currKL = eval3DKLFromLogProbs(trueLogPostValsAtEvalPts, ...
    randEstLogPostValsAtEvalPts, paramSpaceBounds);
  rand_errs(experimentIter, randIter) = currKL;
  fprintf(' Rand Iter: %d, num-pts: %d, KL: %0.4f\n', randIter, ...
    currNumRandPts, currKL);

end

