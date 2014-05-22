% Implements RAND for SDSS

randPts = rand(NUM_AL_ITERS, numDims);
fprintf('Evaluating the Log likelihood at randomly selected points..\n');
[randLogProbs] = evalLogJoint(randPts);

for randIter = 1:numResultsToBeStored

  currNumRandPts = randIter * STORE_RESULTS_EVERY;
  Xtr = randPts(1:currNumRandPts, :);
  Ytr = randLogProbs(1:currNumRandPts);

  % Obtain errors and estimates
  [currKL, randLogJointEst, randProbEst] = evalRegMethodKLProgress(Xtr, Ytr, ...
    gpFitParams, klEvalPts, truePAtEvalPts, evalMCMCParams, optKDEBandWidth);

  % record and report
  rand_errs(experimentIter, randIter) = currKL;
  fprintf(' Rand Iter: %d, #pts: %d, KL: %.4f\n', randIter, currNumRandPts, ...
    currKL);

end

