% Implements RAND for SDSS

randPts = rand(NUM_AL_ITERS, numDims);
fprintf('Evaluating the Log likelihood at randomly selected points..\n');
[randLogProbs] = evalLogJoint(randPts);

for randIter = 1:numResultsToBeStored

  currNumRandPts = STORE_RESULTS_AT(randIter);
  Xtr = randPts(1:currNumRandPts, :);
  Ytr = randLogProbs(1:currNumRandPts);

  % Obtain errors and estimates
  [currErr, randLogJointEst]= evalRegMethodProgress(Xtr, Ytr, ...
    gtPts, gtLogProbs, gpFitParams);

  % record and report
  rand_errs(experimentIter, randIter) = currErr;
  fprintf(' Rand Iter: %d, #pts: %d, Err: %e\n', randIter, ...
    currNumRandPts, currErr);

end

