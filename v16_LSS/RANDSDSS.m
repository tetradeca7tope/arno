% Implements RAND for SDSS

randPts = rand(NUM_AL_ITERS, numDims);
fprintf('Evaluating the Log likelihood at randomly selected points..\n');
[randLogProbs] = evalLogJoint(randPts);

randLogJointEst = regressionWrap(randPts, randLogProbs, noiseLevelGP, ...
  lowestLogliklVal, logLiklRange, cvCostFunc);

% Save results
save(RAND_EXP_FILE, 'randPts', 'randLogProbs', 'noiseLevelGP', ...
'lowestLogliklVal', 'logLiklRange', 'cvCostFunc');

