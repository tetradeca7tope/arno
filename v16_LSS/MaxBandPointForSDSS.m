% Performs Max Band Point for SDSS

% Initialization
initialPts = [];
initialLogProbs = [];

% Set parameters for MBP
phi = @exp; gradPhi = @exp; % Use exponential Transformation
alMbpParams.num_iters = 10;
alMbpParams.init_step_size = 1;
% First obtain the points via MBP
[mbpPts, mbpLogProbs, mbpLipConst] = alMaxBandPoint(evalLogJoint, ...
  initialPts, initialLogProbs, phi, gradPhi, initLipschitzConstant, ...
  paramSpaceBounds, NUM_AL_ITERS, alMbpParams);

% Now perform regression on each of these points to obtain the estimates
mbpLogJointEst = regressionWrap(mbpPts, mbpLogProbs, noiseLevelGP, ...
  lowestLogliklVal, logLiklRange);

% Save results
save(MBP_EXP_FILE, 'mbpPts', 'mbpLogProbs', 'initialPts', 'noiseLevelGP', ...
  'lowestLogliklVal', 'logLiklRange');

