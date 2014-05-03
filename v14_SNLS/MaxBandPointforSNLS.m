% MaxBanPoint for SNLS

initLipschitzConst = 1000; % I think this is an under estimate

% Set params for alMaxBandPoint
phi = @exp; gradPhi = @exp;
alMbpParams.num_grad_desc_init_pts = 300;
alMbpParams.init_step_size = 0.1;
[mbpPts, mbpLogProbs, mbpLipConst] = alMaxBandPoint(evalLogJoint, [], [], ...
  phi, gradPhi, initLipschitzConst, paramSpaceBounds, NUM_AL_ITERS, ...
  alMbpParams);

% Now estimate the regression function
funcHandle = 

