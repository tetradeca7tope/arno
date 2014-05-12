% Script runs MaxBandPoint for Galaxy Coloring

% Initialization
initialPts = [];
initialLogProbs = [];

% Set parameters for MBP
phi = @exp; gradPhi = @exp; % use exponential transformation
alMbpParams.num_iters = 10;
alMbpParams.init_step_size = 1;
% First obtain the points via MBP
[mbpPts, mbpLogProbs, mbpLipConst] = alMaxBandPoint(evalLogJoint, ...
  initialPts, initialLogProbs, phi, gradPhi, INIT_LIPSCHITZ_CONST, ...
  paramSpaceBounds, NUM_AL_ITERS, alMbpParams);

% Now perform regression on each of these points to obtain the estimates.
logJointEst = regressionWrap(mbpPts, mbpLogProbs, noiseLevelGP, ...
  lowestLogliklVal, logLiklRange);

  PLOT_LOCAL_OK = false;
  PLOT_LOCAL_OK = true;
  if PLOT_LOCAL_OK
    % Plot the points picked by MBP
    figure;
    if numDims == 1
      plot(mbpPts, mbpLogProbs, 'rx'); hold on;
      th = linspace(paramSpaceBounds(1), paramSpaceBounds(2), 100)';
      plot(th, logJointEst(th), 'g-');
    elseif numDims == 2
      plot3(mbpPts(:,1), mbpPts(:,2), mbpLogProbs, 'kx', 'MarkerSize', 10);
      hold on;
      plot2DFunction(logJointEst, [paramSpaceBounds(1,:), ...
        paramSpaceBounds(2,:)], 'mesh');
    else
      plot(mbpPts(:,1), mbpPts(:,2), 'rx'); hold on;
      axis([paramSpaceBounds(1,:), paramSpaceBounds(2, :)]);
    end
    title('Pts chosen by MBP');
  end

