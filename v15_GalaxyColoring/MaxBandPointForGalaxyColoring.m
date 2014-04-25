% Script runs MaxBandPoint for Galaxy Coloring

% Parameters for MaxBandPoint
INIT_LIPSCHITZ_CONST = 10;

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
  PARAM_SPACE_BOUNDS, NUM_AL_ITERS, alMbpParams);

% Now perform regression on each of these points to obtain the estimates.
logJointEst = regressionWrap(mbpPts, mbpLogProbs, GP_NOISE_LEVEL, ...
  LOWEST_LOGLIKL_VAL, LOGLIKL_RANGE);

  PLOT_LOCAL_OK = false;
  PLOT_LOCAL_OK = true;
  if PLOT_LOCAL_OK
    % Plot the points picked by MBP
    if NUM_DIMS == 1
      plot(mbpPts, mbpLogProbs, 'rx'); hold on;
      th = linspace(PARAM_SPACE_BOUNDS(1), PARAM_SPACE_BOUNDS(2), 100)';
      plot(th, logJointEst(th), 'g-');
    elseif NUM_DIMS == 2
      figure;
      plot3(mbpPts(:,1), mbpPts(:,2), mbpLogProbs, 'kx', 'MarkerSize', 10);
      hold on;
      plot2DFunction(logJointEst, [PARAM_SPACE_BOUNDS(1,:), ...
        PARAM_SPACE_BOUNDS(2,:)], 'mesh');
    else
      plot(mbpPts(:,1), mbpPts(:,2), 'rx'); hold on;
      axis([PARAM_SPACE_BOUNDS(1,:), PARAM_SPACE_BOUNDS(2, :)]);
    end
    title('Pts chosen by MBP');
  end

