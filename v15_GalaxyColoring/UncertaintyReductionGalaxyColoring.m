% Sets up Uncertainty reduction for Galaxy Coloring

% First obtain the points via Active Learning
ucParams.numALCandidates = numALCandidates;
ucParams.lowestLogliklVal = lowestLogliklVal;
ucParams.alBandwidth = alBandwidth;
ucParams.alScale = alScale;
ucParams.gpNoiseLevel = noiseLevelGP;
[ucPts, ucLogProbs] = alGPUncertaintyReduction(evalLogJoint, [], [], ...
  paramSpaceBounds, NUM_AL_ITERS, ucParams);

% Now perform regression and plot the curve
ucLogJointEst = regressionWrap(ucPts, ucLogProbs, noiseLevelGP, ...
  lowestLogliklVal, logLiklRange, cvCostFunc);

  PLOT_LOCAL_OK = false;
  PLOT_LOCAL_OK = true;
  % Plot points picked by MBP
  if PLOT_LOCAL_OK
    figure;
    if numDims == 1
      plot(ucPts, ucLogProbs, 'rx'); hold on;
      th = linspace(paramSpaceBounds(1), paramSpaceBounds(2), 100)';
      plot(th, ucLogJointEst(th), 'g-');
    elseif numDims == 2
      plot3(ucPts(:,1), ucPts(:,2), ucLogProbs, 'kx', 'MarkerSize', 10);
      hold on;
      plot2DFunction(ucLogJointEst, [paramSpaceBounds(1,:), ...
        paramSpaceBounds(2,:)], 'mesh');
    else
      plot(ucPts(:,1), ucPts(:,2), 'rx'); hold on;
      axis([paramSpaceBounds(1,:), PARAM_SPACE_BOUNDS(2, :)]);
    end
    title('Pts chosen by UC');
  end
