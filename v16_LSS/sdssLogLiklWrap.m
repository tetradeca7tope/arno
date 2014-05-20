function [logLiklVals] = sdssLogLiklWrap(evalAtPts, lowestLogLiklVal)
% evalPts is an numPtsx9 array. Each row is a point at which we want to
% evaluate the log Likelihood.
% These are the cosmological Parameters
% Omega_k     : Spatial Curvature [-1, 0.9]
% Omega_Lambda: Dark Energy Fraction [0, 1]
% omega_C     : Cold Dark Matter Density [0 1.2]
% omega_B     : Baryon Density [0.001 - 0.25]
% n_s         : Scalar Spectral Index [0.5 - 1.7]
% A_s         : Scalar Fluctuation Amplitude = 0.6845
% alpha       : Running of Spectral Index = 0.0
% b           : Galaxy Bias [0.0 3.0]
% Q_nl        : Nonlinear Correction = 30.81

  % Prelims
  numPts = size(evalAtPts, 1);
  numDims = size(evalAtPts, 2); % This should be 9 ?
  OUT_FILE = 'sim/liklOut.txt'; % hard-coded for now

  logLiklVals = zeros(numPts, 1);

  % NOw call the simulator iteratively
  for iter = 1:numPts

    % First create the command
    currEvalPt = evalAtPts(iter, :);
    commandStr = 'cd sim && ./a.out ';
    for k = 1:9
      commandStr = sprintf('%s %f ', commandStr, currEvalPt(k));
    end
    commandStr = sprintf('%s && cd ..', commandStr);
    
    % Execute the command
    system(commandStr);

    % Read from file
    outVal = load(OUT_FILE);
    if isnan(outVal), outVal = -inf; end
    logLiklVals(iter) = max(outVal, lowestLogLiklVal);

  end

end

