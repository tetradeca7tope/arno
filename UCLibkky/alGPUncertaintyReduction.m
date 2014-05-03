function [ucPts, ucVals] = alGPUncertaintyReduction( oracle, initPts, ...
  initVals, bounds, numALIters, params);
% Performs Uncertainty Reduction to Pick points. At each iteration, it fits a GP
% using alBandwidth and alScale and picks the point with the highest
% uncertainty.
% The following are the fields in params.
% - numALCandidates: the number of candidates to choose from each iteration
% - lowestLogliklVal: the lowest value for the log likeliihood
% - alBandWidth: the bandwidth for GP regression
% - alScale: the Scale parameter for GP regression. 
% - gpNoiseLevel: the Noise level for GPs
% - logLiklRange: the range for the likelihood. Something in the ballpart is ok.
% The last is not strictly needed. But if you provide logLiklRange but do not
% provide gpNoiseLevel and/or alScale, the code will automatically pick some pts


  % If Init points are not given, initialize at the centre of the rectangle.
  if isempty(initPts)
    initPts = [(bounds(:,1) + bounds(:,2))/2]';
    initVals = oracle(initPts);
  end

  % Prelims
  numInitPts = size(initPts, 1);
  numDims = size(initPts, 2); 

  % Check for parameters expected in params
  if ~isfield(params, 'numALCandidates')
    params.numALCandidates = 100;
  end
  if ~isfield(params, 'alBandwidth')
    params.alBandwidth = numALIters^(-1/(1.3 + numDims));
    warning('Bandwidth not specified. Using %f\n', params.alBandwidth);
  end
  if isfield(params, 'logLiklRange')
  % then use some heuristics to set the scale and the noise level
    if ~isfield(params, 'alScale')
      params.alScale = params.logLiklRange/2;
    end
    if ~isfield(params, 'gpNoiseLevel')
      params.gpNoiseLevel = params.logLiklRange/100;
    end
  end

  % Define the following before proceeding
  ucPts = initPts;
  ucVals = initVals;

  fprintf('Performing alMaxBandPoint (dim = %d)\n', numDims);
  for ucIter = 1:numALIters

    % prelims
    numUcPts = size(ucPts, 1);

    % First specify the candidates
    ucCandidates = bsxfun( @plus, ...
      bsxfun(@times, rand(params.numALCandidates, numDims), ...
             (bounds(:,2) - bounds(:,1))' ), ...
      bounds(:,1)' );

    % Now, Run GP Regression on the candidates
    gpHyperParams.noise = params.gpNoiseLevel * ones(numUcPts, 1);
    gpHyperParams.meanFunc = @(arg) params.lowestLogliklVal;
    gpHyperParams.sigmaSm = params.alBandwidth;
    gpHyperParams.sigmaPr = params.alScale;
    [ucM, ~, ucK] = GPRegression(ucPts, ucVals, ucCandidates, gpHyperParams);
    
    % Now pick the best point among the candidates.
    ucS = diag(ucK);
    uncert = (exp(ucS) - 1) .* exp(2*ucM + ucS);
    [~, ucMaxIdx] = max(uncert);
    ucMaxPt = ucCandidates(ucMaxIdx, :);

    % Finally evaluate the log likelihood at this point and store it.
    ucMaxPtVal = oracle(ucMaxPt);
    ucPts = [ucPts; ucMaxPt];
    ucVals = [ucVals; ucMaxPtVal];

  end

end
