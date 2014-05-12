classdef SNExperiment < handle

  properties
    data; % the dataset being used
    paramSpaceBounds; % a numDims x 2 matrix giving the lower and upper bounds
                      % for each parameter.
    numObsToUse; % The number of observations to use
    lowestLogLiklVal;
  end

  methods

    % Constructor
    function obj = SNExperiment(dataFile, paramSpaceBounds, numObsToUse, ...
                                lowestLogLiklVal)
      obj.data = load(dataFile);
      obj.paramSpaceBounds = paramSpaceBounds;
      obj.lowestLogLiklVal = lowestLogLiklVal;
      if isempty (numObsToUse), obj.numObsToUse = size(obj.data, 1);
      else, obj.numObsToUse = numObsToUse;
      end
    end

    function normCoords = getNormCoords(obj, trueCoords)
      offset = bsxfun(@minus, trueCoords, obj.paramSpaceBounds(:,1)');
      normCoords = bsxfun(@rdivide, offset, ...
        obj.paramSpaceBounds(:,2)' - obj.paramSpaceBounds(:,1)' ); 
    end

    function trueCoords = getTrueCoords(obj, normCoords)
      scaled = bsxfun(@times, normCoords, ...
        obj.paramSpaceBounds(:,2)' - obj.paramSpaceBounds(:,1)' ); 
      trueCoords = bsxfun(@plus, scaled, obj.paramSpaceBounds(:,1)');
    end

    function [logJointProbs, lumMeans] = normCoordLogJointProbs(obj, evalPts)
      [logJointProbs, lumMeans] = obj.trueCoordLogJointProbs( ...
                                    obj.getTrueCoords(evalPts));
    end

    function [logJointProbs, lumMeans] = trueCoordLogJointProbs(obj, evalPts)
      % Identify points that are outside the domain
      belowDomain = sparse(bsxfun(@le, evalPts, obj.paramSpaceBounds(:,1)') );
      aboveDomain = sparse(bsxfun(@ge, evalPts, obj.paramSpaceBounds(:,2)') );
      outOfDomain = sum(belowDomain + aboveDomain, 2) > 0;
      % Create params for returning
      numPts = size(evalPts, 1);
      logJointProbs = obj.lowestLogLiklVal * ones(numPts, 1);
      lumMeans = -1000 * ones(numPts, obj.numObsToUse);
      % Now compute at the points within the bounds
      [inBoundLJPs, inBoundLMs] = evalSNLSLogLikl(evalPts(~outOfDomain, :), ...
                                                   obj.data, obj.numObsToUse);
      logJointProbs(~outOfDomain, :) = max(inBoundLJPs, obj.lowestLogLiklVal);
      lumMeans(~outOfDomain, :) = inBoundLMs;
    end

  end % methods
end % classdef

