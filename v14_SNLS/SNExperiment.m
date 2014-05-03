classdef SNExperiment < handle

  properties
    data; % the dataset being used
    paramSpaceBounds; % a numDims x 2 matrix giving the lower and upper bounds
                      % for each parameter.
    lowestLogLiklVal = -2000;
  end

  methods

    % Constructor
    function obj = SNExperiment(dataFile, paramSpaceBounds)
      obj.data = load(dataFile);
      obj.paramSpaceBounds = paramSpaceBounds;
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
      [logJointProbs, lumMeans] = evalSNLSLogLikl(evalPts, obj.data);
      logJointProbs = max(logJointProbs, obj.lowestLogLiklVal);
    end

  end % methods
end % classdef

