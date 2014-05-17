classdef GCExperiment < handle

  properties
    problemSpaceBounds; % a numDims x 2 matrix giving the lower and upper bounds
                      % for each parameter.
    lowestLogLiklVal;
  end

  methods

    % Constructor
    function obj = GCExperiment(problemSpaceBounds, lowestLogLiklVal)
      obj.problemSpaceBounds = problemSpaceBounds;
      obj.lowestLogLiklVal = lowestLogLiklVal;
    end

    function normCoords = getNormCoords(obj, trueCoords)
      offset = bsxfun(@minus, trueCoords, obj.problemSpaceBounds(:,1)');
      normCoords = bsxfun(@rdivide, offset, ...
        obj.problemSpaceBounds(:,2)' - obj.problemSpaceBounds(:,1)' ); 
    end

    function trueCoords = getTrueCoords(obj, normCoords)
      scaled = bsxfun(@times, normCoords, ...
        obj.problemSpaceBounds(:,2)' - obj.problemSpaceBounds(:,1)' ); 
      trueCoords = bsxfun(@plus, scaled, obj.problemSpaceBounds(:,1)');
    end

    function [logJointProbs] = normCoordLogJointProbs(obj, evalPts)
      [logJointProbs] = obj.trueCoordLogJointProbs(obj.getTrueCoords(evalPts));
    end

    function [logJointProbs] = trueCoordLogJointProbs(obj, evalPts)
      % Identify points that are outside the domain
      less = @(x,y) x<y;
      great = @(x,y) x>y;
      belowDomain = sparse(bsxfun(less, evalPts,obj.problemSpaceBounds(:,1)'));
      aboveDomain = sparse(bsxfun(great, evalPts,obj.problemSpaceBounds(:,2)'));
      outOfDomain = sum(belowDomain + aboveDomain, 2) > 0;
      % Create params for returning
      numPts = size(evalPts, 1);
      logJointProbs = -inf * ones(numPts, 1);
      % Now compute at the points within the bounds
      [inBoundLJPs] = galaxyLogLikelihoodWrap(evalPts(~outOfDomain, :), ...
                                                obj.lowestLogLiklVal);
      logJointProbs(~outOfDomain, :) = max(inBoundLJPs, obj.lowestLogLiklVal);
    end

  end % methods
end % classdef

