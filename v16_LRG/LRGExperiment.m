classdef LRGExperiment < handle

  properties
    problemSpaceBounds = [-1    0.9 ; ... Omega_k: 0
                          0     1   ; ... Omega_Lambda: 0.762
                          0     1.2 ; ... omega_c: 0.1045
                          0.001 0.25; ... Omega_b: 0.02233
                          0.5   1.7 ; ... n_s: 0.951
                          0.65  0.75; ... A_s: 0.6845
                          -0.1  0.1 ; ... alpha: 0
                          0     3.0 ];  % b: 1.908
    lowestLogLiklVal; % The highest is -9 so this is good enough.
    % Set default values for A_s, alpha and Q_nl
    DFLT_QNL = 30.81;
  end

  methods

    % Constructor
    function obj = LRGExperiment(lowestLogLiklVal)
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
      belowDomain= sparse(bsxfun(less, evalPts, obj.problemSpaceBounds(:,1)'));
      aboveDomain= sparse(bsxfun(great, evalPts, obj.problemSpaceBounds(:,2)'));
      outOfDomain= sum(belowDomain + aboveDomain, 2) > 0;
      % Create params for returning
      numPts = size(evalPts, 1);
      numInBoundPts = numPts - sum(outOfDomain);
      logJointProbs = -5 * abs(obj.lowestLogLiklVal) * ones(numPts, 1);
      % Now add A_s, alpha and Q_nl to evalPts
      augEvalPts = [evalPts(~outOfDomain, 1:8) ...
                    obj.DFLT_QNL * ones(numInBoundPts, 1)];
      % Now compute at the points within the bounds
      inBoundLJPs = lrgLogLiklWrap( augEvalPts, obj.lowestLogLiklVal);
      logJointProbs(~outOfDomain, :) = max(inBoundLJPs, obj.lowestLogLiklVal);
    end

  end % methods
end % classdef

