function [evalProb, logProbs, normConst] = ...
  obtainProbHandle(evalUnNormLogProb, evalPts, bounds)
% Returns
% 1. evalProb: A function handle for the true probability 
% 2. logProbs: the log(prob) values after normalization
% 3. normConst: the Normalizing Constant

  vol = prod( bounds(:,2) - bounds(:,1) );

  % Don't evaluate more than 10000 points at a time
  numPts = size(evalPts, 1);
  unNormLogProbs = zeros(numPts, 1);
  numPtsPerEpoch = 10000;
  for i = 1:ceil(numPts/numPtsPerEpoch)
    currStartIdx = (i - 1) * numPtsPerEpoch + 1;
    currEndIdx = min(i * numPtsPerEpoch, numPts);
    currIdxs = (currStartIdx:currEndIdx)';
    unNormLogProbs(currIdxs) = evalUnNormLogProb(evalPts(currIdxs, :));
  end
    
%   unNormLogProbs = evalUnNormLogProb(evalPts);

  normConst = vol * mean(exp(unNormLogProbs));
  evalProb = @(arg) exp(evalUnNormLogProb(arg)) / normConst;
  logProbs = unNormLogProbs - log(normConst);

end

