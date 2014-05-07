function [evalProb, logProbs, normConst] = ...
  obtainProbHandle(evalUnNormLogProb, evalPts, bounds)
% Returns
% 1. evalProb: A function handle for the true probability 
% 2. logProbs: the log(prob) values after normalization
% 3. normConst: the Normalizing Constant

  vol = prod( bounds(:,2) - bounds(:,1) );
  unNormLogProbs = evalUnNormLogProb(evalPts);
  normConst = vol * mean(exp(unNormLogProbs));
  evalProb = @(arg) exp(evalUnNormLogProb(arg)) / normConst;
  logProbs = unNormLogProbs - log(normConst);

end

