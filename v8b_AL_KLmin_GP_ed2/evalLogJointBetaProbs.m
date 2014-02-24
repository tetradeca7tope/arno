function [log_joint] = evalLogJointBetaProbs(evalPts, sumX, N, ...
  priorParams, noise_level)
% This function obtains some noisy estimates of the joint probability.
% the Points at which the joint should be evaluated should be stored in
% the vector evalPts
% CC is also a parameter of this distribution. 

  alp = priorParams.alp;
  bep = priorParams.bep;

  th = evalPts; % for brevity in the expression below.
  log_joint = (sumX + alp -1)*log(th) + (N - sumX + bep -1)*log(1-th) + ...
              -log(beta(alp, bep));

end
