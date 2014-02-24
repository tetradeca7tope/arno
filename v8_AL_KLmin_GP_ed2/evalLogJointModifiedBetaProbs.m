function [log_joint] = evalLogJointModifiedBetaProbs(evalPts, sumX, N, ...
  priorParams, noise_level)
% This function obtains some noisy estimates of the joint probability.
% the Points at which the joint should be evaluated should be stored in
% the vector evalPts
% CC is also a parameter of this distribution. 

  alp = priorParams.alp;
  bep = priorParams.bep;
  c_param = priorParams.c_param;

  th = evalPts; % for brevity in the expression below.
  p = (th.^2 + (c_param - th).^2);
  log_joint = sumX*log(p) + (N-sumX)*log(1-p) + ...
              (alp-1)*log(th) + (bep-1)*log(1-th) + ...
              -log(beta(alp, bep));

end
