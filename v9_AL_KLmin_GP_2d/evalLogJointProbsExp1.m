function [log_joint] = evalLogJointProbsExp1(evalTh1, evalTh2, sumX, N, ...
  priorParams, noise_level)
% This function obtains some noisy estimates of the joint probability.
% the Points at which the joint should be evaluated should be stored in
% the vector evalPts
% CC is also a parameter of this distribution. 

  alp1 = priorParams(1);
  bep1 = priorParams(2);
  alp2 = priorParams(3);
  bep2 = priorParams(4);

  th1 = evalTh1;
  th2 = evalTh2;
  p = computeP(th1, th2);

  log_joint = sumX*log(p) + (N-sumX)*log(1-p) + ...
              (alp1-1)*log(th1) + (bep1-1)*log(1-th1) + ...
                - log(beta(alp1, bep1)) + ...
              (alp2-1)*log(th2) + (bep2-1)*log(1-th2) + ...
                - log(beta(alp2, bep2));

  % Add some noise
  noise_level = 0;
  log_joint = log_joint + noise_level * randn(size(log_joint));

end
