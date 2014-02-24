function [log_joint, garb1, true_post, marginal_prob, garb2] = ...
  evalLogJointProbsExp2(evalTh1, evalTh2, noise_level)
% This function obtains some noisy estimates of the joint probability.
% the Points at which the joint should be evaluated should be stored in
% the vector evalPts

  P1 = 0.3;
  P2 = 0.5;
  P3 = 1 - P1 - P2;
  M1 = [0.2; 0.2];
  M2 = [0.8; 0.35];
  M3 = [0.6; 0.75];
  % Generate the following fixed covariance matrices
  S1 = 0.1 * [1; 0.5; 0.5; 0.8];
  S2 = 0.2 * [.9; -0.05; -0.05; 1];
  S3 = 0.2 * [1.2; -0.2; -0.2; 0.3];
  % Generate some rando covariance matrices

  th = [Th1(:), Th2(:)]; 
  jp = P1 * mvnpdf(th, M1, S1) + P2 * mvnpdf(th, M2, S2) + ...
       P3 * mvnpdf(th, M3, S3);
  log_joint = log(jp) - 40; % the -40 is just to affect scaling

  if nargout > 1
    joint_probs = exp(log_joint);
    marginal_prob = numerical_2D_integration_wrap(th, joint_probs);
    true_post = joint_probs / marginal_prob;

    garb1 = [0.0; 0.0];
  end

%   noise_level = 0;
  log_joint = log_joint + noise_level * randn(size(log_joint));

end
