function [log_joint, garb1, true_post, marginal_prob, garb2] = ...
  evalLogJointBiNormal(evalPts, noise_level)
% garb1, garb2 are garbage variables.

  P1 = 0.6;
  P2 = 1 - P1;
  M1 = 0.25;
  M2 = 0.65;
  S1 = 0.015;
  S2 = 0.035;

  joint_probs = P1 * normpdf(evalPts, M1, S1) + P2 * normpdf(evalPts, M2, S2);
  log_joint = log(joint_probs) - 25; % the - 25 is just scaling

  % find scaling constant and renormalize
  joint_probs = exp(log_joint);
  marginal_prob = numerical_1D_integration(evalPts, joint_probs);
  true_post = joint_probs / marginal_prob;

  garb1 = 0;
  garb2 = 0;
  % add noise before returning
  log_joint = log_joint + noise_level * 0.5 * randn(size(log_joint));

end
