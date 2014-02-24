function [theta, sumX, true_post, marginal_prob, X] = ...
  sampleExp1(priorParams, evalTh1, evalTh2, N)
% alp and bep gives the prior for drawing theta
% evalPriorAt is the vector of values on which the prior should be evaluated
% (make sure that this is fine enough so as to ensure that the integration is
% accurate.

% Function returns
% theta: the true parameter
% sumX: the number of positive trials in the Bernoulli experiment
% true_post: the true_posterior evaluated at points in evalPriorAt
% X: the actual samlpled values

  alp1 = priorParams(1);
  bep1 = priorParams(2);
  alp2 = priorParams(3);
  bep2 = priorParams(4);

  % cheating here. preferring multimodal distros over balanced ones.
  theta = zeros(2,1);
  theta1 = dirichlet_sample([20, 10]); theta(1) = theta1(1);
  theta2 = dirichlet_sample([20, 10]); theta(2) = theta2(1);
%   theta1 = dirichlet_sample([alp1, bep1]); theta(1) = theta1(1);
%   theta2 = dirichlet_sample([alp2, bep2]); theta(2) = theta2(1);
%   theta1 = 0.2; theta2 = 0.4;

  % brevity
  t1 = theta(1);
  t2 = theta(2);
  p = computeP(t1, t2);

  X = double(rand(N,1) < p);
  sumX = sum(X);
  fprintf('theta = (%.5f,%.5f), p = %f, sumX: %d\n', t1, t2, p, sumX);

  log_joint_probs = evalLogJointProbsExp1 (evalTh1, evalTh2, sumX, N, ...
                      priorParams, 0);
  joint_probs = exp(log_joint_probs);
  marginal_prob = numerical_2D_integration(joint_probs, evalTh1, evalTh2);
  true_post = joint_probs / marginal_prob;

end

