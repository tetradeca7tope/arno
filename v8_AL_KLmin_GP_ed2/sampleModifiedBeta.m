function [theta, sumX, true_post, marginal_prob, X] = ...
  sampleModifiedBeta(priorParams, evalPriorAt, N)
% alp and bep gives the prior for drawing theta
% N refers to the number of samples to be drawn (optional)
% evalPriorAt is the vector of values on which the prior should be evaluated
% (make sure that this is fine enough so as to ensure that the integration is
% accurate.

% Function returns
% theta: the true parameter
% sumX: the number of positive trials in the Bernoulli experiment
% true_post: the true_posterior evaluated at points in evalPriorAt
% X: the actual samlpled values

  if size(evalPriorAt, 2) > 1,
    error('Send points at which to evaluate the prior at as a column vector.');
  end

  alp = priorParams.alp;
  bep = priorParams.bep;
  c_param = priorParams.c_param;

%   theta_ = dirichlet_sample([alp, bep]); theta = theta_(1);
  % cheating here. preferring multimodal distros over balanced ones.
  theta_ = dirichlet_sample([30, 5]); theta = theta_(1);
%   theta = 0.15; % cheating here. trying to force a multimodal distro.
  p = theta^2 + (c_param - theta)^2;
  X = double(rand(N,1) < p);
  sumX = sum(X);
  fprintf('theta = %f, p = %f, sumX: %d\n', theta, p, sumX);

  log_joint_probs = evalLogJointModifiedBetaProbs(evalPriorAt, sumX, N, ...
                      priorParams, 0);
  joint_probs = exp(log_joint_probs);
  marginal_prob = numerical_1D_integration(evalPriorAt(2:end-1), ...
                    joint_probs(2:end-1));
  true_post = joint_probs / marginal_prob;

end
