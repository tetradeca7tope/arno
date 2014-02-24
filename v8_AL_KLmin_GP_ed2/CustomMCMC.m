function [sampled_pts, sampled_log_probs, prop_pts, prop_log_probs] = ...
  CustomMCMC(numSamples, proposal_std, init_pt, init_log_joint, evalLogProb);
% This function generates a bunch of samples using MCMC
% numSamples: number of samples to generate
% proposal_var: variance of the proposal distribution
% init_pt: initial point (leave empty if you want to initialize randomly
% evalLogProb: function handle for computing the probabilty
% The proposal distributions are a set of beta binomial distributions.

BORDER_TOL = 1e-4;

  if isempty(init_pt)
    init_pt = rand(); % pick a random point uniformly.
    init_log_joint = evalLogJointProbs(init_pt);
  end

  sampled_pts = zeros(numSamples+1, 1); % add the init point and then remove
  sampled_probs = zeros(numSamples+1, 1);
%   init_pt,
  sampled_pts(1) = init_pt;
  sampled_log_probs(1) = init_log_joint;
  prop_pts = zeros(numSamples, 1);
  prop_log_probs = zeros(numSamples, 1);

  for i = 2:numSamples+1
    logit_curr_sample = logit(sampled_pts(i-1));
    logit_new_sample = logit_curr_sample + proposal_std * randn();
    new_sample = min( ...
                   max(...
                     logitinv(logit_new_sample), ...
                     BORDER_TOL), ...
                   1-BORDER_TOL );
    new_log_prob = evalLogProb(new_sample);
    prop_pts(i-1) = new_sample;
    prop_log_probs(i-1) = new_log_prob;

    % Compute the ratios and decide whether or not to accept
    log_prob_ratio = new_log_prob - sampled_log_probs(i-1); 
    fprintf('new_sample: %f, log-prob: %f, logp-ratio: %f, ', ...
            new_sample, new_log_prob, log_prob_ratio);
%     new_sample, new_prob, prob_ratio,
    if log(rand()) < log_prob_ratio, % accept the proposal
      sampled_pts(i) = new_sample;
      sampled_log_probs(i) = new_log_prob;
    else % reject the proposal and use the same values
      sampled_pts(i) = sampled_pts(i-1);
      sampled_log_probs(i) = sampled_log_probs(i-1);
    end
  end

  sampled_pts = sampled_pts(2:end);
  sampled_log_probs = sampled_log_probs(2:end);

%   [prop_pts, prop_log_probs],
end
