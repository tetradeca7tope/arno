function [mcmc_pt, mcmc_log_prob, sample, sample_log_prob] = Custom2DMCMC( ...
  proposalStd, currPt, currLogProb, evalLogProb)
% This function generates a new sample using MCMC and returns the sample and
% the sample_log_prob. If the new sample is accepted, mcmc_pt and
% mcmc_log_prob contain the same values.

  BORDER_TOL = 1e-2;
  LAP_SM_CONSTANT = 0.2;
  
  % the proposal induced by the Gaussian and the logit seems to be too
  % slow to converge. Using a beta centred at the current
  % sample instead

%   % Logit transform & Gaussian
%   logit_curr_pt = logit(currPt);
%   logit_sample = logit_curr_pt + proposalStd * randn(1,2);
%   sample = min( ...
%             max( ...
%               logitinv(logit_sample), ...
%               BORDER_TOL), ...
%             1 - BORDER_TOL);
%   proposal_ratio = 1;

  % 2D Beta centred at current point. Now proposalStd gives the inverse
  % precision of the beta
  beta_precision = 1/ proposalStd;
  mod_curr_pt = [( currPt(1) + LAP_SM_CONSTANT ) / (1 + 2*LAP_SM_CONSTANT), ...
                 ( currPt(2) + LAP_SM_CONSTANT ) / (1 + 2*LAP_SM_CONSTANT)];
  prop1 = [mod_curr_pt(1), 1-mod_curr_pt(1)]/proposalStd;
  prop2 = [mod_curr_pt(2), 1-mod_curr_pt(2)]/proposalStd;
  s1_ = dirichlet_sample(prop1); s1 = s1_(1);
  s2_ = dirichlet_sample(prop2); s2 = s2_(1);
  sample = [s1 s2];
  trans1 = s1_ * beta_precision;
  trans2 = s2_ * beta_precision;
  proposal_ratio = evalBetaProb(prop1(1), prop1(2), sample(1)) * ...
                   evalBetaProb(prop2(1), prop2(2), sample(2)) / ...
                   ( evalBetaProb(trans1(1), trans1(2), mod_curr_pt(1)) * ...
                     evalBetaProb(trans2(1), trans2(2), mod_curr_pt(2)) );
  
  % evaluate the joint at the current point
  sample_log_prob = evalLogProb(sample);

  log_prob_ratio = sample_log_prob - currLogProb + log(proposal_ratio); 
  if (log(rand()) < log_prob_ratio), % accept the proposal
    mcmc_pt = sample;
    mcmc_log_prob = sample_log_prob;
  else
    mcmc_pt = currPt;
    mcmc_log_prob = currLogProb;
  end

end
