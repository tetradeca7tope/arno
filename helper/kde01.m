function [est_probs, f, h] = kde01(X, candidate_hs)
% Performs Kernel density estimation on a dataset confined to [0,1]^d by first
% applying the logit transform. Returns the estimated probabilities and
% a function handle
% The candidate_hs need to be in the logit transformed R^d space

  if ~exist('candidate_hs', 'var')
    candidate_hs = [];
  end

  logit_X = logit(X);
  [~, flogit, h] = kde(logit_X, candidate_hs);
  
%   f = @(data) flogit(logitinv(data)) .* ( exp(data) ./ (1 + exp(data).^2) );
  f = @(data) flogit(logit(data)) ./ prod(data .* (1 - data), 2) ;
  est_probs = []; % f(X); doing this for now to speed things up
end
