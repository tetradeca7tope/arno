function [est_probs, f] = kde01(X)
% Performs Kernel density estimation on a dataset confined to [0,1]^d by first
% applying the logit transform. Returns the estimated probabilities and
% a function handle

  logit_X = logit(X);
  [~, flogit] = kde(logit_X);
  
%   f = @(data) flogit(logitinv(data)) .* ( exp(data) ./ (1 + exp(data).^2) );
  f = @(data) flogit(logit(data)) ./ prod(data .* (1 - data), 2) ;
  est_probs = f(X);
end
