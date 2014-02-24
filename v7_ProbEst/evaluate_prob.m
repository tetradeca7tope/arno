function [p] = evaluate_prob(th, Z, w, h)
% th: the points at which you need to evaluate the probability
% Z: the points at which the probs were observed
% w: weights
% h: bandwidth parameter

  data_dim = size(Z,2);
  D = dist2(th, Z);
  K = exp( - D/ (2*h^2));
  gauss_norm_constant = 1/ (2 * pi * h^2) ^ (data_dim/2);
  p = gauss_norm_constant * K * w;

end
