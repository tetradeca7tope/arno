function [p, sumX, true_post, marginal_prob, X] = ...
  sampleBeta(priorParams, evalPriorAt, N)
% priorParams gives the parameters for the prior.

  alp = priorParams.alp;
  bep = priorParams.bep;

  % Sample p from this prior
  theta_ = dirichlet_sample([alp, bep]); p = theta_(1);
  % Now generate N Bernoulli points from this p
  X = double(rand(N,1) < p);
  sumX = sum(X);
  fprintf('p = %f, sumX = %d\n', p, sumX);

  % now compute the true posterior and the marginal
  th = evalPriorAt; % for brevity
  true_post = th.^(alp+sumX-1) .* (1-th).^(bep + N - sumX -1) / ...
              beta(alp+sumX, bep+N-sumX);
  marginal_prob = beta(alp + sumX, bep + N - sumX) / beta(alp, bep);

end

