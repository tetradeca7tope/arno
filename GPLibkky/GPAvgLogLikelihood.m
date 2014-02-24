function [avgloglikl] = GPAvgLogLikelihood(mu, K, observation)
% Computes the average log-likelihood of an observation of a GP.

  n = size(mu,1);

  maxdiag = max(diag(K));
  TOL_RATIO = 1e-3;
  % remove elements whose diagonals are too close to 0.
  idxs_to_keep = ones(n,1);
  for i = 1:n
    if K(i,i)/maxdiag < TOL_RATIO
      idxs_to_keep(i) = 0;
    end
  end
  num_discarded_pts = n - sum(idxs_to_keep);
  idxs_to_keep = idxs_to_keep > 0;
  mu = mu(idxs_to_keep);
  K = K(idxs_to_keep, idxs_to_keep);
  observation = observation(idxs_to_keep);
  if num_discarded_pts > 0
%     fprintf('Discarding %d pts to compute avg-log-likl.\n', ...
%             num_discarded_pts);
  end
  num_pts = n - num_discarded_pts;

%   loglikl = (observation - mu)' * pinv(K) * (observation - mu) / 2;

%   loglikl = (observation - mu)' * (observation - mu) / 2;

  % obtain the positive eigevalues of K
%   eigvals = eig(K);
%   n = size(K,1);
%   eigprod = 1;
%   i = 1;
%   while (i <= n) && (eigvals(i)) > 0
%     eigprod = eigprod * eigvals(i);
%     i = i + 1;
%   end
%   loglikl = -log(2*pi*eigprod) - ...
%             (observation - mu)' * pinv(K) * (observation - mu) / 2;
             
  loglikl = - (n/2)*log(2*pi) - 0.5*sum(log(diag(K))) + ...
            - (observation - mu)' * ( (1./ diag(K)) .* (observation - mu)) / 2;
            
%   loglikl = - (n/2)*log(2*pi) - 0.5*log(det(K)) + ...
%             - (observation - mu)' * ( K \ (observation - mu)) / 2;

  avgloglikl = loglikl/num_pts;

end
