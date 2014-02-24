function [estLogProbs] = estimate_beta_joint_probs(pts, sumX, N, alp, bep, ...
                                                   noise_level, PLOT_OK)

  % This function obtains some noisy estimates of the likelihood function
  % The data is assumed to come from a bernoulli distribution

  % the points at which the likelihood function should be evaluated is
  % assumed to be stored in the vector pts

  % pts refers to the set of points on which the probabilities should be
  % evaluated. N, sumX refer to the total number of observations and the
  % total number of positive observations respectively. alp, bep refer to
  % the parameters of teh prior.

  if ~exist('PLOT_OK')
    PLOT_OK = 1;
  end

  logProbs = sumX*log(pts) + (N - sumX)*log(1 - pts) + ...
               log(pts.^(alp-1) .* (1 - pts).^(bep-1) / ...
                   beta(alp, bep) );

  % Now add some noise to this.
  estLogProbs = logProbs + sqrt(noise_level)*randn(size(logProbs));

  % Plot them out
  if PLOT_OK
    % figure;
    z_th = linspace(0,1);
    log_z_th = sumX * log(z_th) + (N - sumX)*log(1 - z_th) + ...
                 log(z_th.^(alp-1) .* (1 - z_th).^(bep-1) / ...
                     beta(alp, bep) );
    % plot(z_th, log_z_th, 'r--'); hold on,
    plot(pts, logProbs, 'rx');
    plot(pts, estLogProbs, 'bx');
    title('log Probs with noise(b) and true log Joing(r--)');
  end

end
