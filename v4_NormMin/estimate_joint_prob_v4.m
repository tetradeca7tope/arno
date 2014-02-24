% This function obtains some noisy estimates of the likelihood function
% The data is assumed to come from a bernoulli distribution

% the points at which the likelihood function should be evaluated is
% assumed to be stored in the vector cand_pts

% The values eventually returned (stored in the vector estJProbs) is
% renormalized so that max(estJProbs) = 1 so as to avoid any precision
% issues.

logProbs = sumX*log(cand_pts) + (N - sumX)*log(1 - cand_pts) + ...
             log(cand_pts.^(alp-1) .* (1 - cand_pts).^(bep-1) / ...
                 beta(alp, bep) );

% Now add some noise to this.
estLogProbs = logProbs + sqrt(noise_level)*randn(size(logProbs));

% Now deduct the maximum.
estNormLogProbs = estLogProbs - max(estLogProbs);
estNormProbs = exp(estNormLogProbs);

% Plot them out
% figure;
plot(cand_pts, exp(logProbs - max(logProbs)), 'rx-', 'MarkerSize', 7); hold on
plot(cand_pts, estNormProbs, 'bx-');
title('Probs (r) and Probs with noise(b)');

