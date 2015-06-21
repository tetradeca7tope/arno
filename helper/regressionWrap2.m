function logJointEst = regressionWrap2(pts, vals, noise_level, ...
  lowest_loglikl_val, loglikl_range, cv_cost_func, bounds)
% A wrapper which returns a function handle after performing regression on the
% points selected by any algorithm.

  if ~exist('cv_cost_func', 'var')
    cv_cost_func = @(y1, y2) (exp(y1) - exp(y2)).^2;
  end

  % prelims
  num_pts = size(pts, 1);
  num_dims = size(pts, 2);

  if ~exist('bounds', 'var') | isempty(bounds)
    bounds = repmat([-inf inf], num_dims, 1);
  end
  
  % Set up hyper params
  cv_candidates.sigmaSmVals = 0.05*num_pts^(-.1);
  cv_candidates.sigmaPrVals = [32]' * loglikl_range; % 2 1 0.5]' * loglikl_range;
  hyper_params.noise = noise_level * ones(num_pts, 1);
  hyper_params.meanFunc = @(arg) lowest_loglikl_val;
  hyper_params.costFunc = cv_cost_func;

  % Create a dummy point to obtain the hyper params
  dummy_pt = zeros(1, num_dims);
%   [est_log_probs, ~, sigmaSmOpt, sigmaPrOpt] = GPKFoldCV(pts, vals, ...
%     dummy_pt, 2, cv_candidates, hyper_params);
  opt_hyper_params = hyper_params;
  sigmaSmOpt = 0.05 * num_pts^(-.1);
  sigmaPrOpt = 32 * loglikl_range;
  opt_hyper_params.sigmaSm = sigmaSmOpt;
  opt_hyper_params.sigmaPr = sigmaPrOpt;
  runTimeParams.retFunc = true;
  [~, ~, ~, logJointEst] = GPRegression(pts, vals, dummy_pt, ...
    opt_hyper_params, runTimeParams);
  fprintf('sm: %0.4f, pr:%0.4f, ', sigmaSmOpt, sigmaPrOpt);
  fprintf('SmVals: (%.4f, %.4f), PrVals: (%.4f, %.4f)\n', ...
    cv_candidates.sigmaSmVals(1), cv_candidates.sigmaSmVals(end), ...
    cv_candidates.sigmaPrVals(1), cv_candidates.sigmaPrVals(end));

  % Finally modify the estimate to return -inf outside the domain
  logJointEst = @(t) boundedLogJointEst(t, logJointEst, bounds, ...
    lowest_loglikl_val);
  
end

function logJointProbs = boundedLogJointEst(t, logJointEst, bounds, ...
  lowestLogLiklVal)

  less = @(x,y) x<y;
  great = @(x,y) x>y;
  belowDomain = sparse(bsxfun(less, t , bounds(:,1)'));
  aboveDomain = sparse(bsxfun(great, t, bounds(:,2)'));
  outOfDomain = sum(belowDomain + aboveDomain, 2) > 0;
  
  % Create params for returning
  numPts = size(t, 1);
  logJointProbs = -inf * ones(numPts, 1);
  % Now compute at the points within the bounds
  [inBoundLJPs] = logJointEst(t(~outOfDomain, :));
  logJointProbs(~outOfDomain, :) = max(inBoundLJPs, lowestLogLiklVal);

end

