function logJointEst = regressionWrapLRG2(pts, vals, noise_level, ...
  lowest_loglikl_val, loglikl_range, cv_cost_func)
% A wrapper which returns a function handle after performing regression on the
% points selected by any algorithm.

  if ~exist('cv_cost_func', 'var')
    cv_cost_func = @(y1, y2) (exp(y1) - exp(y2)).^2;
  end

  % prelims
  num_pts = size(pts, 1);
  num_dims = size(pts, 2);
  
  % Set up hyper params
  cv_candidates.sigmaSmVals = logspace(-2.6, 1.3, 10)' * (1/num_pts)^(1/10);
  cv_candidates.sigmaPrVals = [0.1 0.2 0.5 1 2 5 10] * loglikl_range;
  hyper_params.noise = noise_level * ones(num_pts, 1);
  hyper_params.meanFunc = @(arg) log( mean(exp(vals))); % lowest_loglikl_val;
  hyper_params.costFunc = cv_cost_func;

  % Create a dummy point to obtain the hyper params
  dummy_pt = zeros(1, num_dims);
  [est_log_probs, ~, sigmaSmOpt, sigmaPrOpt] = GPKFoldCV(pts, vals, ...
    dummy_pt, 10, cv_candidates, hyper_params);
  opt_hyper_params = hyper_params;
  opt_hyper_params.sigmaSm = sigmaSmOpt;
  opt_hyper_params.sigmaPr = sigmaPrOpt;
  runTimeParams.retFunc = true;
  [~, ~, ~, logJointEst] = GPRegression(pts, vals, dummy_pt, ...
    opt_hyper_params, runTimeParams);
  fprintf('sm: %0.4f, pr:%0.4f, ', sigmaSmOpt, sigmaPrOpt);
  fprintf('SmVals: (%.4f, %.4f), PrVals: (%.4f, %.4f)\n', ...
    cv_candidates.sigmaSmVals(1), cv_candidates.sigmaSmVals(end), ...
    cv_candidates.sigmaPrVals(1), cv_candidates.sigmaPrVals(end));

end

