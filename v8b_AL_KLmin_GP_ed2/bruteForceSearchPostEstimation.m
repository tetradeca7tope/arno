% This m-file estimates the posterior using a brute force search
gr_kl_progress = zeros(NUM_AL_ITERS,1);
fprintf('grid iter: ');

for grid_iter = 1:NUM_AL_ITERS
  num_grid_pts = NUM_INITIAL_PTS + grid_iter;
  if ~KNOW_LIKL_AT_BOUNDARY
    % Init uniformly
    grid_pts = linspace(0,1,num_grid_pts+2)';
    grid_pts = grid_pts(2:end-1);
  else
    % Use this initialization
    grid_pts = linspace(BOUNDARY_RIGHT, 1-BOUNDARY_RIGHT, num_grid_pts + 2)';
    grid_pts = grid_pts(2:end-1);
    grid_pts = [linspace(BOUNDARY_LEFT, BOUNDARY_RIGHT, NUM_BOUNDARY_PTS)'; ...
                grid_pts; linspace(1 - BOUNDARY_RIGHT, 1 - BOUNDARY_LEFT, ...
                                   NUM_BOUNDARY_PTS)'];
  end
  %%%%%
    grid_pts = rand(num_grid_pts, 1);
  % evaluate joint probabilities
  grid_joint_log_probs = evalLogJointProbs(grid_pts, sumX);
  % create the hyper-params struct
  bf_hyper_params.noise = NOISE_LEVEL * ones(size(grid_pts,1), 1);
  if strcmp(MEAN_FUNC, 'logbarrier')
    [a, b, c] = fitLogBarrier(grid_pts, grid_joint_log_probs);
    bf_hyper_params.meanFunc = @(arg) (a*log(arg) + b*log(1-arg) + c);
  elseif strcmp(MEAN_FUNC, 'lowest')
    bf_hyper_params.meanFunc = @(arg) LOWEST_LOGLIKL_VAL;
  else
    bf_hyper_params.meanFunc = @(arg) mean(grid_joint_log_probs);
  end
  if USE_OPTIMAL_BANDWIDTH
    candidates.sigmaSmVals = optimal_bandwidth;
    candidates.sigmaPrVals = optimal_scale;
  else
    candidates.sigmaSmVals = gp_sm_param_ratios / sqrt(size(all_pts, 1));
    candidates.sigmaPrVals = gp_scale_params;
  end
%   candidates.sigmaSmVals = gp_smoothness_params;
  [gr_est_joint_log_probs, K, sigmaSmOpt, sigmaPrOpt] = GPKFoldCV( grid_pts, ...
    grid_joint_log_probs, th, NUM_KFOLDCV_PARTITIONS, candidates, ...
    bf_hyper_params);
  gr_est_joint_probs = exp(gr_est_joint_log_probs);
  % 2.3 estimate prob(X)
  gr_est_pX = numerical_1D_integration(th, gr_est_joint_probs);
  gr_est_log_pX = log(gr_est_pX);
  % 2.4 Estimate the posterior
  gr_est_log_post = gr_est_joint_log_probs - gr_est_log_pX;
  gr_est_post = exp(gr_est_log_post);
  % compute the KL divergence between the estimate and the true posterior
  ptwise_KL = true_post .* (log(true_post) - gr_est_log_post);
  true_KL = numerical_1D_integration(th(2:end-1), ptwise_KL(2:end-1));
  gr_kl_progress(grid_iter) = true_KL;
  % Save results
  if SAVE_LOCAL_RESULTS
    gr_ptwise_KL = ptwise_KL;
    save_file_name = sprintf('%s/bf_iter_%d.mat', DEBUG_DIR_NAME, grid_iter);
    save(save_file_name, 'gr_est_joint_log_probs', 'gr_est_post', ...
                         'gr_ptwise_KL');
  end
  fprintf('iter: %d, KL: %f, sigmaSm: %f, sigmaPr: %f\n', grid_iter, ...
          true_KL, sigmaSmOpt, sigmaPrOpt);
end
fprintf('\n');

% figure(fig_klprogress);
% plot(log(gr_kl_progress), 'g+-');

