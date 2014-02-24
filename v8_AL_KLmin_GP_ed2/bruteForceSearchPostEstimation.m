% This m-file estimates the posterior using a brute force search
gr_kl_progress = zeros(NUM_AL_ITERS,1);

for grid_iter = 1:NUM_AL_ITERS
  num_grid_pts = NUM_INITIAL_PTS + grid_iter;
  grid_pts = linspace(0,1,num_grid_pts+2)'; grid_pts = grid_pts(2:end-1);
  grid_joint_log_probs = evalLogJointProbs(grid_pts, sumX);
  % create the hyper-params struct
  bf_hyper_params = hyper_params;
  bf_hyper_params.noise = NOISE_LEVEL * ones(num_grid_pts, 1);
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

end

figure(fig_klprogress);
plot(log(gr_kl_progress), 'g+-');

