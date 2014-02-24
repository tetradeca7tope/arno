% Estimate posterior using a grid search

% NUM_GRID_ITERS = 300; %NUM AL Iters
NUM_GRID_ITERS = NUM_SAMPLES_FOR_OTHER_EXPERIMENTS;

gr_kl_progress = zeros(NUM_GRID_ITERS, 1);

fprintf('Grid Iter: ');
for grid_iter = 1:NUM_GRID_ITERS

  fprintf('%d, ', grid_iter);
  % First select num_grid_pts uniformly from the parameter space.
  % First create a grid of N^2 points (N = ceil(sqrt(num_grid_pts))) and
  % then pick num_grid_pts points randomly from this set.
  num_grid_pts = NUM_INITIAL_PTS_PER_DIM^2 + grid_iter;
  num_pts_per_dim_in_cand_grid = ceil( sqrt(num_grid_pts) );
  cand1 = linspace(0,1, num_pts_per_dim_in_cand_grid + 2);
  cand1 = cand1(2:end-1);
  [C1, C2] = meshgrid(cand1, cand1);
  cand_grid_pts = [C1(:), C2(:)];
  % Now randomly pick a set of num_grid_pts from this set.
  shuffle_order = randperm(num_pts_per_dim_in_cand_grid^2);
  cand_grid_pts = cand_grid_pts(shuffle_order, :);
  grid_pts = cand_grid_pts(1:num_grid_pts, :);
  grid_joint_log_probs = evalLogJointProbs(grid_pts(:,1), grid_pts(:,2), sumX);

  % Create the hyper-params struct
  bf_hyper_params = hyper_params;
  bf_hyper_params.noise = NOISE_LEVEL * ones(num_grid_pts, 1);
  [gr_est_joint_log_probs, K, ~, ~] = GPKFoldCV(grid_pts, ...
    grid_joint_log_probs, th, NUM_KFOLDCV_PARTITIONS, candidates, ...
    bf_hyper_params);
  gr_est_joint_probs = exp(gr_est_joint_log_probs);
  % estimate prob(X)
  gr_est_pX = numerical_2D_integration_wrap(th, gr_est_joint_probs);
  gr_est_log_pX = log(gr_est_pX);
  % Estimate the posterior
  gr_est_log_post = gr_est_joint_log_probs - gr_est_log_pX;
  gr_est_post = exp(gr_est_log_post);
  % compute the KL divergence between the estimate and the true posterior
  gr_kl_progress(grid_iter) = estimate_2D_KL(th, log(true_post), ...
                                             gr_est_log_post);

end

fprintf('\n');
figure(fig_kl_progress);
plot(log(gr_kl_progress), 'g+-');
