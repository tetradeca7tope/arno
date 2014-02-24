% Estimate posterior using a grid search

NUM_GRID_ITERS = NUM_AL_ITERS;

gr_kl_progress = zeros(NUM_GRID_ITERS, 1);
% fig_post = figure;

fprintf('Grid Iter: ');
for grid_iter = 1:NUM_GRID_ITERS

%   fprintf('%d, ', grid_iter);
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

  num_pts_dim_in_base_grid = floor(sqrt(num_grid_pts));
  cand1 = linspace(0,1, num_pts_dim_in_base_grid + 2);
  cand1 = cand1(2:end-1);
  [C1, C2] = meshgrid(cand1, cand1);
  grid_pts = [C1(:), C2(:)];
  num_additional_pts = num_grid_pts - num_pts_dim_in_base_grid^2;
  additional_pts = rand(num_additional_pts, 2);
  additional_pts = BORDER_TOL + (1 - 2*BORDER_TOL) * additional_pts;
  grid_pts = [grid_pts; additional_pts];

  % Just pick them randomly.
  grid_pts = rand(num_grid_pts, 2);

  % Evaluate Log jont prob
  grid_joint_log_probs = evalLogJointProbs(grid_pts(:,1), grid_pts(:,2), sumX);
  % Create the hyper-params struct
  if strcmp(MEAN_FUNC, 'lowest')
    bf_hyper_params.meanFunc = @(arg) min([LOWEST_VAL; grid_joint_log_probs]);
  else
    bf_hyper_params.meanFunc = @(arg) mean(grid_joint_log_probs);
  end
  if SET_PARAMETERS_TO_OPTIMAL
    candidates.sigmaSmVals = optimal_bandwidth;
    candidates.sigmaPrVals = optimal_scale;
  else
    candidates.sigmaSmVals = gp_sm_param_ratios / sqrt(num_grid_pts);
    candidates.sigmaPrVals = gp_scale_params;
  end
  bf_hyper_params.noise = NOISE_LEVEL * ones(num_grid_pts, 1);
  [gr_est_joint_log_probs, K, SmOpt, PrOpt] = GPKFoldCV(grid_pts, ...
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

  % Print the results our
  fprintf('Iter: %d, Sm: %f, Pr: %f KL: %f \n', grid_iter, SmOpt, PrOpt, ...
          gr_kl_progress(grid_iter));
  fprintf('Smvals: %s, Prvals: %s\n', ...
    mat2str([candidates.sigmaSmVals(1), candidates.sigmaSmVals(end)]), ...
    mat2str([candidates.sigmaPrVals(1), candidates.sigmaPrVals(end)]) );

  PLOT_OK_LOCAL = false;
%   PLOT_OK_LOCAL = true;
  if PLOT_OK_LOCAL
    Gr_Est_Post = reshape(gr_est_post, RESOLUTION_PER_DIM, RESOLUTION_PER_DIM);
    figure(fig_post);
    PLOT_FUNC(Th1, Th2, True_post, 'Color', 'g'); hold on,
    PLOT_FUNC(Th1, Th2, Gr_Est_Post, 'Color', 'r');
    % Mark the initial points with a x
    plot(grid_pts(:,1), grid_pts(:,2), 'bx');
    hold off,
  end
end

fprintf('\n');
if PLOT_OK_KL_PROGRESS
  figure(fig_kl_progress);
  plot(log(gr_kl_progress), 'g+-');
end
