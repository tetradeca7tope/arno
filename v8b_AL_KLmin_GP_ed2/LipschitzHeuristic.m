% Uses the Lipshcitz Heuristic to compute 

LIPSCHITZ_CONSTANT = 1000;

lh_curr_pts = initial_pts;
lh_kl_progress = zeros(NUM_AL_ITERS, 1);

for lh_iter = 1:NUM_AL_ITERS

  fprintf('lh_iter: %d, ', lh_iter);

  % prelims
  lh_num_curr_pts = size(lh_curr_pts, 1);

  % Evaluate the joint
  lh_obs_joint_log_probs = evalLogJointProbs(lh_curr_pts, sumX);

  % Active Learning
  % Randomly pick points
  lh_candidates = rand(NUM_AL_CANDIDATES, 1);
  lh_candidates = BORDER_TOL + (1 - 2*BORDER_TOL)*lh_candidates;

  % Now obtain the upper & lower bounds on each of the lh_candidates
  pivots = repmat(lh_obs_joint_log_probs', NUM_AL_CANDIDATES, 1);
  diffs = abs( lh_candidates * ones(1, lh_num_curr_pts) - ...
               ones(NUM_AL_CANDIDATES, 1) * lh_curr_pts');
  upper_bounds = pivots + LIPSCHITZ_CONSTANT * diffs;
  lower_bounds = pivots - LIPSCHITZ_CONSTANT * diffs;
  upper_bounds = min(upper_bounds, [], 2);
  lower_bounds = max(lower_bounds, [], 2);

  % DEBUG
  PLOT_OK_LOCAL = false;
  if PLOT_OK_LOCAL
    figure;
    hold on;
    plot(th, evalLogJointProbs(th, sumX), 'b');
    plot(lh_curr_pts, lh_obs_joint_log_probs, 'rx');
    plot(lh_candidates, upper_bounds, 'ko');
    plot(lh_candidates, lower_bounds, 'ko');
    pause;
    close;
  end

  % Now choose points
  band = exp(upper_bounds) - exp(lower_bounds);
  [~, lh_max_idx] = max(band);
  lh_max_pt = lh_candidates(lh_max_idx, :);
  lh_curr_pts = [lh_curr_pts; lh_max_pt];

  % Perform regression on the estimate and compute the KL
  cv_cands.sigmaSmVals = logspace(-2, 2, 20)' * ...
                           (1/lh_num_curr_pts)^(1/5);
  cv_cands.sigmaPrVals = logspace(-1, 1, 10)' * ...
    (max(lh_obs_joint_log_probs) - min(lh_obs_joint_log_probs) );
  hyper_params.noise = NOISE_LEVEL * ones(lh_num_curr_pts, 1);
  hyper_params.meanFunc = []; %@(arg) mean(lh_obs_joint_log_probs);
  [lh_est_joint_log_probs] = GPKFoldCV( ...
    lh_curr_pts(1:end-1, :), lh_obs_joint_log_probs, th, ...
    NUM_KFOLDCV_PARTITIONS, cv_cands, hyper_params);
  lh_est_joint_probs = exp(lh_est_joint_log_probs);
  lh_est_pX = numerical_1D_integration(th, lh_est_joint_probs);
  lh_est_post = lh_est_joint_probs / lh_est_pX;
  % Compute the KL
  lh_ptwise_kl = true_post .* (log(true_post) - log(lh_est_post));
  lh_true_kl = numerical_1D_integration(th, lh_ptwise_kl);
  lh_kl_progress(lh_iter) = lh_true_kl;

  % Print some results here
  fprintf('Pt-Chosen: %0.4f, KL: %0.4f\n', lh_max_pt, lh_true_kl);
  
end

% Plot the results
PLOT_OK_LOCAL = true;
if PLOT_OK_LOCAL
  fig_lh_post = figure;
  hold on;
  plot(th, true_post, 'g');
  plot(th, lh_est_post, 'r');
  plot(lh_curr_pts(1:NUM_INITIAL_PTS, :), -0.2* ones(NUM_INITIAL_PTS, 1), ...
       'bx');
    for i = 1:NUM_AL_ITERS 
      text(lh_curr_pts(NUM_INITIAL_PTS+i, :), -0.2, num2str(i), 'Color', 'k');
    end
  hold off;

  fig_lh_kl = figure;
  plot(lh_kl_progress);
end
