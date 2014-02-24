% Active Learning in the Uncertainty Reduction Framework

% Initialize Active Learning
uc_pts = initial_pts;
uc_kl_progress = zeros(NUM_AL_ITERS, 1);

% Now run Active Learning
for uc_iter = 1:NUM_AL_ITERS

  % Prelims
  num_uc_pts = size(uc_pts, 1);

  % 1. Evaluate the Joint Probability
  % =================================
  uc_obs_joint_log_probs = evalLogJointProbs(uc_pts(:,1), uc_pts(:,2), sumX);

  % Run GP
  % ======
  % Prepare candidates for K-fold CV
  if SET_PARAMETERS_TO_OPTIMAL
    candidates.sigmaSmVals = optimal_bandwidth;
    candidates.sigmaPrVals = optimal_scale;
  else
    candidates.sigmaSmVals = gp_sm_param_ratios / sqrt(num_uc_pts);
    candidates.sigmaPrVals = gp_scale_params;
  end
  % specify other hyper params
  if strcmp(MEAN_FUNC, 'lowest')
    uc_hyper_params.meanFunc = @(arg) min([LOWEST_VAL; uc_obs_joint_log_probs]);
  else
    uc_hyper_params.meanFunc = @(arg) mean(uc_obs_joint_log_probs);
  end
  uc_hyper_params.noise = NOISE_LEVEL * ones(num_uc_pts, 1);
  % GP Regression / K-Fold CV
  [uc_est_joint_log_probs, ~, ucSigmaSmOpt, ucSigmaPrOpt] = GPKFoldCV( ...
    uc_pts, uc_obs_joint_log_probs, th, NUM_KFOLDCV_PARTITIONS, candidates, ...
    uc_hyper_params);
  uc_est_joint_probs = exp(uc_est_joint_log_probs);
  % Estimate Prob(X)
  uc_est_pX = numerical_2D_integration_wrap(th, uc_est_joint_probs);
  uc_est_post = uc_est_joint_probs / uc_est_pX;
  Uc_Est_Post = reshape(uc_est_post, RESOLUTION_PER_DIM, RESOLUTION_PER_DIM);
  % Compute the KL divergence between the estimate and the true posterior
  true_kl = estimate_2D_KL(th, log(true_post), log(uc_est_post) );
  uc_kl_progress(uc_iter) = true_kl;

  % Now pick the next point
  uc_candidates = rand(NUM_AL_CANDIDATES, 2);
  % constrain them to lie in [delta, 1-delta]^2 
  uc_candidates = BORDER_TOL + (1 - 2*BORDER_TOL) * uc_candidates;
  
  % Now evaluate at GP at these points
  uc_hyper_params.sigmaSm = ucSigmaSmOpt;
  uc_hyper_params.sigmaPr = ucSigmaPrOpt;
  [ucM, ~, uc_cand_K] = GPRegression(uc_pts, uc_obs_joint_log_probs, ...
    uc_candidates, uc_hyper_params);
  ucS = diag(uc_cand_K);
  uncert = (exp(ucS) - 1) .* exp(2*ucM + ucS);
  [~, uc_max_idx] = max(uncert);
  uc_max_pt = uc_candidates(uc_max_idx, :);
  uc_pts = [uc_pts; uc_max_pt];

  % PLOT & PRINT RELEVANT RESULTS
  % =============================
  PLOT_OK_LOCAL = true;
%   PLOT_OK_LOCAL = false;
  if PLOT_OK_LOCAL
    if ~exist('fig_post', 'var'), fig_post = figure; end
    figure(fig_post);
    contour(Th1, Th2, True_post, 'Color', 'g'); hold on,
%     contour(Th1, Th2, Uc_Est_Post, 'Color', 'r');
    % Mark the initial points with a x
    plot(initial_pts(:,1), initial_pts(:,2), 'bx');
    % indicate the points chosen via active learning in order
    for i = 1:uc_iter
      text(uc_pts(NUM_INITIAL_PTS_PER_DIM^2 + i, 1), ...
           uc_pts(NUM_INITIAL_PTS_PER_DIM^2 + i, 2), ...
           num2str(i), 'Color', 'k');
    end
%     title('True posterior(g) vs estimated (r)');
    hold off
  end
  fprintf('uc_iter: %d, KL: %f, max-pt: (%0.3f,%0.3f), Sm: %f, Pr: %f', ...
          uc_iter, true_kl, uc_max_pt(:,1), uc_max_pt(:,2), ucSigmaSmOpt, ...
          ucSigmaPrOpt);
  fprintf('SmVals: (%.4f, %.4f), PrVals: (%.4f, %.4f)\n', ...
    candidates.sigmaSmVals(1), candidates.sigmaSmVals(end), ...
    candidates.sigmaPrVals(1), candidates.sigmaPrVals(end));

end

% Plot the results at the end
if PLOT_OK_KL_PROGRESS
  figure(fig_kl_progress);
  plot(log(uc_kl_progress), 'k-*');
  hold on,
end
