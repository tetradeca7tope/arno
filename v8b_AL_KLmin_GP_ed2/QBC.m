% Use Active Learning - Uncertainty reduction

% Initiate Uncertainy Reduction
uc_curr_pts = initial_pts;
uc_kl_progress = zeros(NUM_AL_ITERS, 1);

% Now run Uncertainty reduction
for uc_iter = 1 : NUM_AL_ITERS

  % prelims
  uc_num_curr_pts = size(uc_curr_pts, 1);

  % First evaluate the joint
  uc_obs_joint_log_probs = evalLogJointProbs(uc_curr_pts, sumX);

  % Run GP
  % ======
  % First prep structure for LOO CV
  if USE_OPTIMAL_BANDWIDTH
    candidates.sigmaSmVals = optimal_bandwidth;
    candidates.sigmaPrVals = optimal_scale;
  else
    candidates.sigmaSmVals = gp_sm_param_ratios / sqrt(uc_num_curr_pts);
    candidates.sigmaPrVals = gp_scale_params;
  end
%   if uc_iter > 15 + randi(5)
%     candidates.sigmaSmVals = gp_sm_param_ratios / sqrt(uc_num_curr_pts);
%     candidates.sigmaPrVals = gp_scale_params;
%   end
  % Specify how to decide mean fn for GP
  if strcmp(MEAN_FUNC, 'logbarrier')
    [a, b, c] = fitLogBarrier(uc_curr_pts, uc_obs_joint_log_probs);
    uc_hyper_params.meanFunc = @(arg) (a*log(arg) + b*log(1 - arg) + c);
  elseif strcmp(MEAN_FUNC, 'lowest')
    uc_hyper_params.meanFunc = @(arg) LOWEST_LOGLIKL_VAL;
  else
    uc_hyper_params.meanFunc = @(arg) mean(uc_obs_joint_log_probs);
  end
  uc_hyper_params.noise = NOISE_LEVEL * ones(uc_num_curr_pts, 1);
  % Run GP with K-fold CV
  [uc_est_joint_log_probs, ~, ucSigmaSmOpt, ucSigmaPrOpt] = GPKFoldCV( ...
    uc_curr_pts, uc_obs_joint_log_probs, th, NUM_KFOLDCV_PARTITIONS, ...
    candidates, uc_hyper_params);
  uc_est_joint_probs = exp(uc_est_joint_log_probs);
  % estimate P(X)
  uc_est_pX = numerical_1D_integration(th, uc_est_joint_probs);
  uc_est_post = uc_est_joint_probs / uc_est_pX;

  % Compute the KL with the true estimate
  uc_ptwise_kl = true_post .* (log(true_post) - log(uc_est_post));
  uc_true_kl = numerical_1D_integration(th, uc_ptwise_kl);
  uc_kl_progress(uc_iter) = uc_true_kl;

  % ACTIVE LEARNING
  % ===============
  % Randomly pick points
  uc_candidates = rand(NUM_AL_CANDIDATES, 1);
  % constrain them to like in [delta, 1-delta]
  uc_candidates = BORDER_TOL + (1 - 2*BORDER_TOL)*uc_candidates;

  % Now pick the point with the highest uncertainty
  % After exponentiating we have a log normal density so the variance should
  % be picked according to the log normal density.
  % 1. Evaluate the predictions at the candidates

%   uc_opt_hps.noise = uc_hyper_params.noise;
%   uc_opt_hps.sigmaSm = ucSigmaSmOpt;
%   uc_opt_hps.sigmaPr = ucSigmaPrOpt;
%   uc_opt_hps.meanFunc = uc_hyper_params.meanFunc;
%   [ucM, ~, uc_cand_K] = GPRegression(uc_curr_pts, uc_obs_joint_log_probs, ...
%     uc_candidates, uc_opt_hps);
%   ucS = diag(uc_cand_K);
%   % Evaluate the objective (the log normal uncertainty)
%   uncert = (exp(ucS) - 1) .* exp(2*ucM + ucS);
%   [~, uc_max_idx] = max(uncert);
%   uc_max_pt = uc_candidates(uc_max_idx);
%   uc_curr_pts = [uc_curr_pts; uc_max_pt];

  uc_opt_hps.noise = uc_hyper_params.noise;
  uc_opt_hps.sigmaSm = ucSigmaSmOpt;
  uc_opt_hps.sigmaPr = ucSigmaPrOpt;
  uc_opt_hps.meanFunc = uc_hyper_params.meanFunc;
  [ucM, ~, uc_cand_K] = GPRegression(uc_curr_pts, uc_obs_joint_log_probs, ...
    uc_candidates, uc_opt_hps);
  curr_gp_samples = GPDrawSamples(ucM, uc_cand_K, NUM_SAMPLE_FUNCTIONS)';
  max_samples = max(curr_gp_samples, [], 2);
  min_samples = min(curr_gp_samples, [], 2);
  diff_samples = exp((max_samples - min_samples).^2) .* exp(2*ucM);
%   diff_samples = exp(max_samples) - exp(min_samples);

  [~, uc_max_idx] = max(diff_samples);
  uc_max_pt = uc_candidates(uc_max_idx);
  uc_curr_pts = [uc_curr_pts; uc_max_pt];
  

  % PLOT & Print Relevant results
  % =============================
  PLOT_OK_LOCAL = false;
  if PLOT_OK_LOCAL
    figure(fig_post);
    plot(th, uc_est_post, 'r--'); hold on,
    plot(th, true_post, 'g');
    plot(uc_candidates, uncert * max(true_post) / max(uncert), 'ko');
    hold off,
    text(uc_max_pt, -1, num2str(uc_iter), 'Color', 'k');
    axis([0,1, -2, 1.05*max([uc_est_post; true_post])]);
    titlestr = sprintf(...
      ['Estimated Post(r--) & True Post (g-)\n', ...
       'Sm: %f, Pr: %.1f, true_KL %0.4f'], ...
      ucSigmaSmOpt, ucSigmaPrOpt, uc_true_kl);
    title(titlestr);
  end
  fprintf('uc_iter: %d, KL: %f, uc-max-pt: %f, Sm: %f, Pr: %f\n', ...
          uc_iter, uc_true_kl, uc_max_pt, ucSigmaSmOpt, ucSigmaPrOpt);
  fprintf('Smvals: %s, Prvals: %s\n', ...
    mat2str([candidates.sigmaSmVals(1), candidates.sigmaSmVals(end)]), ...
    mat2str([candidates.sigmaPrVals(1), candidates.sigmaPrVals(end)]) );

end

% figure(fig_klprogress);
% plot(log(uc_kl_progress), 'k*-'); hold on,
