% Expected Error Reduction Query strategy 

% Initiate Active Learning
curr_pts = initial_pts; % at each iter the joint is evaluated at these pts.
kl_progress = zeros(NUM_AL_ITERS, 1);

% Now run Active Learning
for al_iter = 1: NUM_AL_ITERS

  % We need to do the following at each Active Learning iteration
  % 1. Evaluate the joint
  % 2. Run the GP
  % 3. Generate Samples and compute expected KL
  % 4. Pick the best point

  % prelims
  num_curr_pts = size(curr_pts, 1);

  % prelims
%   gp_smoothness_params = gp_sm_param_ratios / sqrt(num_curr_pts);

  % 1. Evalute the joint
  % ====================
  [obs_joint_log_probs] = evalLogJointProbs(curr_pts, sumX);

  % 2. Run GP
  % =========
  % 2.1 Prep structure for LOO CV
%   candidates.sigmaSmVals = 0.1 * gp_sm_param_ratios * std(curr_pts);
  candidates.sigmaSmVals = gp_sm_param_ratios / sqrt(num_curr_pts);% gp_smoothness_params;
%   candidates.sigmaSmVals = gp_smoothness_params;
  candidates.sigmaSmVals = 0.037 * (30/al_iter)^(0.3);
  candidates.sigmaPrVals = gp_scale_params;
  % Specify how to decide mean function for the GP.
  if strcmp(MEAN_FUNC, 'logbarrier')
    % first fit a log barrier to the data
    [a, b, c] = fitLogBarrier(curr_pts, obs_joint_log_probs);
    hyper_params.meanFunc = @(arg) (a*log(arg) + b*log(1-arg) + c);
  elseif strcmp(MEAN_FUNC, 'lowest')
    hyper_params.meanFunc = @(arg) min(-600, min(obs_joint_log_probs));
  else
    hyper_params.meanFunc = @(arg) mean(obs_joint_log_probs);
  end
  hyper_params.noise = NOISE_LEVEL * ones(num_curr_pts, 1);
  % 2.2 Run GP with K-fold CV
  [est_joint_log_probs, K, sigmaSmOpt, sigmaPrOpt] = GPKFoldCV(curr_pts, ...
    obs_joint_log_probs, th, NUM_KFOLDCV_PARTITIONS, candidates, hyper_params);
  est_joint_probs = exp(est_joint_log_probs);
  % 2.3 estimate prob(X)
  est_pX = numerical_1D_integration(th, est_joint_probs);
  est_log_pX = log(est_pX);
  % 2.4 Estimate the posterior
  est_log_post = est_joint_log_probs - est_log_pX;
  est_post = exp(est_log_post);
  % compute the KL divergence between the estimate and the true posterior
  ptwise_KL = true_post .* (log(true_post) - est_log_post);
  true_KL = numerical_1D_integration(th(2:end-1), ptwise_KL(2:end-1));
  kl_progress(al_iter) = true_KL;

  % 3. Generate Samples and commpute Expected KL at each candidate pt
  % =================================================================
  % Use smaller bandwidth when drawing the samples
  sample_draw_hyper_params = hyper_params;
  if num_curr_pts < 101, divide_bandwidth_by = 2;
  else, divide_bandwidth_by = 1;
  end
  sample_draw_candidates.sigmaSmVals = [sigmaSmOpt/divide_bandwidth_by];
  sample_draw_candidates.sigmaPrVals = candidates.sigmaPrVals; 
  [sample_draw_est_joint_log_probs, sample_draw_K] = GPKFoldCV(curr_pts, ...
    obs_joint_log_probs, th, NUM_KFOLDCV_PARTITIONS, ...
    sample_draw_candidates, hyper_params);
  curr_gp_samples = GPDrawSamples(sample_draw_est_joint_log_probs, ...
    sample_draw_K, NUM_SAMPLE_FUNCTIONS)';
%   curr_gp_samples = GPDrawSamples(est_joint_log_probs, ...
%     K, NUM_SAMPLE_FUNCTIONS)';
    % note the transpose above. now each column is a sample
  sample_log_pX = log(numerical_1D_integration(th, exp(curr_gp_samples)));
  sample_log_post = bsxfun(@minus, curr_gp_samples, sample_log_pX);
  sample_post = exp(sample_log_post);
  pointwise_KL = bsxfun(@times, ...
                        bsxfun(@plus, -sample_log_post, est_log_post), ...
                        est_post);
  sample_KLs = numerical_1D_integration(th, pointwise_KL);
  estim_KL = mean(sample_KLs);

  % Prepare iterations over samples
  % create this struct which we will need later
  % add a 0-noise component to the struct 
  opt_hyper_params.noise = hyper_params.noise;
  opt_hyper_params.sigmaSm = sigmaSmOpt;
  opt_hyper_params.sigmaPr = sigmaPrOpt;
  opt_hyper_params.meanFunc = hyper_params.meanFunc;
  accum_cand_KLs = zeros(NUM_AL_CANDIDATES, 1);

  % Since we cannot perform an exact minimization we generate a bunch of
  % samples and pick the minimum.
  % first create a bunch of samples uniformly from [0, 1] as candidates for
  % the AL algorithm.
%   al_candidates = linspace(0, 1, NUM_AL_CANDIDATES + 2)';
%   al_candidates = al_candidates( 2:(end-1) );
%   al_candidates = al_candidates + ...
%                   (1/NUM_AL_CANDIDATES) * (rand(NUM_AL_CANDIDATES, 1) - 0.5);
  al_candidates = rand(NUM_AL_CANDIDATES, 1);
  % constrain them to [delta, 1-delta] to avoid very large regressands
  al_candidates = BORDER_TOL + (1 - 2*BORDER_TOL)*al_candidates;

%   PLOT_OK_LOCAL = true;
%   if PLOT_OK_LOCAL
%     figure(fig_gpreg);
%     plot(th, est_joint_log_probs); hold on
%     plot(al_candidates, al_cand_vals, 'm--');
%     plot(curr_pts, obs_joint_log_probs, 'rx', 'MarkerSize', 10);
%     title('Regression vs AL cands');
%     hold off,
%   end

  % 4. Now Pick the best point
  % ==========================
  % 4.1 First accumulate the KL divergences
  for cand_iter = 1:NUM_AL_CANDIDATES
    for sample_iter = 1:NUM_SAMPLE_FUNCTIONS

      sample_estimate = interp1(th, curr_gp_samples(:, sample_iter), ...
                                al_candidates(cand_iter));
      % Now create hypothetical regressors and regressands
      reg_y = [sample_estimate; obs_joint_log_probs];
      reg_x = [al_candidates(cand_iter); curr_pts];

      % For the hypothetical regression set the noise param to zero.
      opt_hyper_params.noise = [0; hyper_params.noise]; 
      s_joint_log_probs = GPRegression(reg_x, reg_y, th, opt_hyper_params);
      % Now estimate the joint and then the posterior
      s_pX = numerical_1D_integration(th, exp(s_joint_log_probs));
      s_log_pX = log(s_pX);
      s_log_post = s_joint_log_probs - s_log_pX;
      s_post = exp(s_log_post);
      % Now compute the KL between this sample and the estimate
      s_ptwise_KL = sample_post(:, sample_iter) .* ( ...
        sample_log_post(:, sample_iter) - s_log_post);
      curr_sample_KL = numerical_1D_integration(th(2:end-1), ...
                                                s_ptwise_KL(2:end-1));
      accum_cand_KLs(cand_iter) = accum_cand_KLs(cand_iter) + curr_sample_KL;
    end
  end
  expected_cand_KLs = accum_cand_KLs / NUM_SAMPLE_FUNCTIONS;
  % 4.2 Pick the minimum and add him to the list of current points
  [~, min_idx] = min(expected_cand_KLs);
  min_pt = al_candidates(min_idx);
  old_pts = curr_pts;
  curr_pts = [min_pt; curr_pts];

  curr_gp_samples = GPDrawSamples(ucM, uc_cand_K, NUM_SAMPLE_FUNCTIONS)';
  % Plot & Print Relevant results
  % ==============================
  % Plot the current estimate
  PLOT_OK_LOCAL = true;
  if PLOT_OK_LOCAL
    if ~exist(DEBUG_DIR_NAME, 'dir')
      mkdir(DEBUG_DIR_NAME);
    end
    figure(fig_post);
    plot(th, est_post, 'r--'); hold on,
    plot(th, true_post, 'g-');
    plot(al_candidates, expected_cand_KLs * 15/max(expected_cand_KLs), 'ko');
    hold off,
    text(min_pt, -1, num2str(al_iter), 'Color', 'k');
    axis([0 1 -2 1.05*max([est_post; true_post])]);
    titlestr = sprintf(...
               ['Estimated Post (r--) & True Posterior (g-)\n', ...
                'Prior: Beta(%d, %d, %d), theta: %f, Num initial pts: %d\n', ...
                'Num data points: %d, Num positive points = %d\n', ...
                'Noise level: %f, gp-bandwidth  = %f, gp-scale = %f\n', ...
                'True KL: %f, Estim KL (using GP): %f'], ...
                ALP, BEP, CCP, theta_true(1), NUM_INITIAL_PTS, ...
                NUM_BINOMIAL_SAMPLES, sumX, NOISE_LEVEL, ...
                sigmaSmOpt, sigmaPrOpt, true_KL, estim_KL);
    title(titlestr);
    save_file_name = sprintf('%s/iter%d_results', DEBUG_DIR_NAME, al_iter);
    saveas(fig_post, save_file_name, 'png');
    % Plot the samples drawn.
    figure(fig_samples);
    for i = 1:NUM_SAMPLE_FUNCTIONS
      plot(th, curr_gp_samples(:,i), 'm'); hold on,
    end
    % plot the points evaluated
    plot(old_pts, obs_joint_log_probs, 'kx', 'MarkerSize', 6);
    true_log_joint_probs = log(true_post * marginal_prob);
    plot(th, true_log_joint_probs, 'c');
    hold off,
  end % PLOT_OK

  % PLOT for debugging
  % ==================
  if SAVE_LOCAL_RESULTS 
    if ~exist(DEBUG_DIR_NAME, 'dir')
      mkdir(DEBUG_DIR_NAME);
    end
    log_joint_probs = log(true_post * marginal_prob);
    save_file_name = sprintf('%s/al_iter_%d.mat', DEBUG_DIR_NAME, al_iter);
    save(save_file_name, 'curr_gp_samples', 'est_joint_log_probs', ...
      'log_joint_probs', 'true_post', 'ptwise_KL', 'est_post');
  end

  % also print the results out.
  fprintf('iter: %d, true-KL: %f, min-pt: %f, Sm: %f, Pr:%f\n', ....
          al_iter, true_KL, min_pt, sigmaSmOpt, sigmaPrOpt);
end

% Finally once this is all done plot the final estimate out
fig_ALfinal = figure; plot(th, true_post, 'g-', 'LineWidth', 2); hold on,
plot(th, est_post, 'r--', 'LineWidth', 2);
title('true posterior (g-) vs final estimated posterior (r--)');
plot(log(kl_progress), 'bo-'); hold on,
