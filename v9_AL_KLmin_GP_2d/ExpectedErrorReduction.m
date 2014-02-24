% Active Learning in the Expected Error Reduction Framework

% Initialize Active Learning
curr_pts = initial_pts; % at each iter the joint is evaluated at these points
KL_progress = zeros(NUM_AL_ITERS, 1);

% Now run Active Learning
for al_iter = 1:NUM_AL_ITERS

  % At each iteration perform the following
  % 1. Evaluate the Joint
  % 2. Run GP regression
  % 3. Generate samples and compute expected KL
  % 4. Pick the best point

  % Prelims
  num_curr_pts = size(curr_pts, 1);

  % 1. Evaluate the Joint probability
  % =================================
  obs_joint_log_probs = ...
    evalLogJointProbs(curr_pts(:,1), curr_pts(:,2), sumX);

  % 2. Run GP
  % =========
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
    hyper_params.meanFunc = @(arg) min(LOWEST_VAL, min(obs_joint_log_probs));
  else
    hyper_params.meanFunc = @(arg) mean(obs_joint_log_probs);
  end
  hyper_params.noise = NOISE_LEVEL * ones(num_curr_pts, 1);
  % GP Regression / K-FoldCV 
  [est_joint_log_probs, K, sigmaSmOpt, sigmaPrOpt] = GPKFoldCV(curr_pts, ...
    obs_joint_log_probs, th, NUM_KFOLDCV_PARTITIONS, candidates, hyper_params);
  est_joint_probs = exp(est_joint_log_probs);
  % estimate prob(X);
  est_pX = numerical_2D_integration_wrap(th, est_joint_probs);
  est_log_pX = log(est_pX);
  % Estimate the posterior
  est_log_post = est_joint_log_probs - est_log_pX;
  est_post = exp(est_log_post);
  % Compute the KL divergence between the estimate and the true posterior
  true_KL = estimate_2D_KL(th, log(true_post), est_log_post);
  KL_progress(al_iter) = true_KL;
  fprintf('Sm-cand: %s, Pr-cand:%s\n', ...
    mat2str([min(candidates.sigmaSmVals), max(candidates.sigmaSmVals)]), ...
    mat2str([min(candidates.sigmaPrVals), max(candidates.sigmaPrVals)]) );
  fprintf('Sm: %f, Pr:%f, KL: %f\n', sigmaSmOpt, sigmaPrOpt, true_KL);
  % Plot results
  Est_joint_log_probs = reshape(est_joint_log_probs, ...
    RESOLUTION_PER_DIM, RESOLUTION_PER_DIM);
  Est_post = reshape(est_post, RESOLUTION_PER_DIM, RESOLUTION_PER_DIM);
    % PLOT Intermediate Results
    % =========================
    PLOT_OK_LOCAL = false; 
    if PLOT_OK_LOCAL
      figure;
      PLOT_FUNC(Th1, Th2, Est_post, 'Color', 'm'); xlabel('th1'); ylabel('th2');
      figure;
      mesh(Th1, Th2, Est_joint_log_probs); hold on,
      plot3(curr_pts(:,1), curr_pts(:,2), obs_joint_log_probs, 'rx');
      xlabel('th1'); ylabel('th2');
      pause;
      close;
    end

  % 3. Generate Samples
  % ===================
  curr_gp_samples = GPDrawSamples(est_joint_log_probs,K, NUM_SAMPLE_FUNCTIONS)';
    % note the transpose. Each column is a sample.
  sample_log_pX = zeros(NUM_SAMPLE_FUNCTIONS, 1);
  for gp_sample_iter = 1:NUM_SAMPLE_FUNCTIONS
    sample_joint = exp(curr_gp_samples(:, gp_sample_iter));
    sample_pX = numerical_2D_integration_wrap(th, sample_joint);
    sample_log_pX(gp_sample_iter) = log(sample_pX);
  end
  sample_log_post = bsxfun(@minus, curr_gp_samples, sample_log_pX');
  sample_post = exp(sample_log_post);
%   sample_KLs = zeros(NUM_SAMPLE_FUNCTIONS, 1);
%   for gp_sample_iter = 1:NUM_SAMPLE_FUNCTIONS
%     sample_KLs(gp_sample_iter) = ...
%       estimate_2D_KL(th, sample_log_post(:,gp_sample_iter), est_log_post);
%   end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Compare the samples with the posterior mean
      PLOT_OK_LOCAL = false;
      if PLOT_OK_LOCAL 
        fpost = figure;
        flogjoint = figure;
        for gp_sample_iter = 1: NUM_SAMPLE_FUNCTIONS
          % Compare the estimated posterior and the sample posterior
          figure(fpost);
          Sample_post_curr = reshape(sample_post(:,gp_sample_iter), ...
                              RESOLUTION_PER_DIM, RESOLUTION_PER_DIM);
          PLOT_FUNC(Th1, Th2, Est_post, 'Color', 'm'); hold on,
          PLOT_FUNC(Th1, Th2, Sample_post_curr, 'Color', 'g'); hold off,
          title('True Log Joint (m) and Sample (g)');
          xlabel('th1'); ylabel('th2');
          % Compare the estimated log joint and the sample log joint
          figure(flogjoint);
          Sample_log_joint = reshape(curr_gp_samples(:,gp_sample_iter), ...
                              RESOLUTION_PER_DIM, RESOLUTION_PER_DIM);
          contour(Th1, Th2, Est_joint_log_probs, 'Color', 'm'); hold on,
          contour(Th1, Th2, Sample_log_joint, 'Color', 'g'); hold off,
          title('Estimated Log Joint (m) and Sample (g)');
          xlabel('th1'); ylabel('th2');
          pause;
        end
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Prepare iterations over samples
  % Create this struct we will need later.
  s_hyper_params.noise = [0; hyper_params.noise];
  s_hyper_params.sigmaSm = sigmaSmOpt;
  s_hyper_params.sigmaPr = sigmaPrOpt;
  s_hyper_params.meanFunc = hyper_params.meanFunc;

  % Here pick a set of random points as candidates for active learning
  al_candidates = rand(NUM_AL_CANDIDATES, 2);
  % constrain them to [delta, 1-delta] to avoid very large regressands
  al_candidates = BORDER_TOL + (1 - 2*BORDER_TOL)*al_candidates;
  accum_cand_KLs = zeros(NUM_AL_CANDIDATES, 1);
  
  % print out
  fprintf('Sample-iter: ');

  % 4. Now Pick the best point
  % ==========================
  % 4.1 First compute the expected KL if the likelihood were to be evaluated at
  % the candidate points.
  for sample_iter = 1:NUM_SAMPLE_FUNCTIONS
    fprintf(' %d,', sample_iter);
    % For the current sample, evaluate the value of the estimated log-joint
    % at each of the candidate points. Use bilinear interpolation. 
    Sample_log_joint = reshape(curr_gp_samples(:, sample_iter), ...
                            RESOLUTION_PER_DIM, RESOLUTION_PER_DIM);
    sample_cand_log_joint = interp2(Th1, Th2, Sample_log_joint, ...
                                    al_candidates(:,1), al_candidates(:,2));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compare the interpolated function and true function
        PLOT_OK_LOCAL = false;
        if PLOT_OK_LOCAL | PLOT_OK_GLOBAL
          figure;
          plot3(th(:,1), th(:,2), curr_gp_samples(:, sample_iter), 'mx');
          hold on;
          plot3(al_candidates(:,1), al_candidates(:,2), ...
                sample_cand_log_joint, 'bo');
          mesh(Th1, Th2, Sample_log_joint);
          xlabel('th1'); ylabel('th2');
          pause;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for cand_iter = 1:NUM_AL_CANDIDATES
      % first obtain the estimate for the current point
      sample_estimate = sample_cand_log_joint(cand_iter);
      % Now create hypothetical regressors and regressands
      reg_y = [sample_estimate; obs_joint_log_probs];
      reg_x = [al_candidates(cand_iter, :); curr_pts];
      % Now perform GP regression. Don't use K-fold CV
      s_joint_log_probs = GPRegression(reg_x, reg_y, th, s_hyper_params);

      % Now estimate the joint and then the posterior
      s_pX = numerical_2D_integration_wrap(th, exp(s_joint_log_probs));
      s_log_pX = log(s_pX);
      s_log_post = s_joint_log_probs - s_log_pX;
      s_post = exp(s_log_post);
      % Now compute the KL between this sample and the estimate
      curr_sample_KL = estimate_2D_KL(th, sample_log_post(:, sample_iter), ...
                                      s_log_post);
      accum_cand_KLs(cand_iter) = accum_cand_KLs(cand_iter) + curr_sample_KL;
    end % for cand_iter
  end % for sample_iter

  % 4.2 Now, locate the minimum
  expected_cand_KLs = accum_cand_KLs/ NUM_SAMPLE_FUNCTIONS;
  % Finally pick the minimum and add it to the list. Remove from the candidates
  [min_exp_cand_kl, min_idx] = min(expected_cand_KLs);
  min_pt = al_candidates(min_idx, :);
  old_pts = curr_pts;
  curr_pts = [curr_pts; min_pt];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot out the objective function for picking the point
    PLOT_OK_LOCAL = false;
    if PLOT_OK_LOCAL
      figure;
      plot3(al_candidates(:,1), al_candidates(:,2), expected_cand_KLs, 'bo');
      hold on,
      plot3(al_candidates(min_idx,1), al_candidates(min_idx,2), ...
            min_exp_cand_kl, 'rx');
      pause;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % PLOT & Print Relevant Results
  PLOT_OK_LOCAL = true;
  if PLOT_OK_LOCAL | PLOT_OK_GLOBAL
    figure(fig_post);
    PLOT_FUNC(Th1, Th2, True_post, 'Color', 'g'); hold on,
    PLOT_FUNC(Th1, Th2, Est_post, 'Color', 'r');
    % Mark the initial points with a x
    plot(initial_pts(:,1), initial_pts(:,2), 'bx');
    % indicate the points chosen via active learning in order
    for i = 1:al_iter
      text(curr_pts(NUM_INITIAL_PTS_PER_DIM^2 + i, 1), ...
           curr_pts(NUM_INITIAL_PTS_PER_DIM^2 + i, 2), ...
           num2str(i), 'Color', 'k'); 
    end
    titlestr = sprintf(...
       ['Estimated Post (r), True Posterior(g), \n', ...
        'Prior: (%.2f, %.2f, %.2f, %.2f), theta: (%f, %f), ', ...
        'num-initial-pts: %d^2\n', ...
        'Num-positive-pts/ Num-data-points: %d/%d\n', ...
        'Noise-level: %f, gp-smoothness: %f, gp-scale: %f\n', ...
        'True-KL: %f\n'], ...
        ALP1, BEP1, ALP2, BEP2, theta_true(1), theta_true(2), ...
        NUM_INITIAL_PTS_PER_DIM, ...
        sumX, NUM_DATA_SAMPLES, NOISE_LEVEL, sigmaSmOpt, sigmaPrOpt, ...
        true_KL);
    title(titlestr);
    xlabel('th1'); ylabel('th2');
    hold off
  end
  fprintf('\nAL-iter: %d, true-KL: %f, min-pt: (%f, %f)\n', al_iter, ...
          true_KL, min_pt(1), min_pt(2));
  
end % for al_iter

% Plot the progress of the KL with each AL iteration
figure(fig_kl_progress);
plot(log(KL_progress), 'bo-'); hold on,
