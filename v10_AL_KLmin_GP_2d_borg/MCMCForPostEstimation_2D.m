% A script that tests MCMC for Posterior Estimation for a 2D example
% Script v9.m should be run prior to this

% Set Parameters for MCMC
PROPOSAL_STD = 1/2.5; %2.0;
NUM_MCMC_RUNS = 1;
% NUM_MCMC_SAMPLES = 300; %NUM_AL_ITERS;
NUM_MCMC_SAMPLES = NUM_SAMPLES_FOR_OTHER_EXPERIMENTS;
% Other parameters to get GP regression working on MCMC properly
PARAM_SPACE_LOWER_THRESH = 0.02;
PARAM_SPACE_UPPER_THRESH = 1 - PARAM_SPACE_LOWER_THRESH;

% Figure handles
fig_mcmc_samples = figure;
fig_mcmc_post = figure;
PLOT_FUNC(Th1, Th2, Est_post, 'Color', 'g');

% Prelims
% set the function handle for querying the joint probability
evalLogJoint = @(evalPts) evalLogJointProbs(evalPts(:,1), evalPts(:,2), sumX);

kl_after_each_sample = zeros(NUM_MCMC_SAMPLES, 1);
num_valid_pts = zeros(NUM_MCMC_SAMPLES, 1);

for mcmc_run_iter = 1: NUM_MCMC_RUNS

  fprintf('MCMC_RUN_ITER: %d\n', mcmc_run_iter);
  % As the initial point for MCMC pick the point with the highest probability
  % among initial_pts and initial_probs
  initial_obs_log_joints = evalLogJoint(initial_pts);
  [~, init_idx] = max(initial_obs_log_joints);
  mcmc_init_pt = initial_pts(init_idx, :);
  mcmc_init_log_joint_prob = initial_obs_log_joints(init_idx);

  all_pts = initial_pts;
  all_log_joints = initial_obs_log_joints;
  
  % Initialize the MCMC loop
  mcmc_pt = mcmc_init_pt;
  mcmc_log_joint = mcmc_init_log_joint_prob;
  mcmc_pts = mcmc_pt;

  % Run MCMC for NUM_MCMC_SAMPLES iterations
  for sample_iter = 1:NUM_MCMC_SAMPLES
    fprintf(' iter: %d, ', sample_iter);
    [mcmc_pt, mcmc_log_joint, new_pt, new_log_joint] = ...
      Custom2DMCMC( PROPOSAL_STD, mcmc_pt, mcmc_log_joint, evalLogJoint);
    % add the points to the current collection
    % Discard points too close to the boundary as they screw the regression up.
    if ( (new_pt < PARAM_SPACE_UPPER_THRESH) & ...
         (new_pt > PARAM_SPACE_LOWER_THRESH) )
      all_pts = [all_pts; new_pt];
      all_log_joints = [all_log_joints; new_log_joint];
    end
    mcmc_pts = [mcmc_pts; mcmc_pt];

    % Now Run GP Regressionon the points selected by MCMC
    hyper_params.noise = NOISE_LEVEL * ones(size(all_pts, 1), 1);
    [est_joint_log_probs, K, gpSmOpt, gpPrOpt] = GPKFoldCV(all_pts, ...
      all_log_joints, th, NUM_KFOLDCV_PARTITIONS, candidates, hyper_params);
    est_joint_probs = exp(est_joint_log_probs);
    % estimate prob(X);
    est_pX = numerical_2D_integration_wrap(th, est_joint_probs);
    est_log_pX = log(est_pX);
    % Estimate the posterior
    est_log_post = est_joint_log_probs - est_log_pX;
    est_post = exp(est_log_post);
    % Compute the KL between the estimate and the true posterior
    curr_KL = estimate_2D_KL(th, log(true_post), est_log_post);
    kl_after_each_sample(sample_iter) = kl_after_each_sample(sample_iter) + ...
                                          curr_KL;
    fprintf('mcmc_pt: (%f, %f), new_pt: (%f, %f), KL: %f\n', ...
      mcmc_pt, new_pt, curr_KL);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT Results of MCMC after each iteration
    PLOT_OK_LOCAL = true;
    if PLOT_OK_LOCAL
      figure(fig_mcmc_samples);
      Est_post = reshape(est_post, RESOLUTION_PER_DIM, RESOLUTION_PER_DIM);
      PLOT_FUNC(Th1, Th2, Est_post, 'Color','m'); xlabel('th1'); ylabel('th2');
      hold on,
      PLOT_FUNC(Th1, Th2, True_post, 'Color','g'); xlabel('th1'); ylabel('th2');
      plot(all_pts(:,1), all_pts(:,2), 'ko');
      plot(mcmc_pts(:,1), mcmc_pts(:,2), 'rx');
      plot(initial_pts(:,1), initial_pts(:,2), 'bx');
      plot(mcmc_pt(:,1), mcmc_pt(:,2), 'rx', 'MarkerSize', 10, 'LineWidth', 3);
      plot(new_pt(:,1), new_pt(:,2), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
      titlestr = sprintf('Results of MCMC after %d points', sample_iter);
      title(titlestr);
      hold off,
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
end

% Average the results over the runs
kl_after_each_sample = kl_after_each_sample / NUM_MCMC_RUNS;
figure(fig_kl_progress);
plot(log(kl_after_each_sample), 'r-x'); hold on,
title('KL vs #-pts');

