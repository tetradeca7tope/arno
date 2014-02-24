% This script uses MCMC to generate some data and uses the samples generated
% for comparison
% The script assumes that the script v8 has already been run previously

% Set parameters for MCMC
PROPOSAL_STD = 3.0;
NUM_MCMC_RUNS = 1;
NUM_MCMC_SAMPLES = NUM_AL_ITERS;
% Other parameters
PARAM_SPACE_LOWER_THRESH = 0.02;
PARAM_SPACE_UPPER_THRESH = 1 - PARAM_SPACE_LOWER_THRESH;

% figures
fig_mcmc_samples = figure;
fig_mcmc_post = figure;
plot(th, true_post, 'g-', 'LineWidth', 2);

% Prelims
% set the funciton handle for querying the joint probability 
evalLogJoint = @(evalPts) evalLogJointProbs(evalPts, sumX);

kl_after_each_sample = zeros(NUM_MCMC_SAMPLES, 1);
num_valid_pts = zeros(NUM_MCMC_SAMPLES, 1);

for mcmc_run_iter = 1: NUM_MCMC_RUNS

  fprintf('MCMC-run : %d\n', mcmc_run_iter);
  % As the initial point pick the point with the highest probability among
  % initial_pts and initial_probs
  initial_obs_log_joints = evalLogJoint(initial_pts);
  [~, init_idx] = max(initial_obs_log_joints);
  mcmc_init_pt = initial_pts(init_idx);
  mcmc_init_log_joint_prob = initial_obs_log_joints(init_idx);

  all_pts = initial_pts([1:init_idx-1, init_idx+1:end, init_idx]);
  all_log_joints = initial_obs_log_joints(...
    [1:init_idx-1, init_idx+1:end, init_idx]);

  % Initialize the MCMC loop
  mcmc_pt = mcmc_init_pt;
  mcmc_log_joint = mcmc_init_log_joint_prob;
  mcmc_pts = mcmc_pt;
  % Run MCMC for NUM_MCMC_SAMPLES iterations
  for sample_iter = 1:NUM_MCMC_SAMPLES
    fprintf(' iter: %d, ', sample_iter);
    [mcmc_pt, mcmc_log_joint, new_pt, new_log_joint] = ...
      CustomMCMC(1, PROPOSAL_STD, mcmc_pt, mcmc_log_joint, evalLogJoint);
    if ~(new_pt > PARAM_SPACE_UPPER_THRESH || new_pt < PARAM_SPACE_LOWER_THRESH)
    % don't add pts close to the boundary. screws the regression up for some
    % reason. TODO: work on this
      all_pts = [all_pts; new_pt];
      all_log_joints = [all_log_joints; new_log_joint];
    end
    mcmc_pts = [mcmc_pts; mcmc_pt];

    % Now Run GP Regression on the points selected by MCMC
%     mcmc_candidates.sigmaSmVals = gp_smoothness_params;
%     mcmc_candidates.sigmaPrVals = gp_scale_params;
    hyper_params.noise = NOISE_LEVEL * ones(size(all_pts,1), 1);;
    [est_joint_log_probs, K, gpSmOpt, gpPrOpt] = GPKFoldCV(all_pts, ...
      all_log_joints, th, NUM_KFOLDCV_PARTITIONS, candidates, hyper_params);
    est_joint_probs = exp(est_joint_log_probs);
    est_pX = numerical_1D_integration(th, est_joint_probs);
    est_log_pX = log(est_pX);
    % Estimate the posterior
    est_log_post = est_joint_log_probs - est_log_pX;
    est_post = exp(est_log_post);
    % compute the KL divergence between the estimate and the true posterior
    ptwise_KL = true_post .* (log(true_post) - est_log_post);
    true_KL = numerical_1D_integration(th, ptwise_KL);
    % finally accumulate the KL divergences
    if true_KL < inf
      kl_after_each_sample(sample_iter) = kl_after_each_sample(sample_iter) + ...
                                          true_KL;
      num_valid_pts(sample_iter) = num_valid_pts(sample_iter) + 1;
    end
    fprintf('mcmc_pt: %f, new_pt: %f, KL: %f\n', ...
      mcmc_pt, new_pt, true_KL);

    PLOT_OK = 0;
    if PLOT_OK
      figure(fig_mcmc_samples);
      plot(all_pts, all_log_joints, 'ko'); hold on,
      plot(mcmc_pts, evalLogJointProbs(mcmc_pts, sumX), 'rx', 'MarkerSize', 10);
      plot(th, est_joint_log_probs, 'b'); hold off,
      figure(fig_mcmc_post);
      plot(th, est_post, 'k--'); 
      pause,
    end

  end
end

% Average the KLs out
kl_after_each_sample = kl_after_each_sample ./ num_valid_pts;

figure(fig_klprogress);
plot(log(kl_after_each_sample), 'rx-');
title('KL divergence after each AL iteration');

