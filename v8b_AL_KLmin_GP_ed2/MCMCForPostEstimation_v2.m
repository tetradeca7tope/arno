% File that runs MCMC for estimation of posterior densities
% 2nd version
PLOT_OK = 0;

if PLOT_OK
  fig_mcmc_samples = figure;
  fig_mcmc_post = figure;
end

NUM_MCMC_RUNS = 1;
% NUM_MCMC_SAMPLES = NUM_AL_ITERS;
NUM_MCMC_SAMPLES = NUM_AL_ITERS + NUM_INITIAL_PTS;
% other constants
PROPOSAL_STD = 3;

% declare some required variables
evalLogJointMCMC = @(x) evalLogJointProbs(x, sumX);
kl_after_each_sample = zeros(NUM_MCMC_SAMPLES, 1);
num_valid_kl_accums = zeros(NUM_MCMC_SAMPLES, 1);

for run_iter = 1:NUM_MCMC_RUNS

  % set initial point for MCMC
  initial_joints = evalLogJointMCMC(initial_pts);
  [~, min_initial_idx] = min(initial_joints);
  mcmc_pt = initial_pts(min_initial_idx);
  mcmc_joint = initial_joints(min_initial_idx);

  all_mcmc_samples = initial_pts;
  all_mcmc_joints = initial_joints;
  accepted_points = mcmc_pt;
  accepted_joints = mcmc_joint;

  for sample_iter = 1:NUM_MCMC_SAMPLES
    [mcmc_pt, mcmc_joint, sample, sample_joint] = MCMCPickNextSample( ...
      mcmc_pt, mcmc_joint, evalLogJointMCMC, PROPOSAL_STD);
    if ~ (sample < 0.01 || sample > 0.99)
      all_mcmc_samples = [all_mcmc_samples; sample];
      all_mcmc_joints = [all_mcmc_joints; sample_joint];
      accepted_points = [accepted_points; mcmc_pt];
      accepted_joints = [accepted_joints; mcmc_joint];
    end

    % Now perform GP regression and estimate the posterior
    candidates.sigmaSmVals = .1 * gp_sm_param_ratios * std(all_mcmc_samples);
%     candidates.sigmaPrVals = gp_pr_param_ratios * std(all_mcmc_joints);
    candidates.sigmaPrVals = linspace(100,600,8)';
    hyper_params.noise = NOISE_LEVEL * ones( size(all_mcmc_samples, 1), 1);
    [mcmc_est_log_joints] = GPKFoldCV(all_mcmc_samples, all_mcmc_joints, th, ...
      NUM_KFOLDCV_PARTITIONS, candidates, hyper_params);
    mcmc_est_joints = exp(mcmc_est_log_joints);
    mcmc_probX = numerical_1D_integration(th, mcmc_est_joints);
    mcmc_est_logpost = mcmc_est_log_joints - log(mcmc_probX);
    mcmc_est_post = exp(mcmc_est_logpost);

    % compute the KL divergence
    ptwise_kl = true_post .* ( log(true_post) - mcmc_est_logpost);
    curr_kl = numerical_1D_integration(th, ptwise_kl);
    if abs(curr_kl) < inf
      kl_after_each_sample(sample_iter) = ...
        kl_after_each_sample(sample_iter) + curr_kl;
      num_valid_kl_accums(sample_iter) = num_valid_kl_accums(sample_iter) + 1;
    else
      fprintf('Ivalid KL: %f, run: %d, sample: %d\n', ...
              curr_kl, run_iter, sample_iter);
    end
    fprintf('run %d/ sample: %d:: mcmc_pt: %f, sample; %f, KL: %f\n', ...
      run_iter, sample_iter, mcmc_pt, sample, curr_kl);

  if PLOT_OK
    figure(fig_mcmc_samples);
    plot(all_mcmc_samples, all_mcmc_joints, 'bo'); hold on,
    plot(accepted_points, accepted_joints, 'rx');
    plot(th, mcmc_est_log_joints, 'b-'); hold off,
    figure(fig_mcmc_post);
    plot(th, mcmc_est_post, 'k');
    pause,
  end

  % end for (sample_iter)
  end
% end for (run_iter)
end

% average the KLs out
kl_after_each_sample = kl_after_each_sample ./ num_valid_kl_accums;

figure(fig_klprogress);
plot(log(kl_after_each_sample), 'm-x');
title('KL divergence after each iteration');

