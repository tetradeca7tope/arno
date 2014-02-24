% Use MCMC to generate samples and then perform KDE

PROPOSAL_STD = 2;
NUM_MCMC_QUERIES = S_RATIO*NUM_AL_ITERS + NUM_INITIAL_PTS;
MCMCD_INIT_PT = [0.2];

evalLogJoint = @(evalPts) evalLogJointProbs(evalPts(:,1), sumX);

% fig_mcmc_samples = figure;

mcmcd_collected_samples = MCMCD_INIT_PT;
mcmc_pt = MCMCD_INIT_PT;
mcmc_log_joint = evalLogJoint(mcmc_pt);
mcmcd_kl_progress = zeros(NUM_MCMC_QUERIES, 1);

fprintf('MCMC-iter: ');
for mcmcd_iter = 1:NUM_MCMC_QUERIES

%   fprintf('mcmc-query: %d\n', mcmcd_iter);
  [mcmc_pt, mcmc_log_joint, new_pt, new_log_joint] = ...
      CustomMCMC(1, PROPOSAL_STD, mcmc_pt, mcmc_log_joint, evalLogJoint);
  mcmcd_collected_samples = [mcmcd_collected_samples; mcmc_pt];

  if mod(mcmcd_iter, 100) == 0
    fprintf('%d, ', mcmcd_iter);
  end
  % Now perform a KDE on these points
  num_unique_pts = size(unique(mcmcd_collected_samples(:,1)), 1);
  if mcmcd_iter > NUM_INITIAL_PTS
    if num_unique_pts < 2
      kde_samples = [mcmcd_collected_samples; 0.5];
    else
      kde_samples = mcmcd_collected_samples;
    end
    [~, mcmcd_kde] = kde(kde_samples);
    [mcmcd_est_post] = mcmcd_kde(th);

    ptwise_KL = true_post .* (log(true_post) - log(mcmcd_est_post));
    curr_kl = numerical_1D_integration(th, ptwise_KL);
    mcmcd_kl_progress(mcmcd_iter) = curr_kl; 


    PLOT_OK_LOCAL = false;
%     PLOT_OK_LOCAL = true;
    if PLOT_OK_LOCAL
      fprintf('mcmcd_iter: %d, pt: (%0.4f), # unique pts: %d,  KL: %.5f\n', ...
        mcmcd_iter, mcmc_pt, num_unique_pts, curr_kl);
      figure(fig_mcmc_samples);
      plot(th, true_post, 'g'); hold on,
      plot(th, mcmcd_est_post, 'r'); hold off,
    end
  end

end
fprintf('\n');

% print
% num_unique_pts = size(unique(mcmcd_collected_samples(:,1)), 1);
% fprintf('Number of unque pts: %d\n\n', num_unique_pts);
mcmcd_kl_progress = mcmcd_kl_progress(NUM_INITIAL_PTS +1 : end);

% Plot them out
% figure(fig_kl_progress);
% plot(log(mcmcd_kl_progress), 'm-*');
