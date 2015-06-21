% Use MCMC to generate samples and then perform KDE

PROPOSAL_STD = 1/10;
NUM_MCMC_QUERIES = 1000 + NUM_INITIAL_PTS_PER_DIM^2;
MCMCD_INIT_PT = [0.5, 0.5];

evalLogJoint = @(evalPts) evalLogJointProbs(evalPts(:,1), evalPts(:,2), sumX);

mcmcd_collected_samples = MCMCD_INIT_PT;
mcmc_pt = MCMCD_INIT_PT;
mcmc_log_joint = evalLogJoint(mcmc_pt);
mcmcd_kl_progress = zeros(NUM_MCMC_QUERIES, 1);

for mcmcd_iter = 1:NUM_MCMC_QUERIES

%   fprintf('mcmc-query: %d, ', mcmcd_iter);
  [mcmc_pt, mcmc_log_joint] = Custom2DMCMC(PROPOSAL_STD, mcmc_pt, ...
    mcmc_log_joint, evalLogJoint);
  mcmcd_collected_samples = [mcmcd_collected_samples; mcmc_pt];

  % Now perform a KDE on these points
  if size(size(unique(mcmcd_collected_samples(:,1)), 1) < 2)
    kde_samples = [mcmcd_collected_samples; 0.5, 0.5];
  else
    kde_samples = mcmcd_collected_samples;
  end
  [~, mcmcd_kde] = kde(kde_samples);
  [mcmcd_est_post] = mcmcd_kde(th);
  curr_kl = estimate_2D_KL(th, log(true_post), log(mcmcd_est_post));
  mcmcd_kl_progress(mcmcd_iter) = curr_kl; 
  if mod(mcmcd_iter, 100) == 0
    fprintf('mcmcd_iter: %d, pt: (%0.4f, %0.4f), KL: %.5f\n', ...
      mcmcd_iter, mcmc_pt(1), mcmc_pt(2), curr_kl);
  end

  PLOT_OK_LOCAL = false;% true;
  if PLOT_OK_LOCAL
    figure(fig_mcmc_samples);
    Est_post = reshape(mcmcd_est_post, RESOLUTION_PER_DIM, RESOLUTION_PER_DIM);
    PLOT_FUNC(Th1, Th2, Est_post, 'Color','m'); xlabel('th1'); ylabel('th2');
    hold on,
    PLOT_FUNC(Th1, Th2, True_post, 'Color','g'); xlabel('th1'); ylabel('th2');
    plot(mcmcd_collected_samples(:,1), mcmcd_collected_samples(:,2), 'ko');
    title('MCMC - KDE');
    hold off,
  end

end

mcmcd_kl_progress = mcmcd_kl_progress(NUM_INITIAL_PTS_PER_DIM^2 +1 : end);

% Plot them out
% figure(fig_kl_progress);
% plot(log(mcmcd_kl_progress), 'm-*');

