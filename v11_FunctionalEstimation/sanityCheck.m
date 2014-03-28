% A script to do some sanity checks on the code

close all;
num_samples = 150000;

% Do this for up to 5 randomly selected pairs of dimensions
plot_for = nchoosek(1:NUM_DIMS, 2);
plot_for = plot_for(randperm(size(plot_for, 1)), :);
plot_for = plot_for(1:min(5, size(plot_for, 1)), :);

  % Estimate the posterior
  cv_candidates.sigmaSmVals = logspace(-2, 2, 20)' * ...
                           1/(num_initial_pts + NUM_AL_ITERS)^(1/5);
  cv_candidates.sigmaPrVals = logspace(-1, 1, 10)' * LOGLIKL_RANGE;
  hyper_params.noise = NOISE_LEVEL * ones(num_initial_pts + NUM_AL_ITERS, 1);
  hyper_params.meanFunc = @(arg) LOWEST_LOGLIKL_VAL;
  hyper_params.costFunc = cv_cost_func;
  tic,
  [est_log_probs, ~, sigmaSmOpt, sigmaPrOpt] = GPKFoldCV(mbp_pts, ...
    mbp_log_probs, zeros(1, NUM_DIMS), 20, cv_candidates, hyper_params);
  toc,
  opt_hyper_params = hyper_params;
  opt_hyper_params.sigmaSm = sigmaSmOpt;
  opt_hyper_params.sigmaPr = sigmaPrOpt;
  logJointEst = @(arg) GPRegression(mbp_pts, mbp_log_probs, arg, ...
    opt_hyper_params);
  fprintf('sm: %0.4f, pr:%0.4f, ', sigmaSmOpt, sigmaPrOpt);
  fprintf('SmVals: (%.4f, %.4f), PrVals: (%.4f, %.4f)\n', ...
    cv_candidates.sigmaSmVals(1), cv_candidates.sigmaSmVals(end), ...
    cv_candidates.sigmaPrVals(1), cv_candidates.sigmaPrVals(end));
  % Now perform MCMC
  tic,
  estimate_mcmc_samples = CustomMCMC( ...
    num_samples + NUM_MCMC_BURNIN_FOR_EST, ...
    MCMC_EST_PROPOSAL_STD, MCMC_EST_INIT_PT, logJointEst);
  toc,
  Xmcmc = estimate_mcmc_samples(NUM_MCMC_BURNIN_FOR_EST+1:end, :);
 

for iter = 1:size(plot_for, 1)
  d1 = plot_for(iter, 1);
  d2 = plot_for(iter, 2);
  % First plot the points chosen in the selected order.
  figure(1);
  for j = 1:NUM_AL_ITERS
    text(mbp_pts(num_initial_pts + j, d1), mbp_pts(num_initial_pts + j, d2), ...
          int2str(j), 'Color', 'k');
    hold on;
  end
  plot(mbp_pts(1:num_initial_pts, d1), mbp_pts(1:num_initial_pts, d2), 'cx');
  axis([-1 2 -1 2]);
  xlabel(int2str(d1));
  ylabel(int2str(d2));
  hold off;

  % Plot Xmcmc
  figure(2);
  plot(Xmcmc(:,d1), Xmcmc(:,d2), 'ro');
  axis([-1 2 -1 2]);
  xlabel(int2str(d1));
  ylabel(int2str(d2));

  % Obtain samples from the true distro
  figure(3);
  X = gendata(num_samples, NUM_DIMS, sigma, p1);
  plot(X(:,d1), X(:,d2), 'bx');
  axis([-1 2 -1 2]);
  xlabel(int2str(d1));
  ylabel(int2str(d2));
 
  pause;
end

% Now plot the points out (if in 2D)
if NUM_DIMS == 1
  N = 200;
  th = linspace(PARAM_SPACE_BOUNDS(1), PARAM_SPACE_BOUNDS(2), N);
  ljth = evalLogJoint(th);
  esth = logJointEst(th);
  plot(th, ljth, 'g'); hold on;
  plot(th, esth, 'r');
end
if NUM_DIMS == 2
  figure;
  N = 100;
  t = linspace(PARAM_SPACE_BOUNDS(1), PARAM_SPACE_BOUNDS(2), N);
  [T1, T2] = meshgrid(t, t);
  th = [T1(:), T2(:)];
  ljth = evalLogJoint(th);
  LJTH = reshape(ljth, N, N);
  esth = logJointEst(th);
  ESTH = reshape(esth, N, N);
  contour(T1, T2, LJTH, 'Color', 'g'); hold on;
  contour(T1, T2, ESTH, 'Color', 'r');
end

