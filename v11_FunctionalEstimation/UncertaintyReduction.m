% Uncertainty reduction query strategy for the problem

% Parameter for Active Learning
NUM_AL_CANDIDATES = 200;
NOISE_LEVEL = 0.005;

if INIT_ON_GRID
  if NUM_INIT_PTS_PER_DIM == 1
    init_pts_1dim = 0.5;
  else
    init_pts_1dim = linspace(0, 1, NUM_INIT_PTS_PER_DIM)';
  end
  initial_pts = gengrid(NUM_DIMS, init_pts_1dim);
else
  initial_pts = rand(NUM_INIT_PTS_PER_DIM-2, NUM_DIMS);
  initial_pts = [initial_pts; zeros(1, NUM_DIMS); ones(1, NUM_DIMS)];
end
num_initial_pts = size(initial_pts, 1);

% Initialize Active Learning
uc_pts = initial_pts;
% Error tracker for each functional
uc_err_prog.f1 = zeros(NUM_AL_ITERS, 1);
uc_err_prog.f2 = zeros(NUM_AL_ITERS, 1);
uc_err_prog.f3 = zeros(NUM_AL_ITERS, 1);
uc_err_prog.f4 = zeros(NUM_AL_ITERS, 1);

% Now run Active Learning
for uc_iter = 1:NUM_AL_ITERS

  % prelims
  num_uc_pts = size(uc_pts, 1);

  % 1. Evaluate Log Joint probability
  uc_obs_log_joint_probs = evalLogJoint(uc_pts);

  % al_candidates
  uc_candidates = -1 + 3*rand(NUM_AL_CANDIDATES, NUM_DIMS);

  % 2. Run GP Regression on the al candidates
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 2.1: Hyper parameters
  hyper_params.noise = NOISE_LEVEL * ones(num_uc_pts, 1);
  hyper_params.meanFunc = @(arg) LOWEST_LOGLIKL_VAL; 
  if USE_OPT_PARAMS
    cv_candidates.sigmaSmVals = OPT_BANDWIDTH;  % Use fixed bandwidth
    cv_candidates.sigmaPrVals = OPT_SCALE;  % Use fixed bandwidth
  else
    cv_candidates.sigmaSmVals = logspace(-2, 2, 20)' * ...
                             (num_initial_pts/num_uc_pts)^(1/5);
    cv_candidates.sigmaPrVals = logspace(-1, 1, 10)' * LOGLIKL_RANGE;
  end
%   cv_candidates.sigmaSmVals = OPT_BANDWIDTH;  % Use fixed bandwidth
  % GP Regression
  [ucM, ucK, ~, ~] = GPKFoldCV(uc_pts, ...
    uc_obs_log_joint_probs, uc_candidates, 20, cv_candidates, hyper_params);

  % 3. Now pick the best point
  ucS = diag(ucK);
  uncert = (exp(ucS) - 1) .* exp(2*ucM + ucS);
  [~, uc_max_idx] = max(uncert);
  uc_max_pt = uc_candidates(uc_max_idx, :);
  uc_pts = [uc_pts; uc_max_pt];

  % Print some of the stuff out
  fprintf('ITER: %d, Pt chosen: %s, \n==================\n', uc_iter, ...
          mat2str(uc_max_pt));

  % Estimating Functionals
  % ======================
  fprintf('Estimating Functionals : ...');
  % Now approximate the functionals using a sampling scheme
  % First prepare the function to be passed to the MCMC procedure
  % 1. First redo the GP regression to find the optimal hyper-params
  cv_candidates.sigmaSmVals = logspace(-2, 2, 20)' * ...
                           (1/num_uc_pts)^(1/5);
  cv_candidates.sigmaPrVals = logspace(-1, 1, 10)' * LOGLIKL_RANGE;
  [est_log_probs, ~, sigmaSmOpt, sigmaPrOpt] = GPKFoldCV(uc_pts(1:end-1, :), ...
    uc_obs_log_joint_probs, uc_candidates, 20, cv_candidates, ...
    hyper_params);
  opt_hyper_params = hyper_params;
  opt_hyper_params.sigmaSm = sigmaSmOpt;
  opt_hyper_params.sigmaPr = sigmaPrOpt;
  logJointEst = @(arg) GPRegression(uc_pts(1:end-1, :), ...
    uc_obs_log_joint_probs, arg, opt_hyper_params);
  fprintf('sm: %0.4f, pr:%0.4f, ', sigmaSmOpt, sigmaPrOpt);
  fprintf('SmVals: (%.4f, %.4f), PrVals: (%.4f, %.4f)\n', ...
    cv_candidates.sigmaSmVals(1), cv_candidates.sigmaSmVals(end), ...
    cv_candidates.sigmaPrVals(1), cv_candidates.sigmaPrVals(end));

  % 2. Now, perform MCMC
  estimate_mcmc_samples = CustomMCMC( ...
    NUM_MCMC_ITERS_FOR_EST + NUM_MCMC_BURNIN_FOR_EST, ...
    MCMC_EST_PROPOSAL_STD, MCMC_EST_INIT_PT, logJointEst);
  Xmcmc = estimate_mcmc_samples(NUM_MCMC_BURNIN_FOR_EST+1:end, :);
  % Now evaluate the functionals on Xmcmc
  T1 = f1(Xmcmc); uc_err_prog.f1(uc_iter) = abs(T1 - f_vals.S1)/abs(f_vals.S1);
  T2 = f2(Xmcmc); uc_err_prog.f2(uc_iter) = abs(T2 - f_vals.S2)/abs(f_vals.S2); 
  T3 = f3(Xmcmc); uc_err_prog.f3(uc_iter) = abs(T3 - f_vals.S3)/abs(f_vals.S3);
  T4 = f4(Xmcmc); uc_err_prog.f4(uc_iter) = abs(T4 - f_vals.S4)/abs(f_vals.S4);
  fprintf('emp-errors: f1: %f, f2: %f, f3: %f, f4: %f\n', ...
          uc_err_prog.f1(uc_iter), ...
          uc_err_prog.f2(uc_iter), ...
          uc_err_prog.f3(uc_iter), ...
          uc_err_prog.f4(uc_iter) );

  % Plot out some of the results
  PLOT_OK_LOCAL = false;
  if PLOT_OK_LOCAL
    figure;
    hold on;
      if NUM_DIMS == 1
        for i = 1:uc_iter
          text(uc_pts(NUM_INIT_PTS_PER_DIM^NUM_DIMS + i, 1), -0.1, ...
               num2str(i), 'Color', 'k');
        end
        plot(initial_pts, -0.2*ones(NUM_INIT_PTS_PER_DIM, 1), ...
             'bx', 'MarkerSize', 10);
        th = linspace(-1,2,100)';
        plot(th, exp(logJointEst(th)), 'g');
        plot(th, exp(evalLogJoint(th)), 'r');
        plot(uc_candidates, uncert/max(uncert), 'ko');
      else
        for i = 1:uc_iter
          text(uc_pts(NUM_INIT_PTS_PER_DIM^NUM_DIMS + i, 1), ...
               uc_pts(NUM_INIT_PTS_PER_DIM^NUM_DIMS + i, 2), ...
               num2str(i), 'Color', 'k');
        end
        plot(initial_pts(:,1), initial_pts(:,2), 'bx', 'MarkerSize', 10);
      end
    pause,
    close;
  end

  fprintf('\n');

end

% Plot the errors
% figure;
% plot(uc_err_prog.f1, 'b-o');
% plot(uc_err_prog.f2, 'k-x');
% plot(uc_err_prog.f3, 'g-s');
% plot(uc_err_prog.f4, 'm-+');
% legend('f1', 'f2', 'f3', 'f4');

