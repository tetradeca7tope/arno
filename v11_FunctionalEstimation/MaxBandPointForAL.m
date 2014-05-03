% Script runs MaxBandPoint for Active Learning

% Parameters for MaxBandPoint
INIT_LIPSCHTIZ_CONST = 5; % set this to a small value initially 

% TODO: repeating code in UncertaintyReduction.m. Fix This !
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
initial_log_probs = evalLogJoint(initial_pts);

% Initialie Active Learning
mbp_pts = initial_pts;
% Error tracker for each functional
mbp_error_prog.f1 = zeros(num_results_to_be_stored, 1);
mbp_error_prog.f2 = zeros(num_results_to_be_stored, 1);
mbp_error_prog.f3 = zeros(num_results_to_be_stored, 1);
mbp_error_prog.f4 = zeros(num_results_to_be_stored, 1);

% Set some parameters for mbp
phi = @exp; grad_phi = @exp; % use exponential transformation
al_mbp_params.num_iters = 40;
al_mbp_params.init_step_size = 0.25;
% First obtain the points from Active Learning
[mbp_pts, mbp_log_probs, mbp_lip_const] = alMaxBandPoint(evalLogJoint, ...
  initial_pts, initial_log_probs, phi, grad_phi, INIT_LIPSCHTIZ_CONST, ...
  PARAM_SPACE_BOUNDS, NUM_AL_ITERS, al_mbp_params);

  PLOT_LOCAL_OK = false;
%   PLOT_LOCAL_OK = true;
  if PLOT_LOCAL_OK
    % Plot the points selected by mbp
    if NUM_DIMS == 1
      plot(mbp_pts, mbp_log_probs, 'rx'); hold on,
      th = linspace(PARAM_SPACE_BOUNDS(1), PARAM_SPACE_BOUNDS(2), 100)';
      plot(th, evalLogJoint(th), 'b-');
    else 
      plot(mbp_pts(:,1), mbp_pts(:,2), 'rx'); hold on,
      axis([PARAM_SPACE_BOUNDS PARAM_SPACE_BOUNDS]);
      pause;
      plot2DFunction(evalLogJoint, [PARAM_SPACE_BOUNDS PARAM_SPACE_BOUNDS], ...
                     'mesh');
    end
    title('Pts chosen by MBP');
    pause,
  end

% Now perform Regression on each of these points and obtain the estimates of the
% functionals.
for mbp_iter = 1:num_results_to_be_stored

  curr_num_mbp_pts = mbp_iter * STORE_RESULTS_EVERY + num_initial_pts;
  fprintf('MBP with %d pts: \n', curr_num_mbp_pts);
  % The regressors and regressands for the curent iteration
  Xtr = mbp_pts(1:curr_num_mbp_pts, :);
  Ytr = mbp_log_probs( 1:curr_num_mbp_pts );

  % GP Regression
  logJointEst = regressionWrap(Xtr, Ytr, NOISE_LEVEL, LOWEST_LOGLIKL_VAL, ...
    LOGLIKL_RANGE, cv_cost_func);

  % Now perform MCMC
  estimate_mcmc_samples = CustomMCMC( ...
    NUM_MCMC_ITERS_FOR_EST + NUM_MCMC_BURNIN_FOR_EST, ...
    MCMC_EST_PROPOSAL_STD, MCMC_EST_INIT_PT, logJointEst);
  Xmcmc = estimate_mcmc_samples(NUM_MCMC_BURNIN_FOR_EST+1:end, :);
  % Now evaluate the functionals on Xmcmc
  T1 = f1(Xmcmc); mbp_err_prog.f1(mbp_iter)= abs(T1 - f_vals.S1)/abs(f_vals.S1);
  T2 = f2(Xmcmc); mbp_err_prog.f2(mbp_iter)= abs(T2 - f_vals.S2)/abs(f_vals.S2);
  T3 = f3(Xmcmc); mbp_err_prog.f3(mbp_iter)= abs(T3 - f_vals.S3)/abs(f_vals.S3);
  T4 = f4(Xmcmc); mbp_err_prog.f4(mbp_iter)= abs(T4 - f_vals.S4)/abs(f_vals.S4);
  fprintf('Iter %d emp-errors: f1: %f, f2: %f, f3: %f, f4: %f\n', ...
          mbp_iter, ...
          mbp_err_prog.f1(mbp_iter), ...
          mbp_err_prog.f2(mbp_iter), ...
          mbp_err_prog.f3(mbp_iter), ...
          mbp_err_prog.f4(mbp_iter) );
  
  PLOT_LOCAL_OK = false;
%   PLOT_LOCAL_OK = true;
  if PLOT_LOCAL_OK
    % Plot the points selected by mbp
    figure;
    if NUM_DIMS == 1
      plot(mbp_pts, mbp_log_probs, 'rx'); hold on,
      th = linspace(PARAM_SPACE_BOUNDS(1), PARAM_SPACE_BOUNDS(2), 100)';
      plot(th, evalLogJoint(th), 'b-');
      plot(th, logJointEst(th), 'g-');
    else 
      plot3(mbp_pts(:,1), mbp_pts(:,2), mbp_log_probs, 'kx', 'MarkerSize', ...
            10, 'LineWidth', 3); hold on,
      axis([PARAM_SPACE_BOUNDS PARAM_SPACE_BOUNDS]);
      pause;
      plot2DFunction(evalLogJoint, [PARAM_SPACE_BOUNDS PARAM_SPACE_BOUNDS], ...
                     'mesh', 'b');
      plot2DFunction(logJointEst, [PARAM_SPACE_BOUNDS PARAM_SPACE_BOUNDS], ...
                     'mesh', 'g');
    end
    title('Pts chosen by MBP');
    pause,
    close;
  end
end
