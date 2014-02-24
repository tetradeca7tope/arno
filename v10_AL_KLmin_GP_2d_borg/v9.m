% Tests our Active Learning Algorithm on a 2D example

REUSE_DATA = false;
% REUSE_DATA = true;

if ~REUSE_DATA
%   clear all;
  REUSE_DATA = false;
end
close all;

PLOT_OK_GLOBAL = false;

% Select Experiment
% experiment = 'exp1';
experiment = 'exp2';

% Set the following parameters for the Experiment
ALP1 = 1.1; BEP1 = 1.3; % The prior for the theta1 beta binomial.
ALP2 = 1.2; BEP2 = 1; % The prior for the theta2 beta binomial.
PRIOR_PARAMS = [ALP1, BEP1, ALP2, BEP2]';
NUM_DATA_SAMPLES = 150;
NUM_INITIAL_PTS_PER_DIM = 4; % # points at which to run the initial simulation.
NUM_AL_CANDIDATES_PER_DIM = 20; % # pts per dimension to choose from for AL
% NUM_AL_ITERS = 100; % Number of iterations of active learning
NUM_SAMPLE_FUNCTIONS = 40;
NUM_KFOLDCV_PARTITIONS = 20;
NOISE_LEVEL = 3; % Noise associated with measurements
RESOLUTION_PER_DIM = 40; % resolution per dimension

% Specify figures for the experiment
fig_post = figure;
fig_kl_progress = figure;
% Specify function for plotting. Mesh, contour or contour3
PLOT_FUNC = @contour;

% Create function oracles for the experiment
if strcmp(experiment, 'exp1')
  genExperimentData = @(eval_x, eval_y) sampleExp1(PRIOR_PARAMS, ...
    eval_x, eval_y, NUM_DATA_SAMPLES);
  evalLogJointProbs = @(eval_x, eval_y, sum_x) evalLogJointProbsExp1(...
    eval_x, eval_y, sum_x, NUM_DATA_SAMPLES, PRIOR_PARAMS, NOISE_LEVEL); 
elseif strcmp(experiment, 'exp2')
  
end

% Specify hyper parameters for GP regression. Choose from the following in
% K-fold cross validation
gp_smoothness_params = (0.2: 0.1: 0.6)';
gp_scale_params = linspace(30, 150, 20)';

% Define the following before proceeding
% Specify the parameter space
th1 = linspace(0, 1, RESOLUTION_PER_DIM)';
[Th1, Th2] = meshgrid(th1, th1);
th = [Th1(:), Th2(:)];
% Specify candidates for Active Learning
al = linspace(0, 1, NUM_AL_CANDIDATES_PER_DIM + 2)'; al = al(2: end-1);
[Al1, Al2] = meshgrid(al, al);
al_candidates = [Al1(:), Al2(:)];
starting_al_candidates = al_candidates;
% Specify the points at which the likelihood should be evaluated initially.
ip = linspace(0, 1, NUM_INITIAL_PTS_PER_DIM + 2)'; ip = ip(2: end-1);
[Ip1, Ip2] = meshgrid(ip, ip);
initial_pts = [Ip1(:), Ip2(:)];
% Clear all temporary variables
clear th1 al ip Ip1 Ip2; % Don't clear Al1, Al2. Need it for bilinear-interp.

if ~REUSE_DATA
  fprintf('Generating New Data: ');
  [theta_true, sumX, True_post, marginal_prob] = genExperimentData(Th1, Th2);
  true_post = True_post(:);
else
  fprintf('Reusing old Data: theta = (%0.4f, %0.4f), sumX = %d\n', ...
          theta(1), theta(2), sumX);
end
% Plot the true posterior
figure(fig_post);
PLOT_FUNC(Th1, Th2, True_post, 'Color', 'g'); xlabel('th1'); ylabel('th2');
pause;

% Initialize Active Learning
num_rem_candidates = NUM_AL_CANDIDATES_PER_DIM^2;
rem_cand_indices = (1: NUM_AL_CANDIDATES_PER_DIM^2)';
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
  candidates.sigmaSmVals = gp_smoothness_params;
  candidates.sigmaPrVals = gp_scale_params;
  % specify other hyper params
  hyper_params.noise = NOISE_LEVEL * ones(num_curr_pts, 1);
  hyper_params.meanFunc = [];
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
  fprintf('KL between Estimate and True Posterior: %f\n', true_KL);
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
          xlabel('th1'); ylabel('th2');
          % Compare the estimated log joint and the sample log joint
          figure(flogjoint);
          Sample_log_joint = reshape(curr_gp_samples(:,gp_sample_iter), ...
                              RESOLUTION_PER_DIM, RESOLUTION_PER_DIM);
          mesh(Th1, Th2, Est_joint_log_probs); hold on,
          mesh(Th1, Th2, Sample_log_joint); hold off,
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
  accum_cand_KLs = zeros(num_rem_candidates, 1);
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
    Sample_cand_log_joint = interp2(Th1, Th2, Sample_log_joint, ...
                                    Al1, Al2);
    sample_cand_log_joint = Sample_cand_log_joint(:);
    sample_cand_log_joint = sample_cand_log_joint(rem_cand_indices, :);
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
    % Some error-checking before we proceed
    if (num_rem_candidates ~= size(rem_cand_indices, 1))
      error('Num remaining candidates: indices and counter do not match.');
    end
    if sum( (al_candidates ~= starting_al_candidates(rem_cand_indices, :)) )
      error('al_candidates and indices do not match.');
    end

    for cand_iter = 1:num_rem_candidates
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
  % remove the point from the candidate list and also the log-joint list
  al_candidates = al_candidates([1:min_idx-1, min_idx+1:end], :);
  rem_cand_indices = rem_cand_indices([1:min_idx-1, min_idx+1:end]);
  num_rem_candidates = num_rem_candidates - 1;

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

% Use MCMC to obtain an estimate of the posterior
MCMCForPostEstimation_2D;

% Obtain a baseline performance using a grid search.
bruteForceSearchForPostEstimation_2D;

