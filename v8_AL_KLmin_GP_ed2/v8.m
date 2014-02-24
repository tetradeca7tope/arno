% This script tests Active Learning for efficient posterior estimation in GPs.

REUSE_DATA = false;
% REUSE_DATA = true;

if ~REUSE_DATA
  clear all;
  REUSE_DATA = false;
end
close all;
% warning('off', 'all');
addpath ../GPLibkky/

% SELECT EXPERIMENT
experiment = 'ModifiedBeta1';
% experiment = 'Beta';
% SELECT MEANFUNCTION
% MEAN_FUNC = 'mean';
MEAN_FUNC = 'logbarrier';

% Set the following parameters for the Experiment
ALP = 1.1; BEP = 1; CCP = 1; % The prior for the beta binomial.
NUM_BINOMIAL_SAMPLES = 500;
NUM_INITIAL_PTS = 5; % # points at which to run the simulation initially
NUM_AL_CANDIDATES = 300; % # candidates to choose from for Active Learning
NUM_AL_ITERS = 100; % Number of iterations of active learning
NUM_SAMPLE_FUNCTIONS = 100;
NUM_KFOLDCV_PARTITIONS = 10;
RESOLUTION = 200;
NOISE_LEVEL = 1; % Noise associated with measurements
BORDER_TOL = 1e-2;

% Other preprocessing
priorParams.alp = ALP;
priorParams.bep = BEP;
% Set parameters for each experiment
if strcmp(experiment, 'ModifiedBeta1')
  priorParams.c_param = CCP;
  genExperimentData = @(evalPriorAt) sampleModifiedBeta(priorParams, ...
    evalPriorAt, NUM_BINOMIAL_SAMPLES);
  evalLogJointProbs = @(evalPts, sum_x) evalLogJointModifiedBetaProbs(...
    evalPts, sum_x, NUM_BINOMIAL_SAMPLES, priorParams, NOISE_LEVEL); 

elseif strcmp(experiment, 'Beta')
  genExperimentData = @(evalPriorAt) sampleBeta(priorParams, ...
    evalPriorAt, NUM_BINOMIAL_SAMPLES);
  evalLogJointProbs = @(evalPts, sum_x) evalLogJointBetaProbs(...
    evalPts, sum_x, NUM_BINOMIAL_SAMPLES, priorParams, NOISE_LEVEL); 
end

% Parameters for Running the GP
% K-fold CV parameters
gp_smoothness_params = (0.01:0.01:0.1)';
gp_scale_params = linspace(100,600,8)';
% sometimes setting them to be multiples of the std works better;
gp_sm_param_ratios = logspace(-0.5, 1.2, 12)'; % the constants by which to
  % multiply the standard deviation of X for the bandwidth parameter.
gp_pr_param_ratios = logspace(-0.5, 1.2, 8)'; % the constants by which to to
  % multiply the std-dev of y for the prior scale parameter.

% Set the Active Learning candidates
al_candidates = linspace(0, 1, NUM_AL_CANDIDATES+2)';
al_candidates = al_candidates(2:end-1);

% Instantiate the following figure handles
fig_post = figure;
fig_samples = figure;

% Define the following before proceeding
th = linspace(0,1,RESOLUTION)'; % used for plotting, integration etc.
% candidates for active learning
initial_pts = linspace(0,1, NUM_INITIAL_PTS+2)';
initial_pts = initial_pts(2:end-1);

%% Generate Data for the experiment
if ~REUSE_DATA
  fprintf('Generating New Data: ');
  [theta_true, sumX, true_post, prior_norm_constant] = genExperimentData(th);
else
  fprintf('Reusing prev data: theta_true = %f, sumX = %d\n', theta_true, sumX);
end
% Plot the true posterior.
figure(fig_post); plot(th, true_post, 'g-', 'LineWidth', 2); hold on,

% Initiate Active Learning
num_rem_candidates = NUM_AL_CANDIDATES;
curr_pts = initial_pts; % at each iter the joint is evaluated at these pts.
kl_progress = zeros(NUM_AL_ITERS, 1);

% Now run Active Learning
  for al_iter = 1:NUM_AL_ITERS

  % We need to do the following at each Active Learning iteration
  % 1. Evaluate the joint
  % 2. Run the GP
  % 3. Generate Samples and compute expected KL
  % 4. Pick the best point

  % prelims
  num_curr_pts = size(curr_pts, 1);

  % 1. Evalute the joint
  % ====================
  [obs_joint_log_probs] = evalLogJointProbs(curr_pts, sumX);

  % 2. Run GP
  % =========
  % 2.1 Prep structure for LOO CV
%   candidates.sigmaSmVals = 0.1 * gp_sm_param_ratios * std(curr_pts);
  candidates.sigmaSmVals = gp_smoothness_params;
  candidates.sigmaPrVals = gp_scale_params;
  % Specify how to decide mean function for the GP.
  if strcmp(MEAN_FUNC, 'logbarrier')
    % first fit a log barrier to the data
    [a, b, c] = fitLogBarrier(curr_pts, obs_joint_log_probs);
    hyper_params.meanFunc = @(arg) (a*log(arg) + b*log(1-arg) + c);
  else
    hyper_params.meanFunc = @(arg) mean(obs_joint_log_probs);
  end
  hyper_params.noise = NOISE_LEVEL * ones(num_curr_pts, 1);
  % 2.2 Run GP with K-fold CV
  [est_joint_log_probs, K, sigmaSmOpt, sigmaPrOpt] = GPKFoldCV(curr_pts, ...
    obs_joint_log_probs, th, NUM_KFOLDCV_PARTITIONS, candidates, hyper_params);
  est_joint_probs = exp(est_joint_log_probs);
  % 2.3 estimate prob(X)
  est_pX = numerical_1D_integration(th, est_joint_probs);
  est_log_pX = log(est_pX);
  % 2.4 Estimate the posterior
  est_log_post = est_joint_log_probs - est_log_pX;
  est_post = exp(est_log_post);
  % compute the KL divergence between the estimate and the true posterior
  ptwise_KL = true_post .* (log(true_post) - est_log_post);
  true_KL = numerical_1D_integration(th(2:end-1), ptwise_KL(2:end-1));
  kl_progress(al_iter) = true_KL;

  % 3. Generate Samples and commpute Expected KL at each candidate pt
  % =================================================================
  curr_gp_samples = GPDrawSamples(est_joint_log_probs,K, NUM_SAMPLE_FUNCTIONS)';
    % note the transpose above. now each column is a sample
  sample_log_pX = log(numerical_1D_integration(th, exp(curr_gp_samples)));
  sample_log_post = bsxfun(@minus, curr_gp_samples, sample_log_pX);
  sample_post = exp(sample_log_post);
  pointwise_KL = bsxfun(@times, ...
                        bsxfun(@plus, -sample_log_post, est_log_post), ...
                        est_post);
  sample_KLs = numerical_1D_integration(th, pointwise_KL);
  estim_KL = mean(sample_KLs);

  % Prepare iterations over samples
  % create this struct which we will need later
  % add a 0-noise component to the struct 
  opt_hyper_params.noise = hyper_params.noise;
  opt_hyper_params.sigmaSm = sigmaSmOpt;
  opt_hyper_params.sigmaPr = sigmaPrOpt;
  opt_hyper_params.meanFunc = hyper_params.meanFunc;
  accum_cand_KLs = zeros(num_rem_candidates, 1);

  % 4. Now Pick the best point
  % ======================
  % 4.1 First accumulate the KL divergences
  for cand_iter = 1:num_rem_candidates
    for sample_iter = 1:NUM_SAMPLE_FUNCTIONS

      sample_estimate = interp1(th, curr_gp_samples(:, sample_iter), ...
                                al_candidates(cand_iter));
      % Now create hypothetical regressors and regressands
      reg_y = [sample_estimate; obs_joint_log_probs];
      reg_x = [al_candidates(cand_iter); curr_pts];

      % For the hypothetical regression set the noise param to zero.
      opt_hyper_params.noise = [0; hyper_params.noise]; 
      s_joint_log_probs = GPRegression(reg_x, reg_y, th, opt_hyper_params);
      % Now estimate the joint and then the posterior
      s_pX = numerical_1D_integration(th, exp(s_joint_log_probs));
      s_log_pX = log(s_pX);
      s_log_post = s_joint_log_probs - s_log_pX;
      s_post = exp(s_log_post);
      % Now compute the KL between this sample and the estimate
      s_ptwise_KL = sample_post(:, sample_iter) .* ( ...
        sample_log_post(:, sample_iter) - s_log_post);
      curr_sample_KL = numerical_1D_integration(th(2:end-1), ...
                                                s_ptwise_KL(2:end-1));
      accum_cand_KLs(cand_iter) = accum_cand_KLs(cand_iter) + curr_sample_KL;
    end
  end
  expected_cand_KLs = accum_cand_KLs / NUM_SAMPLE_FUNCTIONS;
  % 4.2 Pick the minimum and add him to the list of current points
  [~, min_idx] = min(expected_cand_KLs);
  min_pt = al_candidates(min_idx);
  old_pts = curr_pts;
  curr_pts = [min_pt; curr_pts];
  al_candidates = al_candidates([1:min_idx-1, min_idx+1:end]);
  num_rem_candidates = num_rem_candidates - 1;

  % Plot & Print Relevant results
  % ==============================
  % Plot the current estimate
  PLOT_OK = 0;
  if PLOT_OK
    figure(fig_post);
    plot(th, est_post, 'r--');
    text(min_pt, -1, num2str(al_iter), 'Color', 'k');
    axis([0 1 -2 1.05*max([est_post; true_post])]);
    titlestr = sprintf(...
               ['Estimated Post (r--) & True Posterior (g-)\n', ...
                'Prior: Beta(%d, %d, %d), theta: %f, Num initial pts: %d\n', ...
                'Num data points: %d, Num positive points = %d\n', ...
                'Noise level: %f, gp-bandwidth  = %f, gp-scale = %f\n', ...
                'True KL: %f, Estim KL (using GP): %f'], ...
                ALP, BEP, CCP, theta_true, NUM_INITIAL_PTS, ...
                NUM_BINOMIAL_SAMPLES, sumX, NOISE_LEVEL, ...
                sigmaSmOpt, sigmaPrOpt, true_KL, estim_KL);
    title(titlestr);
    % Plot the samples drawn.
    figure(fig_samples);
    for i = 1:NUM_SAMPLE_FUNCTIONS
      plot(th, curr_gp_samples(:,i), 'm'); hold on,
    end
    % plot the points evaluated
    plot(old_pts, obs_joint_log_probs, 'kx', 'MarkerSize', 6);
    hold off,
  end % PLOT_OK
  % also print the results out.
  fprintf('iter: %d, true-KL: %f, estim-KL: %f, min-pt: %f\n', ....
          al_iter, true_KL, estim_KL, min_pt);
end

% Finally once this is all done plot the final estimate out
fig_ALfinal = figure; plot(th, true_post, 'g-', 'LineWidth', 2); hold on,
plot(th, est_post, 'r--', 'LineWidth', 2);
title('true posterior (g-) vs final estimated posterior (r--)');
fig_klprogress = figure;
plot(log(kl_progress), 'bo-'); hold on,

% Use MCMC to Obtain an estimation of the posterior
% MCMCForPostEstimation;

% finally obtain a baseline performance using an exhaustive search
bruteForceSearchPostEstimation;

