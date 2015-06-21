% Tests our Active Learning Algorithm on a 2D example

REUSE_DATA = false;
% REUSE_DATA = true;

if ~REUSE_DATA
%   clear all;
  REUSE_DATA = false;
end
close all;

PLOT_OK_GLOBAL = false;

% DEBUG
QUICK_DEBUG_MODE = false;
% QUICK_DEBUG_MODE = true;

% SELECT MEANFUNCTION
MEAN_FUNC = 'mean';
MEAN_FUNC = 'lowest';
% MEAN_FUNC = 'logbarrier';

% Select Experiment
% experiment = 'exp1';
experiment = 'exp2';

% Set the following parameters for the Experiment
ALP1 = 1.1; BEP1 = 1.3; % The prior for the theta1 beta binomial.
ALP2 = 1.2; BEP2 = 1; % The prior for the theta2 beta binomial.
PRIOR_PARAMS = [ALP1, BEP1, ALP2, BEP2]';
NUM_DATA_SAMPLES = 150;
NUM_INITIAL_PTS_PER_DIM = 2; % # points at which to run the initial simulation.
NUM_KFOLDCV_PARTITIONS = 20;
NOISE_LEVEL = 0.3; % Noise associated with measurements
BORDER_TOL = 1e-2;
SET_PARAMETERS_TO_OPTIMAL = true;
PLOT_OK_KL_PROGRESS = false;

if QUICK_DEBUG_MODE
% NUM_AL_ITERS = 2;
  NUM_AL_CANDIDATES = 3;
  NUM_SAMPLE_FUNCTIONS = 2;
  RESOLUTION_PER_DIM = 8;
else
  NUM_AL_ITERS = 100; % Number of iterations of active learning
  NUM_AL_CANDIDATES = 20;
  NUM_SAMPLE_FUNCTIONS = 40;
  RESOLUTION_PER_DIM = 50; % resolution per dimension
end

% Specify figures for the experiment
% fig_post = figure;
if PLOT_OK_KL_PROGRESS
  fig_kl_progress = figure;
  fig_post = figure;
end
% Specify function for plotting. Mesh, contour or contour3
PLOT_FUNC = @contour;

% Create function oracles for the experiment
if strcmp(experiment, 'exp1')
  genExperimentData = @(eval_x, eval_y) sampleExp1(PRIOR_PARAMS, ...
    eval_x, eval_y, NUM_DATA_SAMPLES);
  evalLogJointProbs = @(eval_x, eval_y, sum_x) evalLogJointProbsExp1(...
    eval_x, eval_y, sum_x, NUM_DATA_SAMPLES, PRIOR_PARAMS, NOISE_LEVEL); 
elseif strcmp(experiment, 'exp2')
  genExperimentData = @(eval_x, eval_y) sampleExp2(eval_x, eval_y);
  evalLogJointProbs = @(eval_x, eval_y, garb1) evalLogJointProbsExp2(...
    eval_x, eval_y, NOISE_LEVEL);
end

if strcmp(experiment, 'exp1')
  LOWEST_VAL = -900;
elseif strcmp(experiment, 'exp2')
  LOWEST_VAL = -60;
end
% Specify hyper parameters for GP regression. Choose from the following in
% K-fold cross validation
% gp_smoothness_params = (0.2: 0.1: 0.6)';
% gp_scale_params = linspace(30, 150, 20)';
gp_sm_param_ratios = logspace(log10(0.5), log10(2.0), 20)';
gp_scale_params = logspace(log10(10), log10(300), 20)';
if SET_PARAMETERS_TO_OPTIMAL
  optimal_bandwidth = 0.13;
  optimal_scale = 60;
end

% Define the following before proceeding
% Specify the parameter space
th1 = linspace(0, 1, RESOLUTION_PER_DIM)';
[Th1, Th2] = meshgrid(th1, th1);
th = [Th1(:), Th2(:)];
% Specify the points at which the likelihood should be evaluated initially.
ip = linspace(0, 1, NUM_INITIAL_PTS_PER_DIM + 2)'; ip = ip(2: end-1);
[Ip1, Ip2] = meshgrid(ip, ip);
initial_pts = [Ip1(:), Ip2(:)];
% Clear all temporary variables
clear th1 al ip Ip1 Ip2; % Don't clear Al1, Al2. Need it for bilinear-interp.

if ~REUSE_DATA
  fprintf('Generating New Data: \n');
  [theta_true, sumX, True_post, marginal_prob] = genExperimentData(Th1, Th2);
  true_post = True_post(:);
else
  fprintf('Reusing old Data: theta = (%0.4f, %0.4f), sumX = %d\n', ...
          theta(1), theta(2), sumX);
end
% Plot the true posterior
% figure(fig_post);
% PLOT_FUNC(Th1, Th2, True_post, 'Color', 'g'); xlabel('th1'); ylabel('th2');

% Use EER
% ExpectedErrorReduction;
% QBC;
% KL_progress = uc_kl_progress;

% Use uncertainty reduction
% UncertaintyReduction_2D;

% Use MCMC to obtain an estimate of the posterior
% MCMCForPostEstimation_2D;

% MCMC - Density estimation
MCMCKDE;

% Obtain a baseline performance using a grid search.
% bruteForceSearchForPostEstimation_2D;

% MCMCKDE;

