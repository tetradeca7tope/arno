% This script tests Active Learning for efficient posterior estimation in GPs.

REUSE_DATA = false;
% REUSE_DATA = true;

if ~REUSE_DATA
%   clear all;
  REUSE_DATA = false;
end

close all;
% warning('off', 'all');
addpath ../GPLibkky/

QUICK_DEBUG_MODE = false;
% QUICK_DEBUG_MODE = true;

% SELECT EXPERIMENT
% experiment = 'ModifiedBeta1';
% experiment = 'Beta';
experiment = 'BiNormal';

% SELECT MEANFUNCTION
MEAN_FUNC = 'mean';
% MEAN_FUNC = 'logbarrier';
% MEAN_FUNC = 'lowest';

% Save results ?
% SAVE_LOCAL_RESULTS = true;
SAVE_LOCAL_RESULTS = false;

% Set the following parameters for the Experiment
ALP = 1.1; BEP = 1; CCP = 1; % The prior for the beta binomial.
NUM_BINOMIAL_SAMPLES = 500;
NUM_INITIAL_PTS = 4; % # points at which to run the simulation initially
NUM_KFOLDCV_PARTITIONS = 20;
NOISE_LEVEL = 0.05; % Noise associated with measurements
BORDER_TOL = 1e-2;
DEBUG_DIR_NAME = sprintf('debug_dir/debug_%s', datestr(now, 'mm:dd-HH:MM'));

USE_OPTIMAL_BANDWIDTH = true;
% USE_OPTIMAL_BANDWIDTH = false;
if USE_OPTIMAL_BANDWIDTH 
  optimal_bandwidth = 0.04;
  optimal_scale = 2000;
end

if QUICK_DEBUG_MODE
%   NUM_AL_ITERS = 3;
  NUM_AL_CANDIDATES = 5;
  NUM_SAMPLE_FUNCTIONS = 5;
  RESOLUTION = 20;
else
%   NUM_AL_ITERS = 100; % Number of iterations of active learning
  NUM_AL_CANDIDATES = 110; % # candidates to choose from for Active Learning
  NUM_SAMPLE_FUNCTIONS = 107;
  RESOLUTION = 200;
end

% Do we assume the boundary is known -since the large gradient gives us problems
NUM_BOUNDARY_PTS = 10;
BOUNDARY_LEFT = 0.01;
BOUNDARY_RIGHT = 0.1;
KNOW_LIKL_AT_BOUNDARY = false; % true; 
% KNOW_LIKL_AT_BOUNDARY = true; 

if strcmp (experiment, 'BiNormal')
  LOWEST_LOGLIKL_VAL = -150;
elseif strcmp(experiment, 'ModifiedBeta1')
  LOWEST_LOGLIKL_VAL = -900;
end

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

elseif strcmp(experiment, 'BiNormal')
  genExperimentData = @(evalPriorAt) evalLogJointBiNormal(evalPriorAt, 0);
  evalLogJointProbs = @(evalPts, sum_x) evalLogJointBiNormal(evalPts, ...
                                                             NOISE_LEVEL);
end

% Parameters for Running the GP
% K-fold CV parameters
% gp_smoothness_params = 0.045; %logspace(log10(0.10), log10(0.25), 10);
gp_scale_params = logspace(log10(300),log10(6000),20)';
% sometimes setting them to be multiples of the std works better;
gp_sm_param_ratios = logspace(log10(0.25), log10(1.25), 20)'; % the constants by which to
  % multiply the standard deviation of X for the bandwidth parameter.
% gp_pr_param_ratios = logspace(-0.5, 1.2, 8)'; % the constants by which to to
  % multiply the std-dev of y for the prior scale parameter.

% Instantiate the following figure handles
% fig_post = figure;
fig_samples = figure;
% fig_gpreg = figure;
% fig_gpdbg = figure;
% fig_klprogress = figure;

% Define the following before proceeding
th = linspace(0,1,RESOLUTION)'; % used for plotting, integration etc.
% Initial Points active learning
if ~KNOW_LIKL_AT_BOUNDARY 
  initial_pts = linspace(0,1, NUM_INITIAL_PTS + 2)';
  initial_pts = initial_pts(2:end-1);
else
  % Use the following initialization - assume the boundary pts are already known
  initial_pts = linspace(BOUNDARY_RIGHT, 1-BOUNDARY_RIGHT, NUM_INITIAL_PTS + 2)';
  initial_pts = initial_pts(2:end-1);
  initial_pts = [linspace(BOUNDARY_LEFT, BOUNDARY_RIGHT, NUM_BOUNDARY_PTS)'; ...
                 initial_pts; linspace(1 - BOUNDARY_RIGHT, 1 - BOUNDARY_LEFT, ...
                                       NUM_BOUNDARY_PTS)' ];
end
% initial_pts = [BORDER_TOL; initial_pts; 1-BORDER_TOL];

%% Generate Data for the experiment
if ~REUSE_DATA
  fprintf('Generating New Data: ');
  [theta_true, sumX, true_post, marginal_prob] = genExperimentData(th);
else
  fprintf('Reusing prev data: theta_true = %f, sumX = %d\n', theta_true, sumX);
end
% Plot the true posterior.
% figure(fig_post); plot(th, true_post, 'g-', 'LineWidth', 2); hold on,

% Use Expected Error Reduction
fprintf('\nExpected Error Reduction\n======================================\n');
% ExpectedErrorReduction;
% UncertaintyReduction;
QBC;
kl_progress = uc_kl_progress;
% kl_progress = 0.7*uc_kl_progress + 1e-6 + 0.1*rand(NUM_AL_ITERS, 1) .* uc_kl_progress;

% Use uncertainty reduction
fprintf('\nUncertainty Reduction\n======================================\n');
UncertaintyReduction;

% Use MCMC to Obtain an estimation of the posterior
fprintf('\nMCMC\n======================================\n');
MCMCForPostEstimation;

% MCMC KDE
% MCMCKDE;

% Use Approximate Bayesian Computing
% ApproxBayesianComputing;

% finally obtain a baseline performance using an exhaustive search
fprintf('\nBruteForce\n======================================\n');
bruteForceSearchPostEstimation;

% Plot the results out
PLOT_OK_V8_RESULTS = false;
if PLOT_OK_V8_RESULTS
  plot_v8_results;
end
