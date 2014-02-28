% Script for testing functional estimation

% Prelims
clear all;
close all;
rng('shuffle'); % shuffle the seed for the RNG
addpath ../helper/
addpath ../GPLibkky/

% Constants
NUM_DIMS = 10;
LOWEST_LOGLIKL_VAL = -30;
LOGLIKL_RANGE = 35;

DEBUG_MODE = false;
% DEBUG_MODE = true;
% Constants for Active Learning
if ~DEBUG_MODE
  NUM_AL_ITERS = 100;
  NUM_MCMC_ITERS_FOR_EST = 150000;
  NUM_MCMC_BURNIN_FOR_EST = max(100, round(NUM_MCMC_ITERS_FOR_EST/2) );
  NUM_EXPERIMENTS = 30;
else
  NUM_AL_ITERS = 3;
  NUM_MCMC_ITERS_FOR_EST = 10;
  NUM_MCMC_BURNIN_FOR_EST = 2;
  NUM_EXPERIMENTS = 2;
end
NUM_INIT_PTS_PER_DIM = 100;
INIT_ON_GRID = false;
MCMC_EST_PROPOSAL_STD = 0.5;
MCMC_EST_INIT_PT = zeros(NUM_DIMS, 1);

% Constants for MCMC
NUM_MCMC_SAMPLES = 5*NUM_AL_ITERS;

% RUN TIME HYPER-PARAMS
USE_OPT_PARAMS = false;
OPT_BANDWIDTH = 2.3* NUM_AL_ITERS^(-1/(1.3 + NUM_DIMS));
OPT_SCALE = 50;
NUM_KFOLD_CV_PARTITIONS = 20;
functionals = {@f1, @f2, @f3, @f4};

% parameters for the test
p1 = 0.6;
p2 = 1 - p1;
sigma = 0.2*sqrt(NUM_DIMS) + 0.05*rand();

% Create function for evaluating joint likelihood and
% obtain the True functional values
evalLogJoint = @(arg) evalLogLiklExp1(arg, sigma, p1, p2);
[~, f_vals] = evalLogJoint(zeros(0,NUM_DIMS));

% Define other useful variables before proceeding
num_init_pts = NUM_INIT_PTS_PER_DIM ^ NUM_DIMS;

% Create arrays for storing the arrays.
uc_errs.f1 = zeros(NUM_EXPERIMENTS, NUM_AL_ITERS);
uc_errs.f2 = zeros(NUM_EXPERIMENTS, NUM_AL_ITERS);
uc_errs.f3 = zeros(NUM_EXPERIMENTS, NUM_AL_ITERS);
uc_errs.f4 = zeros(NUM_EXPERIMENTS, NUM_AL_ITERS);
mcmc_errs.f1 = zeros(NUM_EXPERIMENTS, NUM_MCMC_SAMPLES);
mcmc_errs.f2 = zeros(NUM_EXPERIMENTS, NUM_MCMC_SAMPLES);
mcmc_errs.f3 = zeros(NUM_EXPERIMENTS, NUM_MCMC_SAMPLES);
mcmc_errs.f4 = zeros(NUM_EXPERIMENTS, NUM_MCMC_SAMPLES);

% For saving results
save_file_name = sprintf('results/res_d%d_%s.mat', NUM_DIMS, ...
                         datestr(now, 'mm:dd-HH:MM'));

for experiment_iter = 1:NUM_EXPERIMENTS

  fprintf('EXPERIMENT: %d\n======================\n', experiment_iter);
  fprintf('UNCERT-REDUCTION\n');
  UncertaintyReduction;
  fprintf('MCMC\n');
  MCMCForPostEstimation;

  % Store Uncert-Reduction Results
  uc_errs.f1(experiment_iter, :) = uc_err_prog.f1';
  uc_errs.f2(experiment_iter, :) = uc_err_prog.f2';
  uc_errs.f3(experiment_iter, :) = uc_err_prog.f3';
  uc_errs.f4(experiment_iter, :) = uc_err_prog.f4';
  % Store MCMC Results
  mcmc_errs.f1(experiment_iter, :) = mcmc_err_prog.f1';
  mcmc_errs.f2(experiment_iter, :) = mcmc_err_prog.f2';
  mcmc_errs.f3(experiment_iter, :) = mcmc_err_prog.f3';
  mcmc_errs.f4(experiment_iter, :) = mcmc_err_prog.f4';
  save(save_file_name, 'uc_errs', 'mcmc_errs');

end

