% Script for testing functional estimation

% Prelims
clear all;
close all;
rng('shuffle'); % shuffle the seed for the RNG
% addpath ~/software/matlab/libs/kky-matlab/
addpath ../helper/
addpath ../GPLibkky/
addpath ../LipLibkky/

% Constants
NUM_DIMS = 5;
LOWEST_LOGLIKL_VAL = -40;
LOGLIKL_RANGE = 50;

DEBUG_MODE = false;
% DEBUG_MODE = true;
% Constants for Active Learning
if ~DEBUG_MODE
  NUM_AL_ITERS = 500;
  NUM_MCMC_ITERS_FOR_EST = 150000;
  NUM_MCMC_BURNIN_FOR_EST = max(100, round(NUM_MCMC_ITERS_FOR_EST/2) );
  NUM_EXPERIMENTS = 30;
  STORE_RESULTS_EVERY = 50;
else
  NUM_AL_ITERS = 6;
  NUM_MCMC_ITERS_FOR_EST = 10;
  NUM_MCMC_BURNIN_FOR_EST = 2;
  NUM_EXPERIMENTS = 2;
  STORE_RESULTS_EVERY = 2;
end
num_results_to_be_stored = NUM_AL_ITERS / STORE_RESULTS_EVERY;
NUM_INIT_PTS_PER_DIM = 1;
INIT_ON_GRID = false;
MCMC_EST_INIT_PT = zeros(NUM_DIMS, 1);

% Constants for MCMC
NUM_MCMC_SAMPLES = 10*NUM_AL_ITERS;
num_mcmc_results_to_be_stored = NUM_MCMC_SAMPLES / STORE_RESULTS_EVERY;

% RUN TIME HYPER-PARAMS
USE_OPT_PARAMS = false;
OPT_BANDWIDTH = 2.3* NUM_AL_ITERS^(-1/(1.3 + NUM_DIMS));
OPT_SCALE = 50;
NUM_KFOLD_CV_PARTITIONS = 20;
functionals = {@f1, @f2, @f3, @f4};
NOISE_LEVEL = 0.05;

% parameters for the test
p1 = 0.5;
p2 = 1 - p1;
sigma = 0.4*sqrt(NUM_DIMS); %+ 0.05*rand();
PARAM_SPACE_BOUNDS = [-3.5*sigma, 1 + 3.5*sigma];
MCMC_EST_PROPOSAL_STD = sigma;
fprintf('Dims: %d, Bounds: %s\n', NUM_DIMS, mat2str(PARAM_SPACE_BOUNDS));

% Create function for evaluating joint likelihood and
% obtain the True functional values
evalLogJoint = @(arg) evalLogLiklExp1(arg, sigma, p1, p2);
[~, f_vals] = evalLogJoint(zeros(0, NUM_DIMS));

% Define other useful variables before proceeding
num_init_pts = NUM_INIT_PTS_PER_DIM ^ NUM_DIMS;

% Define the Cost function for Cross Validation
cv_cost_func = @(y1, y2) (exp(y1) - exp(y2)).^2;

% Create arrays for storing the arrays.
uc_errs.f1 = zeros(NUM_EXPERIMENTS, num_results_to_be_stored);
uc_errs.f2 = zeros(NUM_EXPERIMENTS, num_results_to_be_stored);
uc_errs.f3 = zeros(NUM_EXPERIMENTS, num_results_to_be_stored);
uc_errs.f4 = zeros(NUM_EXPERIMENTS, num_results_to_be_stored);
mcmc_errs.f1 = zeros(NUM_EXPERIMENTS, num_mcmc_results_to_be_stored);
mcmc_errs.f2 = zeros(NUM_EXPERIMENTS, num_mcmc_results_to_be_stored);
mcmc_errs.f3 = zeros(NUM_EXPERIMENTS, num_mcmc_results_to_be_stored);
mcmc_errs.f4 = zeros(NUM_EXPERIMENTS, num_mcmc_results_to_be_stored);
mcmc_reg_errs.f1 = zeros(NUM_EXPERIMENTS, num_results_to_be_stored);
mcmc_reg_errs.f2 = zeros(NUM_EXPERIMENTS, num_results_to_be_stored);
mcmc_reg_errs.f3 = zeros(NUM_EXPERIMENTS, num_results_to_be_stored);
mcmc_reg_errs.f4 = zeros(NUM_EXPERIMENTS, num_results_to_be_stored);
mbp_errs.f1 = zeros(NUM_EXPERIMENTS, num_results_to_be_stored);
mbp_errs.f2 = zeros(NUM_EXPERIMENTS, num_results_to_be_stored);
mbp_errs.f3 = zeros(NUM_EXPERIMENTS, num_results_to_be_stored);
mbp_errs.f4 = zeros(NUM_EXPERIMENTS, num_results_to_be_stored);

% For saving results
save_file_name = sprintf('results/res_d%d_%s.mat', NUM_DIMS, ...
                         datestr(now, 'mm:dd-HH:MM'));

for experiment_iter = 1:NUM_EXPERIMENTS

  fprintf('EXPERIMENT: %d\n======================\n', experiment_iter);

  fprintf('UNCERT-REDUCTION\n');
  UncertaintyReduction;
  uc_errs.f1(experiment_iter, :) = uc_err_prog.f1';
  uc_errs.f2(experiment_iter, :) = uc_err_prog.f2';
  uc_errs.f3(experiment_iter, :) = uc_err_prog.f3';
  uc_errs.f4(experiment_iter, :) = uc_err_prog.f4';

  fprintf('MCMC\n');
  MCMCForPostEstimation;
  mcmc_errs.f1(experiment_iter, :) = mcmc_err_prog.f1';
  mcmc_errs.f2(experiment_iter, :) = mcmc_err_prog.f2';
  mcmc_errs.f3(experiment_iter, :) = mcmc_err_prog.f3';
  mcmc_errs.f4(experiment_iter, :) = mcmc_err_prog.f4';

  fprintf('MCMC-REG\n');
  mcmc_reg_errs.f1(experiment_iter, :) = mcmc_reg_err_prog.f1';
  mcmc_reg_errs.f2(experiment_iter, :) = mcmc_reg_err_prog.f2';
  mcmc_reg_errs.f3(experiment_iter, :) = mcmc_reg_err_prog.f3';
  mcmc_reg_errs.f4(experiment_iter, :) = mcmc_reg_err_prog.f4';

  fprintf('MAX-BAND-POINT');
  MaxBandPointForAL;
  mbp_errs.f1(experiment_iter, :) = mbp_err_prog.f1';
  mbp_errs.f2(experiment_iter, :) = mbp_err_prog.f2';
  mbp_errs.f3(experiment_iter, :) = mbp_err_prog.f3';
  mbp_errs.f4(experiment_iter, :) = mbp_err_prog.f4';

  % Save the results
  save(save_file_name, 'uc_errs', 'mcmc_errs', 'mcmc_reg_errs', 'mbp_errs');
  % copy to a file named temp
  copyCmd = sprintf('cp %s results/temp.mat', save_file_name); system(copyCmd);

end

