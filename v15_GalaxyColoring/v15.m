% v15: Galaxy Colouring

clear all;
close all;
rng('shuffle'); % Shuffle the seed for Matlab's RNG
addpath ../helper/
addpath ../GPLibkky/
addpath ../LipLibkky/

% Problem Dependent Constants
% NUM_DIMS = 1;
% PARAM_SPACE_BOUNDS = [-1 1];
% NUM_DIMS = 2;
% PARAM_SPACE_BOUNDS = [-1 1; -2 1];
NUM_DIMS = 3;
PARAM_SPACE_BOUNDS = [-1 1; -2 1; -1 0];
% Would be best if all bounds are in [0 1]. Lets work on this later.
INIT_LIPSCHITZ_CONST = 10; % TODO @Ying: Need an estimate for this.
GP_NOISE_LEVEL = 1; % For GP Regression on the selected points
LOWEST_LOGLIKL_VAL = -100; % TODO @Ying Need an estimate for this
LOGLIKL_RANGE = 100; % TODO @Ying: need an estimate for this.

DEBUG_MODE = false;
% DEBUG_MODE = true;
if ~DEBUG_MODE
  NUM_AL_ITERS = 100;
  NUM_EXPERIMENTS = 1;
else
  NUM_AL_ITERS = 5;
  NUM_EXPERIMENTS = 1;
end

% Function Handle
evalLogJoint = @(arg) galaxyLogLikelihoodWrap(arg, LOWEST_LOGLIKL_VAL);

for experiment_iter = 1:NUM_EXPERIMENTS

  fprintf('MaxBandPointForGalaxyColoring\n');
  MaxBandPointForGalaxyColoring;

end

