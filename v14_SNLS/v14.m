% Script for running experiments on the SNLS dataset

% Prelims
clear all;
close all;
rng('shuffle'); % shuffle the seed for the RNG
addpath ../helper/
addpath ../GPLibkky/
addpath ../LipLibkky/
addpath ../MCMCLibkky/
addpath ../ABCLibkky/

% Constants
numDims = 3;
lowestLogLiklVal = - 2000;
logLiklRange = 2000;
paramSpaceBounds = [0 1; 0 1; 0 1];

DEBUG_MODE = false;
if ~DEBUG_MODE
  NUM_AL_ITERS = 200;
  NUM_MCMC_ITERS_FOR_EST = 100000;
  NUM_MCMC_BURNIN_FOR_EST = 100;
  NUM_EXPERIMENTS = 1;
  STORE_RESULTS_EVERY = NUM_AL_ITERS/30;
else
  NUM_AL_ITERS = 6;
end

% Constants for MCMC
NUM_MCMC_SAMPLES = 10*NUM_AL_ITERS;

for i = 1:NUM_EXPERIMENTS

  % Now run each experiment

end
