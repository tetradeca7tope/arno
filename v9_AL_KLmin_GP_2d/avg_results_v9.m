% Runs v8 multiple times and accumulates the runs.

% rng('shuffle');

NUM_EXPERIMENTS_FOR_AVGING = 20;
NUM_AL_ITERS = 100;
NUM_SAMPLES_FOR_OTHER_EXPERIMENTS = 10*NUM_AL_ITERS;

results_al = zeros(NUM_EXPERIMENTS_FOR_AVGING, NUM_AL_ITERS);
results_mcmc = zeros(NUM_EXPERIMENTS_FOR_AVGING, ...
                     NUM_SAMPLES_FOR_OTHER_EXPERIMENTS);
results_bf = zeros(NUM_EXPERIMENTS_FOR_AVGING, ...
                   NUM_SAMPLES_FOR_OTHER_EXPERIMENTS);
results_abc = zeros(NUM_EXPERIMENTS_FOR_AVGING, ...
                   NUM_SAMPLES_FOR_OTHER_EXPERIMENTS);
results_mcmcd = zeros(NUM_EXPERIMENTS_FOR_AVGING, ...
                   NUM_SAMPLES_FOR_OTHER_EXPERIMENTS);
results_uc = zeros(NUM_EXPERIMENTS_FOR_AVGING, ...
                   NUM_AL_ITERS);

tic,
for experiment_iter = 1:NUM_EXPERIMENTS_FOR_AVGING
  fprintf('Experiment: %d / %d\n', experiment_iter, NUM_EXPERIMENTS_FOR_AVGING);
  clearvars -except NUM_EXPERIMENTS_FOR_AVGING NUM_AL_ITERS ...
                    NUM_SAMPLES_FOR_OTHER_EXPERIMENTS ...
                    results_al results_mcmc results_bf results_abc ...
                    results_uc results_mcmcd ...
                    experiment_iter;
  close all;
  v9;
%   results_al(experiment_iter, :) = KL_progress';
%   results_mcmc(experiment_iter, :) = kl_after_each_sample';
%   results_bf(experiment_iter, :) = gr_kl_progress';
  results_mcmcd(experiment_iter, :) = mcmcd_kl_progress';
%   results_uc(experiment_iter, :) = uc_kl_progress';
  % Name the temp file depending on the name of the machine
  [~, hostname] = system('hostname');
  hostname = strsplit(hostname, '.');
  hostname = hostname{1};
  savefile_str = sprintf('experiments/temp_res_%s.mat', hostname);
  savefile_str = 'experiments/mcmcd.mat';
  save(savefile_str, 'results_al', 'results_bf', 'results_mcmc', ...
       'results_abc', 'results_uc', 'results_mcmcd');
end
toc,

figure;
plot(mean(log(results_al)), 'bo-'); hold on,
plot(mean(log(results_mcmc)), 'rx-');
plot(mean(log(results_bf)), 'gs-');
title('Comparison: AL(b-o), MCMC(r-x), Brute-Force(g-s)');

filestr = sprintf('experiments/test-%s.mat', datestr(now, 'yyyy-mm-dd-HH:MM'));
save(filestr, 'results_al', 'results_bf', 'results_mcmc','results_abc', ...
              'results_mcmcd', 'results_uc');

