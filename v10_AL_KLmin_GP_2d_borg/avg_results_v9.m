% Runs v8 multiple times and accumulates the runs.

rng('shuffle');

NUM_EXPERIMENTS_FOR_AVGING = 2;
NUM_AL_ITERS = 100;
NUM_SAMPLES_FOR_OTHER_EXPERIMENTS = 3*NUM_AL_ITERS;

results_al = zeros(NUM_EXPERIMENTS_FOR_AVGING, NUM_AL_ITERS);
results_mcmc = zeros(NUM_EXPERIMENTS_FOR_AVGING, ...
                     NUM_SAMPLES_FOR_OTHER_EXPERIMENTS);
results_bf = zeros(NUM_EXPERIMENTS_FOR_AVGING, ...
                   NUM_SAMPLES_FOR_OTHER_EXPERIMENTS);

tic,
for experiment_iter = 1:NUM_EXPERIMENTS_FOR_AVGING
  fprintf('Experiment: %d / %d\n', experiment_iter, NUM_EXPERIMENTS_FOR_AVGING);
  clearvars -except NUM_EXPERIMENTS_FOR_AVGING NUM_AL_ITERS ...
                    NUM_SAMPLES_FOR_OTHER_EXPERIMENTS ...
                    results_al results_mcmc results_bf experiment_iter;
  close all;
  v9;
  results_al(experiment_iter, :) = KL_progress';
  results_mcmc(experiment_iter, :) = kl_after_each_sample';
  results_bf(experiment_iter, :) = gr_kl_progress';
  [~, hostname] = system('hostname');
  savefile_str = sprintf('experiments/temp_res_%s.mat', hostname);
  save(savefile_str, 'results_al', 'results_bf', 'results_mcmc');
end
toc,

figure;
plot(mean(log(results_al)), 'bo-'); hold on,
plot(mean(log(results_mcmc)), 'rx-');
plot(mean(log(results_bf)), 'gs-');
title('Comparison: AL(b-o), MCMC(r-x), Brute-Force(g-s)');

filestr = sprintf('experiments/test-%s.mat', datestr(now, 'yyyy-mm-dd-HH:MM'));
save(filestr, 'results_al', 'results_bf', 'results_mcmc');

