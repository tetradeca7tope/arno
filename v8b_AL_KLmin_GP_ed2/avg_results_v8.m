% Runs v8 multiple times and accumulates the runs.

NUM_EXPERIMENTS_FOR_AVGING = 30;
NUM_AL_ITERS = 100;
% NUM_EXPERIMENTS_FOR_AVGING = 2;
% NUM_AL_ITERS = 3;

results_al = zeros(NUM_EXPERIMENTS_FOR_AVGING, NUM_AL_ITERS);
results_mcmc = results_al;
results_mcmcd = results_al;
results_bf = results_al;
results_abc = results_al;
results_uc = results_al;
results_lh = results_al;

%% MAIN LOOP
tic,
for experiment_iter = 1:NUM_EXPERIMENTS_FOR_AVGING
  fprintf('Experiment: %d / %d\n', experiment_iter, NUM_EXPERIMENTS_FOR_AVGING);
  clearvars -except NUM_EXPERIMENTS_FOR_AVGING NUM_AL_ITERS ...
                    results_al results_mcmc results_bf results_abc ...
                    results_uc results_mcmcd ...
                    experiment_iter;
  close all;
  v8;
%   results_al(experiment_iter, :) = kl_progress';
  results_uc(experiment_iter, :) = uc_kl_progress';
%   results_mcmc(experiment_iter, :) = kl_after_each_sample';
%   results_mcmcd(experiment_iter, :) = mcmcd_kl_progress';
%   results_bf(experiment_iter, :) = gr_kl_progress'; 
%   results_abc(experiment_iter, :) = abc_kl_progress';
  results_lh(experiment_iter, :) = lh_kl_progress';
  save('experiments/temp_res.mat', 'results_al', 'results_bf', ...
       'results_mcmc', 'results_abc', 'results_uc', 'results_mcmcd', ...
       'results_lh');
%   save('experiments/mcmcd_abc.mat', 'results_abc', 'results_mcmcd');
end
toc,

figure;
plot(mean(log(results_al)), 'bo-'); hold on,
plot(mean(log(results_mcmc)), 'rx-');
plot(mean(log(results_bf)), 'gs-');
plot(mean(log(results_uc)), 'c*-');
title('Comparison: AL(b-o), MCMC(r-x), Brute-Force(g-s), ABC(c-*)');

filestr = sprintf('experiments/test-%s.mat', datestr(now, 'yyyy-mm-dd-HH:MM'));
save(filestr, 'results_al', 'results_bf', ...
     'results_mcmc', 'results_abc', 'results_uc', 'results_mcmcd', ...
     'results_lh');

