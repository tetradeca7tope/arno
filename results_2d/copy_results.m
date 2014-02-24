load results_2d;

load temp;

res_al = [res_al; KL_progress'];
res_bf = [res_bf; gr_kl_progress'];
res_mcmc = [res_mcmc; kl_after_each_sample'];

save('results_2d.mat', 'res_al', 'res_bf', 'res_mcmc');

