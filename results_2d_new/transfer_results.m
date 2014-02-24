% script for transferring results

savefile = 'all_results';
tempfilenames = {'res2d_1.mat', ...
                 'res2d_2.mat', ...
                 'res2d_3.mat' ...
                };
num_al_iters = 50;

r_al = zeros(0, num_al_iters);
r_uc = r_al;
r_bf = r_al;
r_mc = r_al;
r_ab = r_al;
r_md = r_al;

for i = 1:numel(tempfilenames) 
  load(tempfilenames{i});
  r_al = [r_al; results_al];
  r_uc = [r_uc; results_uc];
  r_bf = [r_bf; results_bf];
  r_mc = [r_mc; results_mcmc];
  r_ab = [r_ab; results_abc];
  r_md = [r_md; results_mcmcd];
end

results_al = sort(r_al, 'descend');
results_uc = sort(r_uc, 'descend');
results_bf = sort(r_bf, 'descend');
results_mcmc = sort(r_mc, 'descend');
results_abc = sort(r_ab, 'descend');
results_mcmcd = sort(r_md, 'descend');

num_exps = size(results_al, 1);

load bf_2d;
results_bf = results_bf(1:num_exps, :);

save('final_res2d.mat', 'results_al', 'results_uc', 'results_bf', ...
      'results_mcmc', 'results_abc', 'results_mcmcd');
