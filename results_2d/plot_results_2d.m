PLOT_ERROR_BARS = 0;

load results_2d 
num_pts = size(res_al, 2);
q = size(res_al, 1),

al = res_al(1:q, 1:num_pts);
bf = res_bf(1:q, 1:num_pts);
mc = res_mcmc(1:q, 1:num_pts);

pkm = @(x, s) plot(mean(x,1), s); 
plm = @(x, s) pkm(log(x), s); 
peb = @(x, col) errorbar(1:num_pts, mean(log(x)), std(log(x)), 'Color', col);

% sort the results and remove outliers
start_idx = 1; %floor(q*.25);
end_idx = q; %ceil(q*.75);
alm = sort(al); alm = alm(start_idx:end_idx, :); 
bfm = sort(bf); bfm = bfm(start_idx:end_idx, :); 
mcm = sort(mc); mcm = mcm(start_idx:end_idx, :); 

figure; hold on,
pkm(alm, 'b-o');
pkm(bfm, 'g-s');
pkm(mcm, 'r-x');
title({'KL vs #pts', 'Active Learning(b-o), Grid Search(g-s), MCMC(r-x)'});

figure; hold on,
plm(alm, 'b-o');
plm(bfm, 'g-s');
plm(mcm, 'r-x');
if PLOT_ERROR_BARS
  peb(alm, 'b');
  peb(bfm, 'g');
  peb(mcm, 'r');
end
title({'log(KL) vs #pts', 'Active Learning(b-o), Grid Search(g-s), MCMC(r-x)'});
 
fkl = figure;
fkll = figure;

for i = 1:q

  a = al(i,:);
  b = bf(i,:);
  c = mc(i,:);

  title_str = sprintf('%d',i);
  figure(fkl);
  pkm(a, 'b-+'); hold on,
  pkm(b, 'g-+');
  pkm(c, 'r-+');
  pkm(al, 'b-o');
  pkm(bf, 'g-s');
  pkm(mc, 'r-x');
  title(title_str);

  figure(fkll);
  plm(a, 'b-+'); hold on,
  plm(b, 'g-+');
  plm(c, 'r-+');
  plm(al, 'b-o');
  plm(bf, 'g-s');
  plm(mc, 'r-x'); hold off,
  title(title_str);

  pause,
end
