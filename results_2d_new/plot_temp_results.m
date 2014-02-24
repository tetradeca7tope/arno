PLOT_ERROR_BARS = 0;

% load temp_res
if ~exist('num_pts', 'var')
 num_pts = size(results_al, 2);
end

al = results_al(1:q, 1:num_pts);
bf = results_bf(1:q, 1:num_pts);
mc = results_mcmc(1:q, 1:num_pts);
ab = results_abc(1:q, 1:num_pts);
uc = results_uc(1:q, 1:num_pts);
md = results_mcmcd(1:q, 1:num_pts);

pkm = @(x, s) plot(mean(x,1), s); 
plm = @(x, s) pkm(log(x), s); 
peb = @(x, col) errorbar(1:num_pts, mean(log(x)), std(log(x)), 'Color', col);

% sort the results and remove outliers
start_idx = 1;
end_idx = q;
% start_idx = floor(q*.2);
% end_idx = ceil(q*.6);
alm = sort(al); alm = alm(start_idx:end_idx, :); 
bfm = sort(bf); bfm = bfm(start_idx:end_idx, :); 
mcm = sort(mc); mcm = mcm(start_idx:end_idx, :); 
abm = sort(ab); abm = abm(start_idx:end_idx, :);
ucm = sort(uc); ucm = ucm(start_idx:end_idx, :);
% mdm = sort(md); mdm = mdm(start_idx:end_idx, :);

figure; hold on,
pkm(alm, 'b-o');
pkm(bfm, 'g-s');
pkm(mcm, 'r-x');
pkm(abm, 'c-*');
pkm(ucm, 'k-d');
% pkm(mdm, 'm->');
title({'KL vs #pts', 'Active Learning(b-o), Grid Search(g-s), MCMC(r-x)'});

figure; hold on,
plm(alm, 'b-o');
plm(bfm, 'g-s');
plm(mcm, 'r-x');
plm(abm, 'c-*');
plm(ucm, 'k-d');
% plm(mdm, 'm->');
if PLOT_ERROR_BARS
  peb(alm, 'b');
  peb(bfm, 'g');
  peb(mcm, 'r');
  peb(abm, 'c');
  peb(ucm, 'k');
%   peb(mdm, 'm');
end
title({'log(KL) vs #pts', 'AL(b-o), Grid Search(g-s), MCMC(r-x), ABC(c-*)'});
 
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
