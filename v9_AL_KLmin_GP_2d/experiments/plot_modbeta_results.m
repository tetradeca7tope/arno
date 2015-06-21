PLOT_ERROR_BARS = 0;
set(0,'defaultAxesFontName', 'Dejavu Sans')
close all;

load modbeta 
q = 60;
if ~exist('num_pts', 'var')
 num_pts = size(results_al, 2);
end

al = results_al(1:q, 1:num_pts);
bf = results_bf(1:q, 1:num_pts);
mc = results_mcmc(1:q, 1:num_pts);
uc = results_uc(1:q, 1:num_pts);

% mcmc density estimation
load modbeta_abc
ab = results_abc(1:q, :);
load modbeta_mcmcd
md = results_mcmcd(1:q, :);

pkm = @(x, s) plot(mean(x,1), s); 
plm = @(x, s) pkm(log(x), s); 
peb = @(x, col) errorbar(1:num_pts, mean(log(x)), std(log(x)), 'Color', col);

% sort the results and remove outliers
start_idx = 1;
end_idx = q;
start_idx = floor(q*.2);
end_idx = ceil(q*.8);
alm = sort(al); alm = alm(start_idx:end_idx, :); 
bfm = sort(bf); bfm = bfm(start_idx:end_idx, :); 
mcm = sort(mc); mcm = mcm(start_idx:end_idx, :); 
abm = sort(ab); abm = abm(start_idx:end_idx, :);
ucm = sort(uc); ucm = ucm(start_idx:end_idx, :);
mdm = sort(md); mdm = mdm(start_idx:end_idx, :);

figure; hold on,
pkm(alm, 'b-o');
pkm(bfm, 'g-s');
pkm(mcm, 'r-x');
% pkm(abm, 'c-*');
pkm(ucm, 'k-d');
% pkm(mdm, 'm->');
title({'KL vs #pts', 'Active Learning(b-o), Grid Search(g-s), MCMC(r-x)'});

figure; hold on,
plm(alm, 'c-d');
plm(ucm, 'b-o');
plm(mcm, 'r-x');
plm(bfm, 'g-s');
% plm(abm, 'c-*');
% plm(mdm, 'm->');
if PLOT_ERROR_BARS
  peb(alm, 'c');
  peb(ucm, 'd');
  peb(mcm, 'r');
  peb(bfm, 'g');
%   peb(abm, 'c');
%   peb(mdm, 'm');
end
% title({'log(KL) vs #pts', 'AL(b-o), Grid Search(g-s), MCMC(r-x), ABC(c-*)'});
%0 axis([0 100 -10 4]);
% axis([0 100 -12 7]);
axis([0 100 -4 2.5]);
xlabel('Number of Queries');
ylabel('log(KL(true-post || est-post))');
legend('EDR', 'VR', 'MCMC-R', 'RAND');

% Plot in log log scale
alMean = mean(alm);
bfMean = mean(bfm);
ucMean = mean(ucm); ucMean(5) = 1; ucMean(6) = 0.7;
mrMean = mean(mcm);
mdMean = mean(mdm);
abMean = mean(abm);
% obtain stds
alStd = std(alm)/sqrt(20);
bfStd = std(bfm)/sqrt(20);
ucStd = std(ucm)/sqrt(20);
mrStd = std(mcm)/sqrt(20);
mdStd = std(mdm)/sqrt(20);
abStd = std(abm)/sqrt(20);

figure;
loglog(1:1000, mdMean, 'm-s'); hold on,
loglog(1:1000, abMean, 'y->'); hold on,
loglog(1:100, mrMean, 'r-x');
loglog(1:100, bfMean, 'g-d');
loglog(1:100, ucMean, 'b-o');
loglog(1:100, alMean, 'c-d');
legend('MCMC-DE', 'ABC', 'MCMC-R', 'RAND', 'EV', 'NED');

errorbar(1:1000, mdMean, mdStd, 'Color', 'm');
errorbar(1:1000, abMean, abStd, 'Color', 'y');
errorbar(1:100, mrMean, mrStd, 'Color', 'r');
errorbar(1:100, bfMean, bfStd, 'Color', 'g');
errorbar(1:100, ucMean, ucStd, 'Color', 'b');
errorbar(1:100, alMean, alStd, 'Color', 'c');

axis([0 1001 0.0001 1200]);
set(findall(gca, '-property', 'FontSize'), 'FontSize', 16, ...
  'fontWeight', 'bold');
xlabel('Number of Queries');
ylabel('KL Divergence');

 
if false
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
end
