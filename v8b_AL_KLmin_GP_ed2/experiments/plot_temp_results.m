PLOT_ERROR_BARS = 1;

% load temp_res
% load modbeta1 
q = 60;
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
start_idx = floor(q*.2);
end_idx = ceil(q*.8);
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
plm(alm, 'c-d');
plm(ucm, 'b-o');
plm(mcm, 'r-x');
plm(bfm, 'g-s');
% plm(abm, 'c-*');
% plm(mdm, 'm->');
if PLOT_ERROR_BARS
  peb(alm, 'c');
  peb(ucm, 'b');
  peb(mcm, 'r');
  peb(bfm, 'g');
%   peb(abm, 'c');
%   peb(mdm, 'm');
end
% title({'log(KL) vs #pts', 'AL(b-o), Grid Search(g-s), MCMC(r-x), ABC(c-*)'});
% axis([0 100 -10 3]);
axis([0 100 -12 7]);
% axis([0 100 -4 2.5]);
xlabel('Number of Queries');
ylabel('log(KL(true-post || est-post))');
legend('EDR', 'VR', 'MCMC-R', 'RAND');

 
if FALSE,
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
