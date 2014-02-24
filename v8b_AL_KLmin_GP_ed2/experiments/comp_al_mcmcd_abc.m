% Script to compare AL vs MCMC-D vs ABC

close all;

load modbeta1;
ral = results_al;

load abc;
rabc = results_abc;

load mcmcd;
rmcmcd = results_mcmcd;

n1 = 100;
n2 = 1000;

q = 60;
start_idx = floor(q*.2);
end_idx = ceil(q*.8);

al = sort(ral); al = al(start_idx:end_idx, :); 
mcmcd = sort(rmcmcd); mcmcd = mcmcd(start_idx:end_idx, :);
abc = sort(rabc); abc = abc(start_idx:end_idx, :);

alm = mean(al);
mcmcdm = mean(mcmcd);
abcm = mean(abc);

loglog(1:n1, alm, 'b-o'); hold on,
loglog(1:n2, mcmcdm, 'm-d');
loglog(1:n2, abcm, 'g-s');


errorbar(1:n1, alm, std(al), 'Color', 'b');
errorbar(1:n2, mcmcdm, std(real(mcmcd)), 'Color', 'm');
errorbar(1:n2, abcm, std(real(abc)), 'Color', 'g');
% peb = @(x, col) errorbar(1:num_pts, mean(log(x)), std(log(x)), 'Color', col);

axis([1, 1000, 1e-6, 10^3.5]);

legend('EDR', 'MCMC-D', 'ABC');
xlabel('Number of Queries');
ylabel('KL(true-post, est-post)');

