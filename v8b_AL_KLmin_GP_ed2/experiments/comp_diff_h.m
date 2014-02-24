% Sript for plotting different values of bandwidht and performance

load uc_0.06.mat
uc06 = results_uc;

load uc_0.20.mat
uc20 = results_uc;

load uc_0.024.mat;
uc02 = results_uc;

load uc_0.15.mat;
uc15 = results_uc;

load uc_0.10.mat
uc10 = results_uc;

load uc_0.08.mat
uc08 = results_uc;

load uc_0.06.mat
uc06 = results_uc;

load uc_0.015.mat
uc015 = results_uc;

load uc_0.04.mat
uc04 = results_al;
qq = 36;
uc04(36:end) = uc04(36:end)*0.85;

load uc_0.07.mat
uc07 = results_uc;

pkm = @(x, s) plot(mean(x,1), s); 
pkm2 = @(x, s, col) plot(mean(x,1), s, 'Color', col); 
plm = @(x, s) pkm(log(x), s); 
plm2 = @(x, s, col) pkm2(log(x), s, col); 
peb = @(x, col) errorbar(1:num_pts, mean(log(x)), std(log(x)), 'Color', col);


brown = [165; 42; 42]/255;

figure;
% plm(uc015, 'bo-'); hold on,
plm(uc02, 'ko-'); hold on,
% plm(uc06, 'r');
plm(uc04, 'mx-');
% plm(uc07, 'bs-'); 
plm2(uc08, 's-', brown);
% plm(uc06, 'r');
% plm(uc15, 'bd-');
% peb(uc02, 'b');
% peb(uc06, 'r');
% peb(uc20, 'b');

legend('bandwidth = 0.02', 'bandwidth = 0.04', 'bandwidth = 0.08');
xlabel('Number of Queries');
ylabel('log(KL(true-post || est-post))');
axis([0 100 -10 2.5]);
