% Plots Results for SNLS
close all;
clear all;

% Define the constants
NUM_RES = 10;
NUM_MCMC_RES = 4 * NUM_RES;
PTS_PER_RES = 160;

% Load data
load snls_1600_highbw;
uc2_errs = uc_errs;
load snls_1600_all_2;

q = sum( double( uc_errs(:,1) > 0 ) ),

% remove outliers
start_idx = floor(q*.2);
end_idx = floor(q*.8);

uc_errs = sort(uc_errs); uc_errs = uc_errs(start_idx: end_idx, :);
uc2_errs = sort(uc2_errs); uc2_errs = uc2_errs(start_idx: end_idx, :);
mcmcReg_errs = sort(mcmcReg_errs); mcmcReg_errs = mcmcReg_errs(start_idx: end_idx, :);
rand_errs = sort(rand_errs); rand_errs = rand_errs(start_idx: end_idx, :);
mcmc_errs = sort(mcmc_errs); mcmc_errs = mcmc_errs(start_idx: end_idx, :);
abc_errs = sort(abc_errs); abc_errs = abc_errs(start_idx: end_idx, :);

% Obtain the means and the stds
mean_uc_err = mean(uc_errs);
mean_uc2_err = mean(uc2_errs);
mean_mcmcReg_err = mean(mcmcReg_errs);
mean_rand_err = mean(rand_errs);
mean_mcmc_err = mean(mcmc_errs);
mean_abc_err = mean(abc_errs);
std_uc_err = std(uc_errs);
std_uc2_err = std(uc2_errs);
std_mcmcReg_err = std(mcmcReg_errs);
std_rand_err = std(rand_errs);
std_mcmc_err = std(mcmc_errs);
std_abc_err = std(abc_errs);

% Now plot them out
loglog(PTS_PER_RES * (1:NUM_MCMC_RES), mean_mcmc_err, 'm-s'); hold on,
loglog(PTS_PER_RES * (1:NUM_MCMC_RES), mean_abc_err, 'y->');
loglog(PTS_PER_RES * (1:NUM_RES), mean_mcmcReg_err, 'r-x');
loglog(PTS_PER_RES * (1:NUM_RES), mean_rand_err, 'g-d');
loglog(PTS_PER_RES * (1:NUM_RES), mean_uc_err, 'b-o');
% loglog(PTS_PER_RES * (1:NUM_RES), mean_uc2_err, 'k-o');
legend('MCMC-DE', 'ABC', 'MCMC-R', 'RAND', 'VR');

% Plot errorbars
errorbar(PTS_PER_RES*(1:NUM_MCMC_RES),mean_mcmc_err,std_mcmc_err,'Color','m');
errorbar(PTS_PER_RES*(1:NUM_MCMC_RES),mean_abc_err,std_abc_err,'Color','y');
errorbar(PTS_PER_RES*(1:NUM_RES),mean_mcmcReg_err,std_mcmcReg_err,'Color','r');
errorbar(PTS_PER_RES*(1:NUM_RES),mean_rand_err,std_rand_err,'Color','g');
errorbar(PTS_PER_RES*(1:NUM_RES),mean_uc_err,std_uc_err,'Color','b');
% errorbar(PTS_PER_RES*(1:NUM_RES),mean_uc2_err,std_uc2_err,'Color','k');
axis([0 (NUM_MCMC_RES*PTS_PER_RES + 100) 0.05 1000]);
xlabel('# Queries');
ylabel('KL Divergence');
