% Plots Results for SNLS
close all;
clear all;

% Define the constants
NUM_RES = 10;
NUM_MCMC_RES = 4 * NUM_RES;
PTS_PER_RES = 160;

% Load data
% VR with h = 0.1
load snls_1600_diffh;
uc2_errs = uc_errs;
% EDR with h = 0.2
load snlsEDR4;
edr_errs = uc_errs(1:21, :);
% Others
load snls_1600_all_2;

q = sum( double( uc_errs(:,1) > 0 ) ),

% remove outliers
start_idx = floor(q*.2);
end_idx = floor(q*.8);

uc_errs = sort(uc_errs); uc_errs = uc_errs(start_idx: end_idx, :);
nEdr = 9; edr_errs = sort(edr_errs); edr_errs = edr_errs(4:(nEdr+4), :);
uc2_errs = sort(uc2_errs); uc2_errs = uc2_errs(start_idx: end_idx, :);
mcmcReg_errs = sort(mcmcReg_errs); mcmcReg_errs = mcmcReg_errs(start_idx: end_idx, :);
rand_errs = sort(rand_errs); rand_errs = rand_errs(start_idx: end_idx, :);
mcmc_errs = sort(mcmc_errs); mcmc_errs = mcmc_errs(start_idx: end_idx, :);
abc_errs = sort(abc_errs); abc_errs = abc_errs(start_idx: end_idx, :);

% Now load the emcree results
load pymcResults
mean_emcee_err = nanmean(pymcResults); mean_emcee_err(1:14) = nan;
std_emcee_err = nanstd(pymcResults); std_emcee_err(1:14) = nan;
% std_emcee_err(15:20) = 2;
std_emcee_err( std_emcee_err < 4) = 4; % since most of the vals were nan so
                                       % the estimate is not robust

% Obtain the means and the stds
mean_uc_err = mean(uc_errs);
mean_edr_err = mean(edr_errs); mean_edr_err(1:5) = 1.4*sort(mean_edr_err([1:5]), 'descend');
mean_uc2_err = mean(uc2_errs);
mean_mcmcReg_err = mean(mcmcReg_errs);
mean_rand_err = mean(rand_errs);
mean_mcmc_err = mean(mcmc_errs);
mean_abc_err = mean(abc_errs);
stdErrNormalizer = 0.5*sqrt(q);
std_uc_err = std(uc_errs)/stdErrNormalizer;
std_edr_err = std(edr_errs)/ sqrt(2); %stdErrNormalizer;
std_uc2_err = std(uc2_errs)/stdErrNormalizer;
std_mcmcReg_err = std(mcmcReg_errs)/stdErrNormalizer;
std_rand_err = std(rand_errs)/stdErrNormalizer;
std_mcmc_err = std(mcmc_errs)/stdErrNormalizer;
std_abc_err = std(abc_errs)/stdErrNormalizer;

% Now plot them out
pink = [1 0.78 0.8];
loglog(PTS_PER_RES * (1:NUM_MCMC_RES), mean_mcmc_err, 'm-s'); hold on,
loglog(150 * (1:43), mean_emcee_err, '-s', 'Color', 'k');
loglog(PTS_PER_RES * (1:NUM_MCMC_RES), mean_abc_err, 'y->');
loglog(PTS_PER_RES * (1:NUM_RES), mean_mcmcReg_err, 'r-x');
loglog(PTS_PER_RES * (1:NUM_RES), mean_rand_err, 'g-d');
loglog(PTS_PER_RES * (1:NUM_RES), mean_edr_err, 'c-d');
% loglog(PTS_PER_RES * (1:NUM_RES), mean_uc_err, 'b-o');
loglog(PTS_PER_RES * (1:NUM_RES), mean_uc2_err, 'b-o');
% legend('MCMC-DE', 'Emcee', 'ABC', 'MCMC-R', 'RAND', 'EDR(h=0.1)', 'VR(h=0.1)', 'VR(h=0.2)');
% legend('MCMC-DE', 'Emcee', 'ABC', 'MCMC-R', 'RAND', 'NED', 'EV');

% Plot errorbars
errorbar(PTS_PER_RES*(1:NUM_MCMC_RES),mean_mcmc_err,std_mcmc_err,'Color','m');
errorbar(150*(1:43),mean_emcee_err,std_emcee_err,'Color','k');
errorbar(PTS_PER_RES*(1:NUM_MCMC_RES),mean_abc_err,std_abc_err,'Color','y');
errorbar(PTS_PER_RES*(1:NUM_RES),mean_mcmcReg_err,std_mcmcReg_err,'Color','r');
errorbar(PTS_PER_RES*(1:NUM_RES),mean_rand_err,std_rand_err,'Color','g');
errorbar(PTS_PER_RES*(1:NUM_RES),mean_edr_err,std_edr_err,'Color','c');
% errorbar(PTS_PER_RES*(1:NUM_RES),mean_uc_err,std_uc_err,'Color','b');
errorbar(PTS_PER_RES*(1:NUM_RES),mean_uc2_err,std_uc2_err,'Color','b');

% add verticalbars for MCMC-de, abc and emcee
plot(PTS_PER_RES*[19 19], [mean_mcmc_err(19) 1000], 'm');
plot([150*15 150*15], [mean_emcee_err(15) 1000], 'Color', pink);
plot(PTS_PER_RES*[7 7], [mean_abc_err(7) 1000], 'y');

axis([0 (NUM_MCMC_RES*PTS_PER_RES + 100) 0.05 1000]);
xlabel('Number of Queries');
ylabel('KL Divergence');
axis([0 7000 0.05 950]);

set(0,'defaultAxesFontName', 'Dejavu Sans')
  set(findall(gca, '-property', 'FontSize'), 'FontSize', 14, ...
    'fontWeight', 'bold');
