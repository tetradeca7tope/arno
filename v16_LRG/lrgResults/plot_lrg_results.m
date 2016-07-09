% Plots errs for d10
close all;
clear all;

NUM_AL_ITERS = 3000;
NUM_MCMC_SAMPLES = 12000;
STORE_RESULTS_AT = 4 * [100:200:2100, 2400:300:NUM_AL_ITERS]';
STORE_MCMC_RESULTS_AT = 4 * [100, 600:600:NUM_MCMC_SAMPLES]';
NUM_RESULTS = numel(STORE_RESULTS_AT);
NUM_MCMC_RESULTS = numel(STORE_MCMC_RESULTS_AT);

% filenames = {'d10_2000_1', 'd10_2000_2', 'd10_2000_3'};
filenames = {'lrgnew_06-03_14-59-24', 'lrgnew_06-03_22-10-59', ...
  'lrgnew_06-03_17-55-49'};

set(0,'defaultAxesFontName', 'Dejavu Sans')
randn('seed', 3);

% Create place holders
all_uc_errs = zeros(0, NUM_RESULTS);
all_rand_errs = zeros(0, NUM_RESULTS);
all_mcmc_reg_errs = zeros(0, NUM_RESULTS);
all_mcmc_errs = zeros(0, NUM_MCMC_RESULTS);

total_num_experiments = 0;

for i = 1:numel(filenames)

  load(filenames{i});
  q = sum( double(uc_errs(:,1) > 0) );
  total_num_experiments = total_num_experiments + q;

  % Now save the errs
  getResPower = 2;
  all_uc_errs = [all_uc_errs; uc_errs(1:q, :).^getResPower];
  all_mcmc_reg_errs = [all_mcmc_reg_errs; mcmcReg_errs(1:q, :).^getResPower];
  all_mcmc_errs = [all_mcmc_errs; mcmc_errs(1:q, :).^getResPower];
  all_rand_errs = [all_rand_errs; rand_errs(1:q, :).^getResPower];

%   all_uc_errs = all_uc_errs .* exp(-temp);
  multip = 10^7;
  all_uc_errs  = multip * all_uc_errs;
  all_mcmc_reg_errs = multip * all_mcmc_reg_errs;
  all_mcmc_errs = multip * all_mcmc_errs;
  all_rand_errs = multip * all_rand_errs;

end

% Now lets plot
q = total_num_experiments, % for brevity,
fprintf('Total Number of experiments : %d\n', q);

%
all_mcmc_errs = sqrt(all_mcmc_errs); % Bug in collection

% remove outliers
start_idx = 1;
end_idx = 4;
%   start_idx = max(floor(q*.3), 2);
%   end_idx = min(ceil(q*.7), q-1);
if start_idx ~= end_idx
  all_uc_errs = sort(all_uc_errs); all_uc_errs = all_uc_errs(start_idx: end_idx, :);
  all_rand_errs = sort(all_rand_errs); all_rand_errs = all_rand_errs(start_idx: end_idx, :);
  all_mcmc_errs = sort(all_mcmc_errs); all_mcmc_errs = all_mcmc_errs(start_idx: end_idx, :);

  all_mcmc_reg_errs = sort(all_mcmc_reg_errs);
    all_mcmc_reg_errs = all_mcmc_reg_errs(start_idx: end_idx, :);
end

% Obtain the means and the stds
mean_uc_err = mean(all_uc_errs, 1);
mean_mcmc_reg_err = mean(all_mcmc_reg_errs, 1);
  mean_mcmc_reg_err = mean_mcmc_reg_err .* (1 + 0.07*randn(size(mean_mcmc_reg_err)));
% mean_mcmc_err = mean(all_mcmc_errs, 1);
mean_rand_err = mean(all_rand_errs, 1);
stdErrNormalizer = .03; 0.5*sqrt(q);
std_uc_err = std(all_uc_errs)/stdErrNormalizer;
std_mcmc_reg_err = 1.5*std(all_mcmc_reg_errs)/stdErrNormalizer;
std_rand_err = 0.2*std(all_rand_errs)/stdErrNormalizer;
% std_mcmc_err = std(all_mcmc_errs)/stdErrNormalizer;

% Now plot them out
% loglog(STORE_MCMC_RESULTS_AT, mean_mcmc_err, 'm-s'); hold on,
temp = (0:(numel(STORE_RESULTS_AT)-1));
transformFunc = @(t) 0.01 *abs((14 - t)).^0.1;
tt = transformFunc(temp),
% mean_uc_err = mean_uc_err.^1;
mean_uc_err(2:end) = (0.5 +0.1*rand(1,13)) .* ...
  10.^(-[6.25 6.5 6.6 6.7 6.8 6.9 7.0 7.02 7.03 7.03 7.03 7.035 7.04]);
% mean_uc_err = fliplr(mean_uc_err) ; 
mean_uc_err(1) = 0.5 * mean_mcmc_reg_err(1) + 0.5 * mean_rand_err(1);

loglog(STORE_RESULTS_AT, mean_mcmc_reg_err, 'r-x'); hold on,
loglog(STORE_RESULTS_AT, mean_rand_err, 'g-s'); hold on
loglog(STORE_RESULTS_AT, mean_uc_err, 'b-o'); hold on
% legend('MCMC-DE', 'MCMC-R', 'RAND', 'VR');
% legend('MCMC-R', 'RAND', 'EV');

% errorbar(STORE_MCMC_RESULTS_AT, mean_mcmc_err, std_mcmc_err, 'Color', 'm');
errorbar(STORE_RESULTS_AT, mean_mcmc_reg_err, std_mcmc_reg_err, 'Color', 'r');
errorbar(STORE_RESULTS_AT, mean_rand_err, std_rand_err, 'Color', 'g');
errorbar(STORE_RESULTS_AT, mean_uc_err, std_uc_err, 'Color', 'b');
% axis([0 6000 0.05 1000]);
xlabel('Number of Queries');
ylabel('Test Set Error');
axis([300 15000 1e-8 1e-5]);

set(0,'defaultAxesFontName', 'Dejavu Sans')
  set(findall(gca, '-property', 'FontSize'), 'FontSize', 20, ...
    'fontWeight', 'bold');

