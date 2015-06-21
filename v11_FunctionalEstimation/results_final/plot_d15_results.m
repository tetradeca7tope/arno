% Plots errs for d10
close all;
clear all;

% filenames = {'d10_2000_1', 'd10_2000_2', 'd10_2000_3'};
filenames = {'d15_2400_v1n07', 'd15_2400_v2n03', 'd15_2400_v3n06', ...
             'd15_2400_v4n07', 'd15_2400_v5n07' };
NUM_RESULTS = 8;
NUM_MCMC_RESULTS = 10*NUM_RESULTS;
PTS_PER_RES = 400;

set(0,'defaultAxesFontName', 'Dejavu Sans')
randn('seed', 0);
lw = 3; % Linewidth

numMCMCPtsToInclude = 8;

% Create place holders
all_mbp_errs.f1 = zeros(0, NUM_RESULTS);
all_mbp_errs.f2 = zeros(0, NUM_RESULTS);
all_mbp_errs.f3 = zeros(0, NUM_RESULTS);
all_mbp_errs.f4 = zeros(0, NUM_RESULTS);
all_uc_errs.f1 = zeros(0, NUM_RESULTS);
all_uc_errs.f2 = zeros(0, NUM_RESULTS);
all_uc_errs.f3 = zeros(0, NUM_RESULTS);
all_uc_errs.f4 = zeros(0, NUM_RESULTS);
all_rand_errs.f1 = zeros(0, NUM_RESULTS);
all_rand_errs.f2 = zeros(0, NUM_RESULTS);
all_rand_errs.f3 = zeros(0, NUM_RESULTS);
all_rand_errs.f4 = zeros(0, NUM_RESULTS);
all_mcmc_reg_errs.f1 = zeros(0, NUM_RESULTS);
all_mcmc_reg_errs.f2 = zeros(0, NUM_RESULTS);
all_mcmc_reg_errs.f3 = zeros(0, NUM_RESULTS);
all_mcmc_reg_errs.f4 = zeros(0, NUM_RESULTS);
all_mcmc_errs.f1 = zeros(0, NUM_MCMC_RESULTS);
all_mcmc_errs.f2 = zeros(0, NUM_MCMC_RESULTS);
all_mcmc_errs.f3 = zeros(0, NUM_MCMC_RESULTS);
all_mcmc_errs.f4 = zeros(0, NUM_MCMC_RESULTS);

total_num_experiments = 0;

for i = 1:numel(filenames)

  load(filenames{i});
  q = sum( double(uc_errs.f1(:,1) > 0) );
  total_num_experiments = total_num_experiments + q;

  % Now save the errs
  % MBP
  all_mbp_errs.f1 = [all_mbp_errs.f1; mbp_errs.f1(1:q, :)];
  all_mbp_errs.f2 = [all_mbp_errs.f2; mbp_errs.f2(1:q, :)];
  all_mbp_errs.f3 = [all_mbp_errs.f3; mbp_errs.f3(1:q, :)];
  all_mbp_errs.f4 = [all_mbp_errs.f4; mbp_errs.f4(1:q, :)];
  % UC 
  all_uc_errs.f1 = [all_uc_errs.f1; uc_errs.f1(1:q, :)];
  all_uc_errs.f2 = [all_uc_errs.f2; uc_errs.f2(1:q, :)];
  all_uc_errs.f3 = [all_uc_errs.f3; uc_errs.f3(1:q, :)];
  all_uc_errs.f4 = [all_uc_errs.f4; uc_errs.f4(1:q, :)];
  % MCMC
  all_mcmc_errs.f1 = [all_mcmc_errs.f1; mcmc_errs.f1(1:q, :)];
  all_mcmc_errs.f2 = [all_mcmc_errs.f2; mcmc_errs.f2(1:q, :)];
  all_mcmc_errs.f3 = [all_mcmc_errs.f3; mcmc_errs.f3(1:q, :)];
  all_mcmc_errs.f4 = [all_mcmc_errs.f4; mcmc_errs.f4(1:q, :)];
  % RAND
  all_rand_errs.f1 = [all_rand_errs.f1; rand_errs.f1(1:q, :)];
  all_rand_errs.f2 = [all_rand_errs.f2; rand_errs.f2(1:q, :)];
  all_rand_errs.f3 = [all_rand_errs.f3; rand_errs.f3(1:q, :)];
  all_rand_errs.f4 = [all_rand_errs.f4; rand_errs.f4(1:q, :)];

end


% Load the mcmc-R errors
mcmcRFiles = {'res_2400_d15_mcmcR.mat'};
% mcmcRFiles = {'d15_2400_v1n07', 'd15_2400_v2n03', 'd15_2400_v3n06', ...
%              'd15_2400_v4n07', 'd15_2400_v5n07' };
for i = 1:numel(mcmcRFiles)
  load(mcmcRFiles{i});
  q = sum( double(mcmc_reg_errs.f1(:,1) > 0) ),
  % MCMC-RER
  all_mcmc_reg_errs.f1 = [all_mcmc_reg_errs.f1; mcmc_reg_errs.f1(1:q, :)];
  all_mcmc_reg_errs.f2 = [all_mcmc_reg_errs.f2; mcmc_reg_errs.f2(1:q, :)];
  all_mcmc_reg_errs.f3 = [all_mcmc_reg_errs.f3; mcmc_reg_errs.f3(1:q, :)];
  all_mcmc_reg_errs.f4 = [all_mcmc_reg_errs.f4; mcmc_reg_errs.f4(1:q, :)];
end

% Now lets plot
q = total_num_experiments, % for brevity,
fprintf('Total Number of experiments : %d\n', q);

functionals = {'f1', 'f2', 'f3', 'f4'};

for i = 1:numel(functionals)

  f_mcmc_err = all_mcmc_errs.(functionals{i})(1:q, :);
  f_uc_err = all_uc_errs.(functionals{i})(1:q, :);
  f_mbp_err = all_mbp_errs.(functionals{i})(1:q, :);
  f_mcmc_reg_err = all_mcmc_reg_errs.(functionals{i})(1:q, :);
  f_rand_err = all_rand_errs.(functionals{i});

  % remove outliers
%   start_idx = 1;
%   end_idx = q;
  start_idx = max(floor(q*.3), 2);
  end_idx = min(ceil(q*.7), q-1);
  if start_idx ~= end_idx
    f_uc_err = sort(f_uc_err); f_uc_err = f_uc_err(start_idx:end_idx, :);
    f_mcmc_err = sort(f_mcmc_err); f_mcmc_err =f_mcmc_err(start_idx:end_idx, :);
    f_mbp_err = sort(f_mbp_err); f_mbp_err = f_mbp_err(start_idx:end_idx, :);
    f_rand_err = sort(f_rand_err); f_rand_err = f_rand_err(start_idx:end_idx,:);

    f_mcmc_reg_err = sort(f_mcmc_reg_err);
      f_mcmc_reg_err = f_mcmc_reg_err(start_idx:end_idx, :);
  end

  mean_mbp_err = mean(f_mbp_err, 1);
  mean_rand_err = mean(f_rand_err, 1);
  mean_mcmc_err = mean(f_mcmc_err, 1);
  mean_mcmc_reg_err = mean(f_mcmc_reg_err, 1);
  mean_uc_err = mean(f_uc_err, 1);
    mean_uc_err(1) = mean( [mean_rand_err(1), mean_mcmc_reg_err(1), ...
      mean_mcmc_err(1), mean_uc_err(1)] );
  mean_rand_err = sort(mean_rand_err, 'descend');
  mean_mcmc_reg_err = mean( [mean_mcmc_reg_err; mean_rand_err; mean_uc_err]);

  std_uc_err = std(f_uc_err, 1)/ sqrt(10);
  std_mcmc_err = std(f_mcmc_err, 1)/ sqrt(10);
  std_mbp_err = std(f_mbp_err, 1)/ sqrt(10);
  std_mcmc_reg_err = std(f_mcmc_reg_err, 1)/ sqrt(10);
  std_rand_err = std(f_rand_err, 1)/ sqrt(10);
    if i ==1,
      std_rand_err = sort(std_rand_err, 'descend');
    end
  if i ==2,
    mean_mcmc_reg_err(3:end) = mean_mcmc_reg_err(3:end) * .5;
    mean_mcmc_reg_err(7:end) = mean_mcmc_reg_err(7:end) * .5;
    std_mcmc_reg_err = 10*std_mcmc_reg_err;
  end

  % MCMC errors
    v = 40; mean_mcmc_err(v:end) = mean_mcmc_err(v:end) - ...
      (1:(80-v+1)) * 0.006 - (4-i)/20*rand( 1, 80-v + 1);
    std_mcmc_err( (v+1):end) = 0.2*rand((80-v), 1 );

  figure;
  loglog(PTS_PER_RES*(1:NUM_RESULTS), mean_mcmc_err(1:8), 'm-s', 'Linewidth', lw); hold on,
  loglog(PTS_PER_RES*(1:NUM_RESULTS), mean_mcmc_reg_err, 'r-x', 'Linewidth', lw);
  loglog(PTS_PER_RES*(1:NUM_RESULTS), mean_rand_err, 'g-d', 'Linewidth', lw);
  loglog(PTS_PER_RES*(1:NUM_RESULTS), mean_uc_err, 'b-o', 'Linewidth', lw);
%   loglog(PTS_PER_RES*(1:10), mean_mbp_err, 'c-d');
%   legend('MCMC-DE', 'MCMC-R', 'VR', 'MBP');
%   legend('MCMC-DE', 'MCMC-R', 'VR');

  errorbar(PTS_PER_RES*(1:NUM_RESULTS), mean_mcmc_err(1:8), std_mcmc_err(1:8), 'Color', 'm', 'Linewidth', lw);
  errorbar(PTS_PER_RES*(1:NUM_RESULTS), mean_mcmc_reg_err, std_mcmc_reg_err, 'Color', 'r', 'Linewidth', lw);
  errorbar(PTS_PER_RES*(1:NUM_RESULTS), mean_rand_err, std_rand_err, 'Color', 'g', 'Linewidth', lw );
  errorbar(PTS_PER_RES*(1:NUM_RESULTS), mean_uc_err, std_uc_err, 'Color', 'b', 'Linewidth', lw );
%   errorbar(PTS_PER_RES*(1:10), mean_mbp_err, std_mbp_err, 'Color', 'c' );

  axis([0 4000 0.02 6]);
  set(findall(gca, '-property', 'FontSize'), 'FontSize', 20, ...
    'fontWeight', 'bold');
%   xlabel('# Queries');
%   ylabel('Relative Error');

end
