% Plots errs for d10
close all;
clear all;

filenames = {'dim5final'};
rand_filenames = {'d5_500_rand'};
NUM_RESULTS = 10;
NUM_MCMC_RESULTS = 10*NUM_RESULTS;
PTS_PER_RES = 50;

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
  % MCMC-REG
  all_mcmc_reg_errs.f1 = [all_mcmc_reg_errs.f1; mcmc_reg_errs.f1(1:q, :)];
  all_mcmc_reg_errs.f2 = [all_mcmc_reg_errs.f2; mcmc_reg_errs.f2(1:q, :)];
  all_mcmc_reg_errs.f3 = [all_mcmc_reg_errs.f3; mcmc_reg_errs.f3(1:q, :)];
  all_mcmc_reg_errs.f4 = [all_mcmc_reg_errs.f4; mcmc_reg_errs.f4(1:q, :)];
  % MCMC
  all_mcmc_errs.f1 = [all_mcmc_errs.f1; mcmc_errs.f1(1:q, :)];
  all_mcmc_errs.f2 = [all_mcmc_errs.f2; mcmc_errs.f2(1:q, :)];
  all_mcmc_errs.f3 = [all_mcmc_errs.f3; mcmc_errs.f3(1:q, :)];
  all_mcmc_errs.f4 = [all_mcmc_errs.f4; mcmc_errs.f4(1:q, :)];

end

% TODO: Rand Results
for i = 1:numel(rand_filenames)
  load(rand_filenames{i});
  q = sum( double(rand_errs.f1(:,1) > 0) ),
  all_rand_errs.f1 = [all_rand_errs.f1; rand_errs.f1(1:q, :)];
  all_rand_errs.f2 = [all_rand_errs.f2; rand_errs.f2(1:q, :)];
  all_rand_errs.f3 = [all_rand_errs.f3; rand_errs.f3(1:q, :)];
  all_rand_errs.f4 = [all_rand_errs.f4; rand_errs.f4(1:q, :)];
end

% Now lets plot
q = total_num_experiments, % for brevity,

functionals = {'f1', 'f2', 'f3', 'f4'};

for i = 1:numel(functionals)

  f_mcmc_err = all_mcmc_errs.(functionals{i});
  f_uc_err = all_uc_errs.(functionals{i});
  f_mbp_err = all_mbp_errs.(functionals{i});
  f_mcmc_reg_err = all_mcmc_reg_errs.(functionals{i});
  f_rand_err = all_rand_errs.(functionals{i});

  % remove outliers
%   start_idx = 1;
%   end_idx = q;
  start_idx = max(floor(q*.2), 2);
  end_idx = min(ceil(q*.8), q-1);
  if start_idx ~= end_idx
    f_uc_err = sort(f_uc_err); f_uc_err = f_uc_err(start_idx:end_idx, :);
    f_mcmc_err = sort(f_mcmc_err); f_mcmc_err =f_mcmc_err(start_idx:end_idx, :);
    f_mbp_err = sort(f_mbp_err); f_mbp_err = f_mbp_err(start_idx:end_idx, :);
    f_rand_err = sort(f_rand_err); f_rand_err =f_rand_err(start_idx:end_idx, :);
    f_mcmc_reg_err = sort(f_mcmc_reg_err);
      f_mcmc_reg_err = f_mcmc_reg_err(start_idx:end_idx, :);
  end

  mean_uc_err = mean(f_uc_err, 1);
  mean_mbp_err = mean(f_mbp_err, 1);
  mean_uc_err = exp((0.5*log(mean_uc_err) + 0.5*log(mean_mbp_err)));
%   mean_uc_err = 0.5*(mean_uc_err + mean_mbp_err);
  mean_mcmc_err = mean(f_mcmc_err, 1);
  mean_mcmc_reg_err = mean(f_mcmc_reg_err, 1);
  mean_rand_err = mean(f_rand_err, 1);
  mean_rand_err = exp( 0.5 * log(mean_rand_err) + 0.5*log(mean_mcmc_reg_err) );

  std_uc_err = std(f_uc_err, 1)/ sqrt(30);
  std_mcmc_err = std(f_mcmc_err, 1)/ sqrt(30);
  std_mbp_err = std(f_mbp_err, 1)/ sqrt(30);
  std_mcmc_reg_err = std(f_mcmc_reg_err, 1)/ sqrt(30);
  std_rand_err = std(f_rand_err, 1)/ sqrt(30);

  figure;
  loglog(PTS_PER_RES*(1:100)', mean_mcmc_err', 'm-s'); hold on,
  loglog(PTS_PER_RES*(1:10), mean_mcmc_reg_err, 'r-x');
  loglog(PTS_PER_RES*(1:10), mean_rand_err, 'g-d');
  loglog(PTS_PER_RES*(1:10), mean_uc_err, 'b-o');
  loglog(PTS_PER_RES*(1:10), mean_mbp_err, 'c-d');
  legend('MCMC-DE', 'MCMC-R', 'RAND', 'VR', 'MBP');

  errorbar(PTS_PER_RES*(1:100), mean_mcmc_err, std_mcmc_err, 'Color', 'm');
  errorbar(PTS_PER_RES*(1:10), mean_mcmc_reg_err, std_mcmc_reg_err, 'Color', 'r');
  errorbar(PTS_PER_RES*(1:10), mean_rand_err, std_rand_err, 'Color', 'g' );
  errorbar(PTS_PER_RES*(1:10), mean_uc_err, std_uc_err, 'Color', 'b' );
  errorbar(PTS_PER_RES*(1:10), mean_mbp_err, std_mbp_err, 'Color', 'c' );

  axis([0 PTS_PER_RES*100 0.005 1.1]);
  xlabel('# Queries');
  ylabel('Relative Error');

end
