% Expects that the relevant file has already been loaded

functionals = {'f1', 'f2', 'f3', 'f4'};
func_desc = {'sum-i=1tod Xi', 'sum-i=1tod Xi^2', ...
  'sum-i=1to(d-1) Xi*X(i+1)*X(i+2)', 'sum-i=1to(d-2) Xi*X(i+1)'};

q = sum( double( mcmc_errs.f1(:,1) > 0));

for i = 1:numel(functionals)

  figure;
  f_mcmc_err = mcmc_errs.(functionals{i})(1:q, :);
  f_uc_err = uc_errs.(functionals{i})(1:q, :);
  f_mbp_err = mbp_errs.(functionals{i})(1:q, :);
  f_mcmc_reg_err = mcmc_reg_errs.(functionals{i})(1:q, :);

  % remove some of the outliers
%   start_idx = 1;
%   end_idx = q;
  start_idx = max(floor(q*.2), 2);
  end_idx = min(ceil(q*.8), q-1);
  if start_idx ~= end_idx
    f_uc_err = sort(f_uc_err); f_uc_err = f_uc_err(start_idx:end_idx, :);
    f_mcmc_err = sort(f_mcmc_err); f_mcmc_err =f_mcmc_err(start_idx:end_idx, :);
    f_mbp_err = sort(f_mbp_err); f_mbp_err = f_mbp_err(start_idx:end_idx, :);

    f_mcmc_reg_err = sort(f_mcmc_reg_err);
      f_mcmc_reg_err = f_mcmc_reg_err(start_idx:end_idx, :);
  end
  
  mean_uc_err = mean(f_uc_err, 1);
  mean_mcmc_err = mean(f_mcmc_err, 1);
  mean_mbp_err = mean(f_mbp_err, 1);
  mean_mcmc_reg_err = mean(f_mcmc_reg_err, 1);

  plot(log(mean_mcmc_err), 'b-o'); hold on,
  plot(log(mean_mcmc_reg_err), 'k-<');
  plot(log(mean_uc_err), 'r-x');
  plot(log(mean_mbp_err), 'g-s');
  title_str = sprintf('dimensions = %d\nfunctional = %s',  dim, func_desc{i});
  legend('mcmc', 'mcmc-reg', 'uc', 'mbp');
  title(title_str);

end

