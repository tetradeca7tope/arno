% Expects that the relevant file has already been loaded

q = sum( double( uc_errs(:,1) > 0));

  figure;

  % remove some of the outliers
  start_idx = 1;
  end_idx = q;
%   start_idx = max(floor(q*.2), 2);
%   end_idx = min(ceil(q*.8), q-1);
  if start_idx ~= end_idx
    uc_errs = sort(uc_errs); uc_errs = uc_errs(start_idx:end_idx, :);
    mbp_errs = sort(mbp_errs); mbp_errs = mbp_errs(start_idx:end_idx, :);
    mcmc_errs = sort(mcmc_errs); mcmc_errs = mcmc_errs(start_idx:end_idx, :);
    mcmcReg_errs=sort(mcmcReg_errs); mcmcReg_errs=mcmcReg_errs(start_idx:end_idx,:);
  end
  
  mean_uc_errs = mean(uc_errs, 1);
  mean_mcmc_errs = mean(mcmc_errs, 1);
  mean_mbp_errs = mean(mbp_errs, 1);
  mean_mcmc_reg_errs = mean(mcmcReg_errs, 1);

  plot(log(mean_mcmc_errs), 'b-o'); hold on,
  plot(log(mean_mcmc_reg_errs), 'k-<');
  plot(log(mean_uc_errs), 'r-x');
  plot(log(mean_mbp_errs), 'g-s');
  legend('mcmc', 'mcmc-reg', 'uc', 'mbp');
  title(title_str);

