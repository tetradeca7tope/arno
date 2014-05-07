% Expects that the relevant file has already been loaded

q = sum( double( uc_errs.f1(:,1) > 0));

  figure;

  % remove some of the outliers
  start_idx = 1;
  end_idx = q;
%   start_idx = max(floor(q*.2), 2);
%   end_idx = min(ceil(q*.8), q-1);
  if start_idx ~= end_idx
    uc_err = sort(uc_err); uc_err = uc_err(start_idx:end_idx, :);
    mbp_err = sort(mbp_err); mbp_err = mbp_err(start_idx:end_idx, :);
    mcmc_err = sort(mcmc_err); mcmc_err = mcmc_err(start_idx:end_idx, :);
    mcmcReg_err=sort(mcmcReg_err); mcmcReg_err=mcmcReg_err(start_idx:end_idx,:);
  end
  
  mean_uc_err = mean(uc_err, 1);
  mean_mcmc_err = mean(mcmc_err, 1);
  mean_mbp_err = mean(mbp_err, 1);
  mean_mcmc_reg_err = mean(mcmcReg_err, 1);

  plot(log(mean_mcmc_err), 'b-o'); hold on,
  plot(log(mean_mcmc_reg_err), 'k-<');
  plot(log(mean_uc_err), 'r-x');
  plot(log(mean_mbp_err), 'g-s');
  legend('mcmc', 'mcmc-reg', 'uc', 'mbp');
  title(title_str);

