% MCMC for Posterior Estimation

% Set parameters for MCMC
PROPOSAL_STD = 2.5*sigma;
MCMC_INIT_PT = zeros(NUM_DIMS, 1);

% Error tracker for each functional
mcmc_err_prog.f1 = zeros(num_mcmc_results_to_be_stored, 1);
mcmc_err_prog.f2 = zeros(num_mcmc_results_to_be_stored, 1);
mcmc_err_prog.f3 = zeros(num_mcmc_results_to_be_stored, 1);
mcmc_err_prog.f4 = zeros(num_mcmc_results_to_be_stored, 1);

% Generate the samples from MCMC
[mcmc_samples, mcmc_queries, mcmc_log_probs] = ...
  CustomMCMC(NUM_MCMC_SAMPLES, PROPOSAL_STD, MCMC_INIT_PT, evalLogJoint);

for mcmc_iter = 1: num_mcmc_results_to_be_stored
  T1 = f1(mcmc_samples(1:mcmc_iter * STORE_RESULTS_EVERY, :));
  mcmc_err_prog.f1(mcmc_iter) = abs(T1 - f_vals.S1)/abs(f_vals.S1);
  T2 = f2(mcmc_samples(1:mcmc_iter * STORE_RESULTS_EVERY, :));
  mcmc_err_prog.f2(mcmc_iter) = abs(T2 - f_vals.S2)/abs(f_vals.S2); 
  T3 = f3(mcmc_samples(1:mcmc_iter * STORE_RESULTS_EVERY, :));
  mcmc_err_prog.f3(mcmc_iter) = abs(T3 - f_vals.S3)/abs(f_vals.S3);
  T4 = f4(mcmc_samples(1:mcmc_iter * STORE_RESULTS_EVERY, :));
  mcmc_err_prog.f4(mcmc_iter) = abs(T4 - f_vals.S4)/abs(f_vals.S4);
end

for mcmc_reg_iter = 1:num_results_to_be_stored

  curr_num_mcmc_pts = mcmc_reg_iter * STORE_RESULTS_EVERY;
  fprintf('MCMC-REG with %d points: \n', curr_num_mcmc_pts);
  % Obtain the regressors and the regressands for the current iteration
  Xtr = mcmc_queries(1:curr_num_mcmc_pts, :);
  Ytr = mcmc_log_probs(1:curr_num_mcmc_pts);

  % GP Regression
  mcmcRegLogJointEst = regressionWrap(Xtr, Ytr, NOISE_LEVEL, ...
    LOWEST_LOGLIKL_VAL, LOGLIKL_RANGE, cv_cost_func);

  % Now, do the MCMC on mcmcRegLogJointEst
  mr_estim_mcmc_samples = CustomMCMC( ...
    NUM_MCMC_ITERS_FOR_EST + NUM_MCMC_BURNIN_FOR_EST, MCMC_EST_PROPOSAL_STD, ...
    MCMC_EST_INIT_PT, mcmcRegLogJointEst);
  Xmcmc = mr_estim_mcmc_samples(NUM_MCMC_BURNIN_FOR_EST+1:end, :);
  % Now evaluate the functions on Xmcmc
  T1 = f1(Xmcmc); 
  T2 = f2(Xmcmc); 
  T3 = f3(Xmcmc); 
  T4 = f4(Xmcmc); 
  mcmc_reg_err_prog.f1(mcmc_reg_iter) = abs((T1-f_vals.S1)/f_vals.S1);
  mcmc_reg_err_prog.f2(mcmc_reg_iter) = abs((T2-f_vals.S2)/f_vals.S2);
  mcmc_reg_err_prog.f3(mcmc_reg_iter) = abs((T3-f_vals.S3)/f_vals.S3);
  mcmc_reg_err_prog.f4(mcmc_reg_iter) = abs((T4-f_vals.S4)/f_vals.S4);
  fprintf('Iter %d emp-errors: f1: %f, f2: %f, f3: %f, f4: %f\n', ...
          mcmc_reg_iter, ...
          mcmc_reg_err_prog.f1(mcmc_reg_iter), ...
          mcmc_reg_err_prog.f2(mcmc_reg_iter), ...
          mcmc_reg_err_prog.f3(mcmc_reg_iter), ...
          mcmc_reg_err_prog.f4(mcmc_reg_iter) );

end

