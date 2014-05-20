% Performs A Random search for Posterior estimation

randPts = bsxfun( ...
  @plus, bsxfun(@times, rand(NUM_AL_CANDIDATES, NUM_DIMS), ...
         (PARAM_SPACE_BOUNDS(:,2) - PARAM_SPACE_BOUNDS(:,1))' ), ...
  PARAM_SPACE_BOUNDS(:,1)' );
randLogProbs = evalLogJoint(randPts);

rand_err_prog.f1 = zeros(num_results_to_be_stored, 1);
rand_err_prog.f2 = zeros(num_results_to_be_stored, 1);
rand_err_prog.f3 = zeros(num_results_to_be_stored, 1);
rand_err_prog.f4 = zeros(num_results_to_be_stored, 1);

for randIter = 1:num_results_to_be_stored

  curr_num_rand_pts = randIter * STORE_RESULTS_EVERY;
  fprintf('Rand with %d pts: \n', curr_num_rand_pts);
  % The regressors and regressands for the current iteration
  Xtr = randPts(1:curr_num_rand_pts, :);
  Ytr = randLogProbs(1:curr_num_rand_pts);

  % GP Regression
  logJointEst = regressionWrap(Xtr, Ytr, NOISE_LEVEL, LOWEST_LOGLIKL_VAL, ...
    LOGLIKL_RANGE, cv_cost_func);

  % Now perform MCMC
  estimate_mcmc_samples = CustomMCMC( ...
    NUM_MCMC_ITERS_FOR_EST + NUM_MCMC_BURNIN_FOR_EST, ...
    MCMC_EST_PROPOSAL_STD, MCMC_EST_INIT_PT, logJointEst);
  Xmcmc = estimate_mcmc_samples(NUM_MCMC_BURNIN_FOR_EST+1:end, :);
  % Now evaluate the functionals on Xmcmc
  T1 = f1(Xmcmc); rand_err_prog.f1(randIter)= abs(T1 -f_vals.S1)/abs(f_vals.S1);
  T2 = f2(Xmcmc); rand_err_prog.f2(randIter)= abs(T2 -f_vals.S2)/abs(f_vals.S2);
  T3 = f3(Xmcmc); rand_err_prog.f3(randIter)= abs(T3 -f_vals.S3)/abs(f_vals.S3);
  T4 = f4(Xmcmc); rand_err_prog.f4(randIter)= abs(T4 -f_vals.S4)/abs(f_vals.S4);
  fprintf('Iter %d emp-errors: f1: %f, f2: %f, f3: %f, f4: %f\n', ...
          randIter, ...
          rand_err_prog.f1(randIter), ...
          rand_err_prog.f2(randIter), ...
          rand_err_prog.f3(randIter), ...
          rand_err_prog.f4(randIter) );

end
