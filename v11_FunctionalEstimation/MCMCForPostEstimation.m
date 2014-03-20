% MCMC for Posterior Estimation

% Set parameters for MCMC
PROPOSAL_STD = sigma;
MCMC_INIT_PT = zeros(NUM_DIMS, 1);

% Error tracker for each functional
mcmc_err_prog.f1 = zeros(num_mcmc_results_to_be_stored, 1);
mcmc_err_prog.f2 = zeros(num_mcmc_results_to_be_stored, 1);
mcmc_err_prog.f3 = zeros(num_mcmc_results_to_be_stored, 1);
mcmc_err_prog.f4 = zeros(num_mcmc_results_to_be_stored, 1);

% Generate the samples from MCMC
mcmc_samples = CustomMCMC(NUM_MCMC_SAMPLES, PROPOSAL_STD, MCMC_INIT_PT, ...
  evalLogJoint);

for sample_iter = 1: num_mcmc_results_to_be_stored

  T1 = f1(mcmc_samples(1:sample_iter * STORE_RESULTS_EVERY, :));
  mcmc_err_prog.f1(sample_iter) = abs(T1 - f_vals.S1)/abs(f_vals.S1);
  T2 = f2(mcmc_samples(1:sample_iter * STORE_RESULTS_EVERY, :));
  mcmc_err_prog.f2(sample_iter) = abs(T2 - f_vals.S2)/abs(f_vals.S2); 
  T3 = f3(mcmc_samples(1:sample_iter * STORE_RESULTS_EVERY, :));
  mcmc_err_prog.f3(sample_iter) = abs(T3 - f_vals.S3)/abs(f_vals.S3);
  T4 = f4(mcmc_samples(1:sample_iter * STORE_RESULTS_EVERY, :));
  mcmc_err_prog.f4(sample_iter) = abs(T4 - f_vals.S4)/abs(f_vals.S4);

end

