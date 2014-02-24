% Use the GP as a probability model for the posterior and use it to compute
% the expected gain in KL in choosing a point.

clear all
close all;
warning('off', 'all');
addpath ../helper/
% addpath /home/kirthevasan/HPrograms/Matlab/third_party_libs/gpml/gpml-matlab-v3.4-2013-11-11
% startup;

% Parameters
N = 100; % Number of samples
num_iters = 1;
num_init_pts = 3;
num_candidates = 50;
num_AL_iters = 15;
% Parameters for the GP
sigma_gp_sm_constant = 0.2; % the GP smoothness sigma will be
                            % sigma_gp_sm_constant/num_regression_pts
sigma_gp_pr_constant = 6; % the GP prior std-dev will be sigma_gp_pr_constant/
                          % std-dev(regressands)

% constants
noise_level = 0;
% Define this vector th. It serves purposes of plotting and integration etc.
th = linspace(0,1,200)';
% Define the set of candidates
candidates = linspace(0,1,num_candidates+2)'; candidates = candidates(2:end-1);

%% Generate Data
% First plot the prior
alp = 2; bep = 2; % alpha and beta for the prior
subplot(2,2,1);
sample_beta_distro; % This m-file samples a set of points

% Now obtain some (noisy) estimates for L(theta)
init_pts = linspace(0,1, num_init_pts+2)'; init_pts = init_pts(2:end-1);
% init_pts = linspace(0.01,0.99, num_init_pts)';
subplot(2,2,2);
[obsJLogProbs] = estimate_beta_joint_probs(init_pts, sumX, N, ...
                                        alp, bep, noise_level);
% Set parameters for GP regression
sigma_gp_sm = sigma_gp_sm_constant/num_init_pts;  % smoothness prior for the GP
sigma_gp_pr = sigma_gp_pr_constant*std(obsJLogProbs); % prior std-dev
measurement_noise = (noise_level)*ones(num_init_pts,1);

% Now perform GP regression on the joint log probabilities
[est_log_pXT, ~, K] = GPRegression(init_pts, obsJLogProbs, measurement_noise, ...
                           sigma_gp_sm, sigma_gp_pr, th);
% Use the GPML library to perform GP regression instead of my implementation
% covfunc = @covSEiso;
% likfunc = @likGauss;
% hyp2.cov = [0; 0]; hyp2.lik = log(0.1);
% hyp2 = minimize(hyp2, @gp, -100, @infExact, [], covfunc, likfunc, ...
%                 cand_pts, obsJLogProbs);
% exp(hyp2.lik)
% nlml2 = gp(hyp2, @infExact, [], covfunc, likfunc, cand_pts, obsJLogProbs);
% [~,~,est_log_pXT,~,~,~] = gp(hyp2, @infExact, [], covfunc, likfunc, ...
%                              cand_pts, obsJLogProbs, th);
% K11 = feval(covfunc, hyp2.cov, cand_pts);
% K12 = feval(covfunc, hyp2.cov, cand_pts, th);
% K22 = feval(covfunc, hyp2.cov, th);
% K = K22 - K12' * pinv(K11) * K12;
% subplot(2,2,2);
% plot(th, est_log_pXT, 'y', 'LineWidth', 6); hold on;
% plot(cand_pts, obsJLogProbs, 'k+', 'MarkerSize', 10, 'LineWidth', 12);

% now numerically obtain P(Z)
pX = numerical_1D_integration(th, exp(est_log_pXT));
logPX = log(pX);
fprintf('Estimated P(X) = %f, log(P(X)) = %f\n', pX, logPX);

% now obtain the estimated posterior
est_log_post = est_log_pXT - logPX;
est_post = exp(est_log_post);

% plot the posterior
subplot(2,2,4);
plot(th, est_post, 'k--');
hold on,
% also plot the true posterior
prior = th.^(alp-1) .* (1 - th).^(bep -1) / beta(alp, bep);
true_post = th.^(alp+sumX-1) .* (1 - th).^(bep + N - sumX -1) / ...
            beta(alp+sumX, bep+N-sumX);
% plot(th, prior, 'b');
plot(th, true_post, 'g-', 'LineWidth', 2);
% Finally plot the points we probed on the posterior plots too.
plot(init_pts, exp(obsJLogProbs - logPX), 'kx', 'MarkerSize', 4);

% Perform this last test. We want to compute the KL between the mean of
% posterior GP and the sampling distribution of the true posterior.
% For this we generate random samples from the GP distribution and then
% compute the expected KL between
num_sample_fns_KL = 150; 
U = randn(size(th,1), num_sample_fns_KL); 
V = real(bsxfun(@plus, (sqrtm(K)*U), est_log_pXT));
subplot(2,2,3);
hold on,
for i = 1:num_sample_fns_KL
  plot(th, V(:,i), 'm');
end
% now obtain the posterior probabilities corresponding to these points
sample_log_pX = log(numerical_1D_integration(th, exp(V)));
est_sample_log_post = bsxfun(@minus, V, sample_log_pX);
est_sample_post = exp(est_sample_log_post);
pointwise_KL = bsxfun(@times, ...
                      bsxfun(@plus, -est_sample_log_post, est_log_post), ...
                      est_post);
sample_KLs = numerical_1D_integration(th, pointwise_KL);
estimated_KL = mean(sample_KLs);
fprintf('Estimated KL between posterior mean and samples : %f, std: %f\n', ...
        estimated_KL, std(sample_KLs));

% Now compute the KL between the true posterior and the mean
ptwise_KL = true_post(2:end-1) .* ...
            (log(true_post(2:end-1)) - est_log_post(2:end-1));
true_KL = numerical_1D_integration(th(2:end-1), ptwise_KL);
fprintf('KL divergence between true-post and post-mean of GP: %f\n', true_KL);

titlestr = sprintf(...
           ['Prior (b), Estimated Post (r-) & True Posterior (g--)\n', ...
            'Prior: Beta(%d, %d), p: %f, Num pts in param space: %d\n', ...
            'Num data points: %d, Num positive points = %d\n', ...
            'Noise level: %f, sigma-gp-sm = %f, sigma-gp-pr = %f\n', ...
            'True KL: %f, Estim KL (using GP): %f'], ...
            alp, bep, p, num_init_pts, ...
            N, sumX, noise_level, sigma_gp_sm, sigma_gp_pr, ...
            true_KL, estimated_KL);
title(titlestr);

%% NOW PERFORM ACTIVE LEARNING
num_rem_candidates = num_candidates;
accumulated_cand_KLs = zeros(num_rem_candidates, 1);
for al_iter = 1:num_AL_iters
  fprintf('AL iter %d, ', al_iter);
  for cand_iter = 1 : num_rem_candidates

    for sample_iter = 1:num_sample_fns_KL
  %     subplot(2,2,3); plot(th, V(:,sample_iter), 'k'); pause;
      sample_estimate = V(round(candidates(cand_iter) * size(th,1)), sample_iter);
  %   sample_estimate = estimate_beta_joint_probs(candidates(cand_iter),sumX,N,alp,bep,0,0);
      reg_y = [sample_estimate; obsJLogProbs];
      reg_x = [candidates(cand_iter); init_pts];
      % obtain the GP regression parameters
      gp_pr_param = sigma_gp_pr_constant * std(reg_y);
      gp_sm_param = sigma_gp_sm_constant / size(reg_x,1);
      % Now perform GP
      s_log_pXT = GPRegression(reg_x, reg_y, [], ...
                                 gp_sm_param, gp_pr_param, th);
      s_pX = numerical_1D_integration(th, exp(s_log_pXT));
      s_logPX = log(s_pX);
      % now obtain the estimated posterior
      s_est_log_post = s_log_pXT - s_logPX;
      s_est_post = exp(s_est_log_post);
      % Now compute expected KL
      curr_ptwise_KL = est_sample_post(2:end-1, sample_iter) .* ( ...
        est_sample_log_post(2:end-1, sample_iter) - s_est_log_post(2:end-1));
  %     curr_ptwise_KL = s_est_post(2:end-1) .* ( ...
  %                      s_est_log_post(2:end-1) - est_log_post(2:end-1) );
  %     curr_ptwise_KL = s_est_post(2:end-1) .* ( ...
  %                      s_est_log_post(2:end-1) - log(true_post(2:end-1)) );
      curr_KL = numerical_1D_integration(th(2:end-1), curr_ptwise_KL);
      % finally add this to the curr_KL
      accumulated_cand_KLs(cand_iter) = accumulated_cand_KLs(cand_iter) + ...
                                        curr_KL;
  %     subplot(2,2,3); plot(th, V(:,sample_iter), 'm');
    end  
  end
  expected_cand_KLs = accumulated_cand_KLs / num_sample_fns_KL;
  plot_cand_KLs = accumulated_cand_KLs/max(accumulated_cand_KLs) * ...
                  max(est_post) * 0.5;
  [~, min_idx] = min(expected_cand_KLs);
  subplot(2,2,4);
  plot(candidates, plot_cand_KLs, 'b-o');
%   plot(candidates(min_idx), 0.4, mat2str(al_iter), 'Color', 'r', ...
%       'MarkerSize', 10, 'LineWidth', 3);
  text(candidates(min_idx), 0.4, num2str(al_iter), 'Color', 'r');
  subplot(2,2,3);
  plot(candidates(min_idx), min(est_log_pXT), 'r+', 'MarkerSize', 10, ...
       'LineWidth', 3);

  expected_cand_KLs = accumulated_cand_KLs / num_sample_fns_KL;

  % Finally plot the curve using the new estimate.
  init_pts = [candidates(min_idx); init_pts];
  obsJLogProbs = estimate_beta_joint_probs(init_pts, sumX, N, alp, bep, 0, false);
%   sigma_gp_sm = sigma_gp_sm_constant / size(init_pts,1);
  sigma_gp_pr = sigma_gp_pr_constant * std(obsJLogProbs);
  [est_log_pXT, ~, K] = GPRegression(init_pts, obsJLogProbs, [], sigma_gp_sm, ...
                             sigma_gp_pr, th);
  pX = numerical_1D_integration(th, exp(est_log_pXT));
  logPX = log(pX);
  est_log_post = est_log_pXT - logPX;
  est_post = exp(est_log_post);
  subplot(2,2,4);
  plot(th, est_post, 'c-.');
  % Now obtain samples again
    U = randn(size(th,1), num_sample_fns_KL); 
    V = real(bsxfun(@plus, (sqrtm(K)*U), est_log_pXT));
    subplot(2,2,3);
    hold off,
    for i = 1:num_sample_fns_KL
      plot(th, V(:,i), 'c'); hold on,
    end
    % now obtain the posterior probabilities corresponding to these points
    sample_log_pX = log(numerical_1D_integration(th, exp(V)));
    est_sample_log_post = bsxfun(@minus, V, sample_log_pX);
    est_sample_post = exp(est_sample_log_post);
    pointwise_KL = bsxfun(@times, ...
                          bsxfun(@plus, -est_sample_log_post, est_log_post), ...
                          est_post);
    sample_KLs = numerical_1D_integration(th, pointwise_KL);
    estimated_KL = mean(sample_KLs);
    fprintf('Estimated KL: %f, std: %f. End of iteration.\n', ...
            estimated_KL, std(sample_KLs));
  % remove the point from our current estimate
  candidates = candidates([1:min_idx-1, min_idx+1:end]);
  num_rem_candidates = num_rem_candidates - 1;
  accumulated_cand_KLs = zeros(num_rem_candidates, 1);
%   pause,
end

% finally plot the end result
figure;
plot(th, est_post, 'r--', 'LineWidth', 2); hold on,
plot(th, true_post, 'g-', 'LineWidth', 2)'
title('Final Result');
