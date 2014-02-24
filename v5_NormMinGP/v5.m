% A matlab script to test if the GP can be used as a reliable representation
% for the KL divergence

clear all;
close all;
addpath ../helper/

% parameters
N = 10; %number of samples
num_cand_pts = 5;
noise_level = 0;

%constants
sigma_k = 2/num_cand_pts;

%% Generate Data
alp = 3; bep = 2;
sample_beta_distro; % Generate a set of points in X;

% Now obtain some estimates for L(theta);
cand_pts = linspace(0,1,num_cand_pts+2)'; cand_pts = cand_pts(2:end-1);
subplot(2,2,1);
estimate_joint_prob_v5;

%% TODO
% 1. Write CVX function to estimate the posterior probability.
% 2. Write script to sample from the posterior

% Apply the logit transform to the points we have
zeta = logit(cand_pts);
% First generate the Kernel matrix using the points we have.
D = dist2(zeta, zeta);
K1 = exp(-0.5*D/sigma_k^2);
K = bsxfun(@rdivide, K1, cand_pts .* (1-cand_pts) );

% Now run the following CVX program to obtain the normalized probs
cvx_begin
variable w(num_cand_pts);
minimize( norm(K*w - estNormProbs) );
subject to
  w >= 0;
cvx_end

% finally renormalize the weights.
w = abs(w);
w = w / sum(w);
subplot(2,2,2); bar(w);

% now plot the posterior and the estimated posterior.
subplot (2,2,3);
th = linspace(0,1,227)'; th = th(2:end-1);
probs = zeros(size(th));
for i = 1:num_cand_pts
  probs = probs + w(i) * (1/sqrt(2*pi*sigma_k^2)) * ...
          exp(-(logit(th) - logit(cand_pts(i)) ).^2 / (2*sigma_k^2) ) ;
end
probs = probs ./ ( th .* (1 - th) );
% Also plot the true posterior
truepost = th.^(alp+sumX-1) .* (1 - th).^(bep + N - sumX -1) / ...
           beta(alp+sumX, bep+N-sumX);
plot(th, probs, 'r--', th, truepost, 'g-'); hold on,

% Estimate the KL between the distributions
KL = numerical_1D_integration(th, truepost .* log(truepost ./ probs) );
fprintf('True KL: %f\n', KL);

% GP Regression
%%%%%%%%%%%%%%%

sigma_gp = 100/num_cand_pts;
[est_log_pXT_gp, ~, K] = GPRegression(cand_pts, estLogProbs, [], sigma_gp, th);
%numerically obtain P(Z)
pX_gp = numerical_1D_integration(th, exp(est_log_pXT_gp));
logPX_gp = log(pX_gp);
%Estimate Posterior using GPs
est_log_post_gp = est_log_pXT_gp - logPX_gp;
est_post_gp = exp(est_log_post_gp);
%plot this posterior too
plot(th, est_post_gp, 'b-.');

% Now generate a bunch of samples from the GP
num_sample_fns_gp = 23;
U = randn(size(th,1), num_sample_fns_gp);
V = bsxfun(@plus, (sqrtm(K)*U), est_log_pXT_gp);
subplot(2,2,4); hold on,
for i = 1:num_sample_fns_gp
  plot(th, V(:,i), 'm');
end
title('Sample functions from GP');
% Estimate the expected KL between the estimated posterior (using norm min)
% and the distributions in the GP.
sample_log_pX = log(numerical_1D_integration(th, exp(V)));
est_sample_log_post = bsxfun(@minus, V, sample_log_pX);
est_sample_post = exp(est_sample_log_post);
pointwise_KL = bsxfun(@times, ...
                      bsxfun(@plus, -est_sample_log_post, log(probs)), ...
                      probs);
sample_KLs = numerical_1D_integration(th, pointwise_KL);
KL_est_gp = mean(sample_KLs);
fprintf('Expected KL between estimated posterior and GP : %f, std: %f\n', ...
         KL_est_gp, std(sample_KLs));

% Now compute the KL between the estimated posterior (using norm minimization)
KL_est_gpmean = ...
  numerical_1D_integration(th, est_post_gp .* log(est_post_gp./probs));
fprintf('KL between Estimated posterior and GP mean : %f\n', KL_est_gpmean);

subplot(2,2,3);
titlestr = sprintf(...
           ['Estimated Post (r-) & True Posterior (g--), ', ...
            'GP-mean(b-.), ', ...
            'Prior: Beta(%d, %d), \np: %f ', ...
            'Num data points: %d, Num positive points = %d\n', ...
            'Num pts in parameter space: %d, '  ...
            'sigma(GP) = %f, sigma_k = %f\n', ...
            'KL(true, estim) = %f, KL(GP-mean, estim): %f'...
            ' E[KL(GP, estim)] = %f\n'], ...
            alp, bep, p, N, sumX, num_cand_pts, ...
            sigma_gp, sigma_k, KL, KL_est_gp, KL_est_gpmean);
title(titlestr);
