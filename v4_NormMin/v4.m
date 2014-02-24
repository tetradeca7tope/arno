% A matlab script to test ideas on using a kernel representation for the
% posterior density.

clear all;
close all;
addpath ../helper/

% parameters
N = 10; %number of samples
num_iters = 2;
num_cand_pts_per_iter = 4;
noise_level = 0;
complexity_l1_pen = 0;

%constants
sigma_k = 2/num_cand_pts_per_iter;

%% Generate Data
alp = 3; bep = 2;
sample_beta_distro; % Generate a set of points in X;

% Now obtain some estimates for L(theta);
cand_pts = linspace(0,1,num_cand_pts_per_iter+2)'; cand_pts = cand_pts(2:end-1);
subplot(2,2,1);
estimate_joint_prob_v4;

%% TODO
% 1. Write CVX function to estimate the posterior probability.
% 2. Write script to sample from the posterior

% Initialize the loop.
curr_num_cand_pts = 0;
 
for iter = 1:num_iters

  curr_num_cand_pts = curr_num_cand_pts + num_cand_pts_per_iter;
%   sigma_k = 1/ curr_num_cand_pts; %TODO: use a better rule here.

  % Apply the logit transform to the points we have
  zeta = log(cand_pts ./ (1-cand_pts) );
  % First generate the Kernel matrix using the points we have.
  D = dist2(zeta, zeta);
  K1 = exp(-0.5*D/sigma_k^2);
  K = bsxfun(@rdivide, K1, cand_pts .* (1-cand_pts) );

  % estimate the joint probabilities
  subplot(2,2,1);
  estimate_joint_prob_v4;

  % Now run the following CVX program to obtain the normalized probs
  cvx_begin
  variable w(curr_num_cand_pts);
  minimize( norm(K*w - estNormProbs) + complexity_l1_pen * norm(w,1));
  subject to
    w >= 0;
  cvx_end

  % finally renormalize the weights.
  w = w / sum(w);
  subplot(2,2,2); bar(w);

  % now plot the posterior and the estimated posterior.
  subplot (2,2,3);
  th = linspace(0,1,227)'; th = th(2:end-1);
  probs = zeros(size(th));
  for i = 1:curr_num_cand_pts
    probs = probs + w(i) * (1/sqrt(2*pi*sigma_k^2)) * ...
            exp(-(logit(th) - logit(cand_pts(i)) ).^2 / (2*sigma_k^2) ) ;
  end
  probs = probs ./ ( th .* (1 - th) );
  % Also plot the true posterior
  truepost = th.^(alp+sumX-1) .* (1 - th).^(bep + N - sumX -1) / ...
             beta(alp+sumX, bep+N-sumX);
  plot(th, probs, 'r--', th, truepost, 'g-');

  % Now sample points from here.
  new_cand_pts = zeros(num_cand_pts_per_iter,1);
  num_pts_per_class = mnrnd(num_cand_pts_per_iter, w);
  count = 0;
  for i = 1:curr_num_cand_pts
    for j = 1:num_pts_per_class(i)
      count = count + 1;
      new_cand_pts(count) = logit(cand_pts(i)) + sigma_k*randn(); 
    end
  end
  % now transform them back to the [0 1] space
  new_cand_pts = exp(new_cand_pts) ./ ( 1 + exp(new_cand_pts) );
  hold on,
  plot(new_cand_pts, 0.2*rand(num_cand_pts_per_iter,1), 'mx');
  
  % finally add them to the cand_pts we have.
  cand_pts = [cand_pts; new_cand_pts];

  % Estimate the KL between the distributions
  KL = numerical_1D_integration(th, truepost .* log(truepost ./ probs) );
  fprintf('iter#: %d, Estimated KL: %f\n', iter, KL);
end



