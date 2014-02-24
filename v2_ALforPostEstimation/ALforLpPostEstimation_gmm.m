% This file generates a sample from a GMM and computes an approximate posterior
% The prior for the Multinomial for the GMM is a Dirichlet(2, 3)
% distribution.
clear all;
close all;

% Parameters
N = 10; % Number of samples

%% Generate Data
% First plot the prior and generate theta from the prior
alp = 2; bep = 4; % alpha and beta for the prior
theta = dirichlet_sample([alp, bep]);
th = linspace(0,1);
f1 = figure; plot(th, th.^(alp-1) .* (1 - th).^(bep-1) / beta(alp, bep) );
title('prior');
fprintf('theta = %s\n', mat2str(theta));

% The two distributions are given below.
m1 = 0.2; s1 = 0.06;
m2 = 0.7; s2 = 0.1;
% plot the GMM
x = linspace(0, 1);
p = theta(1)*normpdf(x, m1, s1) + theta(2)*normpdf(x, m2, s2);
figure, plot(x, p); title('density');
X = mvgmmrand([m1, m2]', cat(3, s1^2, s2^2), theta, N);
figure, hist(X),

%% Prepare for Optimization
% Prepare your basis functions
param_base = 2:2:6;
alphas = repmat(param_base, 1, 10)';
betas = repmat(param_base, 10, 1); betas = betas(:);
n = numel(alphas); % the number of functions for basis pursuit
% Obtain the set of candidate points on which to evaluate the pdf
cand_pts = (0.05:0.05:0.95)'; 
m = numel(cand_pts);
% We now evaluate the basis functions at each of the candidate points
% The matrix Q is a mxn matrix
Q = zeros(m, n);
for i = 1:n
  Q(:,i) = cand_pts.^(alphas(i)-1) .* (1 - cand_pts).^(betas(i)-1) / ...
           beta(alphas(i), betas(i));
end
% Now also evaluate P([alpha, beta], X) at each of these points
% compute them in logspace in order to avoid precision issues
logProbs = zeros(m, 1);
for i = 1:m
  t = cand_pts(i);
  lps = t*normpdf(X, m1, s1) + (1 - t)*normpdf(X, m2, s2);
  logProbs(i) = sum(log(lps)) + ...
                log(t.^(alp-1) .* (1 - t).^(bep-1) / beta(alp, bep));
  logProbs(i) = sum(log(lps));
end
% Probs = exp(logProbs);
% compute p(z|theta) * p(theta) for each of these points
% Prep parameters for optimization
l1_pen = 0.5; % penalty for l1 norm of alpha

% % CVX program to solve in Lp norm
% cvx_begin
% variable u(n);
% variable phi;
% % minimize( norm(Q*u - exp(logProbs), 2)  + l1_pen*norm(u,1));
% minimize ( norm( Q*u - exp(log(phi) *logProbs), 2) ) ;
% subject to
%   sum(u) == 1;
%   u >= 0;
% cvx_end
% % u = u / sum(u);

%% CVX program to minimize KL divergence
cvx_begin
variable u(n)
minimize( sum(-entr(Q*u)) - (Q*u)' * logProbs );
subject to
  sum(u) == 1;
  u >= 0;
cvx_end

% plot the resulting posterior 
probs = zeros(size(th));
for i = 1:n
  probs = probs + u(i) * th.^(alphas(i)-1) .* (1 - th).^(betas(i)-1) / ...
          beta(alphas(i), betas(i));
end

figure(f1); hold on,
plot(th, probs, 'r--');
titlestr = sprintf('estimated post(r--) & prior(b-), theta = %s', ...
                   mat2str(theta)); 
title(titlestr);

% Plot u
figure;
plot(u, 'x-');
