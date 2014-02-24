% This file generates a sample from a Bernoulli RV
% clear all;
close all;

% Parameters
N = 10; % Number of samples
num_iters = 1;
num_cand_pts = 10;

%% Generate Data
% First plot the prior
alp = 1; bep = 1; % alpha and beta for the prior
th = linspace(0,1);
f1 = figure; plot(th, th.^(alp-1) .* (1 - th).^(bep-1) / beta(alp, bep) );
title('prior');

% Sample p from this prior
theta = dirichlet_sample([alp, bep]); p = theta(1);
fprintf('p = %f\n', p);
% Now generate N Bernoulli points from this p
X = double(rand(N,1) < p);
sumX = sum(X);
fprintf('Number of positive trials : %d\n', sumX);

%% Prepare Basis functions for Optimization
% Prepare your basis functions
param_base = 1:10;
alphas = repmat(param_base, 1, 10)';
betas = repmat(param_base, 10, 1); betas = betas(:);
coeff1 = N/10;
alphas = coeff1 * [alphas; 10*alphas]; betas = coeff1 * [betas; 10*betas];
n = numel(alphas); % the number of functions for basis pursuit
% Obtain the set of candidate points on which to evaluate the pdf

%% Iterate
cand_pts = linspace(0,1,num_cand_pts+2)'; cand_pts = cand_pts(2:end-1);
for iter = 1:num_iters
  m = numel(cand_pts);
  % We now evaluate the basis functions at each of the candidate points
  % The matrix Q is a mxn matrix
  Q = zeros(m, n);
  for i = 1:n
    Q(:,i) = cand_pts.^(alphas(i)-1) .* (1 - cand_pts).^(betas(i)-1) / ...
             beta(alphas(i), betas(i));
  end
  % Now also evaluate P(X, p) at each of these points
  % compute them in logspace in order to avoid precision issues
  logProbs = sumX*log(cand_pts) + (N - sumX)*log(1 - cand_pts) + ...
             log(cand_pts.^(alp-1) .* (1 - cand_pts).^(bep-1) / ...
                 beta(alp, bep) );
  % logProbs = sum(X)*log(cand_pts) + (N - sum(X))*log(1 - cand_pts);
  % compute p(z|theta) * p(theta) for each of these points
  % Prep parameters for optimization
  l1_pen = 0.5; % penalty for l1 norm of alpha

  %% CVX program to solve in Lp norm
  cvx_begin
  variable u(n);
  variable phi;
  % minimize( norm(Q*u - exp(logProbs), 2)  + l1_pen*norm(u,1));
  minimize ( norm( Q*u - phi*exp(logProbs), 2) ) ;
  subject to
    sum(u) == 1;
    u >= 0;
    phi > 0;
  cvx_end
  
%   %% CVX program to minimize KL divergence
%   cvx_begin
%   variable u(n)
%   minimize( sum(-entr(Q*u)) - (Q*u)' * logProbs );
%   subject to
%     sum(u) == 1;
%     u >= 0;
%   cvx_end

  % plot the resulting posterior 
  probs = zeros(size(th));
  for i = 1:n
    probs = probs + u(i) * th.^(alphas(i)-1) .* (1 - th).^(betas(i)-1) / ...
            beta(alphas(i), betas(i));
  end
  figure(f1); hold on,
  plot(th, probs, 'r--', 'LineWidth', iter);
  
  % Finally sample the new set of candidate points
  
end

% finally obtain the true posterior
truepost = th.^(alp+sumX-1) .* (1 - th).^(bep + N - sumX -1) / ...
           beta(alp+sumX, bep+N-sumX);
figure(f1);
plot(th, truepost, 'g-.');
titlestr = sprintf( ['estimated post(r--) & prior(b-),', ...
             'true post(g-.)\ntheta = %s\nsum(X) = %d, ', ...
             'num basis fns = %d'], ...
             mat2str(theta), sumX, n); 
title(titlestr);

% Plot u
figure;
plot(u, 'x-');
