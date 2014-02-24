%simple 1D example for Active learning for GP regression
clear all;
close all;

% params
N = 100;  % number of candidate points
sigma = 0.1;    % kernel bandwidth
% a = 5.05; b = 10; % parameters to determine variance. a is maximum
%               % variance. a - b is minimum variance.
b = 3; a = 1.01*b; % set a & b this way
C = N/2;  % This is the computational budget

% Generate a uniform grid of points. These are the candidate points.
X = linspace(0,1,N+2)'; X = X(2:end-1);
D = dist2(X,X);
K = exp(-0.5*D/sigma^2);

% Use the following SDP to determine rho
cvx_begin SDP
variable rho(N);
minimize (lambda_max(K + diag(a - b*rho) ) );
subject to
  0 <= rho <= 1;
  sum(rho) == C;
cvx_end

figure;
plot(rho),

% Now simulate GP regression for a simple function.
% we will regress on noisy estimates of 10*x^2 - 5x + 1
y_ = 10*X.^2 - 5*X + 1;
y = y_ + diag(a - b*rho)*randn(N,1); % add noise to it.

% Now obtain plots of the estimate
npts = 100;
z = linspace(0,1, npts)';
[mu, var] = GPRegression(X, y, a - b*rho, sigma, z);
% D11 = dist2(X, X); K11 = exp(-0.5*D11/sigma^2);
% D22 = dist2(z, z); K22 = exp(-0.5*D22/sigma^2);
% D12 = dist2(X, z); K12 = exp(-0.5*D12/sigma^2);
% 
% % obtain outputs
% mu = K12' * ( (K11 + diag(a - b*rho)*eye(N)) \ y);
% K = K22 - K12'* ((K11 + diag(a - b*rho)*eye(N)) \ K12);
% var = 2*sqrt(diag(K));
% 
% figure;
% plot(X, y, 'kx', 'MarkerSize', 10); hold on,
% plot(z, mu)
% plot(z, mu + var, 'g--');
% plot(z, mu - var, 'g--');
%finally also plot the original function
plot(z, 10*z.^2 - 5*z + 1, 'r-.');

