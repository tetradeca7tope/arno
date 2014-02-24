% Unit test for MCMC

close all;
N = 150000;
dims = 10;
init_pt = zeros(dims, 1);
sigma = 0.2*sqrt(dims) + 0.05*rand();
p1 = 0.7;
p2 = 1 - p1;

% MCMC parameters
proposal_std = 0.5;
evalLogLikl = @(arg) evalLogLiklExp1(arg, sigma, p1, p2);

% 1. MCMC
%%%%%%%%%
% Sample via MCMC
num_burnin_samples = max(100, round(N/2));
num_mcmc_samples =  + N;
mcmc_samples = CustomMCMC(num_burnin_samples + N, proposal_std, init_pt, ...
                          evalLogLikl);
Xmcmc = mcmc_samples(num_burnin_samples+1:end, :);
% Evaluate the functionals on the mcmc samples
T1 = f1(Xmcmc);
T2 = f2(Xmcmc);
T3 = f3(Xmcmc);
T4 = f4(Xmcmc);
fprintf('MCMC: Empirical Vals: f1: %f, f2: %f, f3: %f, f4: %f\n', T1, T2, T3, T4);
% Plot the values
plot(Xmcmc(:,1), Xmcmc(:,2), 'x'); hold on,
num_unique_pts = size(unique(Xmcmc(:,1)), 1);
fprintf('Number of Unique points: %d\n', num_unique_pts);

% 2. Rejection Sampling
%%%%%%%%%%%%%%%%%%%%%%%
% Sample via rejection sampling
% boundaries = [-1, 2];
boundaries = [-10, 10];
max_prob = evalLogLikl(zeros(1, dims));
Xrs = RejectionSamplingUniform(N, evalLogLikl, boundaries, max_prob, dims);
% Evaluate the functionals on the samples
T1 = f1(Xrs);
T2 = f2(Xrs);
T3 = f3(Xrs);
T4 = f4(Xrs);
fprintf('RS  : Empirical Vals: f1: %f, f2: %f, f3: %f, f4: %f\n', T1, T2, T3, T4);
% Plot the values
plot(Xrs(:,1), Xrs(:,2), 'mo'); hold on,

% Finally evaluate the functionals
[~, f_vals] = evalLogLikl(zeros(0, dims));
f_vals,

% For comparison plot the contours of the distribution
RES = 100;
t = linspace(-1, 2, RES)';
[T1, T2] = meshgrid(t, t);
th = [T1(:) T2(:)];
[LogLikl] = evalLogLikl(th);
P = reshape(exp(LogLikl), RES, RES);
contour(T1, T2, P);
fprintf('Area: %f\n', numerical_2D_integration(P, T1, T2));

% Baseline performance
Xbase = gendata(N, dims, sigma, p1);
T1 = f1(Xbase);
T2 = f2(Xbase);
T3 = f3(Xbase);
T4 = f4(Xbase);
fprintf('Baseline: Empirical Vals: f1: %f, f2: %f, f3: %f, f4: %f\n', T1, T2, T3, T4);
 
