% A Demo for GP regression.
close all;
clear all;

% Generate data
X = [0.01; 0.25; 0.4; 0.5; 0.6; 0.75;  0.8; 0.99];
m = size(X, 1);
f = @(X)  2.5*(-6*(X - 0.6).^2 + 5*X.*sin(4*X));
y = f(X);
% f = @(X)  X.^2 + X.*sin(4*X);
% m = 8; X = linspace(0,1,m)'; y = f(X) + 0*randn(m,1);
% X = [-4 -3 -1 0 1.5]'; X = (X + 4)/5.5; y = [-2 0 1 2 -1]'; 

sigmaSm = 0.3*std(X);
sigmaPr = 2*std(y);

% Generate test points
N = 100;
Z = linspace(0,1,N)';

% Perform GP-Regression
hyperParams.sigmaSm = sigmaSm;
hyperParams.sigmaPr = sigmaPr;
hyperParams.noise = [];
hyperParams.meanFunc = [];
runtimeParams.plotOn = false;
[post_mean, ~, K] = GPRegression(X, y, Z, hyperParams, runtimeParams);

% Now draw some samples
num_samples = 100;
gp_samples = GPDrawSamples(post_mean, K, num_samples);
% plot the samples
for i = 1:num_samples
  plot(Z, gp_samples(i,:), 'm-'); hold on,
end
% plot the mean on top of it all
% plot(Z, post_mean, 'b', 'LineWidth', 2);
fprintf('Error: %f\n', norm(post_mean - f(Z)));
plot(X, y, 'kx', 'MarkerSize', 10, 'LineWidth', 3);
axis off;

% Now plot the exponentiated samples
figure;
f_samples = exp(gp_samples);
for i = 1:num_samples
  plot(Z, f_samples(i,:), 'm-'); hold on,
end
plot(X, exp(y), 'kx', 'MarkerSize', 10, 'LineWidth', 3);
axis off;
