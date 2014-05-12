% A Unit test in 1 dimension for mbpObjective.m, mbpGradient.m, maxBandPoint.m
% and alMaxBandPoint.m

close all;
clear all;
addpath ../LipLibkky/
addpath ../helper/

% Define phi
phi = @(x) exp(x); gradPhi = @(x) exp(x);
% phi = @(x) x^2; gradPhi = @(x) 2*x;
LIPSCHITZ_CONST = 6;
% For the experiments
num_pts = 10;
resolution = 100;

% Test 1: One Dimension
% ---------------------
fprintf('Test in 1 Dimension:\n');
f = @(t) t.^2 /4 + sin(7*t);
X = rand(num_pts, 1);
% X = [0.1 0.6 0.95]';
y = f(X);
obj = @(t) -mbpObjective(t, X, y, phi, LIPSCHITZ_CONST, [0 1]);
gradObj = @(t) -mbpGradient(t, X, y, phi, gradPhi, LIPSCHITZ_CONST);
% Now determine init_pt;
init_pt = rand();
gd_params.num_iters = 100;
gd_params.init_step_size = 0.1;
% Run Gradient Descent
[fmin, xmin] = gradientDescent(obj, gradObj, init_pt, gd_params);
fprintf('Init Pt: %.4f, Min-Pt: %0.4f\n\n', init_pt, xmin);
% Create some plots
th = linspace(0, 1, resolution)';
figure;
subplot(1,2,1);
  pivots = repmat(y', resolution, 1);
  diffs = sqrt( dist2(th, X) );
  upper_bounds = pivots + LIPSCHITZ_CONST * diffs;
  lower_bounds = pivots - LIPSCHITZ_CONST * diffs;
  upper_bounds = min(upper_bounds, [], 2);
  lower_bounds = max(lower_bounds, [], 2);
  plot(th, upper_bounds, 'b'); hold on,
  plot(th, lower_bounds, 'b-.'); hold on,
subplot(1,2,2);
obj_th = zeros(resolution, 1);
for i = 1:resolution
  obj_th(i) = obj(th(i));
end
plot(th, -obj_th, 'b'); hold on,
plot(init_pt, -obj(init_pt), 'ro');
plot(xmin, -fmin, 'rx');
fprintf('Now Running Max Band Point\n');
mbp_params.bounds = [0 1];
chosen_pt = maxBandPoint(X, y, LIPSCHITZ_CONST, phi, gradPhi, mbp_params);
plot(chosen_pt, -obj(chosen_pt), 'r*', 'MarkerSize', 10);
fprintf('Chosen-Pt: %0.4f\n\n', chosen_pt);

% TEST for alMaxBandPoint
% =======================
figure;
plot(th, f(th), 'b'); hold on,
almbp_params = gd_params;
[mbp_pts, mbp_vals, mbp_lipschitz_const] = alMaxBandPoint( ...
  f, [], [], phi, gradPhi, LIPSCHITZ_CONST, [0 1], ...
  30, almbp_params);
plot(mbp_pts, mbp_vals, 'rx');
title_str = sprintf('L = %f', LIPSCHITZ_CONST);
title(title_str);

