% A Unit test in 1 dimension for mbpObjective.m, mbpGradient.m, maxBandPoint.m
% and alMaxBandPoint.m

close all;
clear all;
addpath ../LipLibkky/

% Define phi
phi = @(x) exp(x); gradPhi = @(x) exp(x);
% phi = @(x) x^2; gradPhi = @(x) 2*x;
LIPSCHITZ_CONST = 1;
% For the experiments
NUM_AL_PTS = 103;
num_pts = 10;
resolution = 100;

% Test 2 : Two Dimensions
% -----------------------
fprintf('Test in 2 Dimensions\n');

% specify the function 
f = @(T) diag(T*T') /4 + sin(2*T(:,1) + 2*T(:, 2));
% Plot the function
t1 = linspace(0,1,resolution); [T1 T2] = meshgrid(t1, t1); t = [T1(:), T2(:)];
ft = f(t); fT = reshape(ft, resolution, resolution);
figure;
mesh(T1, T2, fT);

% Evaluate the function at some randomly chosen points
X = rand(num_pts, 2);
y = f(X);
% Plot the upper and lower bounds
  figure;
  ubound = zeros(size(T1));
  lbound = zeros(size(T1));
  for i = 1:resolution
    for j = 1:resolution
      curr_pt = [T1(i,j); T2(i,j)];
      distances = sqrt( dist2(X, curr_pt') );
      ubound(i,j) = min( y + LIPSCHITZ_CONST*distances );
      lbound(i,j) = max( y - LIPSCHITZ_CONST*distances );
    end
  end
  mesh(T1, T2, ubound); hold on
  mesh(T1, T2, lbound);

% TEST for mbpObjective and mbpGradient
% =====================================
init_pt = rand(2, 1);
% init_pt = [0.5 0.5]';
% Set things up for Gradient Descent
gd_params.num_iters = 200;
gd_params.init_step_size = 0.1;
obj = @(t) -mbpObjective(t, X, y, phi, LIPSCHITZ_CONST, [0 1; 0 1]);
gradObj = @(t) -mbpGradient(t, X, y, phi, gradPhi, LIPSCHITZ_CONST);
[fmin, xmin] = gradientDescent(obj, gradObj, init_pt, gd_params);
fprintf('Init_pt: %s, Min-Pt: %s\n\n', mat2str(init_pt), mat2str(xmin));
% Create some plots
figure;
obj_th = zeros(size(T1));
for i = 1:resolution
  for j = 1:resolution
    obj_th(i,j) = obj( [T1(i,j); T2(i,j) ] );
  end
end
contour(T1, T2, -obj_th); hold on;
plot(init_pt(1), init_pt(2), 'ro', 'MarkerSize', 10);
plot(xmin(1), xmin(2), 'rx', 'MarkerSize', 10);
figure;
mesh(T1, T2, -obj_th); hold on;
plot3(init_pt(1), init_pt(2), -obj(init_pt), 'mo', 'MarkerSize', 10);
plot3(xmin(1), xmin(2), -fmin, 'mx', 'MarkerSize', 10);

% TEST for MaxBandPoint
% =====================
fprintf('Test for Max Band Point\n');
mbp_params.bounds = [0 1; 0 1];
chosen_pt = maxBandPoint(X, y, LIPSCHITZ_CONST, phi, gradPhi, mbp_params);
fprintf('Chosen-Pt: %s\n\n', mat2str(chosen_pt));
plot3(chosen_pt(1), chosen_pt(2), -obj(chosen_pt), 'r*', 'MarkerSize', 10);

% TEST for alMaxBandPoint
% =======================
fprintf('Test for alMaxBandPoint\n');
% Determine the initial values
almbp_params = gd_params;
% al_init_pts = []; %  Algorithm picks the init pt
al_init_pts = [0.9 0.9]; % A bad initialization for the function
al_init_vals = f(al_init_pts);

[mbp_pts, mbp_vals, mbp_lipschitz_const] = alMaxBandPoint( ...
  f, al_init_pts, al_init_vals, phi, gradPhi, LIPSCHITZ_CONST, [0 1; 0 1], ...
  NUM_AL_PTS, almbp_params);
figure;
contour(T1, T2, fT); hold on,
plot(mbp_pts(:,1), mbp_pts(:,2), 'rx');
title_str = sprintf('L = %f, (%d pts), init: %s', ...
  mbp_lipschitz_const, NUM_AL_PTS, mat2str(al_init_pts));
title(title_str);

