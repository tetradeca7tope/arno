% Unit tests for mbpObjective and mbpGradient

close all;

% Define phi
phi = @(x) exp(x); gradPhi = @(x) exp(x);
% phi = @(x) x^2; gradPhi = @(x) 2*x;
LIPSCHITZ_CONST = 9;
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
fprintf('Chosen-Pt: %0.4f\n', chosen_pt);
fprintf('Paused ...\n\n');
pause;
% Finally do alMaxBandPoint
figure;
plot(th, f(th), 'b'); hold on,
[mbp_pts, mbp_vals, mbp_lipschitz_const] = alMaxBandPoint( ...
  f, [], [], phi, gradPhi, LIPSCHITZ_CONST, [0 1], ...
  30, params);
plot(mbp_pts, mbp_vals, 'rx');

% Test 2 : Two Dimensions
% -----------------------
fprintf('Test in 2 Dimensions\n');
f = @(T) diag(T*T') /4 + sin(2*T(:,1) + 2*T(:, 2));
X = rand(num_pts, 2);
y = f(X);
% First of all, Plot the function
t1 = linspace(0,1,resolution); [T1 T2] = meshgrid(t1, t1); t = [T1(:), T2(:)];
ft = f(t); fT = reshape(ft, resolution, resolution);
figure;
mesh(T1, T2, fT);
init_pt = rand(2, 1);
% init_pt = [0.5 0.5]';
% First Plot the upper and lower bounds
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
% Run MaxBand Point
fprintf('Now Running Max Band Point\n');
mbp_params.bounds = [0 1; 0 1];
chosen_pt = maxBandPoint(X, y, LIPSCHITZ_CONST, phi, gradPhi, mbp_params);
fprintf('Chosen-Pt: %s\n', mat2str(chosen_pt));
plot3(chosen_pt(1), chosen_pt(2), -obj(chosen_pt), 'r*', 'MarkerSize', 10);
% Finally do alMaxBandPoint
figure;
% plot(th, f(th), 'b'); hold on,
[mbp_pts, mbp_vals, mbp_lipschitz_const] = alMaxBandPoint( ...
  f, [], [], phi, gradPhi, LIPSCHITZ_CONST, [0 1], ...
  100, params);
plot3(mbp_pts(:,1), mbp_pts(:,2), mbp_vals, 'rx');

