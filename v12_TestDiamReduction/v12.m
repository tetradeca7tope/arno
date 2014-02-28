% Tests the point which minimizes the diameter

GRID_SIZE = 300;

% params for the experiment
L = 5;
y0 = 2;
y1 = 4;
x0 = 0;
x1 = 1;
% Determine the transformation
% phi = @(x) exp(x);
phi = @(x) x.^2;
% phi = @(x) (x);

grid = linspace(0, 1, GRID_SIZE)';

upper_bounds = min(y0 + L*abs(grid - x0), y1 + L*abs(grid - x1));
lower_bounds = max(y0 - L*abs(grid - x0), y1 - L*abs(grid - x1));

plot(grid, lower_bounds, 'b--'); hold on,
plot(grid, upper_bounds, 'b-.');
adverse_points = zeros(GRID_SIZE, 1);
adverse_vals = zeros(GRID_SIZE, 1);

for grid_iter = 1:GRID_SIZE

  possible_vals = linspace(lower_bounds(grid_iter), upper_bounds(grid_iter), ...
                           GRID_SIZE)';
  diameters = zeros(GRID_SIZE, 1);
  for val_iter = 1:GRID_SIZE
    % form upper and lower bounds for each value
    local_ub = possible_vals(val_iter) + L*abs(grid - grid(grid_iter));
    local_lb = possible_vals(val_iter) - L*abs(grid - grid(grid_iter));
    ub =  min(upper_bounds, local_ub);
    lb =  max(lower_bounds, local_lb);

    diffs = phi(ub) - phi(lb);
    diameters(val_iter) = max(diffs);

    if grid_iter == 80 && val_iter == 80
      plot(grid, lb, 'k--'); hold on,
      plot(grid, ub, 'k-.');
      plot(grid, diffs / max(diffs) * 3, 'c');
      diameters(val_iter),
      pause;
    end

  end

  if grid_iter == 80
    f1 = gcf();
    figure;
    plot(diameters); 
    pause
    figure(f1);
  end

  [max_val, max_idx] = max(diameters);
  adverse_points(grid_iter) = possible_vals(max_idx);
  adverse_vals = max_val;
end

plot(grid, adverse_points, 'm');
plot(grid, adverse_vals/ max(adverse_vals) * 3, 'g');
