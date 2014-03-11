function ret = mbpObjective(x, X, y, phi, L, bounds)
% x is the point at which the function needs to be evaluated.
% This (along with mbpGradient) are to be passed to Gradient Descent as function
% handles.
% X, y, L, phi, bounds are as defined in maxBandPoint.m
  % if x lies outside bounds, return -inf
  if sum(double(x < bounds(:,1)) + double(x > bounds(:,2)))
    ret = -inf;
  else
    distances = sqrt( dist2(X, x') );
    upper_bounds = y + L * distances;
    lower_bounds = y - L * distances;
    ret = phi(min(upper_bounds)) - phi(max(lower_bounds));
  end
end

