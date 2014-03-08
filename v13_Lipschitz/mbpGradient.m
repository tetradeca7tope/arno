function ret = mbpGradient(x, X, y, phi, gradPhi, L)
% x is the point at which the function needs to be evaluated.
% X, y, L, phi, gradPhi, bounds are as defined in maxBandPoint.m
  distances = sqrt( dist2(X, x') );
  upper_bounds = y + L * distances;
  lower_bounds = y - L * distances;
  % find the indices corresponding to the lower and upper bounds.
  [ubound, ub_idx] = min(upper_bounds);
  [lbound, lb_idx] = max(lower_bounds);
  % now obtain the gradient
  ret = gradPhi(ubound) * L * sign(x - X(ub_idx, :)') - ...
        gradPhi(lbound) * L * sign(X(lb_idx, :)' - x);
end

