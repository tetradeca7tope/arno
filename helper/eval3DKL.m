function kl = eval3DKL(p, q, bounds)
% Estimates the 3D KL divergence between two distributions.

  if strcmp(class(p), 'function_handle')
    % Evaluate on a grid
    [g1, g2, g3] = gen3DCentrePointGrid(100, bounds);
    pts = [g1(:), g2(:), g3(:)];
    p = p(pts);
    q = q(pts);
  end
  ptwiseKL = p .* log(p./q);
  vol = prod( bounds(:,2) - bounds(:,1) );
  kl = mean(ptwiseKL) * vol;

end

