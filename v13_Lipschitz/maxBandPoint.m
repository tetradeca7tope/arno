function [chosen_pt] = maxBandPoint(X, y, L, phi)
% X and y give the (x,y) pairs at which the function has been evaluated already.
% L is the Lipschitz constraints.
% phi is the transformation we are interested in.

  NUM_GRAD_DESC_TRIALS = 40; % # different initializations for gradient desc.
  GRAD_DESC_TOL = 1e-10; % tolerance for gradient descent.

  % Prelims
  num_dims = size(X, 2);
  m = size(X, 1);

  ax_lb = min(X)';
  ax_ub = max(X)';
  axis_lb = 
  grad_desc_init_pts = 

end
