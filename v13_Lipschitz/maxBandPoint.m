function [chosen_pt] = maxBandPoint(X, y, L, phi, gradPhi, params)
% X and y give the (x,y) pairs at which the function has been evaluated already.
% L is the Lipschitz constraints.
% phi is the transformation we are interested in.
% gradPhi is the gradient of Phi.
% params contains all ancillary information for the algorithm (optional)
% params.bounds is a dx2 matrix giving upper and lower bounds for each dimension
% params.num_grad_desc_init_pts: # pts to init gradient descent.
% params.num_grad_desc_iters: # iterations of gradient descent.
% params.grad_desc_top: tolerance for gradient descent.

  NUM_GRAD_DESC_TRIALS = 40; % # different initializations for gradient desc.
  GRAD_DESC_TOL = 1e-10; % tolerance for gradient descent.

  % Prelims
  num_dims = size(X, 2);
  m = size(X, 1);

  % Create an empty struct if not passed
  if ~exist('params', 'var') params = struct();
  end
  % Determine upper and lower bounds to initialize gradient descent
  if ~isfield(params, 'bounds')
    ax_lb = min(X)';
    ax_ub = max(X)';
    axis_lb = ax_lb - 0.25*(ax_ub - ax_lb);
    axis_ub = ax_ub + 0.25*(ax_ub - ax_lb);
    params.bounds = [axis_lb, axis_ub];
  end
  if ~isfield(params, 'num_grad_desc_init_pts')
    params.num_grad_desc_init_pts = min(10*2^num_dims, 200);
  end
  if ~isfield(params, 'num_grad_desc_iters')
    params.num_grad_desc_iters = 20; 
  end
  if ~isfield(params, 'grad_desc_init_step_size')
    params.grad_desc_init_step_size = 0.1;
  end

  % Set things up for gradient descent
  grad_desc_init_pts = bsxfun( ...
    @plus, bsxfun(@times, rand(params.num_grad_desc_init_pts, num_dims), ...
                  (params.bounds(:,1) - params.bounds(:,2))'), ...
    params.bounds(:, 1)' );
  % Define the objective and the gradient of the objective
  obj = @(t) -mbpObjective(t, X, y, phi, L, [0 1]);
  gradObj = @(t) -mbpGradient(t, X, y, phi, gradPhi, L);
  % Parameters for Gradient Descent
  gd_params.num_iters = params.num_grad_desc_iters;
  gd_params.init_step_size = params.grad_desc_init_step_size;
  
  % For storing the current best result
  curr_min_val = inf;
  curr_min_pt = grad_desc_init_pts(1, :)';
  % Run gradient descent on all the init pts and obtain the maximum.
  for gd_init_iter = 1:params.num_grad_desc_init_pts
    curr_init_pt = grad_desc_init_pts(gd_init_iter, :)'; 
    [fmin, xmin] = gradientDescent(obj, gradObj, curr_init_pt, gd_params);
    if fmin < curr_min_val
      curr_min_val = fmin;
      curr_min_pt = xmin;
    end
  end

  % Finally return the best point
  chosen_pt = curr_min_pt;

end

