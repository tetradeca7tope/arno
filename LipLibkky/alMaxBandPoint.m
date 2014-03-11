function [mbp_pts, mbp_vals, mbp_lipschitz_const] = alMaxBandPoint( ...
  oracle, initPts, initVals, phi, gradPhi, initLipschitzConst, bounds, ...
  numALIters, params);
% Performs Active Learning using the Maximum Band Point procedure. At each
% iteration it invokes maxBandPoint.m which picks the next point.
% oracle is a function handle to the function to be learned. It takes a nxd
% matrix of n pts and returns an nx1 vector.
% In addition it does the follwing.
% 1. Increases the Lipschitz Constant if it detects any violations.
% 2. If init_pts is empty it initializes at the centre of bounds

  % If initPts not given, initalize at the centre of the rectangle.
  if isempty(initPts)
    initPts = [(bounds(:,1) + bounds(:,2))/2]';
    initVals = oracle(initPts); 
  end
  % If phi is not given, use exp(.)
  if isempty(phi)
    phi = @exp;
    gradPhi = @exp;
  end
  % If params.check_for_lipschitz_violations is not given, assume true.
  if ~isfield(params, 'check_for_lipschitz_violations')
    params.check_for_lipschitz_violations = true;
  end
  
  % Define the following before proceeding
  num_init_pts = size(initPts, 1);
  num_dims = size(initPts, 2);
  mbp_pts = initPts;
  mbp_vals = initVals;
  mbp_lipschitz_const = initLipschitzConst;
  params.bounds = bounds; % As it is, I am passing all of params to maxBandPoint

  for al_iter = 1:numALIters
    next_pt = maxBandPoint(mbp_pts, mbp_vals, mbp_lipschitz_const, phi, ...
                           gradPhi, params);
    next_val = oracle(next_pt');
    fprintf('AL iter %d: val: %0.4f, pt: %s\n', ...
            al_iter, next_val, mat2str(next_pt));

    % Check for violations on the Lipschitz assumption
    % Had you set too small a Lipschitz constant, incrase it.
    if params.check_for_lipschitz_violations
      diffs = bsxfun(@minus, mbp_pts, next_pt');
      distances = sqrt(sum(diffs.^2, 2));
      % Violations in Upper bound
      ubnds_on_next_val = mbp_vals + mbp_lipschitz_const * distances;
      [ubnd, u_idx] = min(ubnds_on_next_val);
      if ubnd < next_val
        new_L = 2*(next_val - mbp_vals(u_idx))/ distances(u_idx);
        fprintf('  Lip-Violn: UB: %0.4f, next_val: %0.4f. ', ubnd, next_val);
        fprintf('Changing L to %0.4f from %0.4f\n', new_L, mbp_lipschitz_const);
        mbp_lipschitz_const = new_L;
      end
      % Violations in Lower bound
      lbnds_on_next_val = mbp_vals - mbp_lipschitz_const * distances;
      [lbnd, l_idx] = max(lbnds_on_next_val);
      if lbnd > next_val
        new_L = (mbp_vals(u_idx) - next_val)/ distances(u_idx);
        fprintf('  Lip-Violn: LB: %0.4f, next_val: %0.4f. ', ubnd, next_val);
        fprintf('Changing L to %0.4f from %0.4f\n', new_L, mbp_lipschitz_const);
        mbp_lipschitz_const = new_L;
      end
    end

    % Finally add the new point to the collection
    mbp_pts = [mbp_pts; next_pt'];
    mbp_vals = [mbp_vals; next_val];
  end
  
end
