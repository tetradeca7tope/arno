function [samples, num_samples_collected] = ...
  RejectionSampling(numSamples, evalLogProb, boundaries, max_val, dims)
% Implements Rejection Sampling and returns numSamples samples

  HARD_THRESHOLD = 1000; % Number of iterations

  if ~exist('dims', 'var')
    dims = size(boundaries, 1);
  end

  if size(boundaries, 1) == 1
    b = zeros(dims, 2);
    for i = 1:dims
      b(i, :) = boundaries;
    end
    boundaries = b;
  end

  % Now compute the bound M on f(x) / u(x) where f is evalLogProb and u
  % is the uniform distribution on the recangle described by boundaries
  rect_area = prod(boundaries(:,2) - boundaries(:,1));
  unif_height = 1/rect_area;
  M = max_val / unif_height;

  samples = zeros(numSamples, dims);
  thresh_counter = 0;
  num_samples_collected = 0;
  thresh_iter = 0;

  while (num_samples_collected < numSamples) && ...
        (thresh_iter < HARD_THRESHOLD),

    u = rand(numSamples, 1);
    prop = bsxfun(@plus, ...
      bsxfun(@times, rand(numSamples, dims), ...
             (boundaries(:,2) - boundaries(:,1))' ), ...
      boundaries(:,1)'); 
%     prop(1:5,:),
    propLogProb = evalLogProb(prop);
%     size(u), size(propLogProb), size(M), size(unif_height),
    accepted_idxs = (log(u) <= propLogProb - log(M) - log(unif_height));
    accepted_pts = prop(accepted_idxs, :);
    num_accepted_pts = size(accepted_pts, 1);
    num_pts_to_add = min(num_accepted_pts, numSamples - num_samples_collected); 
    pts_to_add = accepted_pts(1:num_pts_to_add, :);
%     size(pts_to_add), size(samples),
    samples(num_samples_collected+1: ...
            num_samples_collected+num_pts_to_add, :) = pts_to_add;
    num_samples_collected = num_samples_collected + num_pts_to_add;

    thresh_iter = thresh_iter + 1;
%     fprintf('Rejection Sampling iter: %d\n', thresh_iter);
  end
  
  samples = samples(1:num_samples_collected, :);

end
