function post_est = estimatePostProb(Th, obsJointLogProbs)
% Estimates the posterior distribution when you have observed obsJointProbs
% at the points in Th. Th is an num_pts x num_dims matrix and onsJointProbs is
% an num_pts x 1 vector.
% Returns a function handle post_est. post_est takes in an nxnum_dims matrix and
% returns a n-length vector.

  % Set this option to avoid warnings
  options = optimoptions('quadprog','Algorithm','interior-point-convex');

  % Prelims
  NUM_KFCV_PARTS = 10;
  num_data = size(Th, 1);
  num_dims = size(Th, 2);

  % rescale obsJointLogProbs to avoid numerical issues
  sc_log_joint_probs = obsJointLogProbs + max(obsJointLogProbs);
  sc_obs_joint_probs = exp(sc_log_joint_probs);

  % Shuffle the data
  shuffle_order = randperm(num_data);
  Th = Th(shuffle_order, :);
  sc_obs_joint_probs = sc_obs_joint_probs(shuffle_order, :);

  % Determine candidates for cross validation
  silverman_h = 1.06 * std(X) / num_data^(-1/5);
  cv_candidate_hs = logspace(-2,2,10)' * silverman_h;
  num_cands = size(candidate_hs, 1);

  % find the best h via cross validation
  best_err = inf;
  for cand_iter = 1:num_cands
    curr_h = cv_candidate_hs(cand_iter);
    curr_err = KFoldExperiment(Th, sc_obs_joint_probs, curr_h, NUM_KFCV_PARTS);
    if curr_err < best_err
      best_err = curr_err;
      best_h = curr_h;
    end
  end

  % Return the function with the optimal parameter
  post_est = estimatePostFromData(Th, sc_obs_joint_probs, best_h);

end


function error_accum = KFoldExperiment(X, y, h, num_partitions)
% This function splits up the data and performs cross validation to determine
% the optimal bandwidth
  error_accum = 0;
  m = size(X, 1);
  for kfold_iter = 1:num_partitions
    test_start_idx = round( (kfold_iter-1)*m/num_partitions + 1 );
    test_end_idx   = round( kfold_iter*m/num_partitions );
    train_indices = [1:test_start_idx-1, test_end_idx+1:m];
    test_indices = [test_start_idx : test_end_idx];
    Xtr = X(train_indices, :);
    ytr = y(train_indices);
    Xte = X(test_indices,:);
    yte = y(test_indices);

    % Obtain the estimate using Xtr and Ytr
    post_est = estimatePostFromData(Xtr, ytr, h);
    % compute the errors
    curr_error = norm(post_est(Xte) - yte)^2;
    error_accum = error_accum + curr_error;
  end
end


function post_est = estimatePostFromData(Th, sc_obs_joint_probs, h)
% This function does the actual computation of the function by solving a QP.
  num_pts = size(Th, 1);
  num_dims = size(Th, 2);
  K = GaussKernel(h, Th);
  alpha = quadprog(K*K, K* sc_obs_joint_probs, [], [], [], [], zeros(num_dims, 1));
  % rescale the alpha's
  alpha = alpha / sum(alpha);
  % Finally return the function
  post_est = @(T) GaussKernel(h, T, Th) * alpha;
end

