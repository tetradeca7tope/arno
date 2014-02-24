function [est_probs, f, h] = kde(X)
% Performs Kernel density estimation on the data. Picks the bandiwidth using
% 10-fold cross validation. Assumes domain [-inf, inf]^d
% If domain is [0, 1]^d use kde01
% input - X: data in an nxd matrix. 
% output - est_probs: estimated probabilities, f: A function handle to estimate
%   the density on new data. h: optimal bandwidth

  n = size(X, 1);
  d = size(X, 2);
  num_partitions_kfoldcv = min(10, n);

  % If only 1 point, return a gaussian centered at that point. The following
  % procedure is equivalent to that.
  if n == 1
    best_h = 1000;
  end
  % shuffle the data
  X = X(randperm(n), :);

  silverman_h = 1.06 * std(X) / n^(-1/5);
  candidate_hs = logspace(-2,2,10)' * silverman_h;
  num_cands = size(candidate_hs, 1);

  best_likl = -inf;

  for cand_iter = 1:num_cands
    curr_h = candidate_hs(cand_iter);
    cv_loglikl = KFoldExperiment(X, curr_h, num_partitions_kfoldcv);
    if cv_loglikl > best_likl
      best_likl = cv_loglikl;
      best_h = curr_h;
    end
  end

  % finally prep results for returning
  h = best_h;
  f = @(data) ( sum(GaussKernel(h, X, data))'/ n); 
  est_probs = f(X);

end


function [cv_loglikl] = KFoldExperiment(X, h, num_partitions_kfoldcv)

  loglikl_accum = 0;
  n = size(X, 1);

  for kfold_iter = 1:num_partitions_kfoldcv
    % Set the partition up.
    test_start_idx = round( (kfold_iter-1)*n/num_partitions_kfoldcv + 1 );
    test_end_idx   = round( kfold_iter*n/num_partitions_kfoldcv );
    train_indices = [1:test_start_idx-1, test_end_idx+1:n]';
    test_indices = [test_start_idx : test_end_idx]';
    num_test_data = test_end_idx - test_start_idx + 1;
    num_train_data = n - num_test_data;
    Xtr = X(train_indices, :);
    Xte = X(test_indices, :);
    % Compute the log likelihood
    Pte = sum(GaussKernel(h, Xtr, Xte))' / n;
    avg_log_likl = sum(log(Pte)) / num_test_data;
    loglikl_accum = loglikl_accum + avg_log_likl;
  end
  cv_loglikl = loglikl_accum / num_partitions_kfoldcv;
end


