function [w_opt, h_opt, Ztr] = ProbEst(Z, P)
% Given a set of points and the corresponding probability values (up to a
% constant) this function returns returns an estimate of the density of the
% form f(z) = SUM(wi * k(z,zi)) where the kernel depends on h.

  % MACROS
%   DEBUG_MODE = true;
  DEBUG_MODE = false;

  % SET the FOLLOWING
  OBJ_NORM = 2;

  % Sort tests into train and cross validation
  m = size(Z,1); % num data points
  mtr = floor((2/3)*m);
  mcv = m - mtr;
  permute_data = randperm(m);
  Ztr = Z(permute_data(1:mtr), :);
  Zcv = Z(permute_data(mtr+1:end), :);
  Ptr = P(permute_data(1:mtr), :);
  Pcv = P(permute_data(mtr+1:end), :);

  Num_logspace_h_vals = 5;
  Num_linspace_h_vals = 10;
  Num_candidate_h_vals = Num_linspace_h_vals * Num_logspace_h_vals;
  candidate_h_vals = zeros(Num_logspace_h_vals * Num_linspace_h_vals, 1);

  logspace_break_pts = logspace(-2, 2, Num_logspace_h_vals+1);
  for i = 1:Num_logspace_h_vals
    candidate_h_vals( (i-1)*Num_linspace_h_vals + 1 : ...
                      i*Num_linspace_h_vals ) ...
      = linspace(logspace_break_pts(i), logspace_break_pts(i+1), ...
                 Num_linspace_h_vals);
  end

  % set the following
  data_dim = size(Z, 2),
  best_optval = inf;
  optvals = zeros(Num_candidate_h_vals, 1);

  for cand_iter = 1:Num_candidate_h_vals
    
    h = candidate_h_vals(cand_iter);
    Dtr = dist2(Ztr, Ztr);
    Ktr = exp( -Dtr/ (2*h^2) ) / (2*pi*h^2)^(data_dim/2);

    cvx_begin quiet
    variable w(mtr);
%     minimize ( norm( (K*w) ./P  - ones(m,1), OBJ_NORM));
    minimize ( norm( (Ktr*w) - Ptr, OBJ_NORM));
    subject to
      w >= 0;
    cvx_end
    C = sum(w);
    w = w/C;
    Pcv_ = Pcv/C;

    cv_probs = evaluate_prob(Zcv, Ztr, w, h);
    curr_cv_optval = norm( Pcv - cv_probs*C, OBJ_NORM ); % / h^data_dim;
    optvals(cand_iter) = curr_cv_optval;
    if curr_cv_optval < best_optval
      best_optval = curr_cv_optval;
      w_opt = w;
      h_opt = h; 
    end

    fprintf('iteration: %d, h = %f, C = %f, optval = %f\n', ...
            cand_iter, h, C, log(curr_cv_optval));

    if DEBUG_MODE
    % PLOT these out
      th = linspace(-20, 20, 2000)';
      p_th = evaluate_prob(th, Ztr, w, h);
      plot(th, p_th); hold on,
      plot(Ztr, Ptr/C, 'mx'); 
      plot(Zcv, Pcv/C, 'kx'); hold off,
      title_str = sprintf('h = %f, optval: %f, best: %f\ninteg: %.12f', ...
        h, log(curr_cv_optval), log(best_optval), ...
        numerical_1D_integration(th, p_th));
      title(title_str);
      pause;
    end
  end

  % finally perform the optimziation using h_opt over all values.
  D = dist2(Z, Z);
  K = exp(-D/ (2*h_opt^2) ) / (2*pi*h_opt^2)^(data_dim/2);
  cvx_begin quiet
  variable w(m);
  minimize ( norm( (K*w) - P, OBJ_NORM));
  subject to
    w >= 0;
  cvx_end
  C = sum(w);
  w = w/C;
  w_opt = w;

  loglog(candidate_h_vals, optvals);
  title('candidate bandwidths vs optimal values');
  fprintf('Experimented with the following h values:');
  candidate_h_vals,
  w_opt, h_opt,

end
