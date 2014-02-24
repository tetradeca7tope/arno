% This script tests Approximate Bayesian Computing for estimating the posterior

NUM_ABC_PTS = S_RATIO*NUM_AL_ITERS + NUM_INITIAL_PTS;
ABC_EPS_FRACTION = 0.02;

curr_accepted_points = zeros(0, 1);
curr_kl = 1000; %Initialize to a v.large value. Will keep track of the last KL
abc_kl_progress = curr_kl * ones(NUM_ABC_PTS, 1);

% fig_abc = figure;

for abc_iter = 1:NUM_ABC_PTS

  % Sample theta from the prior
%   abc_theta_ = dirichlet_sample([ALP, BEP]);   
%   abc_theta = abc_theta_(1);
% 
%   abc_p = abc_theta^2 + (CCP - abc_theta)^2;
%   abc_X = double(rand(N,1) < abc_p);
%   abc_sumX = sum(abc_X);

  [abc_theta, abc_sumX] = sampleModifiedBeta(priorParams, th, ...
    NUM_BINOMIAL_SAMPLES, true);

%   if ~dontcheat
%     % cheating here. preferring multimodal distros over balanced ones.
%     theta_ = dirichlet_sample([30, 5]); theta = theta_(1);
%   else
%     theta_ = dirichlet_sample([alp, bep]); theta = theta_(1);
%   end
% %   theta = 0.15; % cheating here. trying to force a multimodal distro.
%   p = theta^2 + (c_param - theta)^2;
%   X = double(rand(N,1) < p);
%   sumX = sum(X);

  % First set it to the default value. If you have a good estimate it will be
  % overwritten
  abc_kl_progress(abc_iter) = curr_kl;
  % Now compare it to sumX. If comparable accept abc_theta. Otherwise discard it
  if (abs(sumX - abc_sumX)/sumX < ABC_EPS_FRACTION)
    curr_accepted_points = [curr_accepted_points; abc_theta];
    
    % And now use a KDE to estimate the density
    if size(curr_accepted_points, 1) > 1
      [~, abc_kde] = kde01(curr_accepted_points);
      [abc_est_posterior] = abc_kde(th);
      ptwise_kl = true_post .* (log(true_post) - log(abc_est_posterior));
      curr_kl = numerical_1D_integration(th, ptwise_kl);
      abc_kl_progress(abc_iter) = curr_kl;
      if curr_kl == 0
        pause;
      end
      
      PLOT_OK_LOCAL = false;
%       PLOT_OK_LOCAL = true;
      if PLOT_OK_LOCAL
        fprintf('iter: %d, #samples collected: %d, KL: %f\n', abc_iter, ...
          num_samples_collected, curr_kl);
        figure(fig_abc);
        plot(th, true_post, 'g'); hold on,
        plot(th, abc_est_posterior, 'r'); hold off,
      end
    end
  end
  num_samples_collected = size(curr_accepted_points, 1);
  
end

% Finally get rid of the first NUM_INITIAL_PTS points
abc_kl_progress = abc_kl_progress(NUM_INITIAL_PTS+1 : end);
% 
% figure(fig_klprogress);
% plot(log(abc_kl_progress), 'c*-');
% title('KL divergence after each AL iteration');

