% Script to plot results
for iter = 1:NUM_AL_ITERS
  f1 = figure;
  % We will plot the following here: 1: the actual log-joint prob, 2- the
  % estimated log-joint probs, 3- the brute force estimate, 4- the KL
  % curve for AL estimate and 5- the brute force estimate
  % 6 - the samples
  load_file_al = sprintf('%s/al_iter_%d.mat', DEBUG_DIR_NAME, iter);
  load_file_bf = sprintf('%s/bf_iter_%d.mat', DEBUG_DIR_NAME, iter);
  load(load_file_al);
  load(load_file_bf);

  for i = 1:NUM_SAMPLE_FUNCTIONS
    plot(th, curr_gp_samples(:,i), 'm'); hold on,
  end
  plot(th, gr_est_joint_log_probs, 'g', 'Linewidth', 2); % The BF estimate
  plot(th, est_joint_log_probs, 'b-.', 'Linewidth', 2); % The GP estimate
  plot(th, log_joint_probs, 'k--', 'Linewidth', 2); % True posterior

    save(save_file_name, 'gr_est_joint_log_probs', 'gr_est_post', ...
                         'gr_ptwise_KL');
  title('Comparison of Log Joint: True(k), Est-AL(b), Est-BF(g)');
  
  f2 = figure;
  plot(th, ptwise_KL, 'b-.'); hold on,
  plot(th, gr_ptwise_KL, 'g', 'Linewidth', 2);
  title('Comparison of KL: Est-AL(b), Est-BF(g)');

  f3 = figure;
  plot(th, gr_est_post, 'g'); hold on,
  plot(th, est_post, 'b-.');
  plot(th, true_post, 'k--');

  % Save files
  sf1 = sprintf('%s/iter%d_logjoint', DEBUG_DIR_NAME, iter);
  saveas(f1, sf1, 'png');
  sf2 = sprintf('%s/iter%d_kl', DEBUG_DIR_NAME, iter);
  saveas(f2, sf2, 'png');
  sf3 = sprintf('%s/iter%d_post', DEBUG_DIR_NAME, iter);
  saveas(f3, sf3, 'png');
  close(f1);
  close(f2);
  close(f3);
end

