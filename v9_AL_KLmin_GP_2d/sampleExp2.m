function [theta, garb1, True_Post, marginal_prob, garb2] = sampleExp2(Th1, Th2)

  N11 = size(Th1, 1);
  N12 = size(Th1, 2);

  th = [Th1(:), Th2(:)];
  log_joint_probs = evalLogJointProbsExp2 (Th1(:), Th2(:), 0);
  Log_Joint_Probs = reshape(log_joint_probs, N11, N12);

  joint_probs = exp(log_joint_probs);
  marginal_prob = numerical_2D_integration_wrap(th, joint_probs);
  true_post = joint_probs / marginal_prob;
  True_Post = reshape(true_post, N11, N12);

  % garbage
  theta = zeros(2,1);
  garb1 = 0;
  garb2 = 0;

end

