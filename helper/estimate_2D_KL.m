function [KL] = estimate_2D_KL(th, logP1, logP2)
% Estimates the KL between 2 distributions from their pointwise estimates P1
% and P2. 

  N = sqrt(size(th,1));

  ptwise_KL = exp(logP1) .* (logP1 - logP2);
  ptwise_KL = reshape(ptwise_KL, N, N);
  Th1 = reshape(th(:,1), N, N);
  Th2 = reshape(th(:,2), N, N);
  KL = numerical_2D_integration(ptwise_KL, Th1, Th2);

end
