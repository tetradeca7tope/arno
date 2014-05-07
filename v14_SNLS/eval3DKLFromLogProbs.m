function kl = eval3DKL(logP, logQ, bounds)

  ptwiseKL = exp(logP) .* (logP - logQ);
  vol = prod( bounds(:,2) - bounds(:,1) );
  kl = mean(ptwiseKL) * vol;

end

