function [kl, l2] = estimKLForLRG( klEvalPts, truePAtEvalPts, q)
% klEvalPts are samples from p, truePAtEvalPts is p(X) at klEvalPts
% q is a function handle for the pdf q. The function estimates KL(p,q).

  % Evaluate q
  qAtEvalPts = q(klEvalPts);

  LOW_TOL = 1e-250;

  % Compute the KL
  logPs = log(truePAtEvalPts);
  logQs = log(qAtEvalPts);
  klIdxs = (truePAtEvalPts > LOW_TOL) & (qAtEvalPts > LOW_TOL);
  kl = mean( logPs(klIdxs) - logQs(klIdxs) );
  if isnan(kl)
    fprintf('Warning: KL was nan\n');
    kl = inf;
  end

  % Compute the L2
  l2 = mean( (truePAtEvalPts - qAtEvalPts).^2 ./ truePAtEvalPts );
  
end

