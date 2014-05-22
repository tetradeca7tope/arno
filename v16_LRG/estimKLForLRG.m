function kl = estimKLForLRG( klEvalPts, truePAtEvalPts, q)
% klEvalPts are samples from p, truePAtEvalPts is p(X) at klEvalPts
% q is a function handle for the pdf q. The function estimates KL(p,q).

  LOW_TOL = 1e-250;
  validIdxs = truePAtEvalPts > LOW_TOL;

  kl = mean( log( truePAtEvalPts(validIdxs)) - ...
             log( q(klEvalPts(validIdxs, :)) ) );
end

