function ret = probRatioStatistic(pts, trueLogPAtPts, logPEst, ...
  refPt, refTrueLogP)
% Compute the ratio of the probabilities w.r.t refProb. Ideally, let refPt be
% the MAP. 

  % First compute the true stuff
  Pi = exp( trueLogPAtPts - refTrueLogP );

  % Now the est stuff
  refEstLogQ = logPEst(refPt);
  logQi = logPEst(pts) - refEstLogQ;
  Qi = exp(logQi);

  ret = mean( abs(Pi-Qi) );

end

