function C = construct3DConfSetFromFunc(probDens, res, bounds, alpha)
% constucts a 1-alpha level confidence set by partitioning the space at
% resolution res and allow 

  % Evaluate the densities at the centre points 
  [g1, g2, g3] = get3DCentrePoingGrid(res, bounds);
  centrePts = [g1(:), g2(:), g3(:)];
  probDensVals = probDens(centrePts);

  % Now obtain the probabilities
  vol = (bounds(1,2) - bounds(1,1)) * (bounds(2,2) - bounds(2,1)) * ...
        (bounds(3,2) - bounds(3,1))/ res^3;
  probVals = probDensVals * vol;
  probVals = probVals / sum(probVals);

  % Now obtain the top 95%
  [sortedProbVals, sortedIdxs] = sort(probVals, 'descend');
  sortedCumProbs = cumsum(sortedProbVals);

  % Now tick the top points
  selPoints = sortedCumProbs <= 1 - alpha;
  unrolledC = sparse(selPoints,  )
  

end
