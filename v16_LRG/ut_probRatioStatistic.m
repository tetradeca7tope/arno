% unit test for probRatioStatistic

% First a simple test
X = randn(1000,1);
logPX = log(normpdf(X));
logQEst = @(t) log(normpdf(t));
refPt = 0;
truePAtRefPt = log( normpdf(refPt) );
prs = probRatioStatistic(X, logPX, logQEst, refPt, truePAtRefPt);
fprintf('Computed = %f, True = 0\n\n', prs);

% Do the same but add some noise to logQEst
logQEst = @(t) log(normpdf(t)) + 0.01*randn(size(t,1), 1);
prs = probRatioStatistic(X, logPX, logQEst, refPt, truePAtRefPt);
fprintf('Computed = %f\n\n', prs);


% Now do something more sophisticated
numDims = 1;
numPts = 1000;
mu = 2;
logP = @(t) log( mvnpdf(t, mu * ones(1, numDims), eye(numDims) ) - 20 );
logPAtMu = logP(mu*ones(1, numDims) );
X = mu + randn(numPts, numDims);
  % Obtain a kde
  [~, q1] = kde(X);
  logQ1 = @(t) log( q1(t) );
% X2 = mu + randn(numPts, numDims);
X2 = X;
  [~, q2] = kde(X2 - mu);
  logQ2 = @(t) log( q2(t) );

numRes = 10000;
pts = 12*rand(numRes, numDims) - 6;
logPAtPts = logP( pts);
prs1 = probRatioStatistic(pts, logPAtPts, logQ1, mu*ones(1, numDims), logPAtMu);
prs2 = probRatioStatistic(pts, logPAtPts, logQ2, mu*ones(1, numDims), logPAtMu);
fprintf('prs1 = %0.4f\nprs2 = %0.4f\n', prs1, prs2);

