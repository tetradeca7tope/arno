function [mu, stddev, K] = GPRegression(X, y, Xtest, hyperParams, runtimeParams)
% Function for performing Gaussian Process Regression. Uses a Gaussian kernel.

  sigmaSm = hyperParams.sigmaSm;
  sigmaPr = hyperParams.sigmaPr;
  meanFunc = hyperParams.meanFunc;
  noise = hyperParams.noise;

  if ~exist('runtimeParams', 'var')
    plotOn = false;
  else
    plotOn = runtimeParams.plotOn;
  end

  if isempty(noise)
    noise = zeros(size(X,1), 1);
  end
  if isempty(meanFunc) % if meanFunc is empty use the mean of the y points
    meanFunc = @(arg) mean(y);
  end

  D11 = Dist2GP(X, X);
  K11 = sigmaPr * exp(-0.5*D11/sigmaSm^2);
  D22 = Dist2GP(Xtest, Xtest);
  K22 = sigmaPr * exp(-0.5*D22/sigmaSm^2);
  D12 = Dist2GP(X, Xtest);
  K12 = sigmaPr * exp(-0.5*D12/sigmaSm^2);

  % obtain outputs
  % concatenate these 2 matrices so that we invert the matrix only once.
  QQ = (K11 + diag(noise)) \ [K12 (y - meanFunc(X))];
  mu = meanFunc(Xtest) + K12' * QQ(:, end);
  K = K22 - K12' * QQ(:,1:end-1);
  stddev = real(sqrt(diag(K)));

  if (plotOn && (size(X,2) ==1) )
    plot(X, y, 'kx', 'MarkerSize', 10); hold on,
    plot(Xtest, mu);
    plot(Xtest, mu + 3*stddev, 'g--');
    plot(Xtest, mu - 3*stddev, 'g--');
    title('GP results: original points(kx), estimated vals(b),error(g--)');
  end

end
