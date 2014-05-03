% Unit test for approximate Bayesian Comutation

% its a simple set up. 2 parameter distribution with a bernoulli random variable
% drawn from (theta1 + theta2^2)/2.
t1 = 0.5;
t2 = 0.7;
p = (t1 + t2.^2)/2;

% Generate observations
N = 1000;
X = rand(N, 1) > p;

% Call ABC
likelihoodSampler = @(arg) bsxfun(@ge, rand(size(arg,1), N), ...
  (arg(:,1) + arg(:,2).^2)/2 );
metric = @(arg, obs) abs(sum(arg,2) - sum(obs) )/ N;
paramSamples = rand(10000, 2);
[~, mVals] = approximateBayesianComputing(paramSamples, ...
  likelihoodSampler, metric, X, []);

% Now create some plots
diffThresh = linspace(min(mVals), max(mVals), 7);
diffThresh = diffThresh(2:end);
for i = 1:6
  subplot(2, 3, 3*floor( (i-1)/3 ) + mod(i-1,3) + 1);
  accParams = mVals < diffThresh(i);
  plot(paramSamples(accParams, 1), paramSamples(accParams, 2), 'rx'); 
  titlestr = sprintf('Tol > %f\n', diffThresh(i) );
  title(titlestr);
end

