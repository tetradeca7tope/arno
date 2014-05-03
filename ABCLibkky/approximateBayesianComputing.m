function [acceptedSamples, metricVals] = approximateBayesianComputing( ...
  priorSamples, likelihoodSampler, metric, observation, acceptTol)
% priorSamples is a bunch of parameters sampled from the prior.
% likelihoodSampler is function handle that returns a sample from the likelihood
% given a parameter value and metric is a function handle which takes in
% a bunch of samples and returns the value of the metric to the observation.
% If acceptCriterion is empty then acceptParams should be a struct which has the
% following parameters: observation, metric, acceptTolerance
% Note that all these function handles should be vectorized

  % Simulated observations from the likelihood
  simObs = likelihoodSampler(priorSamples);

  % compute the distance from the simulations to the observation
  metricVals = metric(simObs, observation);

  % Compute which points are accepted
  if ~isempty(acceptTol)
    acceptedIdxs = metricVals < acceptTol;
    acceptedSamples = priorSamples(acceptedIdxs, :);
  else
    acceptedSamples = [];
  end

end

