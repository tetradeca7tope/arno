% getPYMCResults.m

% First set the SNLS Experiment up
% ==============================================================================
% These are the points at which we'll be evaluating the KL for comparison
NUM_EXPERIMENTS = 30;
numResultsToBeStored = 43; %round(10000/160);

% NUM_EXPERIMENTS = 1;
% numResultsToBeStored = 5; %round(10000/160);

% Problem Constants
snlsParamSpaceBounds = [60 80; 0 1; 0 1];
paramSpaceBounds = [0 1; 0 1; 0 1]; % we are normalizing them later on
partitionRes = 100; % For integrating and estimating the KL etc.
[g1, g2, g3] = gen3DCentrePointGrid(partitionRes, paramSpaceBounds);
evalPts = [g1(:), g2(:), g3(:)];
% Obtain the truth
truePostFile = sprintf('truePost_nObs%d_res%d_H0%d-%d.txt', ...
  [], snlsParamSpaceBounds(1,1), snlsParamSpaceBounds(1,2),  ...
  partitionRes);
trueLogPostValsAtEvalPts = load(truePostFile);

% For the results
pymcResults = zeros(NUM_EXPERIMENTS, numResultsToBeStored);
saveResultsFile = 'pymcResults.mat';

for expIter = 1:NUM_EXPERIMENTS
  fprintf('Experiment %d\n', expIter);

  for resIter = 1:numResultsToBeStored

    loadFile = sprintf('emcee/emceeRes/exp%d_stop%d.txt', (expIter-1), ...
      (resIter-1) );
    currMCMCSamples = load(loadFile); 
    currNumMcmcPts = size(currMCMCSamples, 1);
%     % Obtain the current set of points

    % Perform KDE
    [~, mcmcProbEst] = kde01(currMCMCSamples);
    mcmcLogProbEstAtEvalPts = log(mcmcProbEst(evalPts));

    % Evaluate the KL
    currKL = eval3DKLFromLogProbs(trueLogPostValsAtEvalPts, ...
      mcmcLogProbEstAtEvalPts, paramSpaceBounds);
    pymcResults(expIter, resIter) = currKL;
    fprintf(' MCMC Iter: %d, num-points: %d, KL: %0.4f\n', resIter, ...
      currNumMcmcPts, currKL);
  end

  save(saveResultsFile, 'pymcResults');
end

