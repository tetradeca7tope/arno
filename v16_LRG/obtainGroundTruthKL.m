% Obtain Ground Truth for the LRG Experiment

% This is what we'll do here
% First run an MC to collect 2*numSamplesEstEval.
% Use 1/2 of them to construct a KDE and 
% Evalute the KDE on the other half.
% Save the results.
% When we obtain an estimate we will be evaluating on points obtained in the
% second half.

% Perform MCMC to collect samples
fprintf('Collecting Samples\n');
totalNumSamplesEstEval = numBurninSampleEstEval + 2*numSamplesEstEval;
[samplesEstEval, queriesEstEval, logPsEstEval] = CustomMCMC( ...
  totalNumSamplesEstEval, evalMCMCProposalStd, evalMCMCInitPt,evalMCMCLogJoint);
% Remove burnin, and then shuffle
samplesEstEval = logitinv(samplesEstEval);
samplesEstEval = samplesEstEval( (numBurninSampleEstEval+1): end, :);
samplesEstEval = samplesEstEval( randperm(2*numSamplesEstEval), :);

% Split into two
samplesEstDensity = samplesEstEval(1:numSamplesEstEval, :);
klEvalPts = samplesEstEval( (numSamplesEstEval+1):end, :);
% print out some info
numUniqueSamples = size( unique(klEvalPts(:, 1)), 1);
fprintf('#unique-samples/#total-samples: %d/%d = %.3f', ...
  numUniqueSamples, numSamplesEstEval, numUniqueSamples/numSamplesEstEval);

% Now do a KDE
fprintf('Performing KDE\n');
[~, truePostEst, optKDEBandWidth] = kde01(samplesEstDensity);

% Evaluate the estimated density at the other half
fprintf('Evaluating estimate at select pts\n');
truePAtEvalPts = truePostEst(klEvalPts);

% Save results
save(gtFile, 'truePAtEvalPts', 'klEvalPts', 'optKDEBandWidth');
save(allGtFile, 'truePAtEvalPts', 'klEvalPts', 'optKDEBandWidth',  ...
  'numBurninSampleEstEval', ...
  'numSamplesEstEval', 'totalNumSamplesEstEval', 'samplesEstEval', ...
  'queriesEstEval', 'logPsEstEval', 'samplesEstDensity');

