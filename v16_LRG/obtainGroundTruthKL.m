% Obtain Ground Truth for the LRG Experiment

% This is what we'll do here
% First run an MC to collect 2*numSamplesEstEval.
% Use 1/2 of them to construct a KDE and 
% Evalute the KDE on the other half.
% Save the results.
% When we obtain an estimate we will be evaluating on points obtained in the
% second half.

% Don't need to do a KDE if you are just doing mmd
DOING_MMD = true;
% DOING_MMD = false;

% Perform MCMC to collect samples
fprintf('Collecting %d Samples with bw: %0.4f\n', numSamplesEstEval, ...
  evalMCMCProposalStd);

if ~DOING_MMD
  totalNumSamplesEstEval = numBurninSampleEstEval + 2*numSamplesEstEval;
else 
  totalNumSamplesEstEval = numBurninSampleEstEval + numSamplesEstEval;
end
  
% Collect the samples
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[samplesEstEval, queriesEstEval, logPsEstEval] = CustomMCMC( ...
  totalNumSamplesEstEval, evalMCMCProposalStd, evalMCMCInitPt, ...
  evalMCMCLogJoint);
max(max(samplesEstEval)), min(min(samplesEstEval)),
if DOING_LOGIT
  samplesEstEval = logitinv(samplesEstEval);
end
% Remove burnin, and then shuffle
samplesEstEval = samplesEstEval( (numBurninSampleEstEval+1): end, :);
  % print out some info
  numUniqueSamples = size( unique(samplesEstEval(:, 1)), 1);
  fprintf('#unique-samples/#total-samples: %d/%d = %.3f', ...
    numUniqueSamples, numSamplesEstEval, numUniqueSamples/numSamplesEstEval);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~DOING_MMD

  % Shuffle & Split into two
  samplesEstEval = samplesEstEval( randperm(2*numSamplesEstEval), :);
  samplesEstDensity = samplesEstEval(1:numSamplesEstEval, :);
  klEvalPts = samplesEstEval( (numSamplesEstEval+1):end, :);

  % Now do a KDE
  fprintf('Performing KDE\n');
  [~, truePostEst, optKDEBandWidth] = kde01(samplesEstDensity);

  % Evaluate the estimated density at the other half
  fprintf('Evaluating estimate at select pts\n');
  truePAtEvalPts = truePostEst(klEvalPts);

else

  klEvalPts = samplesEstEval;
  truePostEst = [];
  optKDEBandWidth = [];
  truePAtEvalPts = [];
  samplesEstDensity = [];

end

% Finally for probRatioStatistic
prsPts = rand(numALCandidates, numDims);
prsLogProbs = evalLogJoint(prsPts);

% Save results
save(gtFile, 'truePAtEvalPts', 'klEvalPts', 'optKDEBandWidth', ...
  'queriesEstEval', 'logPsEstEval', 'prsPts', 'prsLogProbs');
save(allGtFile, 'truePAtEvalPts', 'klEvalPts', 'optKDEBandWidth',  ...
  'numBurninSampleEstEval', ...
  'numSamplesEstEval', 'totalNumSamplesEstEval', 'samplesEstEval', ...
  'queriesEstEval', 'logPsEstEval', 'samplesEstDensity');

