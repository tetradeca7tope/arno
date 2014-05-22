% Collect samples from true likelihood and estimate for comparison
% Assume that you have already run UC once

numSamples = 25e6;
mcmcProposalStd = 5;
mcmcInitPt = 0.5*ones(numDims, 1); % init at centre point

% First collect true samples
fprintf('True LogL\n');
mcmcEvalLogJoint = @(t) evalLogJoint(logitinv(t));
[logitMCMCSamps, logitMCMCQuers, MCMCLogPs] = CustomMCMC(numSamples, ...
  mcmcProposalStd, mcmcInitPt, mcmcEvalLogJoint);
mcmcSamps = logitinv(logitMCMCSamps);
mcmcQuers = logitinv(logitMCMCQuers);

% Now collect samples from the estimate
fprintf('Est LogL\n');
mcmcEstLogJoint = @(t) ucLogJointEst(logitinv(t));
[logitMCMCEstSamps, logitMCMCEstQuers, MCMCEstLogPs] = CustomMCMC( ...
  numSamples, mcmcProposalStd, mcmcInitPt, mcmcEstLogJoint);
mcmcEstSamps = logitinv(logitMCMCEstSamps);
mcmcEstQuers = logitinv(logitMCMCEstQuers);

filename = sprintf('allCollSamples-%d-%s.mat', ...
  numSamples, datestr(now, 'mm-dd-HH-MM') );
save('allCollectedSamples.mat', ...
  'logitMCMCSamps', 'logitMCMCQuers', 'MCMCLogPs', 'mcmcSamps', 'mcmcQuers', ...
  'logitMCMCEstSamps', 'logitMCMCEstQuers', 'MCMCEstLogPs', ...
  'mcmcEstSamps','mcmcEstQuers');

