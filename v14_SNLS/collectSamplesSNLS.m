% Collect samples from true likelihood and estimate for comparison
% Assume that you have already run UC once

numSamples = 1e7 + 1e5;

% USE_LOGIT = true;
USE_LOGIT = false;

if USE_LOGIT
  mcmcProposalStd = 5;
  mcmcInitPt = logit( [.5 0.24 0.68]'); % init at ML point
else
  mcmcProposalStd = 0.1;
  mcmcInitPt =  [.5 0.24 0.68]'; % init at ML point
end

% First collect true samples
% fprintf('True LogL\n');
% tic,
% mcmcEvalLogJoint = @(t) evalLogJoint(logitinv(t));
% [logitMCMCSamps, logitMCMCQuers, MCMCLogPs] = CustomMCMC(numSamples, ...
%   mcmcProposalStd, mcmcInitPt, mcmcEvalLogJoint);
% toc,
% mcmcSamps = logitinv(logitMCMCSamps);
% mcmcQuers = logitinv(logitMCMCQuers);

% Now collect samples from the estimate
fprintf('Est LogL\n');
tic,
if USE_LOGIT
  mcmcEstLogJoint = @(t) ucLogJointEst(logitinv(t));
  [logitMCMCEstSamps, logitMCMCEstQuers, MCMCEstLogPs] = CustomMCMC( ...
    numSamples, mcmcProposalStd, mcmcInitPt, mcmcEstLogJoint);
  mcmcEstSamps = logitinv(logitMCMCEstSamps);
  mcmcEstQuers = logitinv(logitMCMCEstQuers);
else
  mcmcEstLogJoint = @(t) ucLogJointEst(t);
  [mcmcEstSamps, mcmcEstQuers, MCMCEstLogPs] = CustomMCMC( ...
    numSamples, mcmcProposalStd, mcmcInitPt, mcmcEstLogJoint);
end
toc,
mcmcEstSamps = mcmcEstSamps(100001:end, :);


filename = sprintf('allCollSamples-%d-%s.mat', ...
  numSamples, datestr(now, 'mm-dd-HH-MM') );
save(filename, 'mcmcEstSamps');
% save(filename, ...
%   'logitMCMCSamps', 'logitMCMCQuers', 'MCMCLogPs', 'mcmcSamps', 'mcmcQuers', ...
%   'logitMCMCEstSamps', 'logitMCMCEstQuers', 'MCMCEstLogPs', ...
%   'mcmcEstSamps','mcmcEstQuers');

