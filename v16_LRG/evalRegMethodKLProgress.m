function [kl, logJointEst, probEst] = evalRegMethodKLProgress( Xtr, Ytr, ...
  gpFitParams, klEvalPts, truePAtEvalPts, evalMCMCParams, optKDEBandWidth)
% This is essentially a utility function for VR, MBP, MCMC-REG and RAND.
% It does the following in order
% - First construct an estimate of the log Joint probability using Xtr and Ytr
% - Then collect samples from this estimate via MCMC
% - Perform a KDE on these samples to obtain an estimate of the pdf
% - Evaluate this pdf estimate at klEvalPts
% - Finally use them to estimate the kl between the estimate and the truth

% Inputs:
% Xtr, Ytr: Regressors and regressands for the estimate
% gpFitParams: parameters for fitting the GP
% klEvalPts: Samples from the ground truth
% truePAtEvalPts: Value of the true Density at klEvalPts;

% Returns:
% kl: the estimated KL between the truth and the estimate
% logJointEst: A function handle for the log JOint estimate via GP Regression.
% probEst: A function handle for the pdf estimate obtained via KDE on the
%          samples collected from logJointEst

  % Construct estimate of the log Joint probability 
  logJointEst = regressionWrap(Xtr, Ytr, gpFitParams.noiseLevelGP, ...
    gpFitParams.lowestLogliklVal, gpFitParams.logLiklRange, ...
    gpFitParams.cvCostFunc);

  % Collect Samples via MCMC
  mcmcLogJointEst = @(t) logJointEst(logitinv(t));
  totalNumSamples = evalMCMCParams.numBurninSampleEstEval + ...
                    evalMCMCParams.numSamplesEstEval;
  [logitMcmcSamples] = CustomMCMC( totalNumSamples, ...
    evalMCMCParams.evalMCMCProposalStd, evalMCMCParams.evalMCMCInitPt, ...
    mcmcLogJointEst);
  logitMcmcSamples = ...
    logitMcmcSamples( (evalMCMCParams.numBurninSampleEstEval+1): end, :);
  mcmcSamples = logitinv( logitMcmcSamples );

  % Do the KDE -- use the optimal Bandwidth
  [~, probEst] = kde01(mcmcSamples, optKDEBandWidth);

  % Finally obtain the KL
  kl = estimKLForLRG(klEvalPts, truePAtEvalPts, probEst);
  
end  
