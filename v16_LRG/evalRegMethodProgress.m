function [err, logJointEst] = evalRegMethodKLProgress( Xtr, Ytr, ...
  gtPts, gtLogProbs, gpFitParams)
% This is essentially a utility function for VR, MBP, MCMC-REG and RAND.

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

  errFunc = @(a, b) sqrt(mean( (exp(a) - exp(b)).^2 ));

  % Construct estimate of the log Joint probability 
  numDims = size(Xtr, 2);
  logJointEst = regressionWrapLRG2(Xtr, Ytr, gpFitParams.noiseLevelGP, ...
    gpFitParams.lowestLogliklVal, gpFitParams.logLiklRange, ...
    gpFitParams.cvCostFunc );

  estLogProbs = logJointEst(gtPts);
  err = errFunc(gtLogProbs, estLogProbs);
  
end

