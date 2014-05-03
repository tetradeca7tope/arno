function [logJointProbs, lumMeans] = evalSNLSLogLikl(evalPts, obsData)
% evalPts in an numPts x 3 matrix. The first column is the Hubble constant, the
% second is OmegaM(matter fraction) and the third is
% OmegaL(dark energy fraction)
% Output: 
% logJointProbs is the log likelihood of the observations
% lumMeans is the expected luminosity. It is a numPts x numObs matrix.

  % Define some constants
  LIGHT_SPEED = 299792.458; % in Km/s
  RESOLUTION = 100; % Resolution for 1D integration

  % Prelims
  numPts = size(evalPts, 1);
  numObs = size(obsData, 1);
  % Decompose obsData
  z = obsData(:, 1);
  obs = obsData(:,2);
  obsErr = obsData(:,3);
  
  % Create arrays for storing the outputs
  logJointProbs = zeros(numPts, 1);
  lumMeans = zeros(numPts, numObs);

  % Now iterate through each of the evalPts.
  % I could vectorize this, but that would make the code very less readable.
  % Also, most of the time we will only be querying at just one point.
  for iter = 1:numPts
    
    % Decompose the current observation
    H = evalPts(iter,1);
    OmegaM = evalPts(iter,2);
    OmegaL = evalPts(iter,3);
    currLumMeans = zeros(numObs, 1); % for storing the computed means
    f = @(t) 1 ./ sqrt(OmegaM*(1+t).^3 + OmegaL); % function handle

    % compute the expected observation for each observation
    for obsIter = 1:numObs
      dL = LIGHT_SPEED * (1 + z(obsIter)) / H * integral(f, 0, z(obsIter));
      currLumMeans(obsIter) = 5 * log10(dL) + 25;
    end
    % Now compute the log likelihood for each observation
    obsLikl = normpdf(obs, currLumMeans, obsErr);
    logJointProbs(iter) = sum(log(obsLikl));
    lumMeans(iter, :) = currLumMeans';
  end

end
