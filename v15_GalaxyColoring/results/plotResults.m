
addpath ../../helper/
addpath ../../GPLibkky/

load exp01mbp.mat
load exp01uc.mat

numDims = 2;
paramSpaceBounds = repmat([0 1], numDims, 1);

mbpLogJointEst = regressionWrap(mbpPts, mbpLogProbs, noiseLevelGP, ...
  lowestLogliklVal, logLiklRange);

ucLogJointEst = regressionWrap(mbpPts, mbpLogProbs, noiseLevelGP, ...
  lowestLogliklVal, logLiklRange);

figure;
plot2DFunction(mbpLogJointEst, paramSpaceBounds, 'mesh'); hold on
plot3(mbpPts(:,1), mbpPts(:,2), mbpLogProbs, 'kx');
title('MBP');

figure;
plot2DFunction(ucLogJointEst, paramSpaceBounds, 'mesh');
plot3(ucPts(:,1), ucPts(:,2), ucLogProbs, 'kx');
title('UC');

