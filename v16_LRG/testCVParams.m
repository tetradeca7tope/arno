
numDims = 12;
ucFile = 'gt10_v1_1.mat';

randTestFiles = {'gt10_v1_0.5.mat', 'gt10_v1_0.75.mat', ...
  'gt10_v1_1.5.mat', 'gt10_v1_2.5.mat', 'groundTruth-10v1'};

numDims = 8;
lowestLogliklVal = -200;
logLiklRange = 200;

% Set up problem class
paramSpaceBounds = repmat([0 1], numDims, 1);
lrgExp = LRGExperiment(lowestLogliklVal);
evalLogJoint = @(arg) lrgExp.normCoordLogJointProbs(arg);
MLEPoint = lrgExp.getNormCoords([0 0.762 0.1045 0.02233 0.951 0.6845 0 1.908]);
MLELogP = evalLogJoint(MLEPoint);
% For fitting the GP
noiseLevelGP = logLiklRange/100;
cvCostFunc = @(y1, y2) (exp(y1) - exp(y2)).^2;

gt = load('gt10_v1_1.mat');
gtPts = gt.ucPts;
gtLogProbs = gt.ucLogProbs;
% gtPts = zeros(0, numDims);
% gtLogProbs = zeros(0, 1);
for i = 1:numel(randTestFiles)
  gt = load(randTestFiles{i});
  gtPts = [gtPts; gt.randPts]; gtLogProbs = [gtLogProbs; gt.randLogProbs];
end

size(gtPts),
load(ucFile),

% numPtsToUse = [10:10:100, 1000:1000:5000];
% numPtsToUse = [100:100:2000, 2400:400:4000]; %, 1000:1000:5000];
numPtsToUse = [3000]; %, 1000:1000:5000];
ucErrs = zeros(numel(numPtsToUse), 1);
mrErrs = zeros(numel(numPtsToUse), 1);
randErrs = zeros(numel(numPtsToUse), 1);

% for i = 1:0
for i = 1:numel(numPtsToUse)
  q = numPtsToUse(i);
%   uclj = regressionWrapLRG(ucPts(1:q, :), ucLogProbs(1:q, :), noiseLevelGP, lowestLogliklVal, ...
%     logLiklRange, cvCostFunc, paramSpaceBounds);
  uclj = regressionWrap(ucPts(1:q, :), ucLogProbs(1:q, :), noiseLevelGP, lowestLogliklVal, ...
    logLiklRange, cvCostFunc);
  ucGTEsts = uclj(gtPts);

%   randlj = regressionWrapLRG(randPts(1:q, :), randLogProbs(1:q, :), noiseLevelGP, ...
%     lowestLogliklVal, logLiklRange, cvCostFunc, paramSpaceBounds);
  randlj = regressionWrap(randPts(1:q, :), randLogProbs(1:q, :), noiseLevelGP, ...
    lowestLogliklVal, logLiklRange, cvCostFunc);
  randGTEsts = randlj(gtPts);

%   mrlj = regressionWrapLRG(mrPts(1:q, :), mrLogProbs(1:q, :), noiseLevelGP, lowestLogliklVal, ...
%     logLiklRange, cvCostFunc, paramSpaceBounds);
  mrlj = regressionWrap(mrPts(1:q, :), mrLogProbs(1:q, :), noiseLevelGP, lowestLogliklVal, ...
    logLiklRange, cvCostFunc);
  mrGTEsts = mrlj(gtPts);

  % Now do the evaluation
  errFunc = @(a, b) mean( (exp(a) - exp(b)).^2 );
  ucErr = errFunc(ucGTEsts, gtLogProbs),
  randErr = errFunc(randGTEsts, gtLogProbs),
  mrErr = errFunc(mrGTEsts, gtLogProbs),
  ucErrs(i) = ucErr;
  mrErrs(i) = mrErr;
  randErrs(i) = randErr;
end

[ucErrs, randErrs, mrErrs],
loglog(1:numel(numPtsToUse), ucErrs, 'b'); hold on,
loglog(1:numel(numPtsToUse), randErrs, 'g');
loglog(1:numel(numPtsToUse), mrErrs, 'm');
