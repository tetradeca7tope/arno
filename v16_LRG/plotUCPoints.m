% plotUCPoints.m

load gt10_v1_1
% load gt10_06-04_17-43-01.mat

numDims = 8;
lowestLogliklVal = -200;
logLiklRange = 200;
indices = nchoosek(1:numDims, 2);
NUM_AL_ITERS = 6000;

paramSpaceBounds = repmat([0 1], numDims, 1);
lrgExp = LRGExperiment(lowestLogliklVal);
evalLogJoint = @(arg) lrgExp.normCoordLogJointProbs(arg);
MLEPoint = lrgExp.getNormCoords([0 0.762 0.1045 0.02233 0.951 0.6845 0 1.908]);
MLELogP = evalLogJoint(MLEPoint);

sucp = sort(ucLogProbs, 'descend'); sucp(1:10),
smcmcp = sort(mcmcLogProbs, 'descend'); smcmcp(1:10),
mrp = mcmcLogProbs(1:NUM_AL_ITERS, :); smrp = sort(mrp, 'descend'); smrp(1:10),
srandp = sort(randLogProbs, 'descend'); srandp(1:10),

% Plot the quries
queryAxes = [0, 6000, -100, 0];
figure; plot(ucLogProbs); axis(queryAxes);
figure; plot(randLogProbs); axis(queryAxes);
figure; plot(mrLogProbs); axis(queryAxes);

for i = 1:size(indices, 1)

  figure;
  i1 = indices(i, 1);
  i2 = indices(i, 2);

  threshold = -50;

%   subplot(1,2,1);
%   plot(ucPts(:, i1), ucPts(:,i2), 'kx'); hold on,
%   highLiklIndices = ucLogProbs > threshold;
%   plot(ucPts(highLiklIndices, i1), ucPts(highLiklIndices, i2), 'co');
%   axis([0 1 0 1]);

%   subplot(2,2,1);
  figure;
  plot(mcmcQueries(:, i1), mcmcQueries(:,i2), 'cx'); hold on,
  highMCMCLiklIndices = mcmcLogProbs > threshold;
  plot(mcmcQueries(highMCMCLiklIndices, i1), mcmcQueries(highMCMCLiklIndices, i2), 'ro');
  axis([0 1 0 1]);
%   plot(MLEPoint(:,i1), MLEPoint(:,i2), 'ms', 'MarkerSize', 8, 'LineWidth', 3);

%   subplot(2,2,2);
  figure;
  plot(ucPts(:, i1), ucPts(:,i2), 'cx'); hold on,
  highLiklIndices = ucLogProbs > threshold;
  plot(ucPts(highLiklIndices, i1), ucPts(highLiklIndices, i2), 'ro');
  axis([0 1 0 1]);
%   plot(MLEPoint(:,i1), MLEPoint(:,i2), 'ms', 'MarkerSize', 8, 'LineWidth', 3);

%   subplot(2,2,3);
  figure;
  plot(randPts(:, i1), randPts(:,i2), 'cx'); hold on,
  highRandLiklIndices = randLogProbs > threshold;
  plot(randPts(highRandLiklIndices, i1), randPts(highRandLiklIndices, i2), 'ro');
  axis([0 1 0 1]);
%   plot(MLEPoint(:,i1), MLEPoint(:,i2), 'ms', 'MarkerSize', 8, 'LineWidth', 3);

%   mrPts = mcmcQueries(1:NUM_AL_ITERS, :);
  figure;
  subplot(2,2,4);
  plot(mrPts(:,i1), mrPts(:,i2), 'yx'); hold on,
  highMRPts = mrp > threshold;
  plot(mrPts(highMRPts, i1), mrPts(highMRPts, i2), 'rd');
  axis([0 1 0 1]);
%   plot(MLEPoint(:,i1), MLEPoint(:,i2), 'ms', 'MarkerSize', 8, 'LineWidth', 3);
  pause; close;

end

