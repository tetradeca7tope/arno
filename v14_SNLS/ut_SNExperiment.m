% Unit test for the class SNExperiment

paramSpaceBounds = [60 80; 0 1; 0 1];
davisSN = SNExperiment('davisdata.txt', paramSpaceBounds);

% To test getNormCoords
trueCoords = [70 0.24 0.68; 41 0.5 0.6; 139 0.789 0.21];
normCoords = davisSN.getNormCoords(trueCoords),
trueCoords = davisSN.getTrueCoords(normCoords),

% Now evaluate the log Likelihood at the above point
[davisSN.normCoordLogJointProbs(normCoords), ...
davisSN.trueCoordLogJointProbs(trueCoords)],

% Plot the marginal log likelihood
N = 100;
testPts = [linspace(0,1, N)', 0.268*ones(N,1), 0.683*ones(N,1)];
logProbs = davisSN.normCoordLogJointProbs(testPts);
trueTestCoords = davisSN.getTrueCoords(testPts);
trueTestCoords = trueTestCoords(:,1);
figure; plot(trueTestCoords, logProbs);
axis([min(trueTestCoords), max(trueTestCoords), min(logProbs), max(logProbs)]);
figure; plot(trueTestCoords, exp(logProbs));
axis([min(trueTestCoords), max(trueTestCoords), 0, max(exp(logProbs))]);
