% Unit test for SDSSExperiment and sdssLogLiklWrap.m

fprintf('Unit Test for lssLogLiklWrap.m\n');
evalAtPts =  [ 0 0.762 0.1045 0.02233 0.951 0.6845 0 1.908 30.81; ...
               0 0.684 0.1045 0.02233 0.951 0.6845 0 1.908 30.81; ...
               0 0.462 1.0000 0.24 1.051 0.6845 0 2.908 30.81; ...
               0 0.662 0.2045 0.01 1.651 0.6845 0 0.908 30.81; ...
               0 0.262 0.4045 0.005 0.851 0.6845 0 1.208 30.81; ...
               0 0.962 0.6045 0.15 0.551 0.6845 0 0.408 30.81; ...
               -0.5 0.762 0.1045 0.22233 1.951 0.6845 0 1.908 30.81];
sdssLogLiklWrap(evalAtPts, -100),               

fprintf('\n UT for SDSS6Experiment\n');

sdss6 = SDSS6Experiment();
evalPts = [evalAtPts(:, 1:5) evalAtPts(:, 8)];
normEvalPts = sdss6.getNormCoords(evalPts),

normEvalPts = [normEvalPts; rand(10, 6)];
sdss6.getTrueCoords(normEvalPts),
sdss6.normCoordLogJointProbs(normEvalPts),

normEvalPts = [rand(10000, 6)];
tic,
ljps = sdss6.normCoordLogJointProbs(normEvalPts);
toc,
sortedLjps = sort(ljps, 'descend');
sortedLjps(1:10),

