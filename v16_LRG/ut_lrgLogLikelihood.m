% Unit test for SDSSExperiment and sdssLogLiklWrap.m

fprintf('Unit Test for lssLogLiklWrap.m\n');
evalAtPts =  [ 0 0.762 0.1045 0.02233 0.951 0.6845 0 1.908 30.81; ...
               0 0.762 0.1045 0.02233 0.951 .75 -0.0 1.908 30.81; ...
               0 0.684 0.1045 0.02233 0.951 0.6845 0 1.908 30.81; ...
               0 0.462 1.0000 0.24 1.051 0.6845 0 2.908 30.81; ...
               0 0.662 0.2045 0.01 1.651 0.6845 0 0.908 30.81; ...
               0 0.262 0.4045 0.005 0.851 0.6845 0 1.208 30.81; ...
               0 0.962 0.6045 0.15 0.551 0.6845 0 0.408 30.81; ...
               -0.5 0.762 0.1045 0.22233 1.951 0.6845 0 1.908 30.81];
lrgLogLiklWrap(evalAtPts, -100),               

fprintf('\n UT for SDSS6Experiment\n');

lrg = LRGExperiment();
evalPts = [evalAtPts(:, 1:8)];
normEvalPts = lrg.getNormCoords(evalPts),

normEvalPts = [normEvalPts; rand(10, 8)];
lrg.getTrueCoords(normEvalPts),
lrg.normCoordLogJointProbs(normEvalPts),

normEvalPts = [rand(10000, 8)];
tic,
ljps = lrg.normCoordLogJointProbs(normEvalPts);
toc,
sortedLjps = sort(ljps, 'descend');
sortedLjps(1:10),

