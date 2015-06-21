% 

close all;

numPts = 3750;
subSet = [1: (numPts)];
subSet = subSet(randperm(numPts));

sgtPts = gtPts(subSet, :);
sgtLogProbs = gtLogProbs(subSet) + 0.7*randn(size(gtLogProbs(subSet)));

k = sgtLogProbs > - 20;
fprintf('k = %d\nContinue ? ...\n', sum(k));
% pause;

truth = sgtLogProbs(k);
uclps = uclj( sgtPts(k,:) );
randlps = randlj( sgtPts(k,:) );
mrlps = mrlj( sgtPts(k,:) );

plot(truth, 'k-d'); hold on,
plot(mrlps, 'r-x');
plot(randlps, 'g-s');
plot(uclps, 'b-o');
legend('Truth', 'MCMC-R', 'RAND', 'VR');

axis([0 sum(k) -22 -8]);
% xticklabel('off');

set(0,'defaultAxesFontName', 'Dejavu Sans')
  set(findall(gca, '-property', 'FontSize'), 'FontSize', 20, ...
    'fontWeight', 'bold');
ylabel('Log Likelihood');
set(gca, 'xticklabel', '');
