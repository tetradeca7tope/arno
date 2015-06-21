% plot Smaple results

load allCollectedSamples

indices = [1 2; 1 3; 2 3];

for i = 1:3

  figure;
  i1 = indices(i, 1);
  i2 = indices(i, 2);
  plot(mcmcEstSamps(:, i1), mcmcEstSamps(:,i2), 'c.'); hold on
  plot(ucPts(:,i1), ucPts(:,i2), 'kx');
  



end
