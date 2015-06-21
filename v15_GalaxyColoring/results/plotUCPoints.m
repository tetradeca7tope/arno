% plotUCPoints.m

close all;
% fileName = 'expUC2.mat';
% load(fileName);

k = 2;
k = 3;
indices = nchoosek(1:12, k);

sucp = sort(ucLogProbs, 'descend'); sucp(1:10),
numHighLogPpts = sum(ucLogProbs > -1000);
fprintf('Num high likl ponts : %d\n', numHighLogPpts);

for i = 1:size(indices, 1)

  i1 = indices(i, 1);
  i2 = indices(i, 2);
  if k == 3
    i3 = indices(i, 3);
  end

  threshold = -1000;
  highLiklIndices = ucLogProbs > threshold;

  if k == 2
    plot(ucPts(:, i1), ucPts(:,i2), 'gx'); hold on,
    plot(ucPts(highLiklIndices, i1), ucPts(highLiklIndices, i2), 'ko');
  elseif k == 3
    plot3(ucPts(:, i1), ucPts(:,i2), ucPts(:,i3), 'gx'); hold on,
    plot3(ucPts(highLiklIndices, i1), ucPts(highLiklIndices, i2), ...
      ucPts(highLiklIndices, i3), 'ko');
  end
    
  axis([0 1 0 1 0 1]);

  hold off,
  pause;
%   close;
end
