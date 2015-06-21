% plotUCPoints.m

indices = nchoosek(1:12, 2);
fileName = 's1.mat';
load(fileName); 

numUniquePoints = size(unique(samples(:,1)), 1);
slps = sort(logProbs, 'descend'); slps(1:10, :),
fprintf('number of unique points: %d\n', numUniquePoints);

for i = 1:size(indices, 1)

  figure;
  i1 = indices(i, 1);
  i2 = indices(i, 2);

  plot(samples(:, i1), samples(:, i2), 'kx');
  pause;

end

