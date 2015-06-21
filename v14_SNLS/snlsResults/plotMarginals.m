% plotMarginals.m

addpath ..
close all;
trueLogPostValsAtEvalPts = load('truePost_nObs_res60_H080-100.txt');
truePost = exp(2*trueLogPostValsAtEvalPts);
truePostGrid = reshape(truePost, 100, 100, 100);

indices = [1 2; 1 3; 2 3];
idxLeftOut = [3; 2; 1];
[g1, g2, g3] = gen3DCentrePointGrid(100, repmat([0 1], 3, 1) );
grids = {g1, g2, g3};
t1 = g1(:,:,1);
t2 = g2(:,:,1);

for i = 1:3

  i1 = indices(i, 1);
  i2 = indices(i, 2);

  if i ~= 3 % rescaling for H0
    s1 = 60 + 20*t1;
    pts1 = 60 + 20*ucPts(:,i1);
  else
    s1 = t1;
    pts1 = ucPts(:,i1);
  end

  if i1 == 1, xlabelstr = 'H_0';
  else xlabelstr = '\Omega_M';
  end
  if i2 == 3, ylabelstr = '\Omega_{\Lambda}';
  else ylabelstr = '\Omega_M';
  end


  figure;
%   subplot(1,2,1);
  marginal = sum(truePostGrid, idxLeftOut(i) );
  marginal = reshape(marginal, 100, 100);
  contour(s1, t2, marginal); hold on,
%   xlabel(xlabelstr);
%   ylabel(ylabelstr);

  figure;
%   subplot(1,2,2);
  plot(pts1, ucPts(:,i2), 'kx'); hold on
%   xlabel(xlabelstr);
%   ylabel(ylabelstr);

end
