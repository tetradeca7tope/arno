function [lPtGrid] = genLeftPointGrid(res, bounds)
% Generates the "left" points of a grid generated on the parameter space with
% resolution.

  numDims = size(bounds, 1);
  linT = zeros(numDims, res);
  for d = 1:numDims
    temp = linspace(bounds(d,1), bounds(d,2), res+1);
    linT(d, :) = temp(1:res);
  end

end
