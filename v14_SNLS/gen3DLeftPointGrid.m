function [g1, g2, g3] = gen3DLeftPointGrid(res, bounds)

  t1 = linspace(bounds(1,1), bounds(1,2), res+1); t1 = t1(1:res);
  t2 = linspace(bounds(2,1), bounds(2,2), res+1); t2 = t2(1:res);
  t3 = linspace(bounds(3,1), bounds(3,2), res+1); t3 = t3(1:res);

  [g1, g2, g3] = ndgrid(t1, t2, t3);

end

