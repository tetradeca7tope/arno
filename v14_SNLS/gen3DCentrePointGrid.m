function [g1, g2, g3] = gen3DCentrePointGrid(res, bounds);

  [g1, g2, g3] = gen3DLeftPointGrid(res, bounds);
  g1 = g1 + (bounds(1, 2) - bounds(1,1)) / (2 * res);
  g2 = g2 + (bounds(2, 2) - bounds(2,1)) / (2 * res);
  g3 = g3 + (bounds(3, 2) - bounds(3,1)) / (2 * res);

end
