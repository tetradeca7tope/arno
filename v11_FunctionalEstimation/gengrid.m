function [vertices] = gengrid(num_dims, nodes)
% nodes is a 1D vector specificying at which points along each axis we should
% place a node.

  n = size(nodes, 1);
  num_grid_pts = n^num_dims;

  grid = cell(1, num_dims);
  [grid{:}] = ndgrid(nodes);

  vertices = zeros(num_grid_pts, num_dims);
  for i = 1:num_dims
    vertices(:,i) = grid{i}(:);
  end
end
