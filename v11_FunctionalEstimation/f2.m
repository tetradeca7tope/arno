function T = f2(X)
% estimates E[sum Xi^2]
  T = sum(sum(X.^2)) / size(X, 1);
end
