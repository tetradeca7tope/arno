function T = f1(X)
% Estimates the functional sum(Xi) from the data
  T = sum(sum(X))/size(X, 1);
end
