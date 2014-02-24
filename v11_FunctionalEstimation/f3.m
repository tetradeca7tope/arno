function T = f3(X)

  n = size(X, 1);
  d = size(X, 2);

  sum_terms = zeros(n, 1);
  for i = 1:(d-2)
    sum_terms = sum_terms + X(:,i) .* X(:,i+1) .* X(:,i+2);
  end
  T = sum(sum_terms)/n;
end
