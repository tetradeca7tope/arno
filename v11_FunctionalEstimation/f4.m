function T = f4(X)

  n = size(X, 1);
  d = size(X, 2);

  sum_terms = zeros(n, 1);
  for i = 1:(d-1)
    sum_terms = sum_terms + (X(:,i).^2) .* X(:,i+1);
  end
  T = sum(sum_terms)/n;
  
end
