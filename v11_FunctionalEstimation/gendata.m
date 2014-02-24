function X = gendata(N, dims, sigma, p1)

  B = double( rand(N, 1) < p1);
  X = bsxfun(@times, sigma*randn(N, dims), B) + ...
      bsxfun(@times, 1 + sigma*randn(N, dims), 1-B);

end
