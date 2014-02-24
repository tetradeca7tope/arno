function [loglikl, f_vals] = evalLogLiklExp1(X, sigma, p1, p2)
% Theta is the value at which the function should be evaluated.

  dims = size(X, 2);

  M1 = zeros(1, dims);
  M2 = ones(1, dims);
  S = sigma^2 * ones(1, dims);
  
  likl = p1 * mvnpdf(X, M1, S) + p2*mvnpdf(X, M2, S);
  loglikl = log(likl);
  loglikl = max(loglikl, -30);

  % Now evaluate the functionals.
  P1 = p1/(p1 + p2);
  P2 = 1 - P1;

  % Returns the value of the functionals f1,.. f4  
  f_vals.S1 = P2*dims;
  f_vals.S2 = dims*(P1*sigma^2 + P2*(1 + sigma^2) );
  f_vals.S3 = (dims -2) * P2;
  f_vals.S4 = (dims -1) * (1+sigma^2) * P2; 

end
