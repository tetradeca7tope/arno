% Unit test for estimatePostProb.m

% Test 1: 1D
% ==========
alpha = 6; beta = 2;
fprintf('1 Dimensional Beta(%d, %d) Distribution\n', alpha, beta);
n = 10;
t = rand(n, 1);
p = betapdf(th, alpha, beta);
% Now do the regression
post_est = estim
