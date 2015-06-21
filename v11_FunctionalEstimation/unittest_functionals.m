% Unit test for the functionals f1, f2, f3

% add required functions
addpath ../helper

N = 20000;
dims = 15;

p1 = 0.8;
p2 = 1 - p1;
sigma = 0.1 + 0.05*rand();

B = double( rand(N, 1) < p1);
X = bsxfun(@times, sigma*randn(N, dims), B) + ...
    bsxfun(@times, 1 + sigma*randn(N, dims), 1-B);
close all;
plot(X(:,1), X(:,2), 'bx');
fprintf('Sigma: %f\n', sigma);

T1 = f1(X);
T2 = f2(X);
T3 = f3(X);
T4 = f4(X);

S1 = p2*dims;
S2 = dims*(p1*sigma^2 + p2*(1 + sigma^2) );
S3 = (dims -2) * p2;
S4 = (dims -1) * (1+sigma^2) * p2; 

fprintf('Empirical values: f1: %f, f2: %f, f3: %f, f4: %f\n', T1, T2, T3, T4);
fprintf('True values     : f1: %f, f2: %f, f3: %f, f4: %f\n', S1, S2, S3, S4);

% Test for evalLogLiklExp1.m
fprintf('\nTest for evalLogLiklExp1.m\n');
RES = 100;
t = linspace(-3, 4, RES);
[T1, T2] = meshgrid(t, t);
th = [T1(:), T2(:)];
[loglikl, f_vals] = evalLogLiklExp1(th, sigma, p1, p2);
likl = exp(loglikl);
Likl = reshape(likl, RES, RES);
figure; contour(T1, T2, Likl);
fprintf('Area under curve: %f\n', numerical_2D_integration_wrap(th, likl));
f_vals,

% Another test for evalLogLiklExp1
fprintf('\nAnother test for evalLogLiklExp1.m\n');
dims = 10;
Pts = randn(0,dims);
X = gendata(N, dims, sigma, p1);
[~, f_vals] = evalLogLiklExp1(Pts, sigma, p1, p2);
T1 = f1(X);
T2 = f2(X);
T3 = f3(X);
T4 = f4(X);
fprintf('Empirical values: f1: %f, f2: %f, f3: %f, f4: %f\n', T1, T2, T3, T4);
f_vals,

