% Unit test for function estimKLForLRG.m

numDims = 8;

m1 = zeros(numDims, 1);
m2 = ones(numDims, 1);
S1 = randn(numDims); S1 = S1* S1';
S2 = randn(numDims); S2 = S2* S2';
% m1 = zeros(numDims, 1);
% m2 = ones(numDims, 1);
% S1 = eye(numDims); %randn(numDims); S1 = S1* S1';
% S2 = S1; %randn(numDims); S2 = S2* S2';

% Generate samples
numSamples = 100000;
X = bsxfun(@plus, randn(numSamples, numDims) * chol(S1), m1');
pAtX = mvnpdf(X, m1', S1);
qFuncHandle = @(t) mvnpdf(t, m2', S2);

[kl, l2] = estimKLForLRG(X, pAtX, qFuncHandle);
fprintf('Estimated: KL: %0.6f\n', kl);
trueKL = 0.5* ( log(det(S2)/det(S1)) - numDims + trace(S2\S1) + ...
                (m2- m1)' * ( S1\ (m2 - m1) ) );
fprintf('True     : KL: %0.6f\n', trueKL);

