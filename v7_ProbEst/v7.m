
clear all;
close all;

% First generate some data and the corresponding probabilites
Z = 3*randn(20,1);
pdf = @(lambda) (0.64*normpdf(lambda, -1, 1.2) + 0.36*normpdf(lambda, 2, 0.5));
P = pdf(Z);
q = max(P);
P = P / max(P); % renormalizing

% Now obtain an estimate
[w_opt, h_opt] = ProbEst(Z, P);

% plot the final result
figure;
th = linspace(-6,6,1000)';
p_th = evaluate_prob(th, Z, w_opt, h_opt);
plot(th, p_th, 'b', Z, P*q, 'bx', th, pdf(th), 'r-.');

