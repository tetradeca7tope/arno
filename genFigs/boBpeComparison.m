% Bayesian Posterior Estimation
close all;
clear all;

f = @(t) - 70*(t+0.7).* (t+0.1).* (t-0.3).* (t-0.8);

x1 = -0.8;
x2 = 0.9;
xPts = [-0.4 0 0.45 0.6 0.75];
qPts = [0.2 0.52];

LW = 4;
MS = 15;

c1 = [150 75 0]/255;
c2 = 'r';

x = linspace(x1, x2, 1000);

figure;
plot(x, f(x), 'LineWidth', LW ); hold on
plot(xPts, f(xPts), 'x', 'LineWidth', LW, 'MarkerSize', MS, 'color', c1);
plot(qPts, f(qPts), 'o', 'LineWidth', LW, 'MarkerSize', MS, 'color', c2);
xlim([x1 x2]);
set(gca, 'XTick', [], 'YTick', []);
set(gca,'position',[0.005 0.005 0.99 0.99],'units','normalized');

figure;
plot(x, exp(f(x)), 'LineWidth', LW ); hold on
plot(xPts, exp(f(xPts)), 'x', 'LineWidth', LW, 'MarkerSize', MS, 'color', c1);
plot(qPts, exp(f(qPts)), 'o', 'LineWidth', LW, 'MarkerSize', MS, 'color', c2);
set(gca, 'XTick', [], 'YTick', []);
% text(qPts(1), exp(f(qPts(1))), '(a)', 'color', 'r');
% text(qPts(2), exp(f(qPts(2))), '(b)', 'color', 'r');
% set(gca,'DefaultTextFontSize',14)
xlim([x1 x2]);
ylim([-10 350]);
set(gca,'position',[0.005 0.005 0.99 0.99],'units','normalized');
