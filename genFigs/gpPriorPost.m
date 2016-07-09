% Generates some figures for BO/GPB
close all;
clear all;
addpath ~/libs/kky-matlab/GPLibkky/
addpath ~/libs/kky-matlab/utils/
rng(40);

N = 100;
numSamples = 10;
width = 5;
confRegionMargin = 0.2;

func = @(t) - 70*(t-0).* (t-0.35).* (t+0.55).* (t-0.65).* (t-0.97);

% LW = 3;
% MS = 15;
% LWS = 1.5;
LW = 4;
MS = 25;
LWS = 2;
FS = 25;

yLabPosn = [-0.08 1.1];
xLabPosn = [0.96 -0.85];

c1 = [150 75 0]/255;
c2 = 'r';


% First plot the function
th = linspace(0,1,100)';

% X = [0.1 0.3 0.6 0.77 0.95]';
% X = [0.1 0.4 0.75 0.95]';
X = [0.1 0.4 0.55 0.78 0.95]';
% X = [0.1 0.4 0.55 0.78 0.95 0.05 0.175 0.25 0.325 0.475 0.63 0.7 0.84 0.98]';
% X = [0.1 0.4 0.55 0.78 0.95  0.25 0.48 0.66 0.84 0.98]';
% X = [0.1 0.4 0.55 0.95]';
% X = [0.1 0.4 0.70 0.80 0.95]';
% X = [0.1 0.4 0.6 0.8 0.98]';
Y = func(X);

hyperparams.meanFunc = [];
hyperparams.sigmaSmRange = [];
hyperparams.sigmaPrRange = [];
hyperparams.sigmaSm = 0;
hyperparams.sigmaPr = 0;
hyperparams.noise = 0.0001;

figure;
[mu, K, funcH, bw, scale] = GPMargLikelihood(X, Y, th, hyperparams);
% Create shaded area
top = mu + width*diag(K);
bottom = mu - width*diag(K);
A = [th' fliplr(th')];
B = [top' fliplr(bottom')];

% Just the function
figure;
plot(th, func(th), 'k--', 'LineWidth', LW); hold on,
axis([-0.005 1.005 (-0.6394-confRegionMargin) (1.0911+confRegionMargin)]);
set(gca, 'Xtick', []);
set(gca, 'Ytick', []);
set(gca,'position',[0.14 0.09 0.85 0.9],'units','normalized');
xlabel('$x$', 'Position', xLabPosn, 'FontSize', FS, 'Interpreter', 'Latex');
ylabel('$f(x)$', 'Position', yLabPosn, 'FontSize', FS, ...
  'Interpreter', 'Latex', 'rot', 0);
% xlabel('$\theta$', 'Position', xLabPosn, 'FontSize', FS, 'Interpreter', 'Latex');
% ylabel('$f(\theta)$', 'Position', yLabPosn, 'FontSize', FS, ...
%   'Interpreter', 'Latex', 'rot', 0);

% Function and points
figure;
plot(th, func(th), 'k', 'LineWidth', LW); hold on,
plot(X, Y, 'rx', 'MarkerSize', 1*MS, 'LineWidth', 2*LW);
plot(th, mu, 'r--', 'LineWidth', LW);
axis([-0.005 1.005 (-0.6394-confRegionMargin) (1.0911+confRegionMargin)]);
set(gca, 'Xtick', []);
set(gca, 'Ytick', []);
set(gca,'position',[0.14 0.09 0.85 0.9],'units','normalized');
xlabel('$x$', 'Position', [0.96 -0.75], 'FontSize', 1.5*FS, 'Interpreter', 'Latex');
ylabel('$f_*$', 'Position', [-0.05 0.95], 'FontSize', 1.5*FS, ...
  'Interpreter', 'Latex', 'rot', 0);
text(0.67, 0.7, '$\hat{f}$', 'FontSize', 1.5*FS, 'Color', 'red', 'Interpreter', 'Latex');
pause;
% ylabel('$f(x)$', 'Position', yLabPosn, 'FontSize', FS, ...
%   'Interpreter', 'Latex', 'rot', 0);
% xlabel('$\theta$', 'Position', xLabPosn, 'FontSize', FS, 'Interpreter', 'Latex');
% ylabel('$f(\theta)$', 'Position', yLabPosn, 'FontSize', FS, ...
%   'Interpreter', 'Latex', 'rot', 0);

figure;
h = fill(A, B, [0.9 0.9 0.9]);
set(h, 'EdgeColor', 'None'); hold on,
plot(th, func(th), 'k--', 'LineWidth', LW); hold on,
gpPostSamples = GPDrawSamples(mu, K, numSamples);
plot(th, gpPostSamples, 'LineWidth', LWS);
plot(X, Y, 'kx', 'MarkerSize', MS, 'LineWidth', LW);
axis([-0.005 1.005 (-0.6394-confRegionMargin) (1.0911+confRegionMargin)]);
set(gca, 'Xtick', []);
set(gca, 'Ytick', []);
set(gca,'position',[0.14 0.09 0.85 0.9],'units','normalized');
xlabel('$x$', 'Position', xLabPosn, 'FontSize', FS, 'Interpreter', 'Latex');
ylabel('$f(x)$', 'Position', yLabPosn, 'FontSize', FS, ...
  'Interpreter', 'Latex', 'rot', 0);
% xlabel('$\theta$', 'Position', xLabPosn, 'FontSize', FS, 'Interpreter', 'Latex');
% ylabel('$f(\theta)$', 'Position', yLabPosn, 'FontSize', FS, ...
%   'Interpreter', 'Latex', 'rot', 0);
box on

% Prior
figure;
tunedHPs = hyperparams;
tunedHPs.sigmaSm = bw;
tunedHPs.sigmaPr = scale;
muPrior = ones(size(th,1), 1) * mean(Y);
D = Dist2GP(th, th);
KPrior = scale * exp(-0.5*D/bw^2);
gpPriorSamples = GPDrawSamples(muPrior, KPrior, numSamples);
top = muPrior + width*diag(KPrior);
bottom = muPrior - width*diag(KPrior);
A = [th' fliplr(th')];
B = [top' fliplr(bottom')];

h = fill(A, B, [0.9 0.9 0.9]);
set(h, 'EdgeColor', 'None'); hold on,


plot(th, func(th), 'k--', 'LineWidth', LW); hold on,
plot(th, gpPriorSamples, 'LineWidth', LWS);
axis([-0.005 1.005 (-0.6394-confRegionMargin) (1.0911+confRegionMargin)]);
set(gca, 'Xtick', []);
set(gca, 'Ytick', []);
xlabel('$x$', 'Position', xLabPosn, 'FontSize', FS, 'Interpreter', 'Latex');
ylabel('$f(x)$', 'Position', yLabPosn, 'FontSize', FS, ...
  'Interpreter', 'Latex', 'rot', 0);
% xlabel('$\theta$', 'Position', xLabPosn, 'FontSize', FS, 'Interpreter', 'Latex');
% ylabel('$f(\theta)$', 'Position', yLabPosn, 'FontSize', FS, ...
%   'Interpreter', 'Latex', 'rot', 0);
set(gca,'position',[0.14 0.09 0.85 0.9],'units','normalized');
box on,


