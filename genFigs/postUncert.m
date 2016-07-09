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

yLabPosn = [0.12 1.1];
xLabPosn = [0.96 -0.85];

c1 = [150 75 0]/255;
c2 = 'r';


% First plot the function
th = linspace(0,1,100)';

% X = [0.1 0.3 0.6 0.77 0.95]';
% X = [0.1 0.4 0.75 0.95]';
X = [0.1 0.4 0.55 0.78 0.95]';
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

% figure;
[mu, K, funcH, bw, scale] = GPMargLikelihood(X, Y, th, hyperparams);
% Create shaded area
top = mu + width*diag(K);
bottom = mu - width*diag(K);
A = [th' fliplr(th')];
B = [top' fliplr(bottom')];

% % Just the function
% figure;
% plot(th, func(th), 'k--', 'LineWidth', LW); hold on,
% axis([-0.005 1.005 (-0.6394-confRegionMargin) (1.0911+confRegionMargin)]);
% set(gca, 'Xtick', []);
% set(gca, 'Ytick', []); % set(gca,'position',[0.14 0.09 0.85 0.9],'units','normalized'); % xlabel('$\theta$', 'Position', xLabPosn, 'FontSize', FS, 'Interpreter', 'Latex');
% ylabel('$f(\theta)$', 'Position', yLabPosn, 'FontSize', FS, ...
%   'Interpreter', 'Latex', 'rot', 0);
% 
% % Function and points
% figure;
% plot(th, func(th), 'k--', 'LineWidth', LW); hold on,
% plot(X, Y, 'kx', 'MarkerSize', MS, 'LineWidth', LW);
% axis([-0.005 1.005 (-0.6394-confRegionMargin) (1.0911+confRegionMargin)]);
% set(gca, 'Xtick', []);
% set(gca, 'Ytick', []);
% set(gca,'position',[0.14 0.09 0.85 0.9],'units','normalized');
% xlabel('$\theta$', 'Position', xLabPosn, 'FontSize', FS, 'Interpreter', 'Latex');
% ylabel('$f(\theta)$', 'Position', yLabPosn, 'FontSize', FS, ...
%   'Interpreter', 'Latex', 'rot', 0);

augFactor = 3;
expF = @(t) exp(augFactor*func(t));
normC = mean(expF(th));
postF = @(t) expF(t)/normC;

figure;
h = fill(A, B, [0.9 0.9 0.9]);
% set(h, 'EdgeColor', 'None'); hold on,
plot(th, func(th), 'k--', 'LineWidth', LW); hold on,
gpPostSamples = GPDrawSamples(mu, K, numSamples);
plot(th, gpPostSamples, 'LineWidth', LWS);
plot(X, Y, 'kx', 'MarkerSize', MS, 'LineWidth', LW);
axis([-0.005 1.005 (-0.6394-confRegionMargin) (1.0911+confRegionMargin)]);
set(gca, 'Xtick', []);
set(gca, 'Ytick', []);
% set(gca,'position',[0.14 0.09 0.85 0.9],'units','normalized');
set(gca,'position',[0.02 0.09 0.97 0.9],'units','normalized');
xlabel('$\theta$', 'Position', xLabPosn, 'FontSize', FS, 'Interpreter', 'Latex');
% ylabel('$f_{\bf{X_{obs}}}(\theta)$', 'Position', yLabPosn, 'FontSize', FS, ...
%   'Interpreter', 'Latex', 'rot', 0);
ylabel('Log Joint $GP$', 'Position', yLabPosn + [0.1 0], 'FontSize', FS, ...
  'Interpreter', 'Latex', 'rot', 0);
box on

% Posterior
figure;
plot(th, postF(th), 'k--', 'LineWidth', LW); hold on,
hSamples = exp(augFactor*gpPostSamples);
hNormCs = mean(hSamples, 2);
hSamples = bsxfun(@rdivide, hSamples, hNormCs);
plot(th, hSamples, 'LineWidth', LWS);
set(gca, 'Xtick', []);
set(gca, 'Ytick', []);
% set(gca,'position',[0.14 0.09 0.85 0.9],'units','normalized');
set(gca,'position',[0.02 0.09 0.97 0.9],'units','normalized');
xlabel('$\theta$', 'Position', [xLabPosn(1) -0.03], 'FontSize', FS, 'Interpreter', 'Latex');
% ylabel('$h_{\bf{X_{obs}}}(\theta)$', 'Position',[yLabPosn(1) 4.55], 'FontSize', FS, ...
%   'Interpreter', 'Latex', 'rot', 0);
ylabel('$F_{\theta|X_{obs}}$', 'Position',[yLabPosn(1), 4.35], 'FontSize', 1.3*FS, ...
  'Interpreter', 'Latex', 'rot', 0);
box on
