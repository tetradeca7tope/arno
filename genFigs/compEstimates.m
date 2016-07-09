% compares two estimates

clear all;
close all;
clc;

addpath ~/libs/kky-matlab/utils/

LW = 4;
MS = 25;
LWS = 2;
FS = 26;
normalisedPosn = [0.02 0.1 0.96 0.8];

N = 100;
th = linspace(0,1,N)';

f = @(t) -200 * (t-0.1) .* (t-0.3) .* (t-0.95);
expf = @(t) exp(f(t));
normC = mean( expf(th) );
postf = @(t) expf(t)/normC;


% log joint 
figure; plot(th, f(th), 'k--', 'LineWidth', LW);
set(gca, 'Ytick', []);
set(gca, 'Xtick', []);
title('$\log P(\theta,X_{obs})$', 'Interpreter', 'Latex', 'FontSize', FS);
xlabel('$\Theta$', 'Position', [0.5 -6.5], 'FontSize', FS, 'Interpreter','Latex');
axis([0 1 -6 14]);
set(gca,'position',normalisedPosn,'units','normalized');

% Joint
figure; plot(th, postf(th), 'k--', 'LineWidth', LW);
set(gca, 'Ytick', []);
set(gca, 'Xtick', []);
title('$P(\theta|X_{obs})$', 'Interpreter', 'Latex', 'FontSize', FS);
xlabel('$\Theta$', 'Position', [0.5 -0.2], 'FontSize', FS, 'Interpreter','Latex');
axis([0 1 0 7.5]);
set(gca,'position',normalisedPosn,'units','normalized');


% Now create the estimates
numPts = 8;
order = 1;

% greenPts = linspace(0,1,numPts+2)'; greenPts = greenPts(2:(end-1));
greenPts = linspace(0.01,0.97, numPts)'; 
% greenPts = linspace(0,1, numPts)'; 
greenVals = f(greenPts);
greenEstimates = localPolyKRegressionCV(th, greenPts, greenVals, [], [order]);
greenExp = exp(greenEstimates);
greenNormC = mean(greenExp);
% greenNormC = normC;
greenExp = greenExp/greenNormC;
% Log joint
figure; 
plot(th, f(th), 'k--', 'LineWidth', LW); hold on,
plot(th, greenEstimates, 'g-', 'LineWidth', LW);
plot(greenPts, greenVals, 'kx', 'MarkerSize', MS, 'LineWidth', LW);
set(gca, 'Ytick', []);
set(gca, 'Xtick', []);
title('$\log P(\theta,X_{obs})$', 'Interpreter', 'Latex', 'FontSize', FS);
xlabel('$\Theta$', 'Position', [0.5 -6.5], 'FontSize', FS, 'Interpreter','Latex');
axis([0 1 -6 14]);
set(gca,'position',normalisedPosn,'units','normalized');
% Posterior
figure;
plot(th, postf(th), 'k--', 'LineWidth', LW); hold on,
plot(th, greenExp, 'g-', 'LineWidth', LW);
set(gca, 'Ytick', []);
set(gca, 'Xtick', []);
title('$P(\theta|X_{obs})$', 'Interpreter', 'Latex', 'FontSize', FS);
xlabel('$\Theta$', 'Position', [0.5 -0.2], 'FontSize', FS, 'Interpreter','Latex');
axis([0 1 0 7.5]);
set(gca,'position',normalisedPosn,'units','normalized');

switch numPts

  case 6
    pinkPts = [0.05 0.4 0.90 0.6 0.7 0.8]';

  case 8
    pinkPts = [linspace(0.55, 0.707, 3)  linspace(0.77, 0.85, 2), 0.1  0.25 0.4 ]';

  case 9
%     pinkPts = [linspace(0.55, 0.85, 6) 0.05  0.2 0.4]';
    pinkPts = [linspace(0.55, 0.707, 3)  linspace(0.77, 0.92, 3), 0.05  0.15 0.35 0.45]';

  case 10 
    pinkPts = [linspace(0.55, 0.707, 3)  linspace(0.77, 0.92, 3), 0.05  0.15 0.35 0.45]';

  otherwise 
    error('unknown number of points !');

end
pinkVals = f(pinkPts);
pinkEstimates = localPolyKRegressionCV(th, pinkPts, pinkVals, [], [order]);
pinkExp = exp(pinkEstimates);
pinkNormC = mean(pinkExp);
pinkNormC = .95*normC;
pinkExp = pinkExp/pinkNormC;
% log joint 
figure; 
plot(th, f(th), 'k--', 'LineWidth', LW); hold on,
plot(th, pinkEstimates, 'm-', 'LineWidth', LW);
plot(pinkPts, pinkVals, 'kx', 'MarkerSize', MS, 'LineWidth', LW);
set(gca, 'Ytick', []);
set(gca, 'Xtick', []);
title('$\log P(\theta,X_{obs})$', 'Interpreter', 'Latex', 'FontSize', FS);
xlabel('$\Theta$', 'Position', [0.5 -6.5], 'FontSize', FS, 'Interpreter','Latex');
axis([0 1 -6 14]);
set(gca,'position',normalisedPosn,'units','normalized');
% Posterior
figure;
plot(th, postf(th), 'k--', 'LineWidth', LW); hold on,
plot(th, pinkExp, 'm-', 'LineWidth', LW);
set(gca, 'Ytick', []);
set(gca, 'Xtick', []);
title('$P(\theta|X_{obs})$', 'Interpreter', 'Latex', 'FontSize', FS);
xlabel('$\Theta$', 'Position', [0.5 -0.2], 'FontSize', FS, 'Interpreter','Latex');
axis([0 1 0 7.5]);
set(gca,'position',normalisedPosn,'units','normalized');

