
addpath ../v9_AL_KLmin_GP_2d/

% close all;

NUM_PTS = 50;
NUM_INIT = 25;

pts_chosen = curr_pts;
% pts_chosen = uc_pts;

RES = 100;
t = linspace(0,1,RES);
[T1, T2] = meshgrid(t,t);
th = [T1(:), T2(:)];

[~, ~, true_post] = sampleExp2(th(:,1), th(:,2));
TruePost = reshape(true_post, RES, RES);

% First the initial points
init_pts = pts_chosen(1:NUM_INIT, :);
plot(init_pts(:,1), init_pts(:,2), 'bx');
hold on,

contour(T1, T2, TruePost, 'Color', 'g');

for i = 1:NUM_PTS
  idx = NUM_INIT + i;
  text(pts_chosen(idx,1), pts_chosen(idx, 2), num2str(i), 'Color', 'k');
end

axis([0,1,0,1]);
% axis square;
