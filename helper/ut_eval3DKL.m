% unit test for 3D KL

p1 = @(t) 2 * double (2 - 2*t(:,1));
p2 = @(t) ones(size(t,1), 1);
KL = eval3DKL(p1, p2, repmat([0 1], 3, 1)),

% 2 Gaussians
p1 = @(t) mvnpdf(t);
p2 = @(t) mvnpdf(t, 1);
KL = eval3DKL(p1, p2, repmat([-4 5], 3, 1)),

