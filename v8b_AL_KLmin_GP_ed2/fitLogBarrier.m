function [a, b, c] = fitLogBarrier(x, y)

% fits a log barrier function to the given data
% i.e. y ~= alog(x) + blog(1-x) + c

n = size(x,1);

P = [log(x) log(1-x) ones(n,1)];

% create the following convex program
cvx_begin quiet
variables params(3);
minimize( norm(P*params - y,2) );
subject to
  params(1) >= 1;
  params(2) >= 1;
cvx_end

a = params(1);
b = params(2);
c = params(3);

end

