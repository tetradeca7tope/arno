function p = ComputeP(th1, th2)
% function for computing p

%   % Creates a curvy posterior
%   p = (th1).^2 .* (th2).^2 + ...
%       (1 - th1).^2 .* (th2).^2 + ...
%       (1 - th1).^2 .* (1 - th2).^2;
 
  % Creates a bimodal posterior with modes close to (0,0), and (1,1);
  p = (th1).^2 .* (th2).^2 + ...
      (1 - th1).^2 .* (1 - th2).^2;

%   % Creates a posterior with modes around centres;
%   s = 0.08;
%   centres = [0.2, 0.4; 0.8, 0.2; 0.7, 0.8];
%   x = linspace(0,1)'; y = x; [X, Y] = meshgrid(x,y);
%   Z = zeros(size(X));
%   p = 0;
%   for i = 1:size(centres,1)
%     Z = Z + exp( -0.5/(s^2) * ((X - centres(i,1)).^2 + (Y - centres(i,2)).^2));
%     p = p + exp( -0.5/(s^2) * ((th1 - centres(i,1)).^2 + ...
%                                (th2 - centres(i,2)).^2));
%   end
%   norm_constant = max(max(Z));
%   p = p / norm_constant;
  
end
