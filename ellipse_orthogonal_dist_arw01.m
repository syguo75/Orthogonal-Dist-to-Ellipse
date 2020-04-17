function [xp,yp,d] = ellipse_orthogonal_dist_arw01(x,y,ep,tol)
% [xp,yp,dist] = ellipse_orthogonal_dist_arw01(x,y,ep)
%   Computes using the algorithm in (Ahn 2001) the orthogonal distances 
% from points (x,y) to the ellipse ep, and the orthogonal contacting 
% points of (x,y) on ep.
%
%   x and y are two same-sized arrays of the coordinates 
% of the points under concern, preferably two column vectors.
%
%   ep is a vector with a length of 5, which gives the
% parameters of the ellipse under concern in the form of
% (Xc, Yc, a, b, theta), where (Xc,Yc) being the center of
% the ellipse, a and b the major and minor axis length respectively,
% and theta the orientation of the major axis (in rads).
%
%   tol is the tolerance for the Newton's method while solving 
% the nonlinear equations regarding xp and yp. It actually has
% the pixel unit, and can be thus set according to the accuracy
% requirement.
%
%   xp and yp are two column vectors giving the coordinates of the
% orthogonal contacting points, and dist is a column vector of
% the orthogonal distances.
%
%   The iteration of equation solving in arw01_mex may diverge. In
% arw01_mex, a maximum of 50 iterations is set. If the iteration
% doesn't produce a converged result when this limit is reached,
% this iteration is then treated as a divergent one, and the results
% are given as NaNs.
%

xc = ep(1);
yc = ep(2);
a = ep(3);
b = ep(4);
theta = ep(5);

x = x(:) - xc;
y = y(:) - yc;
q = [x y] * [cos(theta) -sin(theta); sin(theta) cos(theta)];
x = q(:,1);
y = q(:,2);

[xe,ye] = arw01_mex(x,y,a,b,tol);

d = sqrt((x-xe).^2 + (y-ye).^2);

q = [xe ye]*[cos(theta) sin(theta); -sin(theta) cos(theta)];
xp = q(:,1) + xc;
yp = q(:,2) + yc;
