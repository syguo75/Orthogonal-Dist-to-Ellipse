function [xp,yp,d] = ellipse_orthogonal_dist_ci(x,y,ep,tol)
% [xp,yp,dist] = ellipse_orthogonal_dist_ci(x,y,ep,tol)
%   Computes using the convergent iterative algorithm the orthogonal 
% distances from points (x,y) to the ellipse ep, and the orthogonal 
% contacting points of (x,y) on ep.
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
% the nonlinear equation regarding s = sin(phi). Note that the tolerance
% is given to directly control the accuracy of xp and yp, and the
% MEX function ci_mex will convert it into the tolerance of s.
%
%   xp and yp are two column vectors giving the coordinates of the
% orthogonal contacting points, and dist is a column vector of
% the orthogonal distances.
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

sx = ones(length(x),1);
sy = sx;
sx(x<0) = -1;
sy(y<0) = -1;

sinv = ci_mex(abs(x),abs(y),a,b,tol);
cosv = sqrt(1-sinv.^2);

xp = a*cosv.*sx;
yp = b*sinv.*sy;
d = sqrt((x-xp).^2 + (y-yp).^2);

q = [xp yp]*[cos(theta) sin(theta); -sin(theta) cos(theta)];

xp = q(:,1) + xc;
yp = q(:,2) + yc;