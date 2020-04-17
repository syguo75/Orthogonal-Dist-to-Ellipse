function [xp,yp,dist] = ellipse_orthogonal_dist_exact(x,y,ep)
% [xp,yp,dist] = ellipse_orthogonal_dist_exact(x,y,ep)
%   Computes using the exact algorithm the orthogonal distances 
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
%   xp and yp are two column vectors giving the coordinates of the
% orthogonal contacting points, and dist is a column vector of
% the orthogonal distances.
%

a = ep(3);
b = ep(4);
d = a^2-b^2;

% Transform the points to the intrinsic coordinate system,
% i.e., Xc = Yc = 0, theta = 0.
x = x(:) - ep(1);
y = y(:) - ep(2);
q = [x y] * [cos(ep(5)) -sin(ep(5)); sin(ep(5)) cos(ep(5))];
x = q(:,1);
y = q(:,2);

% Reflect the points into the first quadrant.
xq = abs(x);
yq = abs(y);

a3 = 2*yq*b/d;
a2 = (a^2*xq.^2 + b^2*yq.^2)/d^2 - 1;
a1 = -a3;
a0 = -b^2*yq.^2/d^2;

% Solve the quartic equation. If no valid root is found, set xp,
% yp and dist as Inf.
sinv = real(solve_quartic_eq(a0,a1,a2,a3));
cosv = sqrt(1-sinv.^2);
i = sinv<0 | sinv>1;
sinv(i) = inf;
cosv(i) = inf;
xp = a*cosv;
yp = b*sinv;

% Find the root that minimizes the distance.
dist2 = (repmat(xq,1,size(xp,2))-xp).^2 + (repmat(yq,1,size(yp,2))-yp).^2;
[~,j] = min(dist2,[],2);
ind = sub2ind(size(dist2),(1:size(dist2,1))',j);
sinv = sinv(ind);

% Deal with the degrading cases.
i = yq==0 & xq>=d/a;
sinv(i) = 0;
i = yq==0 & xq<d/a;
sinv(i) = sqrt(1 - a^2*xq(i).^2/d^2);

% Find the orthogonal contacting points in the first quadrant.
cosv = sqrt(1-sinv.^2);
xp = a*cosv;
yp = b*sinv;
dist = sqrt((xq-xp).^2 + (yq-yp).^2);

% Reflect the orthogonal contacting points back to the original
% quadrant in the intrinsic system.
i = x<0;
xp(i) = -xp(i);
i = y<0;
yp(i) = -yp(i);

% Transform the orthogonal contacting point back to the original
% coordinate system.
q = [xp yp]*[cos(ep(5)) sin(ep(5)); -sin(ep(5)) cos(ep(5))];

xp = q(:,1) + ep(1);
yp = q(:,2) + ep(2);

function s = solve_quartic_eq(a0,a1,a2,a3)

a0 = a0(:);
a1 = a1(:);
a2 = a2(:);
a3 = a3(:);

b2 = -a2;
b1 = a1.*a3 - 4*a0;
b0 = 4*a2.*a0 - a1.^2 -a3.^2.*a0;

y = solve_cubic_eq(b0,b1,b2);

a0 = repmat(a0, 1, size(y,2));
a1 = repmat(a1, 1, size(y,2));
a2 = repmat(a2, 1, size(y,2));
a3 = repmat(a3, 1, size(y,2));

R = sqrt(a3.^2/4 - a2 + y);
S = zeros(size(R));
T = zeros(size(R));

i = R==0;
if any(any(i))
    S(i) = 3*a3(i).^2/4 - 2*a2(i);
    T(i) = 2*sqrt(y(i).^2 - 4*a0(i));
end

i = ~i;
if any(any(i))
    S(i) = 3*a3(i).^2/4 - R(i).^2 - 2*a2(i);
    T(i) = (4*a3(i).*a2(i) - 8*a1(i) - a3(i).^3)./(4*R(i));
end

D = sqrt(S+T);
E = sqrt(S-T);
    
s = [-a3/2+R+D -a3/2+R-D -a3/2-R+E -a3/2-R-E]/2;

function z = solve_cubic_eq(a0,a1,a2)

a0 = a0(:);
a1 = a1(:);
a2 = a2(:);

Q = (3*a1 - a2.^2)/9;
R = (9*a2.*a1 - 27*a0 - 2*a2.^3)/54;
D = Q.^3 + R.^2;
SD = sqrt(D);
S = zeros(size(R));
T = zeros(size(R));

i = find(D>=0);
if ~isempty(i)
    S(i) = nthroot(R(i) + SD(i), 3);
    T(i) = nthroot(R(i) - SD(i), 3);
end

i = find(D<0);
if ~isempty(i)
    S(i) = (R(i) + SD(i)).^(1/3);
    T(i) = (R(i) - SD(i)).^(1/3);
end
z = [-a2/3+(S+T)/2 -a2/3-(S+T)/2+sqrt(-3)*(S-T)/2 -a2/3-(S+T)/2-sqrt(-3)*(S-T)/2];

% function z = real_root_of_cubic_eq(a0,a1,a2)
% 
% a0 = a0(:);
% a1 = a1(:);
% a2 = a2(:);
% 
% z = zeros(length(a0),1);
% 
% Q = (3*a1 - a2.^2)/9;
% R = (9*a2.*a1 - 27*a0 - 2*a2.^3)/54;
% D = Q.^3 + R.^2;
% SD = sqrt(D);
% 
% i = find(D>=0);
% if ~isempty(i)
%     S = nthroot(R(i) + SD(i), 3);
%     T = nthroot(R(i) - SD(i), 3);
%     z(i) = -a2(i)/3 + S + T;
% end
% 
% i = find(D<0);
% if ~isempty(i)
%     S = (R(i) + SD(i)).^(1/3);
%     T = (R(i) - SD(i)).^(1/3);
%     tz = [-a2(i)/3+(S+T)/2 -a2(i)/3-(S+T)/2+sqrt(-3)*(S-T)/2 -a2(i)/3-(S+T)/2-sqrt(-3)*(S-T)/2];
%     c = abs(cos(mod(angle(tz),pi)));
%     [~,j] = max(c,[],2);
%     ind = sub2ind(size(tz),(1:length(S))',j);
% %     [~,j] = max(abs(tz));
% %     z(i) = real(tz(j))
%     z(i) = real(tz(ind));
% end