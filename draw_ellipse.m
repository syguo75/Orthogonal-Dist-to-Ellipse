function draw_ellipse(p,S)
% draw_ellipse(p)
%   根据椭圆参数p绘制椭圆。
%
% p是长度为5的矢量，p(1)和p(2)分别为xc和yc；p(3)和p(4)
% 分别为长轴a和短轴b；p(5)为长轴朝向角，以弧度为单位。
%
% draw_ellipse将在当前窗口中绘制椭圆。若需要叠加，请事先
% 设置hold on。
%

t = 0:0.1:360;
x = p(3)*cosd(t);
y = p(4)*sind(t);
xy = [cos(p(5)) -sin(p(5)); sin(p(5)) cos(p(5))]*[x;y];
plot(xy(1,:)+p(1), xy(2,:)+p(2), S);