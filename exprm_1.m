function exprm_1()

ep = [0 0 1 .5 pi/6];
p = [0 0
    -0.3 0
    0.8 0
    0 -0.6
    -0.5 sqrt(3)/4
    0.6 -0.6
    0.01 0.02
    ];
p = p*[cos(pi/6) sin(pi/6);-sin(pi/6) cos(pi/6)];
x = p(:,1);
y = p(:,2);

[xe1,ye1,d1] = ellipse_orthogonal_dist_arw01(x,y,ep,1e-3);
[xe2,ye2,d2] = ellipse_orthogonal_dist_exact(x,y,ep);
[xe3,ye3,d3] = ellipse_orthogonal_dist_ci(x,y,ep,1e-6);
fprintf('ARW01:\n');
fprintf('  Xe = %f, Ye = %f, d = %f\n', [xe1,ye1,d1]');
fprintf('EA:\n');
fprintf('  Xe = %f, Ye = %f, d = %f\n', [xe2,ye2,d2]');
fprintf('CA:\n');
fprintf('  Xe = %f, Ye = %f, d = %f\n', [xe3,ye3,d3]');

figure, draw_ellipse(ep,'b-'), hold on;
plot(x,y,'k.');
plot(xe3,ye3,'ro');
plot([x';xe3'],[y';ye3'],'g--');
axis equal;