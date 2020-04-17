function exprm_2()

ep = [0 0 1 .5 pi/6];
[x,y] = meshgrid(-10:0.01:10);
p = [x(:) y(:)]*[cos(pi/6) sin(pi/6);-sin(pi/6) cos(pi/6)];
x = p(:,1);
y = p(:,2);

profile on;
[~,~,d1] = ellipse_orthogonal_dist_arw01(x,y,ep,1e-3);
[~,~,d2] = ellipse_orthogonal_dist_exact(x,y,ep);
[~,~,d3] = ellipse_orthogonal_dist_ci(x,y,ep,1e-10);
profile report;

fprintf('Failure statistics:\n');
fprintf('  ARW01: #failures = %d. rAF = %0.2f\n', sum(isnan(d1)), sum(isnan(d1))/length(d1));
fprintf('  EA:    #failures = %d. rAF = %0.2f\n', sum(isinf(d2)), sum(isinf(d2))/length(d2));
fprintf('\n');

md = min([d1 d2 d3],[],2);
fprintf('Detection statistics:\n');
fprintf('  ARW01: rSD = %0.2f\n', sum(d1==md)/length(md));
fprintf('  EA   : rSD = %0.2f\n', sum(d2==md)/length(md));
fprintf('  CA   : rSD = %0.2f\n', sum(d3==md)/length(md));
fprintf('\n');

fprintf('Distance Overhead:\n');
i = d1>md;
fprintf('  ARW01: dMO = %0.1e\n', max(d1(i)-md(i)));
i = d2>md;
fprintf('  EA   : dMO = %0.1e\n', max(d2(i)-md(i)));
i = d3>md;
fprintf('  CA   : dMO = %0.1e\n', max(d3(i)-md(i)));