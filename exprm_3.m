function exprm_3()

rng('default');

center = (rand(100,2)-0.5)*200;
a = rand(100,1)*49+1;
b = (rand(100,1)*0.8+0.1).*a;
phi = (rand(100,1)-0.5)*pi;
np = randi([100 1000],100,1);
d = (rand(100,1)*0.09) + 0.01;

ep = zeros(100,5);
dataset = cell(100,1);
for i = 1:100
    ep(i,:) = [center(i,:) a(i) b(i) phi(i)];
    theta = randi([0 365],np(i),1);
    p = [cosd(theta)*a(i) sind(theta)*b(i)] + (rand(np(i),2)-0.5)*2*d(i);
    p = p*[cos(phi(i)) sin(phi(i));-sin(phi(i)) cos(phi(i))];
    dataset{i} = [p(:,1)+center(i,1) p(:,2)+center(i,2)];
end

d1 = [];
d2 = [];
d3 = [];
profile on;
for i = 1:100
    [~,~,cd1] = ellipse_orthogonal_dist_arw01(dataset{i}(:,1),dataset{i}(:,2),ep(i,:),1e-3);
    [~,~,cd2] = ellipse_orthogonal_dist_exact(dataset{i}(:,1),dataset{i}(:,2),ep(i,:));
    [~,~,cd3] = ellipse_orthogonal_dist_ci(dataset{i}(:,1),dataset{i}(:,2),ep(i,:),1e-6);
    d1 = [d1;cd1];
    d2 = [d2;cd2];
    d3 = [d3;cd3];
end
profile report;

fprintf('Failure statistics:\n');
fprintf('  ARW01: # failures = %d, rAF = %0.2f\n', sum(isnan(d1)), sum(isnan(d1))/length(d1));
fprintf('  EA   : # failures = %d, rAF = %0.2f\n', sum(isinf(d2)), sum(isinf(d2))/length(d2));
fprintf('\n');

fprintf('Detection statistics:\n');
md = min([d1 d2 d3],[],2);
fprintf('  ARW01: rSD = %0.2f\n',sum(d1==md)/length(md));
fprintf('  EA   : rSD = %0.2f\n',sum(d2==md)/length(md));
fprintf('  CA   : rSD = %0.2f\n',sum(d3==md)/length(md));
fprintf('\n');

fprintf('Distance overhead:\n');
i = d1>md;
fprintf('  ARW01: dMO = %0.1e\n', max(d1(i)-md(i)));
i = d2>md;
fprintf('  EA   : dMO = %0.1e\n', max(d2(i)-md(i)));
i = d3>md;
fprintf('  CA   : dMO = %0.1e\n', max(d3(i)-md(i)));