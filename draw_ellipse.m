function draw_ellipse(p,S)
% draw_ellipse(p)
%   ������Բ����p������Բ��
%
% p�ǳ���Ϊ5��ʸ����p(1)��p(2)�ֱ�Ϊxc��yc��p(3)��p(4)
% �ֱ�Ϊ����a�Ͷ���b��p(5)Ϊ���ᳯ��ǣ��Ի���Ϊ��λ��
%
% draw_ellipse���ڵ�ǰ�����л�����Բ������Ҫ���ӣ�������
% ����hold on��
%

t = 0:0.1:360;
x = p(3)*cosd(t);
y = p(4)*sind(t);
xy = [cos(p(5)) -sin(p(5)); sin(p(5)) cos(p(5))]*[x;y];
plot(xy(1,:)+p(1), xy(2,:)+p(2), S);