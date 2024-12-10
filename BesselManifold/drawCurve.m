function drawCurve(curve)

global phi;

% 椭圆的参数方程
t = linspace(0, 1, 50); % 参数theta从0到2*pi
pos=zeros(2,size(t,2));
for i=1:1:size(t,2)
    U=[1 t(i) t(i)^2 t(i)^3];
    pos(:,i)=curve.p'*phi'*U';
end

% 绘制椭圆
plot(pos(1,:), pos(2,:));
axis equal; % 保持x和y轴的比例一致
grid on;
title('椭圆的绘制');
xlabel('x');
ylabel('y');
