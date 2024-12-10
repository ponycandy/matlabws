function drawEclipse(Eclipse)

x_c=Eclipse.xc;
rot=Eclipse.theta;

a=Eclipse.a;
b=Eclipse.b;
sumup=a+b;
c=max(a,b);
a=c;
b=sumup-c;

% 椭圆的参数方程
theta = linspace(0, 2*pi, 100); % 参数theta从0到2*pi
pos=zeros(2,size(theta,2));
for i=1:1:size(theta,2)
    pos(:,i)=x_c+[cos(rot),-sin(rot);sin(rot),cos(rot)]*[a*cos(theta(i));b*sin(theta(i))];
end

% 绘制椭圆
plot(pos(1,:), pos(2,:));
axis equal; % 保持x和y轴的比例一致
grid on;
title('椭圆的绘制');
xlabel('x');
ylabel('y');
