initall();
%测试椭圆的绘制
global currentCurve;
% 定义长轴和短轴的一半
a =2; % 长轴的一半
b =0.5; % 短轴的一半

% 定义旋转角度
theta =pi/4; % 旋转角度

% 构造对角矩阵D
D = [1/a^2, 0; 0, 1/b^2];

% 构造旋转矩阵R
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];

% 计算矩阵P
P = R * D * R';

Eclipse=struct('xc',[2;4],'P',P,'a',a,'b',b,'theta',theta,'R',R);
%接下来的例子检测绘制曲线函数
curve=struct('p',[0,3;1,2;3,3;5,5]);
currentCurve=curve;
Obs={Eclipse};
NearbyCur={};
Nowcur=curve;
opts=OptimizeOver(Obs,NearbyCur,Nowcur);

