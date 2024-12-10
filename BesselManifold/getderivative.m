syms t a b c c1 c2 rot
P=[a b;b c];%椭圆P矩阵
xc=[c1;c2];%椭圆中心
%首先获取椭圆点的参数化表达式

pos=xc+[cos(rot),-sin(rot);sin(rot),cos(rot)]*[a*cos(t);b*sin(t)];
%获取垂线的表达式
%投影矩阵
dfx=(pos-xc)'*P;
PI_x=eye(2)-dfx'*(dfx*dfx')^(-1)*dfx;


% 定义符号变量
syms x

% 定义符号表达式
f = x^2 + 3*x + 2;

% 求一阶导数
df = diff(f, x);

% 求二阶导数
d2f = diff(f, x, 2);

% 显示导数表达式
disp('一阶导数：');
pretty(df)
disp('二阶导数：');
pretty(d2f)

% 计算一阶导数在 x = 1 处的值
df_at_1 = vpa(subs(df, x, 1));
disp('一阶导数在 x = 1 处的值：');
disp(df_at_1)