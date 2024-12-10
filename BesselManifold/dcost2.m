function df=dcost2(vec)
global local_vec;
local_vec=vec;

% 定义被积函数
f = @(x) pathlength(x);

% 计算0-1上的定积分
result = integral(f, 0, 1);

% 显示结果
df=result;

end

function val=pathlength(x)
global local_vec;
global currentCurve;
global phi;

vec=local_vec;
curve=currentCurve;
val=x;

for i=1:size(x,2)
    dU=[0 1 2*x(i) 3*x(i)^2];
    val(i)=dU*phi*(vec*curve.p'+curve.p*vec')*phi'*dU'/(2*(dU*phi*curve.p*curve.p'*phi'*dU')^0.5);
end

end