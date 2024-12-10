function results=findminpoints(Eclipse,curve)



%第一步，给出初始点
thetalist = linspace(0, 2*pi, 6); % 参数theta从0到2*pi
value=zeros(1,6);
res={};
for i = 1:length(thetalist)
    results=startiter(Eclipse,curve,thetalist(i));
    res{i}=results;
    value(i)=results.d;
end
%然后找最小值
[~, minIndex] = min(value);
results=res{minIndex};
end



function res=startiter(Eclipse,curve,theta)

global phi;
P=Eclipse.P;
x_c=Eclipse.xc;
rot=Eclipse.theta;
a=Eclipse.a;
b=Eclipse.b;
iter=1;
while 1%开始循环
    %参数化获取theta对应的x坐标
    x=x_c+[cos(rot),-sin(rot);sin(rot),cos(rot)]*[a*cos(theta);b*sin(theta)];
    %求解曲线上的对应点
    %step1求解投影矩阵PIx
    innerP=(P*(x-x_c))'*P*(x-x_c);
    PIx=eye(2)-P*(x-x_c)*(x-x_c)'*P/innerP;
    %PIx*curve.p'*phi'*U==PIx*x;%两行方程，所有方程的公共解就是落点
    A=PIx*curve.p'*phi';
    B=PIx*x;
    % 定义多项式系数
    p3= [A(1,4),A(1,3), A(1,2), A(1,1)-B(1)];
    % 求解多项式的根
    r = roots(p3);
    %检验每一个根：1.非虚数。3.01之间
    rvalid=[];
    % 遍历根的向量
    k=1;
    for i = 1:length(r)
        % 检测当前根是否是实数
        if ~isreal(r(i))
            % 如果不是实数，那么它是虚数
            %             fprintf('根 %d 是虚数：%.2f + %.2fi\n', i, real(r(i)), imag(r(i)));
        else
            % 如果是实数，那么它不是虚数
            %继续检测有效性
            if(r(i)<0)
                r(i)=0;
            end
            if(r(i)>1)
                r(i)=1;
            end
            rvalid(k)=r(i);
            k=k+1;
        end
    end
    %寻找比较小的那一个作为根，如果最小根在端点上，说明超限了，此时使用
    %另一套方法求导数
    dis=[];
    k=1;
    for i = 1:length(rvalid)
        ui=rvalid(i);
        U=[1 ui ui^2 ui^3];
        g=curve.p'*phi'*U';
        dis(k)=norm(g-x);k=k+1;
    end
    [~, minIndex] = min(dis);
    
    %接下来求解应对向量V的距离函数的导数:
    V_x=[cos(rot),-sin(rot);sin(rot),cos(rot)]*[-a*sin(theta);b*cos(theta)];
    %第一种情况，超限
    if(rvalid(minIndex)==0 || rvalid(minIndex)==1)
        %这个时候g是常数
        Ddis=-2*g'*V_x+2*x'*V_x;
    end
    %第二种情况，正常
    if(rvalid(minIndex)~=0 && rvalid(minIndex)~=1)
        ui=rvalid(minIndex);
        U=[1 ui ui^2 ui^3];
        g=curve.p'*phi'*U';
        
        dPI_x=-(P*V_x*(x-x_c)'*P+P*(x-x_c)*V_x'*P)/innerP-P*(x-x_c)*(x-x_c)'*P*(-2*(x-x_c)'*P^2*V_x)/innerP^2;
        k1=2*ui;k2=3*ui*ui;%k1 k2必然大于0
        B=V_x+dPI_x*x-dPI_x*curve.p'*phi'*U';%方程右侧
        A=PIx*curve.p'*phi';%方程左侧
        %求解方程
        du=B(1)/(A(1,2)+A(1,3)*k1+A(1,4)*k2);
        u=[0,du,k1*du,k2*du];%由椭圆切空间计算U切空间
        dg=curve.p'*phi'*u';%由U切空间计算贝塞尔切空间
        %计算距离对x导数
        Ddis=2*dg'*g-2*(dg'*x+g'*V_x)+2*x'*V_x;
    end
    %theta的运动方向定为-alpha*Ddis
    alpha=0.2;
    theta=theta-alpha*Ddis;
    %theta=theta-0.1*sign(Ddis);
    if(abs(Ddis)<0.01 || iter>30)
        res=struct('Ct',rvalid(minIndex),'Cp',g,'Ep',x,'d',norm(g-x));
        return
    end
    %
%     clf;
%     drawEclipse(Eclipse);
%     hold on;
%     drawCurve(curve);
%     hold on;
%     plot(g(1),g(2),'o');
%     hold on
%     plot(x(1),x(2),'o');
%     hold on
%     plot([x(1) g(1)],[x(2) g(2)]);
%     hold on
%     vscs=P*(x-x_c)/norm(P*(x-x_c));
%     plot([x(1) x(1)+vscs(1)],[x(2) x(2)+vscs(2)]);
%     hold on
%     plot([x(1) x(1)-vscs(2)],[x(2) x(2)+vscs(1)]);
%     hold on
%     axis equal
%     pause(0.1);
    %然后一直循环
    iter=iter+1;
    
end

end