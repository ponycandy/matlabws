function results=findwithin(Eclipse,curve)



%第一步，给出初始点
t = linspace(0, 1, 3); % 参数theta从0到2*pi
value=zeros(1,3);
res={};
for i = 1:length(t)
    results=startiter(Eclipse,curve,t(i));
    res{i}=results;
    value(i)=results.d;
end
%然后找最小值
[~, minIndex] = min(value);
results=res{minIndex};
end



function res=startiter(Eclipse,curve,t)

global phi;
P=Eclipse.P;
x_c=Eclipse.xc;
rot=Eclipse.theta;
a=Eclipse.a;
b=Eclipse.b;
iter=1;
while 1%开始循环
    x=x_c;
    ui=t;
    U=[1 ui ui^2 ui^3];
    dU=[0 1 2*ui 3*ui^2];
    %单位化
    dU=dU/norm(dU);
    
    g=curve.p'*phi'*U';
    dfdt=2*(g-x_c)'*P*curve.p'*phi'*dU';
    %注意，这个值并不是
    %到圆心的最大值，是需要按照椭圆扁度进行扭曲变形的
    %所以动画没有错，这就是目标值
    
    %切一下梯度:
    cliff=10;
    if(abs(dfdt)>cliff)
        dfdt=sign(dfdt)*cliff;
    end
    alpha=0.002;
    t=t-alpha*dfdt;
    if(abs(dfdt)<0.01 || t<0 || t>1 || iter>30)
        if(t<0)
            t=0;
        end
        if(t>1)
            t=1;
        end
        res=struct('in',(g-x_c)'*P*(g-x_c)<1,'Ct',t,'Cp',g,'Ep',x,'d',(g-x_c)'*P*(g-x_c));
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