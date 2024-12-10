function df=dcost1(vec)
%计算第一个损失函数对决策变量P的导数
%只计算一层
global phi;
global ObsList;
global currentCurve;

global cost_1_var_flag;
global cost_2_var_flag;

global cost_1_varchache;
global cost_2_varchache;

% global nearbyCurveList;
%计算最小值
%首先判断曲线与椭圆的相交性，相交的时候
%使用简单代价函数

df=0;
for i = 1:length(ObsList)
    Eclipse=ObsList{i};
    
    results=findwithin(Eclipse,currentCurve);
    %在椭圆外面用圆滚一圈，得到新的椭圆，设圆半径为r
    %新椭圆的长轴变为a+r，短轴变为b+r
    a=Eclipse.a;
    b=Eclipse.b;
    r=0.5;
    K=[(a/(a+r))^2,0;0,(b/(b+r))^2];
    P_new=Eclipse.P*Eclipse.R*K*Eclipse.R';
    
    u=results.Ct;
    U=[1 u u^2 u^3];
    x=Eclipse.xc;%最大化到圆心的距离    
    g=currentCurve.p'*phi'*U';
    
    if((g-x)'*P_new*(g-x)<1)
        if(results.in)
            u=results.Ct;
            p=currentCurve.p;
            U=[1 u u^2 u^3];
            x=Eclipse.xc;%最大化到圆心的距离
        else
            results=findminpoints(Eclipse,currentCurve);
            u=results.Ct;
            p=currentCurve.p;
            U=[1 u u^2 u^3];
            x=results.Ep;%最大化到最近点的距离
        end
        df=df-(U*phi*(p*vec'+vec*p')*phi'*U'-2*U*phi*vec*x);
    else
        df=df+0;
    end
end
return