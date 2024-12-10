function curve=OptimizeOver(Obs,NearbyCur,Nowcur)

global currentCurve;
global phi;
global ObsList;
global Nearby;
global cost_1_var_flag;
global cost_2_var_flag;

global cost_1_varchache;
global cost_2_varchache;

cost_1_var_flag=0;
cost_2_var_flag=0;
iterationnum=10;

currentCurve=Nowcur;
ObsList=Obs;
Nearby=NearbyCur;
alpha=0.1;
dfdt2=1;
cache=0;
while 1
    [grad,dfdt]=findgradient(Nowcur);
    
    cliff=10;
    if(abs(dfdt)>cliff)
        dfdt=sign(dfdt)*cliff;
    end
    
    %更新，在这里需要捕捉目标曲线的局部循环状态
    %或者，随着捕获的渐进，逐渐减小参数alpha
    Nowcur.p=Nowcur.p+alpha*grad;%使用一下数值二阶导
    currentCurve=Nowcur;
    clf
    drawCurve(currentCurve);
    hold on
    drawEclipse(ObsList{1});
    axis equal
%     gd=dcost2(grad);
  %  uud=norm(grad-cache);
    dfdt2=(dfdt-cache)*0.2;
   disp([dfdt2,dfdt])
    
    if(abs(dfdt)<0.001 )%优化问题实际上是一个3玩家Nash均衡问题，所以终止条件捕获为三玩家的Nash均衡即可
        curve=struct('p',Nowcur.p);
        break
    end
    pause(0.1);
    cache=dfdt;
end