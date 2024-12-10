function [grad,df]=findgradient(curve)
global currentCurve;
global normbase;


currentCurve=curve;



fp=@(vec) -1*dcost1(vec)-0.2*dcost2(vec);%+dcost3(vec)*vec;
grad=fp(normbase{1})*normbase{1}+fp(normbase{2})*normbase{2}+fp(normbase{3})*normbase{3}+fp(normbase{4})*normbase{4};
df=fp(grad);