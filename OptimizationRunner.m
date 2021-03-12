clear; clc; close all

x0 = [80, 150e3];
str = "Case 2";

fun = @(x)modelOpt(x,str);
options = optimset('Display','iter','PlotFcns',@optimplotfval,'TolX',1e2    ,'TolFun',1e7);	
[x,fval,exitflag,output] = fminsearch(fun,x0,options) 

save("optOutput_"+ str,'x','str'); 