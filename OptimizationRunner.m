clear; clc; close all

x0 = [1,1];
str = "Mixed";

fun = @(x)modelOpt(x,str);
options = optimset('Display','iter','PlotFcns',@optimplotfval,'TolX',1e-3    ,'TolFun',1e6);	
[x,fval,exitflag,output] = fminsearch(fun,x0,options) 

save("optOutput_"+ str,'x','str'); 