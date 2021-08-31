clear; clc; close all

x0 = [1.2];
str = "ISSM_center";

fun = @(x)modelOpt(x,str);
options = optimset('Display','iter','PlotFcns',@optimplotfval,'TolX',1e-2    ,'TolFun',1e6);	
[x,fval,exitflag,output] = fminsearch(fun,x0,options) 

save("optOutput_"+ str,'x','str'); 