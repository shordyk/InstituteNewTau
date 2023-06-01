clear; clc; close all

x0 = [60e3,120e3];
str = "Depth_center";

fun = @(x)modelOpt(x,str);
options = optimset('Display','iter','PlotFcns',@optimplotfval,'TolX',1e3    ,'TolFun',1e6);	
[x,fval,exitflag,output] = fminsearch(fun,x0,options) 

save("optOutput_"+ str); 