function [lambda] = calcAdvection(T,u,v,xy,dx,rho_i,C_p)
% calcAdvection(T,u,v,xy,dx,rho_i,C_p) 
% Calculates advection term lambda, returns interpolation function

%% Make Grids
    xi = min(xy(:,1))-dx:dx:max(xy(:,1))+dx;
    yi = (min(xy(:,2))-dx:dx:max(xy(:,2))+dx)';
    [Xi,Yi] = ndgrid(xi,yi);

    us = scatteredInterpolant(xy(:,1),xy(:,2),u,'linear','nearest');
    vs = scatteredInterpolant(xy(:,1),xy(:,2),v,'linear','nearest');
    Ts = scatteredInterpolant(xy(:,1),xy(:,2),T,'linear','linear');
    
%% grid to find derivatives
    ug = us(Xi,Yi); % [m/s]
    vg = vs(Xi,Yi); % [m/s]
    Tg = Ts(Xi,Yi); % [m/s]
    
%% Smooth
    ug = imgaussfilt(ug,2);
    vg = imgaussfilt(vg,2);
    Tg = imgaussfilt(Tg,5);
    
%% Calc Lambda
    [Tx, Ty] = gradient(Tg,dx,dx);
    lam = rho_i*C_p*(Tx.*ug + Ty.*vg)./2; %divide by 2 for dept int factor
    lambda = griddedInterpolant(Xi,Yi,lam,'linear','nearest');
    