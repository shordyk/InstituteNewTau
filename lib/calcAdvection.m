function [lambda] = calcAdvection(T,u,v,xy,dx,rho_i,C_p)
% calcAdvection(T,u,v,xy,dx,rho_i,C_p) 
% Calculates advection term lambda, returns interpolation function

%% Make Grids
% *s is a scattered variable, *g is a gridded variable, *x or *y is a derivative
% of the variable (except for dx which is grid sizing)
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
    ug = imgaussfilt(ug,5e3/dx);
    vg = imgaussfilt(vg,5e3/dx);
    Tg = imgaussfilt(Tg,10e3/dx);
    
%% Calc Lambda
    [Tx, Ty] = gradient(Tg,dx,dx);
    lam = rho_i*C_p*(Tx.*ug + Ty.*vg)./2; 
    % Divide by 2 assumes bed is temperate, and so the depth integrated
    % temperature gradient is 1/2 that of the surface temp gradient.
    lambda = griddedInterpolant(Xi,Yi,lam,'linear','nearest');
    