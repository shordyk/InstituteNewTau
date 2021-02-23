clear
clc

addpath('lib') 

% %Institute
% overgrab = 0;
% xmax = -7e5;
% xmin = -10e5;
% ymax =  5e5;
% ymin =  1e5;
% xx = -9.5e5:dx*26/42:-9e5;
% yy = 2.2e5:dx:3e5;

%Thwaites
overgrab = 0;
xmin = -15.3e5;
xmax = -13.5e5;
ymin =  -5.6e5;
ymax =  -2.5e5;
start = [-14.5e5,-3.2e5];
stop = [-14.4e5,-4e5];
dist = stop-start;
xx = start(1):dist(1)/30:stop(1);
yy = start(2):dist(2)/30:stop(2);

dx = 1e3;
smth = 10e3;
xi = xmin-dx*overgrab:dx:xmax+dx*overgrab;
yi = ymin-dx*overgrab:dx:ymax+dx*overgrab;
[Xi,Yi] = ndgrid(xi,yi);
xy = [Xi(1:end)',Yi(1:end)'];

%% Physical parameters
a     = 2e-26^(-1/3);       % a:     flow parameter pre-factor [Pa s^1/3] @-35C from cuffey 
nn    = 3;                  % Glens law power
p     = 4/3;                % p:     flow parameter power [ ]
g     = 9.81;               % g:     acceleration due to gravity [m/s^2]
nu    = .4;                 % Thermal relaxation parameter [ ]
rho   = 917;                % rho:   density of ice [kg/m^3]
rho_w = 1000;               % rho_w: density of water [kg/m^3]
C_p   = 2050;               % specific heat of ice [J/Kg/K]
K     = 2.1;                % thermal conductivity of ice [W/m/K]
A_m   = 2.4e-24;            % Meyer's prefactor [Pa^-3 s^-1]
T_m   = 273;                % Ice melting point [k] 


dz = .1;  %vertical resolution of thermal depth profiles (frac of H) [ ]

%% Import data and smooth 
[Acc, T_s] = loadALBMAP();
T = T_s(xy(:,1),xy(:,2));
bm_b =  bedmachine_interp('bed',Xi,Yi);
bm_s =  bedmachine_interp('surface',Xi,Yi);
[u,v] = measures_interp('velocity',xy(:,1),xy(:,2)); %[m/yr]
[spd] = measures_interp('speed',Xi,Yi)'; %[m/yr]
u = u/3.154e7;
v = v/3.154e7;
%% Smoothing
% Numerator is the window we're smoothing over in [m], spacing of these grids

smoothbed = imgaussfilt(bm_b,smth/dx);
smoothsurf = imgaussfilt(bm_s,smth/dx);
% rectify rock above ice issue, force that ice is non-zero thickness everywhere
smoothsurf(smoothbed > smoothsurf) = smoothbed(smoothbed > smoothsurf) + 1; %Pe and Br = 0 result in NaNs

h_b_init = griddedInterpolant(Xi,Yi,smoothbed);
h_s_init = griddedInterpolant(Xi,Yi,smoothsurf);
h_init = griddedInterpolant(Xi,Yi,smoothsurf-smoothbed);

%% Meyer Model
        ep_dot = calcTrigridStrain(u,v,xy,dx); %returns intperolation object
        lambda  = calcAdvection(T,u,v,xy,dx/4,rho,C_p); %TODO better derivatives, analytic?

        % Brinkman number [ ]
        Br =@(x,y) 2*subplus(h_s_init(x,y)-h_b_init(x,y)).^2./(K*(T_m-T_s(x,y))).*((subplus(ep_dot(x,y)).^(nn+1))/A_m).^(1/nn);

        % Peclet number  [ ]
        Pe =@(x,y) rho*C_p.*Acc(x,y).*subplus(h_s_init(x,y)-h_b_init(x,y))./K;

        % Vertical Peclet number  [ ]
%         La =@(x,y) lambda(x,y).*subplus(h_s_init(x,y)-h_b_init(x,y)).^2./(K*(T_m-T_s(x,y)));
        La =@(x,y) zeros(size(x)); %no advect for now
        % Critical Strain [s^-1]
        ep_star =@(x,y) ((La(x,y)/2 + ((Pe(x,y).^2)/2)./(Pe(x,y)-1+exp(-Pe(x,y))))).^(nn/(nn+1))...
        .*(K*(T_m-T_s(x,y))./(A_m.^(-1/nn).*(subplus(h_s_init(x,y)-h_b_init(x,y))).^2)).^(nn/(nn+1));

        % Temp profile at xy [K]
        t_z =@(x,y) tempProfile(ep_dot(x,y),ep_star(x,y),Pe(x,y),Br(x,y),La(x,y),T_s(x,y),T_m,dz); 

        % Enhancement Factor []
        E_t =@(x,y) depthIntEnhancement(t_z(x,y),a.^(-3),dz);

        % Mean Temp [K]
        T_bar = @(x,y) trapz(t_z(x,y),2)*dz;

%% Plotting, allow cross section        

        
% Plotting temp depth averaged
figure(1)
clf
p = surf(reshape(xy(:,1),size(Xi)),reshape(xy(:,2),size(Xi)),reshape(T_bar(xy(:,1),xy(:,2))-273.15,size(Xi)));
hold on
plot(xx,yy,'r*-','linewidth',2)
contour(xi,yi,spd, [10, 10] , 'k:','HandleVisibility','off')
contour(xi,yi,spd, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd, [1000, 1000] , 'k-','LineWidth',2)
title('Temp Avg')
set(p, 'edgecolor', 'none');
colorbar
caxis([-40 0])
view(2)

figure(2)
clf
tempM = t_z(xx',yy')'-273.15;
height = zeros(size(tempM));
x_along = zeros(size(tempM));
for i = 1:length(xx)
    height(:,i) = h_b_init(xx(i),yy(i)) + (0:dz:1)'*h_init(xx(i),yy(i));
    x_along(:,i) = sqrt((xx(1) - xx(i))^2 + (yy(1) - yy(i))^2);
end   
p = surf(x_along, height, tempM);
set(p, 'edgecolor', 'none');
view(2)
colorbar;
% caxis([-40 0])