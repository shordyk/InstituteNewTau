% clear
% clc
% close all

addpath('lib') 


%For files must download layer data and data (mvdr) from https://data.cresis.ku.edu/data/rds/
% file = "";
% %Institute
% xmax = -8.7e5;
% xmin = -9.8e5;
% ymax =  3.3e5;
% ymin =  2e5;
% start = [-9.5e5,2.2e5];
% stop = [-9e5,3e5];


%Thwaites
% xmin = -15.3e5;
% xmax = -13.5e5;
% ymin =  -5.6e5;
% ymax =  -2.5e5;
% start = [-14.5e5,-3.2e5];
% stop = [-14.4e5,-4e5];

%Dragon Margin 
% xmin = -8e5;
% xmax = -5e5;
% ymin =  -10e5;
% ymax =  -7e5;
% % start = [-6.4e5,-8.2e5]; %Dragon Margin
% % stop = [-6.6e5,-8.6e5];
% start = [-5.6e5,-8.4e5]; %East Margin of B1 (10 dB?)
% stop = [-5.3e5,-8.2e5];

% Siple Coast IRMCR1B_20131126_01_031
% xmin = -8e5;
% xmax = -5e5;
% ymin =  -9.75e5;
% ymax =  -6.75e5;
% start = [-6.55e5,-7.48e5]; 
% stop = [-6.42e5,-7.02e5];
% file = "radarData/53908034/Data_20131126_01_031";

% Siple Coast IRMCR1B_20131126_01_044
% xmin = -8e5;
% xmax = -5e5;
% ymin =  -8.75e5;
% ymax =  -5.75e5;
% start = [-6.55e5,-7.45e5]; 
% stop = [-6.42e5,-7.02e5];
% file = "radarData/53908058/Data_20131126_01_044";


% IRMCR1B_20141115_06_004 () Institute
% xmin = -975e3;
% xmax = -865e3;
% ymin =  215e3;
% ymax =  355e3;
% start = [-935e3,282e3]; 
% stop = [-895e3,288e3];
% file = 'radarData/120358122/Data_20141115_06_004';


% IRMCR1B_20141115_06_006
% xmin = -910e3;
% xmax = -860e3;
% ymin =  235e3;
% ymax =  310e3;
% start = [-860e3,280e3]; 
% stop = [-900e3,268e3];
% file = 'radarData/120358087/Data_20141115_06_006';

% file = 'radarData_good/Data_20131126_01_035';


if(file ~= "")
    load(file + ".mat");
    load(file + "_Layers.mat")
    [start(1),start(2)] = ll2ps(Latitude(1),Longitude(1));
    [stop(1),stop(2)]   = ll2ps(Latitude(end),Longitude(end));
    xmax = max([start(1) stop(1)])+50e3;
    xmin = min([start(1) stop(1)])-50e3;
    ymax = max([start(2) stop(2)])+50e3;
    ymin = min([start(2) stop(2)])-50e3;
%     bottom    = ncread(file,'Bottom');       %[s]
%     surface   = ncread(file,'Surface');      %[s]
%     Time  = ncread(file,'fasttime');     %[µs]
%     Data = ncread(file,'amplitude');    %[dB]
%     latitude  = ncread(file,'lat');
%     longitude =ncread(file,'lon');
end

dist  = stop-start;
% xx    = start(1):dist(1)/200:stop(1);
% yy    = start(2):dist(2)/200:stop(2); 
[xx,yy] = ll2ps(Latitude,Longitude);
% Grid to plot around transect
dx = 1e3;
smth = 10e3;
xi = xmin-dx:dx:xmax+dx;
yi = ymin-dx:dx:ymax+dx;
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
K_dif = 34.4/3.154e7;       % Thermal diffusivity of ice [m^2/s]
A_m   = 2.4e-24;            % Meyer's prefactor [Pa^-3 s^-1]
T_m   = 273;                % Ice melting point [k] 
sig_p = 6.6e-6;             % ice conditivity, [S m^{-1}]  Siemens/meter
E_p   = 0.55;               % Ice E, [eV]
T_r   = 251;                % ref temp [K]
k_b     = 8.617e-5;           % boltzmann constant [eV K^{-1}]

dz = .025;  %vertical resolution of thermal depth profiles (frac of H) [ ]

%% Import data and smooth 
[Acc, T_s] = loadALBMAP(); %accumulation rate and surface temp [m/s] [K]
Geo = loadGEO(); %geothermal heat flux from Shen [W/m^2]
T = T_s(xy(:,1),xy(:,2));
bm_b =  bedmachine_interp('bed',Xi,Yi);
bm_s =  bedmachine_interp('surface',Xi,Yi);
[u,v] = measures_interp('velocity',xy(:,1),xy(:,2)); %[m/yr]
[spd] = measures_interp('speed',Xi,Yi)'; %[m/yr]
u = u/3.154e7;
v = v/3.154e7;
%% Smoothing
% Numerator is the window we're smoothing over in [m], spacing of these grids

smoothbed = bm_b;%imgaussfilt(bm_b,smth/dx);
smoothsurf = bm_s;%imgaussfilt(bm_s,smth/dx);
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
        
        % Geothermal heat flux
%         Geo =@(x,y) 50e-3; % geothermal heat flux, W/m^2
        
        % Robin temp profile from Rezvanbehbahani2019. I've capped at 0C.
        % (z = 0 is bed)
        q =@(x,y) sqrt(Acc(x,y)./(2.*K_dif.*h_init(x,y)));
        t_robin =@(x,y) min(T_s(x,y) - Geo(x,y)*sqrt(pi)./(2*K.*q(x,y)).*(erf((0:dz:1).*h_init(x,y).*q(x,y))-erf(h_init(x,y).*q(x,y))),273.15);
        
        % sigma(t)
        sig =@(t) sig_p .* exp(E_p./k_b.*(1./T_r - 1./t)); %Pure ice only for now
        
        % Attenuation
        c2a = 0.912e6; %conversion factor from sigma [S/m] to antenuation [dB/m], see Macgregor et al 2007 eq (10)
        atten =@(x,y) (trapz(sig(t_z(x,y)).*c2a,2).*dz.*h_init(x,y)./1e3);
        
        % Attenuation
        c2a = 0.912e6; %conversion factor from sigma [S/m] to antenuation [dB/m], see Macgregor et al 2007 eq (10)
        atten_robin =@(x,y) (trapz(sig(t_robin(x,y)).*c2a,2).*dz.*h_init(x,y)./1e3);
        
        % Combo Temp
        t_combo =@(x,y) max(t_robin(x,y),t_z(x,y));
        
        % Combo Attenuation
        c2a = 0.912e6; %conversion factor from sigma [S/m] to antenuation [dB/m], see Macgregor et al 2007 eq (10)
        atten_combo =@(x,y) (trapz(sig(t_combo(x,y)).*c2a,2).*dz.*h_init(x,y)./1e3);
        
        % Enhancement Factor []
        E_t =@(x,y) depthIntEnhancement(t_z(x,y),a.^(-3),dz);
        
        % Mean Temp [K]
        T_bar = @(x,y) trapz(t_z(x,y),2)*dz;

%% Plotting, allow cross section        

        
% Plotting temp depth averaged
% figure(1)
% clf
% % p = surf(reshape(xy(:,1),size(Xi)),reshape(xy(:,2),size(Xi)),reshape(T_bar(xy(:,1),xy(:,2))-273.15,size(Xi)));
% p = surf(reshape(xy(:,1),size(Xi)),reshape(xy(:,2),size(Xi)),zeros(size(spd')),spd');
% hold on
% plot(xx,yy,'r*-','linewidth',2)
% scatter(xx(1),yy(1),100,'kp')
% contour(xi,yi,spd, [10, 10] , 'k:','HandleVisibility','off')
% contour(xi,yi,spd, [30, 30] , 'k--','HandleVisibility','off')
% contour(xi,yi,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
% contour(xi,yi,spd, [1000, 1000] , 'k-','LineWidth',2)
% % title('Temp Avg')
% set(p, 'edgecolor', 'none');
% colorbar
% view(2)
% % axis equal
% 
% figure(2)
% clf
% subplot(231)
tempM = t_z(xx',yy')'-273.15;
height = zeros(size(tempM));
x_along = zeros(size(tempM));
for i = 1:length(xx)
    height(:,i) = h_b_init(xx(i),yy(i)) + (0:dz:1)'*h_init(xx(i),yy(i));
    x_along(:,i) = sqrt((xx(1) - xx(i))^2 + (yy(1) - yy(i))^2);
end   
% p = surf(x_along, height, tempM);
% set(p, 'edgecolor', 'none');
% view(2)
% hold on
% plot(x_along(1,:),measures_interp('speed',xx,yy),'LineWidth', 3)
% scatter(0,h_b_init(xx(1),yy(1))-100,100,'kp')
% colorbar;
% caxis([min(T_s(xx,yy))-273.15, 0])
% title('Temp Profile Meyer')
% 
% subplot(232)
% tempM = t_robin(xx',yy')'-273.15;
% height = zeros(size(tempM));
% x_along = zeros(size(tempM));
% for i = 1:length(xx)
%     height(:,i) = h_b_init(xx(i),yy(i)) + (0:dz:1)'*h_init(xx(i),yy(i));
%     x_along(:,i) = sqrt((xx(1) - xx(i))^2 + (yy(1) - yy(i))^2);
% end   
% p = surf(x_along, height, tempM);
% set(p, 'edgecolor', 'none');
% view(2)
% hold on
% plot(x_along(1,:),measures_interp('speed',xx,yy),'LineWidth', 3)
% scatter(0,h_b_init(xx(1),yy(1))-100,100,'kp')
% colorbar;
% caxis([min(T_s(xx,yy))-273.15, 0])
% title('Temp Profile Robin')
% 
% subplot(233)
% tempM = t_combo(xx',yy')'-273.15;
% height = zeros(size(tempM));
% x_along = zeros(size(tempM));
% for i = 1:length(xx)
%     height(:,i) = h_b_init(xx(i),yy(i)) + (0:dz:1)'*h_init(xx(i),yy(i));
%     x_along(:,i) = sqrt((xx(1) - xx(i))^2 + (yy(1) - yy(i))^2);
% end   
% p = surf(x_along, height, tempM);
% set(p, 'edgecolor', 'none');
% view(2)
% hold on
% plot(x_along(1,:),measures_interp('speed',xx,yy),'LineWidth', 3)
% scatter(0,h_b_init(xx(1),yy(1))-100,100,'kp')
% colorbar;
% caxis([min(T_s(xx,yy))-273.15, 0])
% title('Temp Profile Combo')
% 
% subplot(212)
% plot(x_along(1,:),atten(xx',yy'),'LineWidth', 3)
% hold on
% plot(x_along(1,:),atten_robin(xx',yy'),'LineWidth', 3)
% plot(x_along(1,:),atten_combo(xx',yy'),'LineWidth', 3)
% scatter(0,atten(xx(1),yy(1))*.95,100,'kp','HandleVisibility','off')
% title("Expected Thermal Attenuation")
% ylabel("Attenuation [dB]")
% xlabel("Along profile [m]")
% legend("Meyer et al","Robin","combo");
% 
% figure(3)
% clf
% p = surf(reshape(xy(:,1),size(Xi)),reshape(xy(:,2),size(Xi)),-1*reshape(atten(xy(:,1),xy(:,2)),size(Xi)));
% hold on
% plot(xx,yy,'r*-','linewidth',2)
% scatter(xx(1),yy(1),100,'kp')
% contour(xi,yi,spd, [10, 10] , 'k:','HandleVisibility','off')
% contour(xi,yi,spd, [30, 30] , 'k--','HandleVisibility','off')
% contour(xi,yi,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
% contour(xi,yi,spd, [1000, 1000] , 'k-','LineWidth',2)
% title('Round trip thermal attentuation')
% set(p, 'edgecolor', 'none');
% colorbar
% caxis([min(-1*atten(xy(:,1),xy(:,2))),max(-1*atten(xy(:,1),xy(:,2)))])
% view(2)

%%
% figure(4)
% clf
% bedmachine_profile(xx,yy)

%%
if(file ~= "")
    Surface_layer = layerData{1}.value{2}.data;
    Bottom_layer = layerData{2}.value{2}.data;
 
    slowtime = 1:numel(Bottom);
    bedPower = 10*log10(interp2(slowtime,Time,Data,slowtime,Bottom_layer));
    
    range_u = 30;
    range_d = 40;
    dt          = Time(2)-Time(1);
    bed_i       = floor(Bottom_layer/dt);
    Bottom_2    = zeros(size(Bottom_layer));
    bedPower_2  = zeros(size(Bottom_layer));
    for in = 1:numel(bed_i)
    [bp,bed_2_i] = max(Data((bed_i(in)-range_u):(bed_i(in)+range_d) ,in)); % Find max in neighborhood
    Bottom_2(in) = Time(bed_2_i+bed_i(in)-range_u-1);              % Find Time of max
    bedPower_2(in) = 10*log10(bp);              % convert power to dB for plotting
    end
    
    figure
    clf
    
    subplot(421)
        imagesc(slowtime,Time,10*log10(Data))
        hold on
        plot(slowtime,Surface_layer,'b')
        plot(slowtime,Bottom_layer,'-','color',rgb('dark gray'),'linewidth',1)
        plot(slowtime,Bottom_2,'-','color',rgb('red'),'linewidth',1)
        hold off
        colorbar
        ylim([0 6]*1e-5)
        title('Radargram with bedpick')
    
    subplot(423)
        plot(slowtime,bedPower,'color',rgb('light gray'))
        hold on
        plot(slowtime,bedPower_2,'color',rgb('light red'))
        plot(slowtime,movmean(bedPower,15),'--','linewidth',2,'color',rgb('dark gray'))
        plot(slowtime,movmean(bedPower_2,15),'--','linewidth',2,'color',rgb('dark red'))
        title('Bed Echo Power')
    
    subplot(425)
        plot(slowtime,atten_combo(xx',yy'),'LineWidth', 3)
        title('Expected Thermal Attenuation')
    
    subplot(427)
        plot(slowtime,Bottom_layer - Surface_layer)
        title('Ice thickness along profile')
    
    
    subplot(222)
        p = surf(reshape(xy(:,1),size(Xi)),reshape(xy(:,2),size(Xi)),zeros(size(spd')),spd');
        hold on
        plot(xx,yy,'r*-','linewidth',2)
        scatter(xx(1),yy(1),100,'kp')
        contour(xi,yi,spd, [10, 10] , 'k:','HandleVisibility','off')
        contour(xi,yi,spd, [30, 30] , 'k--','HandleVisibility','off')
        contour(xi,yi,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
        contour(xi,yi,spd, [1000, 1000] , 'k-','LineWidth',2)
        % title('Temp Avg')
        set(p, 'edgecolor', 'none');
        colorbar
        view(2)
        title('Map view of transect with surface velocity')
    
    subplot(224)
        tempM = t_combo(xx',yy')'-273.15;
        height = zeros(size(tempM));
        x_along = zeros(size(tempM));
        for i = 1:length(xx)
            height(:,i) = h_b_init(xx(i),yy(i)) + (0:dz:1)'*h_init(xx(i),yy(i));
            x_along(:,i) = sqrt((xx(1) - xx(i))^2 + (yy(1) - yy(i))^2);
        end   
        p = surf(x_along, height, tempM);
        set(p, 'edgecolor', 'none');
        view(2)
        hold on
        plot(x_along(1,:),measures_interp('speed',xx,yy),'LineWidth', 3)
        scatter(0,h_b_init(xx(1),yy(1))-100,100,'kp')
        colorbar;
        caxis([min(T_s(xx,yy))-273.15, 0])
        title('Modeled Temp Profile along transect with surface V overlay')
    
        drawnow
end