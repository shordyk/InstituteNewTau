% Requires download and install the following software: 
% Distmesh https://popersson.github.io/distmesh/index.html
% CVX http://cvxr.com/cvx/download/
% MEaSUREs matlab plug in + data files
% https://www.mathworks.com/matlabcentral/fileexchange/47329-measures
% BedMachine matlab plug in + data files
% https://www.mathworks.com/matlabcentral/fileexchange/69159-bedmachine

% clc
% clear
% close all

%% Initialization
% Scenario to run if running one at a time comment out below, run this file
% directly

% str = 'Mixed';
% mapFile = 'gridInstitute5000.mat';

% Comment so we know what's happening, thats always nice.
disp("Running " + str + " now...");

% Load input files
initializeInputs();

% Build model bed, geometery, functions
initializeModel();

%% Define \tau
% Define map of basal strength according to named scenario. 
tau_c = defineTau(str);

%% Build System
buildSystem();


%% Thermomechanical coupling loop
for t_i = 1:100  
    % Thermocouple fields to update everyloop
    % Strain rate [s^-1]
        ep_dot = calcTrigridStrain(u,v,xy,dx); %returns intperolation object
        
%         if(true)  adjust if Advection calculation self-consistent
%             T_calc = T;
%         else
%             T_calc = T_bar(xy(:,1),xy(:,2));
%         end
%         lambda  = calcAdvection(T_calc,u,v,xy,dx,rho,C_p); 

        % Brinkman number [ ]
        Br =@(x,y) 2*subplus(h_s_init(x,y)-h_b_init(x,y)).^2./(K*(T_m-T_s(x,y))).*((subplus(ep_dot(x,y)).^(nn+1))/A_m).^(1/nn);

        % Peclet number  [ ]
        Pe =@(x,y) rho*C_p.*Acc(x,y).*subplus(h_s_init(x,y)-h_b_init(x,y))./K;

        % Horizontal Peclet number  [ ]
%         La =@(x,y) lambda(x,y).*subplus(h_s_init(x,y)-h_b_init(x,y)).^2./(K*(T_m-T_s(x,y)));
        La =@(x,y) zeros(size(x)); %exclude advection
        
        % Critical Strain [s^-1]
        ep_star =@(x,y) ((La(x,y)/2 + ((Pe(x,y).^2)/2)./(Pe(x,y)-1+exp(-Pe(x,y))))).^(nn/(nn+1))...
        .*(K*(T_m-T_s(x,y))./(A_m.^(-1/nn).*(subplus(h_s_init(x,y)-h_b_init(x,y))).^2)).^(nn/(nn+1));

        % Temp profile at xy [K]
        t_z =@(x,y) tempProfile(ep_dot(x,y),ep_star(x,y),Pe(x,y),Br(x,y),La(x,y),T_s(x,y),T_m,dz); 

        % Enhancement Factor []
        E_t =@(x,y) depthIntEnhancement(t_z(x,y),a.^(-3),dz);

        % Mean Temp [K]
        T_bar = @(x,y) trapz(t_z(x,y),2)*dz;

    % Calc Enhancement Factors, relax into solution. Have max value for
    % stabilization
    cap = 30^(-1/nn); %stability cap on enhancement TODO remove
    e_new = (E_t(xy_c(:,1),xy_c(:,2))).^(-1/nn);
    e_new(e_new < cap) = cap;  % max enhancement is a min viscosity
    
    if(t_i == 1)
        enhance = e_new;
    end

    if(t_i ~= 1) %first step we don't relax, we use E = 1 everywhere (zero strain is also an options)
        enhance = (1-nu) * enhance + nu*e_new;
        res = norm(e_new - enhance) / norm(e_new);
        disp("Residual: " + res);
        
        if (res < 1e-3) %check for thermal stabilization
            break; 
        end
    end

    %% Solve
    % Unused BCs
    %       v(xy(:,2) > ymax - dx/2) == 0;
    %       v(xy(:,2) < ymin - dx/2) == 0;
    
    % Options include     cvx_precision low, cvx_begin quiet
    % CVX may throw a warning about non-empty problems here, that is OK.
    cvx_begin quiet
        variables u(nN) v(nN)
        obj = 2.*a./p.*sum(enhance.*h_av.*tau_area.*pow_pos(norms([A*u,B*v,1/2*(B*u+A*v)],2,2),p)) + ...
              F*tau_c(xy(:,1),xy(:,2),u,v) + ...
              rho*g*sum(h_av.*((A*h_s).*(D*u) + (B*h_s).*(D*v)));
        subject to
            u(nw_bound) == spd_BC_u./3.154E7;
            v(nw_bound) == spd_BC_v./3.154E7;
            u(se_bound) == spd_BC_u2./3.154E7;
            v(se_bound) == spd_BC_v2./3.154E7;
            u(sw_bound) == spd_BC_uL./3.154E7;
            v(sw_bound) == spd_BC_vL./3.154E7;
        minimize(obj)
    cvx_end
    if(~strcmp(cvx_status,"Solved"))
        disp("CVX has issues: " + cvx_status);
        if(~contains(cvx_status,"Solved"))
            break;
        end
    end
    % u and v are [m/s]    
    %% Visualization in loop 
    %(uncomment to see avg temp, enhancement, and Pe, Lambda, Br every loop
    inLoopPlotting;
end
%% Save data to data file
clear fg1 fg2
mpClean = erase(mapFile, [".mat","workingGrid_"]);
save("data/data_" + mpClean + str + "noAdvect.mat");

%% Vis out of loop
spd2 = measures_interp('speed',xy(:,1),xy(:,2)); %[m/yr]

figure('Position', [0 0 1200 600]);
clf
sgtitle(str);
subplot(141)

trisurf(t,xy(:,1),xy(:,2),zeros(size(spd2)),(spd2),...
       'edgecolor','none')
hold on
title('Speed of Measures')
xlabel('X')
ylabel('Y')
f = gca;
f.ColorScale = 'log';
view(2)
colorbar
view(2)
axis equal

subplot(142)
trisurf(t,xy(:,1),xy(:,2),h_s_init(xy(:,1),xy(:,2)),(sqrt(u.^2 + v.^2)*3.154E7),...
       'edgecolor','none')   
caxis([0.3323  381.5379])
title('Speed')
xlabel('X')
ylabel('Y')
colorbar
f = gca;
f.ColorScale = 'log';
view(2)
axis equal

subplot(143)
trisurf(t,xy(:,1),xy(:,2),tau_c(xy(:,1),xy(:,2),u,v)./norms([u,v],2,2),...
       'edgecolor','none')
% hold on
% trisurf(t,xy(:,1),xy(:,2),h_b_init(xy(:,1),xy(:,2)),...
%        'edgecolor','black','facecolor','none')
colorbar
caxis([0e3 100e3]);
colormap(gca, Cmap/255.0)
title('Basal \tau')
xlabel('X')
ylabel('Y')
view(2)
axis equal

subplot(144)
trisurf(t_c,xy_c(:,1),xy_c(:,2),df,...
    'edgecolor','none');
title('Driving force')
xlabel('X')
ylabel('Y')
colorbar
caxis([0e3 150e3]);
colormap(gca, Cmap/255.0)
view(2)
axis equal

% load institute_antiflow/vel_profiles_paul.mat
 load vel_profiles_paul_gl_str_2022.mat

[um,vm] = measures_interp('velocity',xy(:,1),xy(:,2));

figure

AF1 = load("data/AntiFlow2_4.mat");
[antiflow_x, antiflow_y] = ll2ps(profile_lat(:,2),profile_lon(:,2));

subplot(221)
trisurf(t,xy(:,1),xy(:,2),zeros(size(xy(:,1))),(sqrt(u.^2 + v.^2)*3.154E7),...
       'edgecolor','none')   
caxis([0.3323  381.5379])
hold on
plot(antiflow_x,antiflow_y,'r','linewidth',3)
title('Speed Model')
xlabel('X')
ylabel('Y')
colorbar
f = gca;
f.ColorScale = 'log';
view(2)
axis equal

subplot(222)
trisurf(t,xy(:,1),xy(:,2),zeros(size(xy(:,1))),(sqrt(um.^2 + vm.^2)),...
       'edgecolor','none')   
caxis([0.3323  381.5379])
hold on
plot(antiflow_x,antiflow_y,'r','linewidth',3)
title('Speed MEaSUREs')
xlabel('X')
ylabel('Y')
colorbar
f = gca;
f.ColorScale = 'log';
view(2)
axis equal

spd_interp = scatteredInterpolant(xy(:,1),xy(:,2),(sqrt(u.^2 + v.^2)*3.154E7));

subplot(212)
plot(profile_path(:,2)-30.5E3,profile_cross(:,2),'LineWidth',3)
hold on
plot(profile_path(:,2)-30.5E3,spd_interp(antiflow_x,antiflow_y),'LineWidth',3)
plot(AF1.xy(AF1.xy(:,2) > 1500-AF1.dx/10,1),AF1.u(AF1.xy(:,2) > 1500-AF1.dx/10),'LineWidth',3)
legend('MEaSUREs','MapView Model','AntiFlow Model','Location','NorthWest')
xlim([-30E3,80E3])
ylabel('[m/yr]')
xlabel('[m]')
figure

AF1 = load("data/AntiFlow2022_3_4.mat");
[antiflow_x, antiflow_y] = ll2ps(profile_lat(:,3),profile_lon(:,3));
subplot(221)
trisurf(t,xy(:,1),xy(:,2),zeros(size(xy(:,1))),(sqrt(u.^2 + v.^2)*3.154E7),...
       'edgecolor','none')   
caxis([0.3323  381.5379])
hold on
plot(antiflow_x,antiflow_y,'r','linewidth',3)
title('Speed Model')
xlabel('X')
ylabel('Y')
colorbar
f = gca;
f.ColorScale = 'log';
view(2)
axis equal

subplot(222)
trisurf(t,xy(:,1),xy(:,2),zeros(size(xy(:,1))),(sqrt(um.^2 + vm.^2)),...
       'edgecolor','none')   
caxis([0.3323  381.5379])
hold on
plot(antiflow_x,antiflow_y,'r','linewidth',3)
title('Speed MEaSUREs')
xlabel('X')
ylabel('Y')
colorbar
f = gca;
f.ColorScale = 'log';
view(2)
axis equal

spd_interp = scatteredInterpolant(xy(:,1),xy(:,2),(sqrt(u.^2 + v.^2)*3.154E7));

subplot(212)
plot(profile_path(:,3)-52.5E3,profile_cross(:,3),'LineWidth',3)
hold on
plot(profile_path(:,3)-52.5E3,spd_interp(antiflow_x,antiflow_y),'LineWidth',3)
plot(AF1.xy(AF1.xy(:,2) > 1500-AF1.dx/10,1),AF1.u(AF1.xy(:,2) > 1500-AF1.dx/10),'LineWidth',3)
legend('MEaSUREs','MapView Model','AntiFlow Model','Location','NorthWest')
xlim([-30E3,80E3])
ylabel('[m/yr]')
xlabel('[m]')

[um,vm] = measures_interp('velocity',xy(:,1),xy(:,2));
figure
trisurf(t,xy(:,1),xy(:,2),tau_c(xy(:,1),xy(:,2),u,v)./norms([u,v],2,2),...
       'edgecolor','none')
title('Basal Strength');
caxis([0 1e5])
view(2)
colorbar
