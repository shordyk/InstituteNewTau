% File used to test/sandbox out functions on the grid. Mostly for work importing from external datasets and make plots
%% Make Grid
clear
close all

if(true) %gate to actually make a new grid
    dx = 4e3;    %dx: nominal grid spacing [m]

    xmax = -8.75e5;
    xmin = -9.75e5;
    ymax =  3.5e5;
    ymin =  1.5e5;

    x_range = xmax-xmin; 
    y_range = ymax-ymin;

    theta = pi/3; %[rad] angle of rotation
    % Mesh Generation
    clf
    pv = [xmin,ymax;xmin,ymin;xmax,ymin;xmax,ymax];
    pv = makePerimeter(pv,dx); %Places points on the perimeter for distMesh to use
    figure(2)
    [xy,t] = distmesh2d(@dpoly,@huniform,dx/1e5,[xmin,ymin;xmax,ymax]/1e5,pv/1e5,pv/1e5);
    % scale to real grid size
    xy = xy*1e5;
    %grab boundry points before rotating
    dwnSt_bound = xy(:,1) < xmin + dx*1/3;
    upSt_bound  = xy(:,1) > xmax - dx*1/3;
    lfSt_bound  = xy(:,2) < ymin + dx*1/3;
    rtSt_bound  = xy(:,2) > ymax - dx*1/3;
    %Rotate to align with flow (anti flow really)
    xy = ((xy - mean(xy)) * [cos(theta) -sin(theta); sin(theta) cos(theta)]) + mean(xy);
    save("gridInstitute" + dx + ".mat");
end
% load('workingGrid_x12_2.mat')
icey = cbrewer('div','BrBG',48);


dz = .02;

%% Physical parameters
a     = 2.1E8;       %a:     flow parameter pre-factor [Pa s^1/3]
p     = 4/3;         %p:     flow parameter power []
g     = 10;          %g:     acceleration due to gravity [m/s^2]
rho   = 917;         %rho:   density of ice [kg/m^3]
rho_w = 1000;        %rho_w: density of water [kg/m^3]
C_p   = 2050;        %specific heat of ice [J/Kg/K]
K     = 2.1;         %thermal conductivity of ice [W/m/K]
A_m   = a^(-3); %prefactor [Pa^-3 s^-1]
nn    = 3;           %Glens law power
T_m   = 273;         %Ice melting point [k] 
 
%% Get observed Data
overgrab = 20;
xi = xmin-dx*overgrab:dx/3:xmax+dx*overgrab;
yi = ymin-dx*overgrab:dx/3:ymax+dx*overgrab;
[Xi,Yi] = meshgrid(xi,yi);
bdmp_b =  bedmap2_interp(Xi,Yi,'bed');
bdmp_s =  bedmap2_interp(Xi,Yi,'surface');
bdmc_b =  bedmachine_interp('bed',Xi,Yi);
bdmc_s =  bedmachine_interp('surface',Xi,Yi);
spd = measures_interp('speed',Xi,Yi);
[u,v] = measures_interp('velocity',Xi,Yi);
[sdx, sdy] = gradient(bdmp_s);
ds = sqrt(sdx.^2+sdy.^2);
n=16;
xmid = (xmin+xmax)/2;
ymid = (ymin+ymax)/2;
startX = xmax*ones(1,n)+2e4;
startY = ymin:(ymax-ymin)/(n-1):ymax;
xybar = mean(xy);

%Clean up, smooth, clip fields
filter     = 1;
% smoothbed  = bdmc_b;%imgaussfilt(bdmc_b,filter);
% smoothsurf = bdmc_s;%imgaussfilt(bdmc_s,filter);

smoothbed = sgolayfilt(bdmc_b,2,3);
smoothsurf = sgolayfilt(bdmc_s,2,2*floor(10e3/dx)+1);

smoothspd  = imgaussfilt(spd,filter);
smoothu    = imgaussfilt(u,filter);
smoothv    = imgaussfilt(v,filter);




%% HydroPotential
% phi = phi_init(Xi,Yi);
% phiS = imgaussfilt(phi,4);
% [phix,phiy] = gradient(-1000*phiS);
% phiN = 120;
% phiX0 = [ones(1,phiN)*(xmax+x_bar)/2,ones(1,phiN)*xmax];
% phiY0 = [linspace(ymin,ymax,phiN),linspace(ymin,ymax,phiN)];

% Antiflow lines
% load vel_profilesAll.mat
% [antiFlowx, antiFlowy] = ll2ps(profile_lat(:,1:end-3),profile_lon(:,1:end-3));
% [st2x,st2y] = ll2ps(-76.4085,-103.4856);

% Run to get streams
% StreamRouting

%% Plot

% figure
%     subplot(121)
%         trisurf(t,xy(:,1),xy(:,2),xi(xy(:,1),xy(:,2)),'edgecolor','none');
%         colorbar
%         view(2)
%         title('\xi')
% %         axis equal
%     subplot(122)
%         trisurf(t,xy(:,1),xy(:,2),ep_dot(xy(:,1),xy(:,2))./ep_star(xy(:,1),xy(:,2)),'edgecolor','none');
%         colorbar
%         view(2)
%         title('Strain Rate /Critical Strain Rate')
%         caxis([0 3])
% %         axis equal
% 
% figure
%     trisurf(t,xy(:,1),xy(:,2),E_t(xy(:,1),xy(:,2)),'edgecolor','none');
%     colorbar
%     view(2)
%     title('Enhancement')
%         axis equal


% figure
%     subplot(121)
%         trisurf(t,xy(:,1),xy(:,2),Pe(xy(:,1),xy(:,2)),'edgecolor','none');
%         colorbar
%         caxis([0 50])
%         view(2)
%         title('Peclet')
% %         axis equal
%     subplot(122)
%         trisurf(t,xy(:,1),xy(:,2),Br(xy(:,1),xy(:,2)),'edgecolor','none');
%         colorbar
%         caxis([0 10])
%         view(2)
%         title('Brinkman')
% %         axis equal
        
%% 
figure
    clf
    h = pcolor(Xi,Yi,spd);
    hold on;
    set(h, 'EdgeColor', 'none');
    
%     scatter(xy(:,1),xy(:,2),'wo')
    title('Bed Height');
    alpha(h,1);
%     scatter(xy(lfSt_bound,1),xy(lfSt_bound,2),'ko')
%     scatter(xy(rtSt_bound,1),xy(rtSt_bound,2),'ko')
    trisurf(t,xy(:,1),xy(:,2),zeros(size(xy(:,1))),'edgecolor','w','facecolor','none')
    plot(xy(upSt_bound,1),xy(upSt_bound,2),'k','Linewidth',4)
    plot(xy(dwnSt_bound,1),xy(dwnSt_bound,2),'k','Linewidth',4)
    plot(xy(lfSt_bound,1),xy(lfSt_bound,2),'k','Linewidth',4)
    plot(xy(rtSt_bound,1),xy(rtSt_bound,2),'k','Linewidth',4)
    scatter(xy(upSt_bound,1),xy(upSt_bound,2),[],rgb('light red'),'filled')
    scatter(xy(dwnSt_bound,1),xy(dwnSt_bound,2),[],rgb('dark red'),'filled')
%     p = plot(S,'-','linewidth',3,'color',rgb('blue'));
%     scatter(st2x,st2y,200,'mp','filled');
%     plot(antiFlowx(:,[2,5,10]),antiFlowy(:,[2,5,10]),'k','linewidth',3)
    % slphi = streamline(stream2(Xi,Yi,phix,phiy,phiX0,phiY0,[.1]));
    % set(slphi,'Color','red')
%     slVel = streamline(stream2(Xi,Yi,u,v,startX([1,8,end]),startY([1,8,end])));
%     set(slVel,'Color','blue','linewidth',3)
%     contour(xi,yi,log10(ds), [1.8:.2:2.2] , 'r-','HandleVisibility','off');
    contour(xi,yi,spd, [30, 30] , 'k--','HandleVisibility','off');
    contour(xi,yi,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
    contour(xi,yi,spd, [1000, 1000] , 'k-','LineWidth',2)
    view(2);
    colorbar
%     caxis([-.6 3.6]); 
%     legend('Log_{10} Speed','Model Domain','TIME Site 2','Antiflow lines','Flow Lines'...,'Vel Streamlines','Phi Steamlines','Meltwater Routing',);
%         );
%     axis off
    axis equal
%% 
% figure(2)
%     clf
%     h = pcolor(Xi,Yi,bdmc_s-bdmc_b);
%     hold on
%     contour(xi,yi,spd, [30, 30] , 'k--','HandleVisibility','off')
%     contour(xi,yi,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
%     contour(xi,yi,spd, [1000, 1000] , 'k-','LineWidth',2)
%     colorbar
%     title('Depth of bed [m]');
%     set(h, 'EdgeColor', 'none');
%     alpha(h,1);
%     
%     hold on
%     plot(ex1_x_clean,ex1_y_clean,'ko-')
%     plot(ex2_x_clean,ex2_y_clean,'ko-')
%     plot(pv(:,1),pv(:,2),'r*-')
%     axis equal
% 
% figure(2)
%     clf
%     h = pcolor(Xi,Yi,bdmp_b);
%     colorbar
%     title('Depth of bedMap [m]');
%     set(h, 'EdgeColor', 'none');
%     alpha(h,1);
% 
% figure(3)
%     clf
%     set(gca,'FontSize',24)
%     trisurf(t,xy(:,1),xy(:,2),h_s_init(xy(:,1),xy(:,2)),...
%            'edgecolor','none')
%     % caxis([0 max(-u)*3.154E7]);
%     hold on
%     colormap(icey);
%     trisurf(t,xy(:,1),xy(:,2),h_b_init(xy(:,1),xy(:,2)),...
%            'edgecolor','none')%,'facecolor','interp')
%     title('Bed and Surface')
%     
%     xlabel('X')
%     ylabel('Y')
%     zlabel('Elevation [m]')
%     colorbar

% figure(4)
%     clf
%     trisurf(t,xy(:,1),xy(:,2),h_s_init(xy(:,1),xy(:,2))-h_b_init(xy(:,1),xy(:,2)),...
%         'edgecolor','none')%,'facecolor','interp');
%     colorbar
%     colormap(icey);
%     title('H [m]');
% 
% figure(5)
%     clf
%     hist(reshape(bdmp_b-bdmc_b,1,numel(bdmp_b)),50);
%     title('Distribution of \Delta Map vs Machine');