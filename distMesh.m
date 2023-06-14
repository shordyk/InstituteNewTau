clear
close all
if(~ismac && ~ispc)
    %Start-up business on sherlock is hard 
    addpath('lib') 
    addpath ../MATLAB
    run startup.m
end
% Institute box
xbox = [-10.00  -9.25   -8.25   -9.00   -10.00]*1e5;
ybox = [  2.33   3.40    2.80    1.7      2.33]*1e5;

dx = 2e3;
overgrab = 10;
xmax =  max(xbox);
xmin = min(xbox);
ymax =  max(ybox);
ymin =  min(ybox);

pv = [xbox; ybox]';
pv = makePerimeter(pv,dx);
xi = xmin-dx*overgrab:dx/2:xmax+dx*overgrab;
yi = ymin-dx*overgrab:dx/2:ymax+dx*overgrab;

%% Making flow lines, optional 
[Xi,Yi] = meshgrid(xi,yi);

% [ui,vi] = measures_interp('velocity',Xi,Yi);
% spd = measures_interp('speed',Xi,Yi);

% figure
% quiver(Xi,Yi,ui,vi)
% hold on
% contour(xi,yi,spd, [10, 10] , 'k:','HandleVisibility','off')
% contour(xi,yi,spd, [30, 30] , 'k--','HandleVisibility','off')
% contour(xi,yi,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
% contour(xi,yi,spd, [1000, 1000] , 'k-','HandleVisibility','off','LineWidth',2)
% title('Quiver')
% plot(pv(:,1),pv(:,2),'-*b')
% h1 = stream2(xi,yi,ui,vi,pv(1,1)-1e4,pv(1,2));
% h2 = stream2(xi,yi,ui,vi,pv(2,1),pv(2,2));
% line1 = streamline(h1);
% line2 = streamline(h2);
% set(line1,'Color','red');
% set(line2,'Color','red');
% plot(h1{1,1}(1:100:2700,1),h1{1,1}(1:100:2700,2),'*-','color',rgb('Forest Green'))
% plot(h2{1,1}(1:100:1300,1),h2{1,1}(1:100:1300,2),'*-','color',rgb('Forest Green'))
% pv_new = [(h1{1,1}(1:10:2700,1)-40e3),h1{1,1}(1:10:2700,2)-10e3;flipud(h2{1,1}(1:10:1300,1)+35e3),flipud(h2{1,1}(1:10:1300,2))+10e3];
% pv_new = [pv_new ; pv_new(1,:)];
% plot(pv_new(:,1),pv_new(:,2),'c-*')

% pv = pv_new;

%%
[xy,t]=distmesh2d(@dpoly,@huniform,.05,[xmin,ymin;xmax,ymax]/1e5,[],pv/1e5);
xy = xy*1e5;
disp("inital mesh complete")

[u,v] = measures_interp('velocity',xy(:,1),xy(:,2));

[xx,yy] = meshgrid(xmin:1e3:xmax,ymin:1e3:ymax);
spd = measures_interp('speed',xx,yy);

if(ismac)
    figure
    surf(xx,yy,zeros(size(spd)),log10(spd),'edgecolor','none')
    hold on 
    view(2)
    plot(pv(:,1),pv(:,2),'r-')
    axis equal
end

%%
if(sum(isnan(u)) + sum(isnan(v)) > 0)
    uFill = scatteredInterpolant(xy(~isnan(u),1),xy(~isnan(u),2),u(~isnan(u)));
    vFill = scatteredInterpolant(xy(~isnan(v),1),xy(~isnan(v),2),v(~isnan(v)));
    u(isnan(u)) = uFill(xy(isnan(u),1),xy(isnan(u),2));
    v(isnan(v)) = vFill(xy(isnan(v),1),xy(isnan(v),2));
    clear uFill vFill
end

if(ismac)
	figure
	trisurf(t,xy(:,1),xy(:,2),zeros(size(xy(:,1))),sqrt(u.^2+v.^2),...
	       'edgecolor','none')  
	colorbar
	view(2)
end
% ep  = calcTrigridStrain(u,v,xy,3e3);
    dx =2e3;
    xi = min(xy(:,1))-dx:dx:max(xy(:,1))+dx;
    yi = (min(xy(:,2))-dx:dx:max(xy(:,2))+dx)';
    [Xi,Yi] = ndgrid(xi,yi);
    us = scatteredInterpolant(xy(:,1),xy(:,2),u);
    vs = scatteredInterpolant(xy(:,1),xy(:,2),v);

    % grid to find strain rates
    ug = us(Xi,Yi); % [m/s]
    vg = vs(Xi,Yi); % [m/s]

    % Smooth
%     ug = imgaussfilt(ug,1);
%     vg = imgaussfilt(vg,1);
    
    % Calc Strain on grid
    [n,m] = size(ug);
    strain = zeros(n-1,m-1,2,2);
    strain(:,:,1,1) = diff(vg(1:end-1,:),1,2)/dx;
    strain(:,:,2,2) = diff(ug(:,1:end-1))/dx;
    strain(:,:,1,2) = .5 .* (diff(ug(1:end-1,:),1,2) + diff(vg(:,1:end-1)))/dx;
    strain(:,:,2,1) =  strain(:,:,1,2);
%   Total  Effective Strain
    EffStrain = sqrt(.5*(strain(:,:,1,1).^2 + strain(:,:,2,2).^2 +...
            (strain(:,:,1,1) + strain(:,:,2)).^2) + strain(:,:,1,2).^2);
    EffStrain(isnan(EffStrain)) = 1e-14; % this doesn't happen anymore Yay!
    EffStrain = imgaussfilt(EffStrain,5);
    ep = griddedInterpolant(Xi(2:end,2:end)-dx/2, Yi(2:end,2:end)-dx/2, EffStrain,'linear','nearest');




fun =@(x,y) min(subplus(log(ep(x*1e5,y*1e5))+9),5);

if(ismac)
	figure
	scatter(xy(:,1),xy(:,2),[],fun(xy(:,1)/1e5,xy(:,2)/1e5),'filled')
	colorbar
end

% figure(2)
% scatter(xx(:),yy(:),[],lgSpd(:),'filled')
% colorbar

%%
fd=@(p) dpoly(p,pv/1e5);
% fh=@(p) 10-log(measures_interp('speed',p(:,1)*1e5,p(:,2)*1e5)+1);
fh=@(p) ones(size(p,1),1) - fun(p(:,1),p(:,2))/6;

edgeLength = .035;

[xy,t]=distmesh2d(fd,fh,edgeLength,[xmin,ymin;xmax,ymax]/1e5,[]);
disp("Successfully meshed at " + edgeLength);
xy = xy*1e5;

%%
% Flowline boundaries
% figure
% clf
% scatter(xy(:,1),xy(:,2),'k.')
% hold on
% idxSw = knnsearch(xy,[h1{1,1}(1:10:2700,1)-40e3,h1{1,1}(1:10:2700,2)-10e3]);
% idxNe = knnsearch(xy,[h2{1,1}(1:10:1300,1)+35e3,h2{1,1}(1:10:1300,2)+10e3]);
% idxNw = knnsearch(xy,[linspace(pv(end-1,1),pv(end,1),100)',linspace(pv(end-1,2),pv(end,2),100)']);
% pvPick = 270;
% idxSe = knnsearch(xy,[linspace(pv(pvPick,1),pv(pvPick+1,1),100)',linspace(pv(pvPick,2),pv(pvPick+1,2),100)']);
% 
% sw_bound = false(size(xy(:,1)));
% ne_bound = false(size(xy(:,1)));
% nw_bound = false(size(xy(:,1)));
% se_bound = false(size(xy(:,1)));
% sw_bound(idxSw) = 1;
% ne_bound(idxNe) = 1;
% nw_bound(idxNw) = 1;
% se_bound(idxSe) = 1;


%% linear boundaries
se_bound = (ybox(4)-ybox(3))/(xbox(4)-xbox(3))*xy(:,1) - xy(:,2)  > (ybox(4)-ybox(3))/(xbox(4)-xbox(3))*xbox(3) - ybox(3)-dx/3;
ne_bound = (ybox(3)-ybox(2))/(xbox(3)-xbox(2))*xy(:,1) - xy(:,2)  < (ybox(3)-ybox(2))/(xbox(3)-xbox(2))*xbox(2) - ybox(2)+dx/3;
nw_bound = (ybox(2)-ybox(1))/(xbox(2)-xbox(1))*xy(:,1) - xy(:,2)  < (ybox(2)-ybox(1))/(xbox(2)-xbox(1))*xbox(1) - ybox(1)+dx/3;
sw_bound = (ybox(1)-ybox(4))/(xbox(1)-xbox(4))*xy(:,1) - xy(:,2)  > (ybox(1)-ybox(4))/(xbox(1)-xbox(4))*xbox(4) - ybox(4)-dx/3;


if(ismac)
    figure(2)
    clf
    scatter(xy(:,1),xy(:,2),'k.')
    hold on
%     scatter(xy(sw_bound2,1),xy(sw_bound2,2),'r')
%     scatter(xy(ne_bound2,1),xy(ne_bound2,2),'r')
%     scatter(xy(nw_bound2,1),xy(nw_bound2,2),'r')
%     scatter(xy(se_bound2,1),xy(se_bound2,2),'r')
    scatter(xy(sw_bound,1),xy(sw_bound,2),'g*')
    scatter(xy(ne_bound,1),xy(ne_bound,2),'r*')
    scatter(xy(nw_bound,1),xy(nw_bound,2),'c*')
    scatter(xy(se_bound,1),xy(se_bound,2),'b*')
end
%%
if(ismac)
	figure
	clear u
	u = measures_interp('speed',xy(:,1),xy(:,2));
	uFill = scatteredInterpolant(xy(~isnan(u),1),xy(~isnan(u),2),u(~isnan(u)));
	u(isnan(u)) = uFill(xy(isnan(u),1),xy(isnan(u),1));
	   
	trisurf(t,xy(:,1),xy(:,2),zeros(size(xy(:,1))),u,...
	       'edgecolor','none')   
	view(2)
	colorbar   
end

try
    save("grid/strainMesh" + strrep(string(edgeLength),"0.","") + ".mat");
catch
    %% The grid folder might not exist first time running this
    warning("Error saving, grid directory might be missing. Making one now.")
    mkdir grid
    save("grid/strainMesh" + strrep(string(edgeLength),"0.","") + ".mat");
end

disp("Successfully Saved")
