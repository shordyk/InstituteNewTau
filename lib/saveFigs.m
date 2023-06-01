% Making figures for fig file. No figures are save automatically, but it is
% easy to do so. 
clear; clc; close all
cd data/
files = dir;
directoryNames = {files.name};
directoryNames = directoryNames(~ismember(directoryNames,{'.','..','.DS_Store'}));
cd ../

ftsize = 18;
for ii = 1:length(directoryNames)
load("data/"+directoryNames{ii});
str = directoryNames{ii};
str_clean = erase(directoryNames{ii}, [".mat","data_"])

figure(ii*2-1)
trisurf(t,xy(:,1),xy(:,2),tau_c(xy(:,1),xy(:,2),u,v)./norms([u,v],2,2)/1e3,...
       'edgecolor','none','facecolor','interp')
hold on
trisurf(t,xy(:,1),xy(:,2),-10e2*ones(size(xy(:,1))),...
       'edgecolor','black','facecolor','none')
c = colorbar;
c.Label.String = 'Basal Strength [kPa]';
caxis([50 150]);
colormap(gca, Cmap/255.0)
title("Basal Strength " + str_clean, 'FontSize',ftsize+2)
xlabel('Easting')
ylabel('Northing')
view(2)
axis equal
axis off
f = gca;
f.XAxis.FontSize = ftsize;
f.YAxis.FontSize = ftsize;
c.FontSize = ftsize;

figure(ii*2)
trisurf(t,xy(:,1),xy(:,2),h_s,sqrt(u.^2 + v.^2)*3.154E7,...
       'edgecolor','none','facecolor','interp')
hold on
trisurf(t,xy(:,1),xy(:,2),zeros(size(xy(:,1))),...
       'edgecolor','black','facecolor','none')
caxis([10^1 10^3.6]);   
title("Speed " + str_clean,  'FontSize',ftsize+2)
xlabel('Easting')
ylabel('Northing')
c = colorbar;
c.Label.String = 'Ice Speed [m/yr]';
view(2)
axis equal
axis off
f = gca;
f.XAxis.FontSize = ftsize;
f.YAxis.FontSize = ftsize;
f.ColorScale = 'log';
c.FontSize = ftsize;


clearvars -except ii directoryNames files ftsize
end

%% Bed and Surface
load("data/"+directoryNames{1});

figure(ii*2+1)
set(gcf,'Position',[100 100 1000 800])
clf
trisurf(t,xy(:,1),xy(:,2),h_b_init(xy(:,1),xy(:,2)),...
        'edgecolor','none','facecolor','interp');
hold on
[C,H] = tricontour(t,xy(:,1),xy(:,2),h_b_init(xy(:,1),xy(:,2)),[-1500:200:500]);
for i = 1:length(H)
    H(i).LineWidth = 1;
    H(i).EdgeColor = 'k';
    H(i).LineStyle = '--';
    H(i).ZData = round((h_b_init(H(i).XData(1),H(i).YData(1)) * ones(size(H(i).XData))-100)/200,0)*200+100;
end
colormap(icey);
caxis([-1500 1500])
cb = colorbar;
xlabel('Easting')
ylabel('Northing')
zlabel('Elevation [m]')
cb.Label.String = 'Elevation [m]';
f = gca;
f.XAxis.FontSize = ftsize;
f.YAxis.FontSize = ftsize;
f.ZAxis.FontSize = ftsize;
cb.FontSize = ftsize;
trisurf(t,xy(:,1),xy(:,2),h_s_init(xy(:,1),xy(:,2)),...
       'edgecolor','none','facecolor','interp','facealpha',.9)
hold on
colormap(icey);
daspect([1 1 1/50])
cb = colorbar;
caxis([-1500 1500])
xlabel('Easting')
ylabel('Northing')
zlabel('Elevation [m]')
cb.Label.String = 'Elevation [m]';
f = gca;
f.XAxis.FontSize = ftsize;
f.YAxis.FontSize = ftsize;
f.ZAxis.FontSize = ftsize;
cb.FontSize = ftsize;
sgtitle('Bed and Surface','FontSize',ftsize+2)
set(gcf,'Renderer','Painter')
view(-148.1119,22.5274)
