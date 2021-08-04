clear; close all
load("x2dThwaites.mat",'x2dThwaites','y2dThwaites','BasalDragThwaites');
load ../workingGrid_ex5.mat
uB = scatteredInterpolant(x2dThwaites,y2dThwaites,BasalDragThwaites);

figure
scatter(x2dThwaites,y2dThwaites,[],BasalDragThwaites)
hold on
scatter(xy(:,1),xy(:,2),'k.')
colorbar

figure
trisurf(t,xy(:,1),xy(:,2),uB(xy(:,1),xy(:,2)),'edgecolor','none')
colorbar
view(2)
axis equal