clc
clear
%% Initialization
% Physical parameters
p = 4/3;
rho = 900;
g = 9.8;
alpha = 2.4E-3;
f = rho*g*sin(alpha);

% Domain parameters
dx = 880;                     %dx: nominal grid spacing
vert_scale = 0.375;
dy = vert_scale*dx;
% Numerical experiment #1 (uniform strength)
% beta = 4.1E6;
% rock_trans = 15.1E3;
% sed_const = 23.95E3;
% sed_var = 0E3;
% sed_trans_l = -23E3;
% sed_trans_r = 72E3;
% ch_str = 0;
% ch_decay = 7E3;
% ch_loc = 0E3;
% Numerical experiment #2 (linear strength)
% beta = 4E6;
% rock_trans = 13.5E3;
% sed_const = 16.55E3;
% sed_var = 16E3;
% sed_trans_l = -23E3;
% sed_trans_r = 71E3;
% ch_str = 0;
% ch_decay = 7E3;
% ch_loc = 0E3;
% Numerical experiment #3 (all plastic)
% beta = 0;
% rock_trans = -23E3;
% sed_const = 17E3;
% sed_var = 15E3;
% sed_trans_l = -23E3;
% sed_trans_r = 71E3;
% ch_str = 500;
% ch_decay = 20E3;
% ch_loc = 13E3;
% Numerical experiment #4 (rock/sediment with channel)
beta_var = 5E6;
beta = 2.9E6;
rock_trans = 10E3;
sed_const = 19.4E3;
sed_var = 8E3;
sed_trans_l = -23E3;
sed_trans_r = 71E3;
ch_str = 24.5E3;
ch_decay = 5.9E3;
ch_loc = 10E3;
% Numerical experiment #4 (all plastic with channel)
% beta = 4.1E6;
% rock_trans = 0E3;
% sed_const = 23.95E3;
% sed_var = 0E3;
% sed_trans_l = -23E3;
% sed_trans_r = 72E3;
% ch_str = 500E3;
% ch_decay = 10E3;
% ch_loc = 13E3;

tau_c =@(x,y,u) (...
                 heaviside(sed_trans_l - x).*1E8.*abs(u) + ...
                 heaviside(x - sed_trans_l).*heaviside(rock_trans - x).*(beta_var*(x/80E3) + beta).^2.*pow_abs(u,5/2) + ...
                 heaviside(x - rock_trans).*(sed_var*((80E3-x)/80E3) + ...
                                             ch_str*(exp(-abs(x-ch_loc)/ch_decay)) + ...
                                             sed_const + ...
                                             1E8*smearedHeavi(x - sed_trans_r,0)).*abs(u)... 
                 ).*heaviside((450+dy/2) - y);


%% Mesh Generation
pv = 40E3.*[0,1;2,1;2,.93;.3,.9;0,.925;-.75,.925;-.75,1;0,1];
[xy,t] = distmesh2d(@dpoly,@huniform,dx,40E3.*[-.85,.75;2.1,1.1],pv,pv);
xy(:,2) = vert_scale*((xy(:,2)-40E3) - min(xy(:,2)-40E3));

nN = size(xy,1);                     %nN: number of nodes
nE = size(t,1);                      %nE: number of elements
b = unique(boundedges(xy,t));        %b:  boundary node numbers
nB = size(b,1);                      %nB: number of boundary nodes

e = [t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];
e = sort(e,2);
[foo,~,ifoo] = unique(e,'rows');
eB = foo(accumarray(ifoo,1) == 1,:); %eB: boundary edges

T = 25/1500*(1500-xy(:,2)) + 248;
T_old = 25/1500*(1500-xy(:,2)) + 248;
res = 1;
while(res > 1E-3)
T = .3*T + (1-.3)*T_old;

A = 3.5e-25; %A: preexponential constant[1/s*Pa^3]
E = 1; %E: enhancement factor[]
Q_h = 0*xy(:,1); %Q_h: melting point adjusted activation energy[J/mol]
R = 8.314; %R: gas constant[J/K*mol]
Tstar = 263.15; %Tstar: activation threshold[K] (-10C)
p0 = 7e-8; %p0: pressure heating coef[K/Pa]
P = rho*g*(1500-xy(:,2)); %P: hydrostatic pressure[Pa]
T_h = T + p0*P; %T_h: melting point adjusted temperature[K]
Tstar_h = Tstar*(0*xy(:,1) + 1) + p0*P; %Tstar_h: adjusted activation threshold[K]
Q2 = 60e3; %Q2: lower activation energy[J/mol]
Q3 = 115e3; %Q3: higher activation energy[J/mol]
Q_h(heaviside(T_h-Tstar_h)==0) = Q2;
Q_h(heaviside(T_h-Tstar_h)==1) = Q3;
a = 1/2*(A*E*exp((-Q_h./R).*((1./T_h)-(1./Tstar_h)))).^(-1/3); %a = A^(-1/n)

%% Build System
A = sparse(nE,nN);
B = sparse(nE,nN);
D = zeros(nE,nN);
F = zeros(1,nN);

tau_area = zeros(nE,1);
for E = 1:nE  % integration over each element
  nodes = t(E,:);
  xyE = [ones(3,1),xy(nodes,:)];
  Area = abs(det(xyE))/2;
  tau_area(E) = Area;
  C = inv(xyE);
  A_E = C(2,:);
  B_E = C(3,:);
  D_E = 1/3*ones(1,3);
  F_E = f*Area/3*ones(1,3);
  A(E,nodes) = A(E,nodes) + A_E;
  B(E,nodes) = B(E,nodes) + B_E;
  D(E,nodes) = D(E,nodes) + D_E;
  F(nodes) = F(nodes) + F_E;
end

b_dx = zeros(size(b,1),1);
for N = 1:size(b,1) % integration over each boundary node
    edges = eB(sum(eB == b(N),2) == 1,:);
    for i = 1:size(edges,1)
        edge = edges(i,:);
        b_dx(N) = b_dx(N) + .5*sqrt(sum(diff(xy(edge,:),1).^2,2));
    end
end


%% Solve

cvx_begin
cvx_quiet true
    variables u(nN)
    obj = (1/p)*sum((D*a).*tau_area.*pow_pos(norms([A*u,B*u],2,2),p)) - F*u + sum(b_dx.*tau_c(xy(b,1),xy(b,2),u(b)));
    subject to
        u >= 0;
        u(xy(:,1) > 80E3 - dx/2) == 0;
        u(xy(:,1) < -30E3 + dx/2) == 0;
    minimize(obj)
cvx_end

u_x = A*u;
u_y = B*u;
mu = (D*a).*(sqrt((u_x).^2 + (u_y).^2)).^(1/2);
tau_E = sqrt((mu.*u_x).^2 + (mu.*u_y).^2); %tau_E: effective stress[Pa]
epsilon_E = sqrt((u_x/2).^2 + (u_y/2).^2); %epsilon_E: effective strain rate[1/s]
f_therm = 2*tau_E.*epsilon_E;

k1 = 9.828; %k1: conductivity preexponential[W/m*K]
k2 = 5.7; %k2: conductivity postexponential[1/K]
k = k1*exp(-k2*1e-3.*T);
tau_T =@(y) (55./(1500-y.*heaviside((450+dy/2) - y))).*heaviside((450+dy/2) - y);

F_therm = zeros(1,nN);
for E = 1:nE  % integration over each element
  nodes = t(E,:);
  xyE = [ones(3,1),xy(nodes,:)];
  Area = abs(det(xyE))/2;
  F_therm_E = f_therm(E)*Area/3*ones(1,3);
  F_therm(nodes) = F_therm(nodes) + F_therm_E;
end

T_old = T;
cvx_begin
cvx_quiet true
    variables T(nN)
    obj = sum((D*k).*tau_area.*pow_pos(norms([A*T,B*T],2,2),2)) - 2E10*F_therm*T - sum(b_dx.*k(b).*tau_T(xy(b,2)).*T(b));
    subject to
        T <= 273;
        T(xy(:,2) > 1500-dx/10) == 248;
    minimize(obj)
cvx_end

res = norm(T-T_old)/norm(T);
disp(res)

end
%% Visualization
u = u*pi*1E7;
T = T - 273;

trisurf(t,xy(:,1),xy(:,2),u,u,...
       'edgecolor','none','facecolor','interp');
hold on
%tricontour(t,xy(:,1),xy(:,2),u,...
%           [0,.001,.002,.003,.004,.005,.006,.007,.008,.009,.010,.011,.012,.013,.014,.015,.016]);
plot3(xy(b(u(b) < 1e-6),1),xy(b(u(b) < 1e-6),2),max(u)+0*b(u(b) < 1e-6),'r.','MarkerSize',20)
plot3(xy(b(b_dx==0),1),xy(b(b_dx==0),2),max(u)+0*b(b_dx==0),'b.','MarkerSize',20)
view(2), colorbar

load vel_profile_full.mat
figure('Position',[10 10 500 400])

subplot(2,1,1)
plot(profile_path-30.5E3,profile_cross,'LineWidth',3)
hold on
plot(xy(xy(:,2) > 1500-dx/10,1),u(xy(:,2) > 1500-dx/10),'LineWidth',3)
axis([-30E3,80E3,0,400])

subplot(2,1,2)
plot(profile_path-30.5E3,gradient(profile_cross),'LineWidth',3)
hold on
plot(xy(xy(:,2) > 1500-dx/10,1),gradient(u(xy(:,2) > 1500-dx/10)),'LineWidth',3)
axis([-30E3,80E3,-80,50])
setFontSize(16)

figure('Position',[10 10 1000 200])
load iceColorMap
colormap(iceColorMap)
% tricontf(xy(:,1),xy(:,2),t,T,...
%          [0,-5,-10,-15,-20,-25]);
trisurf(t,xy(:,1),xy(:,2),T,T,...
       'edgecolor','none','facecolor','interp');
axis off
view(2)

figure('Position',[10 10 500 150])
yyaxis left
plot(sed_trans_l:1E2:rock_trans,beta+beta_var*(sed_trans_l:1E2:rock_trans)/80E3,'LineWidth',3)
axis([-30E3,80E3,1E6,4E6])
yyaxis right
plot(rock_trans:1E2:sed_trans_r,sed_var*((80E3-(rock_trans:1E2:sed_trans_r))/80E3) + ch_str*(exp(-abs((rock_trans:1E2:sed_trans_r)-ch_loc)/ch_decay)) + sed_const,'LineWidth',3)
%axis([-30E3,80E3,15E3,32E3])
axis([-30E3,80E3,15E3,55E3])
setFontSize(16)

%figure 
%bedmap2_profile(profile_lat(profile_path < 110E3),profile_lon(profile_path < 110E3))