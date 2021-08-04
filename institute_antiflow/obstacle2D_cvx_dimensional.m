clc
clear
%% Initialization
dx = 1E3;                     %dx: nominal grid spacing

% Physical parameters
a = 1.4286E24;
f = 29.8E3;
tau_c =@(x,y,u) (...
                 heaviside(20E3 - x).*1E6.*pow_abs(u,5/2) + ...
                 heaviside(x - 20E3).*(4E3 + 1E25*heaviside(x - 60E3)).*abs(u)...
                 ).*heaviside(-1150 - y);
% tau_c =@(x,y,u) (...
%                  heaviside(20E3 - x).*4E6.*pow_abs(u,5/2) + ...
%                  heaviside(x - 20E3).*(5E3 + 35E3*heaviside(x - 60E3)).*abs(u)...
%                  ).*heaviside(-1150 - y);

%% Mesh Generation
pv = [0,0;65E3,0;65E3,-1200;10E3,-1700;0,-1200;0,0];
[xy,t] = distmesh2d(@dpoly,@huniform,dx,[-1E3,-2700;66E3,1000],pv,pv);
nN = size(xy,1);                     %nN: number of nodes
nE = size(t,1);                      %nE: number of elements
b = unique(boundedges(xy,t));        %b:  boundary node numbers
nB = size(b,1);                      %nB: number of boundary nodes

e = [t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];
e = sort(e,2);
[foo,~,ifoo] = unique(e,'rows');
eB = foo(accumarray(ifoo,1) == 1,:); %eB: boundary edges


%% Build System
A = sparse(nE,nN);
B = sparse(nE,nN);
D = zeros(1,nB);
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
  F_E = Area/3*ones(1,3);
  A(E,nodes) = A(E,nodes) + A_E;
  B(E,nodes) = B(E,nodes) + B_E;
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
    variables u(nN)
    obj = a*3/4*sum(tau_area.*pow_pos(norms([A*u,B*u],2,2),4/3)) - f*F*u + sum(b_dx.*tau_c(xy(b,1),xy(b,2),u(b)));
    subject to
        u >= 0;
    minimize(obj)
cvx_end

%% Visualization
u = 3.154E7*u;

trisurf(t,xy(:,1),xy(:,2),u,u,...
       'edgecolor','none','facecolor','interp');
hold on
%tricontour(t,xy(:,1),xy(:,2),u,...
%           [0,.001,.002,.003,.004,.005,.006,.007,.008,.009,.010,.011,.012,.013,.014,.015,.016]);
plot3(xy(b(u(b) < 1e-6),1),xy(b(u(b) < 1e-6),2),max(u)+0*b(u(b) < 1e-6),'r.','MarkerSize',20)
plot3(xy(b(b_dx==0),1),xy(b(b_dx==0),2),max(u)+0*b(b_dx==0),'b.','MarkerSize',20)
view(2), axis equal, colorbar

figure
plot(xy(xy(:,2) > -50,1),u(xy(:,2) > -50),'LineWidth',3)
hold on
load vel_profiles.mat
j = 6;
plot(profile_path(:,j),profile_cross(:,j),'LineWidth',3)

figure
plot(xy(xy(:,2) > -50,1),gradient(u(xy(:,2) > -50)))
hold on
plot(profile_path(:,j),gradient(profile_cross(:,j)),'LineWidth',3)