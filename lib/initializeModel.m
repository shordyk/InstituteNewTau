% Initialize Model parameters, must be done after loading inputs

%% Model Vars
nN = size(xy,1);                     %nN: number of nodes
nE = size(t,1);                      %nE: number of elements
b = unique(boundedges(xy,t));        %b:  boundary node numbers
nB = size(b,1);                      %nB: number of boundary nodes
e = [t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];
e = sort(e,2);
[foo,~,ifoo] = unique(e,'rows');
eB = foo(accumarray(ifoo,1) == 1,:); %eB: boundary edges
res = 1; %start residual arb high number

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
enhance = ones(nE,1);       % Thermal enhancement factor [ ]
E_man = 1; %(1/6)^(-1/3);       % Manual Adjustment to thermocoupling enhancement

% initialize to zero velocity case [m/s] 
u = zeros(size(xy(:,1))); 
v = u;

dz = .1;  %vertical resolution of thermal depth profiles (frac of H) [ ]

%% Import bedMachine data and smooth 
% import data on higher resolution square grid that is larger than the model grid by
% 'overgrab'
overgrab = 20;
xi = xmin-dx*overgrab:dx/2:xmax+dx*overgrab;
yi = ymin-dx*overgrab:dx/2:ymax+dx*overgrab;
[Xi,Yi] = meshgrid(xi,yi);

% Raw fields
bm_b =  bedmachine_interp('bed',Xi,Yi);
bm_s =  bedmachine_interp('surface',Xi,Yi);

% Smoothing
% Numerator is the window we're smoothing over in [m], spacing of these grids
% is actually dx/2 for bm_X grids hence the extra "*2".

smoothbed = imgaussfilt(bm_b,2.5e3*2/dx);
smoothsurf = imgaussfilt(bm_s,10e3*2/dx);

% smoothbed = sgolayfilt(bm_b,2,2*floor(10e3/dx)+1);
% smoothsurf = sgolayfilt(bm_s,2,2*floor(10e3/dx)+1);

% rectify rock above ice issue, force that ice is non-zero thickness everywhere
smoothsurf(smoothbed > smoothsurf) = smoothbed(smoothbed > smoothsurf) + 1; %Pe and Br = 0 result in NaNs


%% Build bed and surf, correct for thinning and floatation
h_real =@(x,y) interp2(xi,yi,bm_s-bm_b,x,y);
% rock_mask =@(x,y) interp2(xi,yi,rock,x,y,'nearest');
h_b_init =@(x,y) interp2(xi,yi,smoothbed,x,y);
h_s_init =@(x,y) interp2(xi,yi,smoothsurf,x,y);
% phi_init =@(x,y) rho/rho_w*h_s_init(x,y) + (rho_w-rho)/rho_w*h_b_init(x,y); %hydropotential per unit water weight
clear bm_b bm_s;
h = subplus(h_s_init(xy(:,1),xy(:,2)) - h_b_init(xy(:,1),xy(:,2))); %h: smoothed ice thickness [m]
h_re = subplus(h_real(xy(:,1),xy(:,2))); % non-smoothed ice thickness [m]

% %% Define a few globals vars
% phi_max = max(max(phi_init(xy(:,1),xy(:,2))));
% phi_min = min(min(phi_init(xy(:,1),xy(:,2))));

%% Mask Between rock/sed (0 is rock, 1 is sed)

rockSed = double(Yi > -(Xi + 5e5)/1.7); %divide between rock/sed midstream
% rockSed2 = heaviside(Yi + (Xi + 5e5)/1.7);
rockSedMask = griddedInterpolant(Xi',Yi',rockSed','nearest');   
clear rockSed;

%% Create vectors of bed/surface for numerical solving
h_s = h_s_init(xy(:,1),xy(:,2));
h_b = h_b_init(xy(:,1),xy(:,2));