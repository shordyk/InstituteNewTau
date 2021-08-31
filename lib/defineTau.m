function [tau_c] = defineTau(str,x0)
% [tau_c] = defineTau(str,h_s_init,h_b_init,phi_init,phi_max,phi_min,[x0])
% returns tau and string for requested string-name of tau scenario. 

    opt = false;
    if(nargin == 2)
        opt = true;
    end
    if(str == "Uniform") % Uniform Plastic Bed, for reference, not in paper
     if(opt)
        yield_base = x0;
    else
        yield_base = 90.559e3;
     end
    
    tau_c = @(x,y,u,v) norms([u,v],2,2)... 
       .*yield_base; 
    elseif(str == "Mixed")  % also from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
%     load("ISSM/x2dThwaites.mat",'x2dThwaites','y2dThwaites','BasalDragThwaites');
%     uB = scatteredInterpolant(x2dThwaites,y2dThwaites,BasalDragThwaites,'linear','linear');
    if(opt)
        scale_p = abs(x0(1));
        scale_l = abs(x0(2));
    else 
        scale_p = .20;
        scale_l = .49;
    end
    xi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
    yi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
    tau  = ncread("ISSM/JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
    [xx,yy] = ndgrid(xi - 3072000,yi - 3072000);
    
    uB = griddedInterpolant(xx,yy,tau(:,:,21).*3.154e7./(measures_interp('speed',xx,yy)));
    uB_p = griddedInterpolant(xx,yy,tau(:,:,21));
    tau_c = @(x,y,u,v,grounded) ...
    scale_p .* ( heaviside(y + (x + 5e5)/1.7))  .* norms([u,v],2,2) .* subplus(uB_p(x,y)) + ...
    scale_l .* ( 1 - heaviside(y + (x + 5e5)/1.7))  .* pow_pos(norms([u,v],2,2),2) .* subplus(uB(x,y))+...
    1e8*heaviside(-(x-(-8.1e5)).^2 - (y-(3.85e5)).^2 + 1.4e5^2).* norms([u,v],2,2);

%     tau_c = @(x,y,u,v,grounded) ...
%         scale_p .* ( rockSedMask(x,y))  .* norms([u,v],2,2) .* subplus(uB_p(x,y)) + ...
%         scale_l .* ( 1 - rockSedMask(x,y))  .* pow_pos(norms([u,v],2,2),2) .* subplus(uB(x,y));
elseif(str == "Mixed_east")  % also from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
%     load("ISSM/x2dThwaites.mat",'x2dThwaites','y2dThwaites','BasalDragThwaites');
%     uB = scatteredInterpolant(x2dThwaites,y2dThwaites,BasalDragThwaites,'linear','linear');
    if(opt)
        scale_p = abs(x0(1));
        scale_l = abs(x0(2));
    else 
        scale_p = .50;
        scale_l = .7;
    end
    xi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
    yi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
    tau  = ncread("ISSM/JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
    [xx,yy] = ndgrid(xi - 3072000,yi - 3072000);
    
    uB = griddedInterpolant(xx,yy,tau(:,:,21).*3.154e7./(measures_interp('speed',xx,yy)));
    uB_p = griddedInterpolant(xx,yy,tau(:,:,21));
    tau_c = @(x,y,u,v) ...
    scale_p .* ( heaviside(y + (x + 5e5)/1.7))  .* norms([u,v],2,2) .* subplus(uB_p(x,y)) + ...
    scale_l .* ( 1 - heaviside(y + (x + 5e5)/1.7))  .* pow_pos(norms([u,v],2,2),2) .* subplus(uB(x,y))+...
    1e8*heaviside(-(x-(-8.1e5)).^2 - (y-(3.85e5)).^2 + 1.3e5^2).* norms([u,v],2,2);
    elseif(str == "Mixed_center")  % also from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
%     load("ISSM/x2dThwaites.mat",'x2dThwaites','y2dThwaites','BasalDragThwaites');
%     uB = scatteredInterpolant(x2dThwaites,y2dThwaites,BasalDragThwaites,'linear','linear');
    if(opt)
        scale_p = abs(x0(1));
        scale_l = abs(x0(2));
    else 
        scale_p = .50;
        scale_l = .7;
    end
    xi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
    yi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
    tau  = ncread("ISSM/JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
    [xx,yy] = ndgrid(xi - 3072000,yi - 3072000);
    
    uB = griddedInterpolant(xx,yy,tau(:,:,21).*3.154e7./(measures_interp('speed',xx,yy)));
    uB_p = griddedInterpolant(xx,yy,tau(:,:,21));
    tau_c = @(x,y,u,v,grounded) ...
    scale_p .* ( heaviside(y + (x + 5e5)/1.7))  .* norms([u,v],2,2) .* subplus(uB_p(x,y)) + ...
    scale_l .* ( 1 - heaviside(y + (x + 5e5)/1.7))  .* pow_pos(norms([u,v],2,2),2) .* subplus(uB(x,y))+...
    1e8*heaviside(-(x-(-8.1e5)).^2 - (y-(3.85e5)).^2 + 1.4e5^2).* norms([u,v],2,2);
    elseif(str == "ISSM")  % from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
    if(opt)
        scale = x0(1);
    else 
        scale = 1.1963;
    end
    xi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
    yi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
    tau  = ncread("ISSM/JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
    [xx,yy] = ndgrid(xi - 3072000,yi - 3072000);
    uB = griddedInterpolant(xx,yy,tau(:,:,21));
    
    tau_c = @(x,y,u,v) norms([u,v],2,2) .*... %Plastic
        subplus(uB(x,y))*scale;
    elseif(str == "ISSM_center")  % from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
    if(opt)
        scale = x0(1);
    else 
        scale = 1.125;
    end
    xi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
    yi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
    tau  = ncread("ISSM/JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
    [xx,yy] = ndgrid(xi - 3072000,yi - 3072000);
    uB = griddedInterpolant(xx,yy,tau(:,:,21));
    
    tau_c = @(x,y,u,v) norms([u,v],2,2) .*... %Plastic
        subplus(uB(x,y))*scale+...
        1e8*heaviside(-(x-(-8.1e5)).^2 - (y-(3.85e5)).^2 + 1.4e5^2).* norms([u,v],2,2);
    elseif(str == "ISSM_east")  % from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
    if(opt)
        scale = x0(1);
    else 
        scale = 1.125;
    end
    xi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
    yi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
    tau  = ncread("ISSM/JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
    [xx,yy] = ndgrid(xi - 3072000,yi - 3072000);
    uB = griddedInterpolant(xx,yy,tau(:,:,21));
    
    tau_c = @(x,y,u,v,grounded) norms([u,v],2,2) .*... %Plastic
        subplus(uB(x,y))*scale+...
        1e8*heaviside(-(x-(-8.1e5)).^2 - (y-(3.85e5)).^2 + 1.3e5^2).* norms([u,v],2,2);
    elseif(str == "Depth_center")  % from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
    if(opt)
        sed = x0(1);
        rock = x0(1);
    else 
        sed = 20e3;
        rock = 115e3;
    end
%         grid = load('gridInstitute3000.mat','xy','dx');
%         [xx,yy] = ndgrid(min(grid.xy(:,1)):dx:max(grid.xy(:,1)),min(grid.xy(:,2)):dx:max(grid.xy(:,2)));
%         bed
        bed = @(xx,yy) bedmachine_interp('bed',xx,yy);  
        tau_c = @(x,y,u,v) ... %Plastic
        sed * norms([u,v],2,2) +...
        (rock-sed) * (1-heaviside(-1000 - bed(x,y))).* norms([u,v],2,2) +...
        1e8  * heaviside(-(x-(-8.1e5)).^2 - (y-(3.85e5)).^2 + 1.4e5^2).* norms([u,v],2,2);
    elseif(str == "ISSM_Linear")  % also from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
%     load("ISSM/x2dThwaites.mat",'x2dThwaites','y2dThwaites','BasalDragThwaites');
%     uB = scatteredInterpolant(x2dThwaites,y2dThwaites,BasalDragThwaites,'linear','linear');
    xi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
    yi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
    tau  = ncread("ISSM/JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
    [xx,yy] = ndgrid(xi - 3072000,yi - 3072000);
    uB = griddedInterpolant(xx,yy,tau(:,:,21).*3.154e7./(measures_interp('speed',xx,yy)));
    
    tau_c = @(x,y,u,v) pow_pos(norms([u,v],2,2),2) .*... %linear sliding law
        subplus(uB(x,y))*1;
    else
        error("Invalid Tau_c Scenario String")
    end
end
