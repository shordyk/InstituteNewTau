function [tau_c] = defineTau(str,rockSedMask,x0)
% [tau_c] = defineTau(str,h_s_init,h_b_init,phi_init,phi_max,phi_min,[x0])
% returns tau and string for requested string-name of tau scenario. 

    opt = false;
    if(nargin == 3)
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
        scale = x0(1);
    else 
        scale = 1.1963;
    end
    xi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
    yi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
    tau  = ncread("ISSM/JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
    [xx,yy] = ndgrid(xi - 3072000,yi - 3072000);
    
    uB = griddedInterpolant(xx,yy,tau(:,:,21).*3.154e7./(measures_interp('speed',xx,yy)));
    uB_p = griddedInterpolant(xx,yy,tau(:,:,21));
    
    tau_c = @(x,y,u,v,grounded) ...
        ( 1 + rockSedMask(x,y))  .* pow_pos(norms([u,v],2,2),2) .* subplus(uB(x,y));
    
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
    
    tau_c = @(x,y,u,v,grounded) norms([u,v],2,2) .*... %Plastic
        subplus(uB(x,y))*scale;
    
    elseif(str == "ISSM_Linear")  % also from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
%     load("ISSM/x2dThwaites.mat",'x2dThwaites','y2dThwaites','BasalDragThwaites');
%     uB = scatteredInterpolant(x2dThwaites,y2dThwaites,BasalDragThwaites,'linear','linear');
    xi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
    yi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
    tau  = ncread("ISSM/JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
    [xx,yy] = ndgrid(xi - 3072000,yi - 3072000);
    uB = griddedInterpolant(xx,yy,tau(:,:,21).*3.154e7./(measures_interp('speed',xx,yy)));
    
    tau_c = @(x,y,u,v,grounded) pow_pos(norms([u,v],2,2),2) .*... %linear sliding law
        subplus(uB(x,y))*1;
    else
        error("Invalid Tau_c Scenario String")
    end
end
