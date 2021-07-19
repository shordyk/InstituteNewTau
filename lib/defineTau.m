function [tau_c] = defineTau(str,h_s_init,h_b_init,phi_init,phi_max,phi_min,x0)
% [tau_c] = defineTau(str,h_s_init,h_b_init,phi_init,phi_max,phi_min,[x0])
% returns tau and string for requested string-name of tau scenario. 

    opt = false;
    if(nargin == 7)
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
   
    elseif(str == "Case 1")  % Overburden Plastic Bed, Case 1
    if(opt)
        scale = abs(x0(1));
        base = abs(x0(2));
    else
        scale = 43;
        base = 40e3;
    end
    tau_c = @(x,y,u,v) (subplus((h_s_init(x,y)-h_b_init(x,y))...% Overburden
       .*subplus(scale+base)))...   %tau values cannot be negative, give floor of 0
       .*norms([u,v],2,2); 
    
    elseif(str == "Case 2")  % Depth weakening, Case 2
    if(opt)
        scale = abs(x0(1));
        base = x0(2);
    else
        scale = 31;
        base = 125e3;
    end
    tau_c = @(x,y,u,v) (h_b_init(x,y)...% Basal elevation
       .*scale+base)... 
       .*norms([u,v],2,2);
   
    elseif(str == "Case 2.1")  % Depth weakening, with strengthening, Case 2.1
    if(opt)
        scale = abs(x0(1)); 
        base = x0(2); 
        rock_hi = 300e3;
        cutoff = -200;
    else
        scale = 31;
        base = 125e3;
        rock_hi = 250e3;
        cutoff = -200;
    end
    tau_c = @(x,y,u,v) ((h_b_init(x,y)...% Basal elevation
       .*scale + base) .* (1-heaviside(-cutoff + h_b_init(x,y))) ...
       + rock_hi*heaviside(-cutoff + h_b_init(x,y)))...
       .*norms([u,v],2,2); 
   
    elseif(str == "Case 2.2")  % Depth weakening w/ channels, Case 2.2
    load('Streams_BedMachine.mat','drain_dist'); %Comes from StreamRouting.m
    if(opt)
        
        scale = abs(x0(1));
        base = x0(2);
        yield_ramp = 6e3;
    else
        scale = 13;
        base = 96e3;
        yield_ramp = 6e3;
    end
    tau_c = @(x,y,u,v) (h_b_init(x,y)...% Basal elevation
       .*scale + base + yield_ramp.*drain_dist(x,y,3e3))...  
       .*(1 - inPoly(x,y,S))... % 0 strength in lake regions
       .*norms([u,v],2,2);
 
    elseif(str == "Case 3.2") % Case 3.2
    load('Streams_BedMachine.mat','drain_dist'); %Comes from StreamRouting.m
    if(opt)
        yield_ramp = 4e3; 
        yield_base = x0(1);
        phi_scale = x0(2);
    else
        yield_ramp = 4.3e3; 
        yield_base = 27e3;
        phi_scale = 103e3;
    end
    tau_c = @(x,y,u,v) norms([u,v],2,2)... 
       .*(1 - inPoly(x,y,S))... % 0 strength in lake regions
       .*(yield_ramp.*drain_dist(x,y,3e3) +...
         phi_scale.*(1-(phi_init(x,y)-phi_min)/(phi_max-phi_min)) + yield_base); 
    
    elseif(str == "Case 3") % Phi based Bed Strength, case 3
    if(opt)
        scale = x0(1);
        base =  x0(2);
    else
        scale = 101.5e3;
        base =  36.6e3;
    end
    tau_c = @(x,y,u,v) norms([u,v],2,2)... 
        .*(scale.*(1-(phi_init(x,y)-phi_min)/(phi_max-phi_min)) + base);
    
    elseif(str == "Case 3.1") % Case 3.1
    if(opt)
        scale = x0(1);
        base =  x0(2);
    else
        scale = 101.5e3;
        base =  36.6e3;
        rock_hi = 250e3;
        cutoff = -200;
    end
    tau_c = @(x,y,u,v) norms([u,v],2,2)... 
        .*((scale.*(1-(phi_init(x,y)-phi_min)/(phi_max-phi_min)) + base).* (1-heaviside(-cutoff + h_b_init(x,y))) +...;
        + rock_hi*heaviside(-cutoff + h_b_init(x,y)));
    elseif(str == "ISSM")  % also from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
%     load("ISSM/x2dThwaites.mat",'x2dThwaites','y2dThwaites','BasalDragThwaites');
%     uB = scatteredInterpolant(x2dThwaites,y2dThwaites,BasalDragThwaites,'linear','linear');
    xi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
    yi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
    tau  = ncread("ISSM/JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
    [xx,yy] = ndgrid(xi - 3072000,yi - 3072000);
    uB = griddedInterpolant(xx,yy,tau(:,:,21));
    
    tau_c = @(x,y,u,v,grounded) norms([u,v],2,2).*... %Plastic
        subplus(uB(x,y))*1;
    elseif(str == "ISSM_Linear")  % also from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
%     load("ISSM/x2dThwaites.mat",'x2dThwaites','y2dThwaites','BasalDragThwaites');
%     uB = scatteredInterpolant(x2dThwaites,y2dThwaites,BasalDragThwaites,'linear','linear');
    xi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
    yi   = ncread("ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
    tau  = ncread("ISSM/JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
    [xx,yy] = ndgrid(xi - 3072000,yi - 3072000);
    uB = griddedInterpolant(xx,yy,tau(:,:,21).*3.154e7./(measures_interp('speed',xx,yy)));
    
    tau_c = @(x,y,u,v,grounded) pow_pos(norms([u,v],2,2),2) .*... %Plastic
        subplus(uB(x,y))*1;
    else
        error("Invalid Tau_c Scenario String")
    end
end
