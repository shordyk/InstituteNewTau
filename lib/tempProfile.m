function [t_z] = tempProfile(ep_dot, ep_star, Pe, Br, La, T_s, T_m, dz)
% [t_z] = tempProfile(ep_dot, ep_star, Pe, Br, La, T_s, T_m, dz)
% returns temp profile for region. Height dz in faction of H
    z = (0:dz:1).*ones(size(ep_dot,1),1.0/dz+1);
    T_s = T_s.*ones(size(z));
    t_z = T_m*ones(size(z));
    Pe = Pe.*ones(size(z));
    Br = Br.*ones(size(z));
    La = La.*ones(size(z));
    xi = temperateTransition(ep_dot, ep_star, Pe, Br, La).*ones(size(ep_dot,1),1.0/dz+1);
    La = La.*(1-xi);       %reduce magnitude of advected term as no temp can advect in thermal zone.  
    t_z(z >= xi) = T_s(z >= xi) + (T_m-T_s(z >= xi)).*((Br(z >= xi) - La(z >= xi))./Pe(z >= xi)).*(1 - z(z >= xi) + 1./Pe(z >= xi).*...
        exp(Pe(z >= xi).*(xi(z >= xi)-1))-1./Pe(z >= xi).*exp(Pe(z >= xi).*(xi(z >= xi)-z(z >= xi))));
end