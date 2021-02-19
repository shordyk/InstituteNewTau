function [xi_h] = temperateTransition(ep_dot, ep_star, Pe, Br, La)
% [xi_h] = temperateTransition(ep_dot, ep_star, Pe, Br) 
% Finds xi as a fraction of H given ep, Pe, Br
    xi_h = zeros(size(ep_dot));
    boo = ep_dot > ep_star;
    xi_h(boo) = ... 
        1 - Pe(boo)./(Br(boo)-La(boo)) - 1./Pe(boo).*(1 + lambertw(-exp(-Pe(boo).^2./(Br(boo)-La(boo))-1)));
end