function [E] = depthIntEnhancement(T_z,A_0,dz)
%[A] = depthIntEnhancement(T_z, H, dz)
%returns depth averaged value of enhancement, following Minchew et al
%no lubrication of crystals is included, only temp viscous dependence.
    E = trapz(calcAfromT(T_z),2)/A_0*dz;
end