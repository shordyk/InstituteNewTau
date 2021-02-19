function [A] = calcAfromT(T)
% [A] = calcAfromT(T)
% Calculates prefactor A from temp input. From Cuffey & Patterson
% Works for scalar, vector, matrix inputs
    T_star = 263;    %[k] transition temp
    A_star = 3.5e-25; %[Pa^-3 s^-1] base prefactor
    E = 1;           %Enhancement
    Q1 = 6e4;        %[J mol^-1] low temp activation energy
    Q2 = 1.15e5;     %[J mol^-1] high temp activation energy
    R  = 8.314;      %[J mol^-1 K^-1] Gas constant
    
    A = zeros(size(T));
    A(T <= T_star) = A_star.*exp(-Q1/R.*(1./T(T <= T_star)-1/T_star));
    A(T >  T_star) = A_star.*exp(-Q2/R.*(1./T(T >  T_star)-1/T_star));
end